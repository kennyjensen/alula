// SPDX-License-Identifier: GPL-2.0-or-later
// Sorenson, NASA TM-81198, Eq. (18): boundary-normal second derivatives
// from a prescribed first derivative and two interior points. We extend
// the cubic-Hermite formula to nonuniform mass coordinates and prescribe
// only the angle. Normal speed is determined by Laplacian(eta)=0, rather
// than imposing both angle and stand-off using a transverse Poisson source.
// This is a derived specialization, not original GRAPE or MSET source.
const dot = (a, b) => a.x * b.x + a.y * b.y;
const finitePoint = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);

// Sorenson Eq. (17): damp and limit the control increment. The sensitivity
// is with respect to the requested value while the previous iterate is fixed.
export function relaxBoundaryStretch({ requested, previous, relaxation = .3, changeLimit = 1 }) {
  if (![requested, previous, relaxation, changeLimit].every(Number.isFinite) || !(relaxation > 0 && relaxation <= 1 && changeLimit > 0))
    throw new Error('Invalid boundary-control update.');
  const difference = requested - previous, relaxed = relaxation * Math.abs(difference), limit = changeLimit * Math.max(Math.abs(previous), 1);
  if (![difference, relaxed, limit].every(Number.isFinite)) throw new Error('Unresolved boundary-control increment.');
  return { value: previous + Math.sign(difference) * Math.min(relaxed, limit), sensitivity: relaxed < limit ? relaxation : 0,
    limited: relaxed >= limit };
}

export function reconstructOrthogonalBoundary({ point, tangent, secondTangent, firstInterior, secondInterior,
  firstDistance: a, secondDistance: b, normalSign = 1 }) {
  if (![point, tangent, secondTangent, firstInterior, secondInterior].every(finitePoint)
    || ![a, b].every(Number.isFinite) || !(a > 0 && b > a) || ![1, -1].includes(normalSign))
    throw new Error('Orthogonal boundary reconstruction requires finite derivatives, two ordered inward distances and an oriented normal.');
  const gamma = dot(tangent, tangent), speed = Math.sqrt(gamma);
  if (!(speed > 0) || !Number.isFinite(gamma)) throw new Error('Singular boundary tangent.');
  const normal = { x: -normalSign * tangent.y / speed, y: normalSign * tangent.x / speed };
  // For r(rho)=r0+rho*rhoDerivative+rho²*rhoSecond/2+rho³*c,
  // K - m*rhoDerivative equals rhoSecond exactly, including unequal steps.
  const A = 2 * b / (a * a * (b - a)), B = -2 * a / (b * b * (b - a)), m = 2 * (1 / a + 1 / b);
  const K = Object.fromEntries(['x', 'y'].map(key => [key,
    A * (firstInterior[key] - point[key]) + B * (secondInterior[key] - point[key])]));
  const curvature = dot(secondTangent, normal) / gamma, Kn = dot(K, normal);
  // The normal projection with beta=0, Q=0 is
  // curvature*lambda² - m*lambda + Kn = 0.
  // Select the positive root continuous with the flat-boundary limit.
  const discriminantRatio = 1 - 4 * curvature * (Kn / m) / m;
  if (!Number.isFinite(Kn) || !Number.isFinite(m) || !Number.isFinite(discriminantRatio))
    throw new Error('No resolved positive harmonic-mass normal-speed branch at this boundary.');
  if (!(Kn > 0) || !(discriminantRatio > 64 * Number.EPSILON))
    throw Object.assign(new Error('No resolved positive harmonic-mass normal-speed branch at this boundary.'), {
      code: 'orthogonal-normal-branch', diagnostics: { Kn, m, curvature, discriminantRatio },
    });
  const normalSpeed = (2 * Kn / m) / (1 + Math.sqrt(discriminantRatio));
  const alpha = normalSpeed * normalSpeed;
  const stretch = -dot(secondTangent, tangent) / gamma - dot(K, tangent) / alpha;
  const normalSecond = Object.fromEntries(['x', 'y'].map(key => [key, K[key] - m * normalSpeed * normal[key]]));
  if (!(normalSpeed > 0) || ![normalSpeed, alpha, stretch].every(Number.isFinite) || !finitePoint(normalSecond))
    throw new Error('Unresolved orthogonal boundary control.');
  // Exact derivative of the selected root, including its normal component.
  // These two interior-point sensitivities permit implicit line updates;
  // explicitly lagging a coefficient proportional to 1/a² can be unstable.
  const rootSlope = m * Math.sqrt(discriminantRatio), tangentialK = dot(K, tangent);
  const gradient = Object.fromEntries(['x', 'y'].map(key => [key,
    -tangent[key] / alpha + 2 * (tangentialK / alpha) * normal[key] / normalSpeed / rootSlope]));
  const interiorDerivatives = [A, B].map(coefficient => Object.fromEntries(['x', 'y'].map(key => [key, coefficient * gradient[key]])));
  if (interiorDerivatives.some(d => !finitePoint(d))) throw new Error('Unresolved orthogonal control sensitivity.');
  return { stretch, normalSpeed, normal, normalSecond, curvature, gamma, discriminantRatio, interiorDerivatives };
}

const sampleBoundary = (nodes, xi, eta, discretization, side, i) => {
  const sign = side === 'lower' ? 1 : -1, j = sign === 1 ? 0 : eta.length - 1;
  const row = nodes[i], a = xi[i] - xi[i - 1], b = xi[i + 1] - xi[i];
  const derivatives = ['x', 'y'].map(key => {
    const left = row[j][key] - nodes[i - 1][j][key], right = nodes[i + 1][j][key] - row[j][key];
    return { key, first: discretization === 'giles-1985' ? (left + right) / (a + b) : (b * left / a + a * right / b) / (a + b),
      second: 2 * (right / b - left / a) / (a + b) };
  });
  try {
    return reconstructOrthogonalBoundary({ point: row[j], tangent: Object.fromEntries(derivatives.map(d => [d.key, d.first])),
      secondTangent: Object.fromEntries(derivatives.map(d => [d.key, d.second])), firstInterior: row[j + sign], secondInterior: row[j + 2 * sign],
      firstDistance: sign * (eta[j + sign] - eta[j]), secondDistance: sign * (eta[j + 2 * sign] - eta[j]), normalSign: sign });
  } catch (error) {
    throw Object.assign(new Error(`${side} boundary station ${i}: ${error.message}`, { cause: error }),
      error.code === 'orthogonal-normal-branch' ? { code: error.code, diagnostics: { ...error.diagnostics, side, station: i } } : {});
  }
};

const cornerStations = (corners, sides, nx) => {
  if (!corners || typeof corners !== 'object' || Array.isArray(corners) || Object.keys(corners).some(side => !sides.includes(side))
    || Object.values(corners).some(row => !Array.isArray(row) || row.some((i, k) => !Number.isInteger(i) || i < 2 || i > nx - 2
      || k && i - row[k - 1] < 2)))
    throw new Error('Corner controls require ordered isolated stations with one smooth interior station on each side.');
  return Object.fromEntries(sides.map(side => [side, new Set(corners[side] ?? [])]));
};

const controlUpdateState = (previous, sides, nx, relaxation, changeLimit) => {
  if (previous === undefined) return null;
  relaxBoundaryStretch({ requested: 0, previous: 0, relaxation, changeLimit });
  if (!previous || typeof previous !== 'object' || Array.isArray(previous) || Object.keys(previous).some(side => !sides.includes(side))
    || sides.some(side => !Array.isArray(previous[side]) || previous[side].length !== nx + 1 || !previous[side].every(Number.isFinite)))
    throw new Error('Previous boundary controls must match every selected side and station.');
  return { previous: Object.fromEntries(sides.map(side => [side, previous[side].slice()])), relaxation, changeLimit };
};

// Angle conditions apply only to explicitly selected boundary stations.
// Inactive stations retain the prescribed background control. A cut or a
// body endpoint can thus be distinguished from a smooth solid-wall segment.
const controlStations = (activeStations, sides, nx, corners) => {
  if (activeStations !== undefined && (!activeStations || typeof activeStations !== 'object' || Array.isArray(activeStations)
    || Object.keys(activeStations).some(side => !sides.includes(side))
    || Object.values(activeStations).some(row => !Array.isArray(row) || row.length !== nx + 1
      || Array.from(row).some(value => typeof value !== 'boolean'))))
    throw new Error('Active boundary stations must be boolean arrays matching selected sides and streamwise nodes.');
  const selected = Object.fromEntries(sides.map(side => [side,
    Array.from({ length: nx + 1 }, (_, i) => Boolean(i && i < nx && (activeStations?.[side]?.[i] ?? true)))]));
  for (const side of sides) for (const i of corners[side])
    if (selected[side][i] && (!selected[side][i - 1] || !selected[side][i + 1]))
      throw new Error('An active corner control requires both neighboring boundary stations active.');
  return selected;
};

// Sorenson, printed p.8: replace a sharp-corner control by the average of
// its two neighboring controls. Our Q=0 specialization averages boundary
// F or P, rather than the original program's two P/Q amplitudes. The caller identifies
// topological corners explicitly; no angle threshold or grid repair is used.
const controlledBoundary = (nodes, xi, eta, discretization, side, i, corners, update, sourceForm) => {
  const field = sourceForm === 'poisson' ? 'poisson' : 'stretch';
  const sample = station => {
    let control = sampleBoundary(nodes, xi, eta, discretization, side, station);
    if (sourceForm === 'poisson') {
      // On the orthogonal boundary J²=alpha*gamma, so
      // P = |grad xi|² F = F/gamma. The fixed boundary metric has
      // no derivative with respect to either interior point.
      const poisson = control.stretch / control.gamma;
      const interiorDerivatives = control.interiorDerivatives.map(d => ({ x: d.x / control.gamma, y: d.y / control.gamma }));
      if (!Number.isFinite(poisson) || interiorDerivatives.some(d => !finitePoint(d))) throw new Error('Unresolved boundary Poisson source.');
      control = { ...control, poisson, interiorDerivatives, derivativeOf: 'poisson' };
    }
    if (!update) return control;
    const result = relaxBoundaryStretch({ ...update, requested: control[field], previous: update.previous[side][station] });
    return { ...control, [sourceForm === 'poisson' ? 'requestedPoisson' : 'requestedStretch']: control[field], [field]: result.value, limited: result.limited,
      interiorDerivatives: control.interiorDerivatives.map(d => ({ x: result.sensitivity * d.x, y: result.sensitivity * d.y })) };
  };
  if (!corners.has(i)) return sample(i);
  // Average AFTER the individual Eq. (17) updates, as in the source's
  // sharp-corner treatment. Limiting and averaging do not commute.
  const neighbors = [i - 1, i + 1].map(k => ({ station: k, control: sample(k) }));
  return { [field]: .5 * (neighbors[0].control[field] + neighbors[1].control[field]), corner: true, neighbors };
};

export function createOrthogonalBoundaryControl({ nodes, xi, eta, discretization = 'giles-1985', sides = ['lower', 'upper'], corners = {},
  activeStations, previous, relaxation = .3, changeLimit = 1, sourceForm = 'metric-stretch' }) {
  const ordered = (a, min) => Array.isArray(a) && a.length >= min && a[0] === 0 && a.at(-1) === 1
    && a.every((v, i) => Number.isFinite(v) && (!i || v > a[i - 1]));
  if (!ordered(xi, 3) || !ordered(eta, 4) || !['giles-1985', 'quadratic'].includes(discretization) || !['metric-stretch', 'poisson'].includes(sourceForm)
    || !Array.isArray(sides) || !sides.length || sides.some(side => !['lower', 'upper'].includes(side))
    || new Set(sides).size !== sides.length || !Array.isArray(nodes) || nodes.length !== xi.length
    || nodes.some(row => !Array.isArray(row) || row.length !== eta.length || !row.every(finitePoint)))
    throw new Error('Boundary-angle control requires a finite grid, ordered normalized coordinates and two interior rows.');
  const selectedCorners = cornerStations(corners, sides, xi.length - 1);
  const active = controlStations(activeStations, sides, xi.length - 1, selectedCorners);
  const update = controlUpdateState(previous, sides, xi.length - 1, relaxation, changeLimit);
  const result = {};
  for (const side of sides) {
    result[side] = nodes.map((row, i) => {
      // Ends have no elliptic PDE. Marked corners average adjacent smooth
      // controls without evaluating an undefined tangent at the corner.
      if (!active[side][i]) return null;
      return controlledBoundary(nodes, xi, eta, discretization, side, i, selectedCorners[side], update, sourceForm);
    });
  }
  return { ...result, coordinateEquations: { eta: 'Laplace', boundaryAngle: Math.PI / 2, sourceForm },
    corners: Object.fromEntries(sides.map(side => [side, [...selectedCorners[side]]])), activeStations: active,
    source: 'Sorenson Eq. (18), nonuniform cubic-Hermite reconstruction with Q=0 and angle only', exactMsetSpacingLaw: false };
}

// Extend the change from a fixed background in mass fraction, using either
// linear weights or two matched exponential tails.
// The value is evaluated from the current grid, not a lagged control state.
// By default derivatives refer to the current line. The optional fourth
// evaluate argument selects another interior row for a coupled-line solve.
export function createOrthogonalBoundaryFeedback({ nodes, xi, eta, background, sides = ['lower', 'upper'], discretization = 'giles-1985', decay, corners = {},
  activeStations, previous, relaxation = .3, changeLimit = 1, sourceForm = 'metric-stretch' }) {
  const active = createOrthogonalBoundaryControl({ nodes, xi, eta, sides, discretization, corners, activeStations,
    previous, relaxation, changeLimit, sourceForm }).activeStations;
  if (!Array.isArray(background) || background.length !== xi.length
    || background.some(row => !Array.isArray(row) || row.length !== eta.length || !row.every(Number.isFinite)))
    throw new Error('Orthogonal control requires a finite background matching the grid.');
  if (decay !== undefined && (!decay || !['lower', 'upper'].every(side => Number.isFinite(decay[side]) && decay[side] > 0)))
    throw new Error('Boundary decay rates must be finite and positive in normalized mass coordinates.');
  const prescribed = background.map(row => row.slice()), selected = sides.slice(), x = xi.slice(), e = eta.slice(), nt = eta.length - 1;
  const selectedCorners = cornerStations(corners, selected, xi.length - 1), hasCorners = Object.values(selectedCorners).some(row => row.size);
  const update = controlUpdateState(previous, selected, xi.length - 1, relaxation, changeLimit);
  const field = sourceForm === 'poisson' ? 'poisson' : 'stretch';
  // Sorenson's two exponential boundary tails. Solve their 2x2 amplitude
  // matching system, so one boundary's tail cannot alter the other's value.
  // Stable differences also recover linear weights as both rates tend to 0.
  const weights = e.map(v => decay ? {
    lower: Math.exp(-decay.lower * v) * Math.expm1(-(decay.lower + decay.upper) * (1 - v)) / Math.expm1(-decay.lower - decay.upper),
    upper: Math.exp(-decay.upper * (1 - v)) * Math.expm1(-(decay.lower + decay.upper) * v) / Math.expm1(-decay.lower - decay.upper)
  } : { lower: 1 - v, upper: v });
  // A boundary control depends on five points, not on the receiving eta
  // line or derivative row. Reuse its reconstruction across those queries.
  // Compare coordinates rather than grid identity: SLOR updates grids in
  // place, and rejected line-search trials must not leave stale controls.
  const cached = Object.fromEntries(selected.map(side => [side, []]));
  const boundaryControl = (grid, side, i) => {
    let entry = cached[side][i];
    if (!entry) {
      const boundary = side === 'lower' ? 0 : nt, sign = side === 'lower' ? 1 : -1, points = [];
      for (const station of selectedCorners[side].has(i) ? [i - 1, i + 1] : [i]) {
        for (const k of [station - 1, station, station + 1]) points.push([k, boundary]);
        for (const distance of [1, 2]) points.push([station, boundary + sign * distance]);
      }
      entry = cached[side][i] = { points, coordinates: new Float64Array(2 * points.length) };
    }
    let unchanged = entry.control !== undefined;
    for (let k = 0; unchanged && k < entry.points.length; k++) {
      const [station, j] = entry.points[k], p = grid[station][j];
      unchanged = Object.is(p.x, entry.coordinates[2 * k]) && Object.is(p.y, entry.coordinates[2 * k + 1]);
    }
    if (unchanged) return entry.control;
    const control = controlledBoundary(grid, x, e, discretization, side, i, selectedCorners[side], update, sourceForm);
    for (let k = 0; k < entry.points.length; k++) {
      const [station, j] = entry.points[k], p = grid[station][j];
      entry.coordinates[2 * k] = p.x; entry.coordinates[2 * k + 1] = p.y;
    }
    entry.control = control;
    return control;
  };
  return { evaluate: (grid, i, j, derivativeRow = j) => {
    if (!Number.isInteger(derivativeRow) || derivativeRow < 0 || derivativeRow > nt)
      throw new Error('Boundary feedback derivative row is outside the grid.');
    let value = prescribed[i][j]; const derivative = { x: 0, y: 0 };
    const neighborDerivatives = hasCorners ? [{ x: 0, y: 0 }, { x: 0, y: 0 }] : null;
    for (const side of selected) {
      if (!active[side][i]) continue;
      const lower = side === 'lower', boundary = lower ? 0 : nt, weight = weights[j][side];
      const control = boundaryControl(grid, side, i);
      value += weight * (control[field] - prescribed[i][boundary]);
      const distance = lower ? derivativeRow : nt - derivativeRow;
      if (distance === 1 || distance === 2) {
        if (control.corner) control.neighbors.forEach((neighbor, k) => {
          for (const key of ['x', 'y']) neighborDerivatives[k][key] += .5 * weight * neighbor.control.interiorDerivatives[distance - 1][key];
        });
        else for (const key of ['x', 'y']) derivative[key] += weight * control.interiorDerivatives[distance - 1][key];
      }
    }
    return { [field]: value, derivative, ...(hasCorners ? { neighborDerivatives } : {}) };
  } };
}
