// SPDX-License-Identifier: GPL-2.0-or-later
// Drela (1986), Appendix E, E.3--E.5: inverse Laplace grid equations.
// Mass-fraction eta and prescribed uniform/nonuniform xi. The separately
// named giles-1985 mode follows the original ELLIP secants and line updates.
// Default boundary coordinates are fixed. Giles's fixed-indexed-y/adjacent-x
// rule (printed pp.186--189) also permits NORLIN's sloped outer rows.
// Optional fixed P adds the standard inverse-Poisson streamwise control:
// Laplacian(xi)=P, Laplacian(eta)=0. Sorenson, NASA TM-81198, Eq. (3).
// This is control-term infrastructure, not a recovered MSET X-spacing law.
import { createWasmSmoother } from './wasm-smoothing.js';
import { createNormalGraphBoundary } from './normal-graph-boundary.js';
import { createOrthogonalBoundaryFeedback } from './orthogonal-boundary-control.js';
import { solveBlockTridiagonal2 } from '../numerics/block-tridiagonal-2.js';
import { solveBlockTridiagonal } from '../numerics/block-tridiagonal.js';
import { grapePoissonDrift } from './grape-poisson-drift.js';
import { slorNumericalError, slorErrorTermination, tagSlorTermination } from './slor-termination.js';

const copy = nodes => nodes.map(row => row.map(p => ({ x: p.x, y: p.y })));

// Giles (1985), printed p.179, OUTLIN: delta-xi is proportional to the
// fourth root of boundary segment length. Its common DXAVG factor cancels
// when normalizing xi. This is the ancestor's rule, not modern MSET's X.
export function gilesStreamwiseCoordinates(boundary) {
  if (!Array.isArray(boundary) || boundary.length < 3 || boundary.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))
    throw new Error('Giles station coordinates require a finite boundary path.');
  const xi = [0];
  for (let i = 1; i < boundary.length; i++) {
    const length = Math.hypot(boundary[i].x - boundary[i - 1].x, boundary[i].y - boundary[i - 1].y);
    if (!(length > 0) || !Number.isFinite(length)) throw new Error('Giles station coordinates require nonzero finite boundary segments.');
    xi.push(xi.at(-1) + Math.sqrt(Math.sqrt(length)));
  }
  const total = xi.at(-1);
  if (!Number.isFinite(total)) throw new Error('Unresolved Giles station coordinate extent.');
  return xi.map(v => v / total);
}

export function createEllipticStreamtubeGrid({ nodes, massFlows, streamwiseSource, streamwiseStretch, harmonicMapControl, orthogonalBoundaryControl, streamwiseCoordinates,
  discretization = 'quadratic', streamwiseSourceDiscretization = 'centered', lineLinearization = 'frozen-metrics', lineGrouping = 'single', lineSearch = 'none', boundaryConditions = {}, boundaryCurves = {} }) {
  const nx = nodes?.length - 1, nt = massFlows?.length;
  if (!Array.isArray(nodes) || nx < 2 || !Array.isArray(massFlows) || nt < 2
    || !massFlows.every(m => Number.isFinite(m) && m > 0)
    || !nodes.every(row => Array.isArray(row) && row.length === nt + 1
      && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))))
    throw new Error('Invalid elliptic streamtube grid.');
  if (!['quadratic', 'giles-1985'].includes(discretization)) throw new Error('Unknown elliptic SLOR discretization.');
  if (!['centered', 'grape-1980'].includes(streamwiseSourceDiscretization)) throw new Error('Unknown elliptic source discretization.');
  if (!['frozen-metrics', 'full-metrics'].includes(lineLinearization)
    || lineLinearization === 'full-metrics' && (discretization !== 'giles-1985' || orthogonalBoundaryControl === undefined))
    throw new Error('Full line metrics currently require Giles secants and boundary-angle control.');
  if (!['single', 'boundary-pairs'].includes(lineGrouping) || lineGrouping === 'boundary-pairs' && lineLinearization !== 'full-metrics')
    throw new Error('Boundary row pairs require full line metrics.');
  if (!['none', 'armijo'].includes(lineSearch) || lineSearch === 'armijo' && lineGrouping !== 'boundary-pairs')
    throw new Error('Newton line backtracking requires paired boundary rows.');
  if (!boundaryConditions || typeof boundaryConditions !== 'object' || Array.isArray(boundaryConditions)
    || Object.keys(boundaryConditions).some(key => !['lower', 'upper'].includes(key)))
    throw new Error('Invalid elliptic boundary conditions.');
  const boundaries = Object.freeze(Object.fromEntries(['lower', 'upper'].map(side => {
    const mode = boundaryConditions[side] === undefined ? 'fixed' : boundaryConditions[side];
    if (!['fixed', 'giles-vertical', 'giles-indexed-y', 'normal-curve'].includes(mode)) throw new Error('Unknown elliptic boundary condition.');
    const j = side === 'lower' ? 0 : nt;
    // The horizontal-only mode preserves its physical line. The separate
    // indexed-y mode follows ELLIP for nonconstant NORLIN y values too; its
    // changing x values do not preserve a fixed physical boundary curve.
    if (mode === 'giles-vertical' && nodes.some(row => row[j].y !== nodes[0][j].y))
      throw new Error('Giles vertical-copy boundary conditions require a horizontal boundary.');
    return [side, mode];
  })));
  if (!boundaryCurves || typeof boundaryCurves !== 'object' || Array.isArray(boundaryCurves)
    || Object.keys(boundaryCurves).some(side => !['lower', 'upper'].includes(side) || boundaries[side] !== 'normal-curve'))
    throw new Error('Invalid elliptic boundary curves.');
  const curves = Object.fromEntries(['lower', 'upper'].filter(side => boundaries[side] === 'normal-curve')
    .map(side => [side, createNormalGraphBoundary(boundaryCurves[side] ?? {})]));
  const sideAt = j => j === 0 ? 'lower' : 'upper';
  const freeBoundary = j => (j === 0 || j === nt) && boundaries[sideAt(j)] !== 'fixed';
  const freeBoundaryX = j => (j === 0 || j === nt)
    && ['giles-vertical', 'giles-indexed-y'].includes(boundaries[sideAt(j)]);
  if (streamwiseCoordinates !== undefined && (!Array.isArray(streamwiseCoordinates) || streamwiseCoordinates.length !== nx + 1
    || streamwiseCoordinates[0] !== 0 || streamwiseCoordinates[nx] !== 1
    || streamwiseCoordinates.some((v, i) => !Number.isFinite(v) || (i && !(v > streamwiseCoordinates[i - 1])))))
    throw new Error('Streamwise coordinates must increase strictly from zero to one and match the grid.');
  const xi = Object.freeze(streamwiseCoordinates ? [...streamwiseCoordinates] : Array.from({ length: nx + 1 }, (_, i) => i / nx));
  if (orthogonalBoundaryControl !== undefined && discretization !== 'giles-1985')
    throw new Error('Boundary-angle feedback currently requires the Giles discretization and metrics updated between lines.');
  if (orthogonalBoundaryControl !== undefined && (!orthogonalBoundaryControl || typeof orthogonalBoundaryControl !== 'object'
    || Array.isArray(orthogonalBoundaryControl) || Object.keys(orthogonalBoundaryControl).some(k => !['background', 'sides', 'decay', 'corners', 'activeStations', 'previous', 'relaxation', 'changeLimit', 'sourceForm'].includes(k))
    || streamwiseSource !== undefined || streamwiseStretch !== undefined || harmonicMapControl !== undefined
    || Object.values(boundaries).some(mode => mode !== 'fixed')))
    throw new Error('Boundary-angle feedback requires fixed boundaries and cannot be combined with other Poisson controls.');
  if (streamwiseSource !== undefined && (!Array.isArray(streamwiseSource) || streamwiseSource.length !== nx + 1
    || streamwiseSource.some(row => !Array.isArray(row) || row.length !== nt + 1 || !row.every(Number.isFinite))))
    throw new Error('Streamwise Poisson source must be a finite array matching the grid.');
  if (streamwiseStretch !== undefined && (streamwiseSource !== undefined || harmonicMapControl !== undefined
    || !Array.isArray(streamwiseStretch) || streamwiseStretch.length !== nx + 1
    || streamwiseStretch.some(row => !Array.isArray(row) || row.length !== nt + 1 || !row.every(Number.isFinite))))
    throw new Error('Streamwise stretch must match the grid and cannot be combined with other Poisson controls.');
  if (harmonicMapControl !== undefined && (streamwiseSource !== undefined || !Array.isArray(harmonicMapControl) || harmonicMapControl.length !== nx + 1
    || harmonicMapControl.some(row => !Array.isArray(row) || row.length !== nt + 1 || row.some(c => !c
      || !['xiXi', 'xiEta', 'etaEta'].every(key => Number.isFinite(c[key]))))))
    throw new Error('Harmonic-map control must match the grid and cannot be combined with a fixed Poisson source.');
  // P is specified at logical nodes, has physical units length^-2, and is
  // immutable throughout relaxation. Q is zero so eta remains a physical
  // harmonic mass coordinate; a general geometric Q would change that flow.
  const poissonP = streamwiseSource?.map(row => Object.freeze([...row]));
  if (poissonP) Object.freeze(poissonP);
  // Δxi = |grad xi|² F, Δeta = 0 gives the inverse drift alpha*F.
  // F is dimensionless for normalized xi, and is fixed in computational
  // space. Its physical Poisson source changes with the grid metrics.
  const stretch = streamwiseStretch?.map(row => Object.freeze([...row]));
  if (stretch) Object.freeze(stretch);
  const mapControl = harmonicMapControl?.map(row => Object.freeze(row.map(c => Object.freeze({ xiXi: c.xiXi, xiEta: c.xiEta, etaEta: c.etaEta }))));
  if (mapControl) Object.freeze(mapControl);
  const mapped = mapControl?.some(row => row.some(c => c.xiXi !== 0 || c.xiEta !== 0 || c.etaEta !== 0)) ?? false;
  const stretched = stretch?.some(row => row.some(f => f !== 0)) ?? false;
  const controlled = orthogonalBoundaryControl !== undefined || mapped || stretched || (poissonP?.some(row => row.some(p => p !== 0)) ?? false);
  const initial = copy(nodes), totalMass = massFlows.reduce((s, m) => s + m, 0), eta = [0];
  for (const m of massFlows) eta.push(eta.at(-1) + m / totalMass);
  eta[nt] = 1;
  if (!Number.isFinite(totalMass) || eta.some((v, i) => i && !(v > eta[i - 1])))
    throw new Error('Unresolved elliptic streamtube mass coordinates.');
  const angleFeedback = orthogonalBoundaryControl === undefined ? null : createOrthogonalBoundaryFeedback({
    nodes: initial, xi, eta, discretization, ...orthogonalBoundaryControl });
  const physicalAngleSource = angleFeedback && orthogonalBoundaryControl.sourceForm === 'poisson';
  const origin = initial[0][0]; let lengthScale = 0;
  for (const row of initial) for (const p of row)
    lengthScale = Math.max(lengthScale, Math.hypot(p.x - origin.x, p.y - origin.y));
  if (!(lengthScale > 0) || !Number.isFinite(lengthScale)) throw new Error('Degenerate elliptic grid.');
  const differences = (a, b) => ({
    first: discretization === 'giles-1985' ? [-1 / (a + b), 0, 1 / (a + b)]
      : [-b / (a * (a + b)), (b - a) / (a * b), a / (b * (a + b))],
    second: [2 / (a * (a + b)), -2 / (a * b), 2 / (b * (a + b))] });
  const derivatives = eta.map((v, j) => !j || j === nt ? null : differences(v - eta[j - 1], eta[j + 1] - v));
  const xiDerivatives = xi.map((v, i) => !i || i === nx ? null : streamwiseCoordinates
    ? differences(v - xi[i - 1], xi[i + 1] - v)
    : { first: [-nx / 2, 0, nx / 2], second: [nx * nx, -2 * nx * nx, nx * nx] });
  const validate = grid => {
    if (!Array.isArray(grid) || grid.length !== nx + 1 || !grid.every(row => Array.isArray(row) && row.length === nt + 1
      && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))) throw new Error('Invalid elliptic grid state.');
    for (const j of [0, nt]) if (curves[sideAt(j)]) for (let i = 0; i <= nx; i++) {
      const point = grid[i][j], expected = curves[sideAt(j)].evaluate(point.x).point;
      if (Math.abs(point.y - expected.y) > 64 * Number.EPSILON * Math.max(lengthScale, Math.abs(point.y)))
        throw new Error('Elliptic farfield node left its prescribed curve.');
    }
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) if (!i || i === nx || !j || j === nt) {
      if (i && i !== nx && curves[sideAt(j)]) {
        continue;
      }
      if (grid[i][j].y !== initial[i][j].y
        || ((!i || i === nx || !freeBoundaryX(j)) && grid[i][j].x !== initial[i][j].x))
        throw new Error('Prescribed elliptic grid boundary coordinates must remain fixed.');
    }
    for (const j of [0, nt]) if (curves[sideAt(j)])
      for (let i = 1; i <= nx; i++) if (!(grid[i][j].x > grid[i - 1][j].x))
        throw new Error('Elliptic farfield stations are not strictly ordered.');
  };
  validate(initial);
  // Difference form avoids cancellation under translation. The existing
  // stencil uses quadratic-exact nonuniform differences. Giles's printed
  // ELLIP instead uses centered secants and updates metrics between lines.
  const etaFirst = (grid, i, j, key) => {
    const [a, , c] = derivatives[j].first;
    return a * (grid[i][j - 1][key] - grid[i][j][key]) + c * (grid[i][j + 1][key] - grid[i][j][key]);
  };
  const xiFirst = (grid, i, j, key) => {
    const [a, , c] = xiDerivatives[i].first;
    return a * (grid[i - 1][j][key] - grid[i][j][key]) + c * (grid[i + 1][j][key] - grid[i][j][key]);
  };
  const sourceStencil = (i, drift) => streamwiseSourceDiscretization === 'grape-1980'
    ? grapePoissonDrift({ drift, leftDistance: xi[i] - xi[i - 1], rightDistance: xi[i + 1] - xi[i] })
    : { coefficients: xiDerivatives[i].first.map(w => drift * w), first: xiDerivatives[i].first };
  const sourceFirst = (grid, i, j, key, drift) => {
    if (streamwiseSourceDiscretization === 'centered') return xiFirst(grid, i, j, key);
    const [a, , c] = sourceStencil(i, drift).first;
    return a * (grid[i - 1][j][key] - grid[i][j][key]) + c * (grid[i + 1][j][key] - grid[i][j][key]);
  };
  const metricAt = (grid, i, j, withJacobian = true) => {
      const dx = xiFirst(grid, i, j, 'x'), dy = xiFirst(grid, i, j, 'y');
      const ex = etaFirst(grid, i, j, 'x'), ey = etaFirst(grid, i, j, 'y');
      const alpha = ex * ex + ey * ey, beta = dx * ex + dy * ey, gamma = dx * dx + dy * dy;
      const jacobian = dx * ey - dy * ex;
      if (![alpha, beta, gamma, jacobian].every(Number.isFinite) || !(alpha > 0 && gamma > 0) || jacobian === 0)
        throw new Error(`Singular elliptic metric at station ${i}, streamline ${j}.`);
      const control = mapControl?.[i][j];
      const feedback = angleFeedback?.evaluate(grid, i, j);
      const feedbackScale = physicalAngleSource ? jacobian * jacobian : alpha;
      const streamwiseDrift = feedback ? feedbackScale * (physicalAngleSource ? feedback.poisson : feedback.stretch) : mapped ? alpha * control.xiXi - 2 * beta * control.xiEta + gamma * control.etaEta
        : stretched ? alpha * stretch[i][j] : poissonP?.[i][j] ? jacobian * jacobian * poissonP[i][j] : 0;
      if (!Number.isFinite(streamwiseDrift)) throw new Error('Nonfinite elliptic Poisson control coefficient.');
      const sx = withJacobian && feedback ? sourceFirst(grid, i, j, 'x', streamwiseDrift) : dx;
      const sy = withJacobian && feedback ? sourceFirst(grid, i, j, 'y', streamwiseDrift) : dy;
      const feedbackMatrix = d => [feedbackScale * sx * d.x, feedbackScale * sx * d.y, feedbackScale * sy * d.x, feedbackScale * sy * d.y];
      let metricLineJacobian;
      if (withJacobian && lineLinearization === 'full-metrics') {
        // With Giles secants, r_eta has no dependence on this eta line.
        // Only g=r_xi changes its metrics: d(alpha)=0, d(beta)=r_eta.dg,
        // d(gamma)=2*g.dg, d(J)=ey*dg.x-ex*dg.y. Include d(J²*P) as
        // well as the separate boundary-control dP already assembled below.
        // g depends on the two xi neighbors, preserving block tridiagonality.
        const [a, , c] = derivatives[j].second, [ax, , cx] = xiDerivatives[i].first;
        metricLineJacobian = ['x', 'y'].flatMap((key, k) => {
          const center = grid[i][j][key], de = etaFirst(grid, i, j, key);
          const etaSecond = a * (grid[i][j - 1][key] - center) + c * (grid[i][j + 1][key] - center);
          const mixed = ax * (etaFirst(grid, i - 1, j, key) - de) + cx * (etaFirst(grid, i + 1, j, key) - de);
          const jacobianFactor = physicalAngleSource ? 2 * streamwiseDrift / jacobian * (k ? sy : sx) : 0;
          return [-2 * ex * mixed + 2 * dx * etaSecond + jacobianFactor * ey,
            -2 * ey * mixed + 2 * dy * etaSecond - jacobianFactor * ex];
        });
      }
      return { alpha, beta, gamma, jacobian, streamwiseDrift,
        ...(metricLineJacobian ? { metricLineJacobian } : {}),
        ...(withJacobian && feedback ? { streamwiseJacobian: feedbackMatrix(feedback.derivative),
          ...(feedback.neighborDerivatives ? { neighboringStreamwiseJacobians: feedback.neighborDerivatives.map(feedbackMatrix) } : {}) } : {}) };
  };
  const metrics = (grid, withJacobian = true) => {
    validate(grid);
    return Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
      if (!i || i === nx || !j || j === nt) return null;
      return metricAt(grid, i, j, withJacobian);
    }));
  };
  const operator = (grid, i, j, key, { alpha, beta, gamma, streamwiseDrift = 0 }) => {
    const center = grid[i][j][key], [a, , c] = derivatives[j].second;
    const [ax, , cx] = xiDerivatives[i].second, [a1, , c1] = xiDerivatives[i].first;
    const dxx = ax * (grid[i - 1][j][key] - center) + cx * (grid[i + 1][j][key] - center);
    const dyy = a * (grid[i][j - 1][key] - center) + c * (grid[i][j + 1][key] - center);
    const mid = etaFirst(grid, i, j, key);
    const dxy = a1 * (etaFirst(grid, i - 1, j, key) - mid) + c1 * (etaFirst(grid, i + 1, j, key) - mid);
    return alpha * dxx - 2 * beta * dxy + gamma * dyy + streamwiseDrift * sourceFirst(grid, i, j, key, streamwiseDrift);
  };
  const residuals = (grid, coefficients = metrics(grid, false)) => {
    let residual = 0; const rows = [];
    for (let i = 1; i < nx; i++) for (let j = 1; j < nt; j++) {
      const { alpha, gamma } = coefficients[i][j];
      const row = { i, j };
      for (const key of ['x', 'y']) {
        row[key] = operator(grid, i, j, key, coefficients[i][j]) / ((alpha + gamma) * lengthScale);
        if (!Number.isFinite(row[key])) throw new Error('Nonfinite elliptic grid residual.');
        residual = Math.max(residual, Math.abs(row[key]));
      }
      rows.push(row);
    }
    const interiorResidual = residual, boundaryRows = []; let boundaryResidual = 0;
    for (const j of [0, nt]) if (freeBoundary(j)) {
      const adjacent = j === 0 ? 1 : nt - 1, width = Math.abs(eta[j] - eta[adjacent]);
      for (let i = 1; i < nx; i++) {
        // Giles modes measure x_eta, including when indexed y varies: that
        // is not a normal condition on a sloped physical boundary. Curve
        // mode instead measures tangential separation. Both use the same
        // mass-coordinate width and global physical length normalization.
        const curve = curves[sideAt(j)], point = grid[i][j], neighbor = grid[i][adjacent];
        const tangent = curve ? curve.evaluate(point.x).derivative : { x: 1, y: 0 };
        const x = ((point.x - neighbor.x) * tangent.x + (point.y - neighbor.y) * tangent.y)
          / Math.hypot(tangent.x, tangent.y) / width / lengthScale;
        if (!Number.isFinite(x)) throw new Error('Nonfinite elliptic boundary residual.');
        boundaryResidual = Math.max(boundaryResidual, Math.abs(x));
        boundaryRows.push({ i, j, x });
      }
    }
    return { residual: Math.max(interiorResidual, boundaryResidual), interiorResidual, boundaryResidual, rows, boundaryRows };
  };
  const quality = (grid, cellSines, firstChanged = 0, lastChanged = nt - 1) => {
    let minCornerSine = Infinity; const invalidCells = [];
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
      let cellMinimum = cellSines?.[i * nt + j];
      if (!cellSines || j >= firstChanged && j <= lastChanged) {
        const points = [grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]];
        cellMinimum = Infinity;
        for (let k = 0; k < 4; k++) {
          const p = points[k], q = points[(k + 1) % 4], r = points[(k + 2) % 4];
          const ax = q.x - p.x, ay = q.y - p.y, bx = r.x - q.x, by = r.y - q.y;
          const sine = (ax * by - ay * bx) / (Math.hypot(ax, ay) * Math.hypot(bx, by));
          cellMinimum = Math.min(cellMinimum, Number.isFinite(sine) ? sine : -Infinity);
        }
        if (cellSines) cellSines[i * nt + j] = cellMinimum;
      }
      minCornerSine = Math.min(minCornerSine, cellMinimum);
      if (!(cellMinimum > 1e-12)) invalidCells.push({ i, j });
    }
    return { valid: !invalidCells.length, minCornerSine, invalidCells };
  };
  const rowGroups = Array.from({ length: nt - 1 }, (_, j) => [j + 1]);
  if (lineGrouping === 'boundary-pairs') {
    // Merge the first two interior rows at each controlled boundary. If
    // those pairs overlap on a small grid, solve their union once.
    for (const side of orthogonalBoundaryControl.sides ?? ['lower', 'upper']) {
      const pair = side === 'lower' ? [1, 2] : [nt - 2, nt - 1];
      const hits = rowGroups.map((rows, k) => rows.some(j => pair.includes(j)) ? k : -1).filter(k => k >= 0);
      rowGroups.splice(hits[0], hits.length, rowGroups.slice(hits[0], hits.at(-1) + 1).flat());
    }
  }
  const lineWorkspaces = new Map();
  const linearizeRows = (grid, rows, reuseWorkspace = false) => {
    if (!angleFeedback || discretization !== 'giles-1985' || !Array.isArray(rows) || !rows.length
      || rows.some((j, k) => !Number.isInteger(j) || j < 1 || j >= nt || k && j !== rows[k - 1] + 1))
      throw new Error('Coupled rows require boundary-angle equations and consecutive interior rows.');
    validate(grid);
    const size = 2 * rows.length;
    let workspace = reuseWorkspace && lineWorkspaces.get(size);
    if (!workspace) {
      const buffer = new Float64Array((nx - 1) * (3 * size ** 2 + 2 * size));
      let offset = 0;
      const blocks = width => Array.from({ length: nx - 1 }, () => {
        const view = buffer.subarray(offset, offset + width); offset += width; return view;
      });
      workspace = { lower: blocks(size ** 2), diagonal: blocks(size ** 2), upper: blocks(size ** 2),
        rhs: blocks(size), scales: blocks(size) };
      if (reuseWorkspace) lineWorkspaces.set(size, workspace);
    }
    // Every entry is overwritten below. Public linearizations retain ownership;
    // only the synchronous internal sweep borrows these buffers.
    const { lower, diagonal, upper, rhs, scales } = workspace;
    for (let i = 1; i < nx; i++) for (const [jr, j] of rows.entries()) {
      const c = metricAt(grid, i, j, false), g = ['x', 'y'].map(key => xiFirst(grid, i, j, key)), h = ['x', 'y'].map(key => etaFirst(grid, i, j, key));
      const wx = xiDerivatives[i], we = derivatives[j], source = sourceStencil(i, c.streamwiseDrift);
      const second = (axis, key) => {
        const p = grid[i][j][key], [a, , b] = (axis === 'xi' ? wx : we).second;
        return axis === 'xi' ? a * (grid[i - 1][j][key] - p) + b * (grid[i + 1][j][key] - p)
          : a * (grid[i][j - 1][key] - p) + b * (grid[i][j + 1][key] - p);
      };
      const xixi = ['x', 'y'].map(key => second('xi', key)), etaeta = ['x', 'y'].map(key => second('eta', key));
      const xieta = ['x', 'y'].map((key, k) => wx.first[0] * (etaFirst(grid, i - 1, j, key) - h[k])
        + wx.first[2] * (etaFirst(grid, i + 1, j, key) - h[k]));
      const driftFirst = ['x', 'y'].map(key => sourceFirst(grid, i, j, key, c.streamwiseDrift));
      const fieldScale = physicalAngleSource ? c.jacobian ** 2 : c.alpha;
      for (let r = 0; r < 2; r++) {
        rhs[i - 1][2 * jr + r] = -operator(grid, i, j, r ? 'y' : 'x', c);
        scales[i - 1][2 * jr + r] = (c.alpha + c.gamma) * lengthScale;
      }
      for (const [kr, k] of rows.entries()) {
        const feedback = angleFeedback.evaluate(grid, i, j, k), de = k - j, we1 = we.first[de + 1] ?? 0, we2 = we.second[de + 1] ?? 0;
        for (let di = -1; di <= 1; di++) {
          const block = [lower, diagonal, upper][di + 1][i - 1];
          const scalar = c.alpha * wx.second[di + 1] * Number(k === j) - 2 * c.beta * wx.first[di + 1] * we1
            + c.gamma * we2 * Number(di === 0) + source.coefficients[di + 1] * Number(k === j);
          const dControl = di === 0 ? feedback.derivative : feedback.neighborDerivatives?.[di < 0 ? 0 : 1] ?? { x: 0, y: 0 };
          for (let v = 0; v < 2; v++) {
            const dg = k === j ? wx.first[di + 1] : 0, dh = di === 0 ? we1 : 0;
            const da = 2 * h[v] * dh, db = h[v] * dg + g[v] * dh, dc = 2 * g[v] * dg;
            const dJ = v === 0 ? h[1] * dg - g[1] * dh : -h[0] * dg + g[0] * dh;
            const dDrift = (physicalAngleSource ? 2 * c.streamwiseDrift / c.jacobian * dJ : c.streamwiseDrift / c.alpha * da)
              + fieldScale * dControl[v ? 'y' : 'x'];
            for (let r = 0; r < 2; r++) block[(2 * jr + r) * size + 2 * kr + v] = (r === v ? scalar : 0)
              + da * xixi[r] - 2 * db * xieta[r] + dc * etaeta[r] + dDrift * driftFirst[r];
          }
        }
      }
    }
    return { lower, diagonal, upper, rhs, scales };
  };
  const sweep = (grid, omega = 1.3, { lineTolerance = 0 } = {}) => {
    if (!Number.isFinite(omega) || !(omega > 0 && omega < 2)) throw new Error('SLOR omega must lie between zero and two.');
    if (!Number.isFinite(lineTolerance) || lineTolerance < 0) throw new Error('Invalid Newton line tolerance.');
    if (lineGrouping === 'boundary-pairs') {
      let next = copy(grid), maxUpdate = 0; const rowSteps = [];
      // The accepted grid's metrics are valid everywhere. Interior row updates
      // only change a three-row stencil. Boundary feedback reads the first two
      // interior rows, so changes there deliberately invalidate every metric.
      const acceptedMetrics = metrics(next, false);
      const cellSines = new Float64Array(nx * nt);
      quality(next, cellSines);
      const controlledRows = new Set((orthogonalBoundaryControl.sides ?? ['lower', 'upper'])
        .flatMap(side => side === 'lower' ? [1, 2] : [nt - 2, nt - 1]));
      for (const rows of rowGroups) {
        const matrix = linearizeRows(next, rows, true);
        const lineResidual = Math.max(...matrix.rhs.flatMap((row, i) => Array.from(row, (r, k) => Math.abs(r / matrix.scales[i][k]))));
        if (lineSearch === 'armijo' && lineResidual <= lineTolerance) {
          rowSteps.push({ rows: rows.slice(), skipped: true, residual: lineResidual, tolerance: lineTolerance }); continue;
        }
        const correction = solveBlockTridiagonal(matrix);
        const baseMerit = .5 * matrix.rhs.reduce((sum, row, i) => sum + row.reduce((s, r, k) => s + (r / matrix.scales[i][k]) ** 2, 0), 0);
        let accepted = false; const trials = [];
        const firstMetricRow = rows.some(j => controlledRows.has(j)) ? 1 : Math.max(1, rows[0] - 1);
        const lastMetricRow = rows.some(j => controlledRows.has(j)) ? nt - 1 : Math.min(nt - 1, rows.at(-1) + 1);
        // Share untouched points read-only; trial updates exclusively own the
        // selected rows. Rejected trials are overwritten from the accepted state.
        const trial = next.map(column => {
          const result = column.slice();
          for (const j of rows) result[j] = { ...column[j] };
          return result;
        });
        for (let halving = 0; halving <= (lineSearch === 'armijo' ? 20 : 0); halving++) {
          const fraction = omega * 2 ** -halving; let movement = 0;
          for (let i = 1; i < nx; i++) for (const [r, j] of rows.entries()) for (const [k, key] of ['x', 'y'].entries()) {
            const delta = fraction * correction[i - 1][2 * r + k]; trial[i][j][key] = next[i][j][key] + delta;
            movement = Math.max(movement, Math.abs(delta) / lengthScale);
          }
          if (lineSearch === 'none') { next = trial; maxUpdate = Math.max(maxUpdate, movement); accepted = true; break; }
          const q = quality(trial, cellSines, rows[0] - 1, rows.at(-1)); let merit = Infinity, rejection = q.valid ? null : 'folded grid';
          if (q.valid) try {
            validate(trial);
            // Retain checks outside the solved rows wherever either the local
            // metric stencil or the global boundary control can have changed.
            for (let i = 1; i < nx; i++) for (let j = firstMetricRow; j <= lastMetricRow; j++)
              acceptedMetrics[i][j] = metricAt(trial, i, j, false);
            const trialMetrics = acceptedMetrics; merit = 0;
            for (let i = 1; i < nx; i++) for (const [r, j] of rows.entries()) for (const [k, key] of ['x', 'y'].entries())
              merit += .5 * (operator(trial, i, j, key, trialMetrics[i][j]) / matrix.scales[i - 1][2 * r + k]) ** 2;
          } catch (error) { rejection = error.message; }
          // W is fixed at this line's initial metric. J*d=-R implies
          // d(0.5*||W R||²)/d(step)=-||W R||² exactly at step zero.
          const sufficient = !rejection && merit <= (1 - 2e-4 * fraction) * baseMerit;
          trials.push({ fraction, merit, rejected: rejection, accepted: sufficient });
          if (sufficient) {
            next = trial; maxUpdate = Math.max(maxUpdate, movement); accepted = true;
            rowSteps.push({ rows: rows.slice(), baseMerit, merit, fraction, minCornerSine: q.minCornerSine, trials }); break;
          }
        }
        if (!accepted) throw slorNumericalError(`No admissible decreasing Newton step for eta rows ${rows.join(', ')}.`, 'no-admissible-decreasing-line-step');
      }
      validate(next); return { nodes: next, maxUpdate, ...(lineSearch === 'armijo' ? { rowSteps } : {}) };
    }
    const coefficients = metrics(grid), next = copy(grid); let maxUpdate = 0;
    const copyBoundary = j => {
      if (!freeBoundary(j)) return;
      const adjacent = j === 0 ? 1 : nt - 1;
      for (let i = 1; i < nx; i++) {
        const curve = curves[sideAt(j)], point = curve ? curve.projectNormal(next[i][adjacent])
          : { x: next[i][adjacent].x, y: next[i][j].y };
        if (curve) {
          const tangent = curve.evaluate(point.x).derivative, neighbor = next[i][adjacent];
          const inward = (neighbor.y - point.y) * tangent.x - (neighbor.x - point.x) * tangent.y;
          if (!(j === 0 ? inward > 0 : inward < 0))
            throw new Error('Adjacent elliptic node lies outside its normal farfield boundary.');
        }
        maxUpdate = Math.max(maxUpdate, Math.abs(point.x - next[i][j].x) / lengthScale,
          Math.abs(point.y - next[i][j].y) / lengthScale);
        next[i][j] = point;
      }
    };
    // ELLIP encounters the lower physical farfield before its adjacent
    // interior line, and the upper farfield after that line. Direct copies
    // are not over-relaxed; endpoints remain prescribed in both cases.
    copyBoundary(0);
    // Frozen metric; Gauss--Seidel between eta lines. The tridiagonal solve
    // includes the same-line part of the nonuniform mixed derivative.
    for (let j = 1; j < nt; j++) {
      const lower = new Float64Array(nx - 1), diagonal = new Float64Array(nx - 1), upper = new Float64Array(nx - 1);
      const rhs = { x: new Float64Array(nx - 1), y: new Float64Array(nx - 1) };
      const [, b1] = derivatives[j].first, [, b2] = derivatives[j].second;
      // Original ELLIP (Giles pp.187–188) evaluates metrics on the current
      // grid while assembling each line, after preceding lines were updated.
      if (discretization === 'giles-1985') for (let i = 1; i < nx; i++) coefficients[i][j] = metricAt(next, i, j);
      for (let i = 1; i < nx; i++) {
        const k = i - 1, { alpha, beta, gamma, streamwiseDrift } = coefficients[i][j];
        const [ax1, bx1, cx1] = xiDerivatives[i].first, [ax2, bx2, cx2] = xiDerivatives[i].second;
        const [sourceLower, sourceDiagonal, sourceUpper] = sourceStencil(i, streamwiseDrift).coefficients;
        lower[k] = alpha * ax2 - 2 * beta * ax1 * b1 + sourceLower;
        upper[k] = alpha * cx2 - 2 * beta * cx1 * b1 + sourceUpper;
        diagonal[k] = alpha * bx2 - 2 * beta * bx1 * b1 + gamma * b2 + sourceDiagonal;
        // Solve for a correction, keeping boundary corrections exactly zero.
        // Evaluating the residual with updated neighboring lines gives the
        // same frozen-coefficient SLOR system without absolute-coordinate RHS.
        for (const key of ['x', 'y']) {
          rhs[key][k] = -operator(next, i, j, key, coefficients[i][j]);
        }
      }
      if (angleFeedback) {
        // Eliminate the boundary F or P analytically. Its derivative is
        // always included. Optional full-metrics adds the same-line metric
        // derivatives; the original frozen-metrics iteration is unchanged.
        // x/y coupling enters diagonal blocks and, at averaged corners,
        // neighboring blocks. Thomas retains linear work and line order.
        const neighborBlock = (a, k, side) => (coefficients[k + 1][j].neighboringStreamwiseJacobians?.[side] ?? [0, 0, 0, 0])
          .map((v, q) => v + (q === 0 || q === 3 ? a : 0)
            + (coefficients[k + 1][j].metricLineJacobian?.[q] ?? 0) * xiDerivatives[k + 1].first[side ? 2 : 0]);
        const correction = solveBlockTridiagonal2({ lower: Array.from(lower, (a, k) => neighborBlock(a, k, 0)),
          upper: Array.from(upper, (a, k) => neighborBlock(a, k, 1)),
          diagonal: Array.from(diagonal, (d, k) => coefficients[k + 1][j].streamwiseJacobian.map((v, q) => v + (q === 0 || q === 3 ? d : 0))),
          rhs: Array.from(rhs.x, (v, k) => [v, rhs.y[k]]) });
        for (let k = 0; k < nx - 1; k++) for (const [c, key] of ['x', 'y'].entries()) {
          const delta = omega * correction[k][c]; next[k + 1][j][key] += delta;
          maxUpdate = Math.max(maxUpdate, Math.abs(delta) / lengthScale);
        }
        continue;
      }
      for (let k = 0; k < nx - 1; k++) {
        const scale = Math.abs(lower[k]) + Math.abs(diagonal[k]) + Math.abs(upper[k]);
        if (!(Math.abs(diagonal[k]) > 32 * Number.EPSILON * scale)) throw new Error('Singular elliptic SLOR line.');
        if (k + 1 < nx - 1) {
          const factor = lower[k + 1] / diagonal[k]; diagonal[k + 1] -= factor * upper[k];
          for (const key of ['x', 'y']) rhs[key][k + 1] -= factor * rhs[key][k];
        }
      }
      for (let k = nx - 2; k >= 0; k--) for (const key of ['x', 'y']) {
        rhs[key][k] = (rhs[key][k] - (k === nx - 2 ? 0 : upper[k] * rhs[key][k + 1])) / diagonal[k];
        const delta = omega * rhs[key][k]; next[k + 1][j][key] += delta;
        maxUpdate = Math.max(maxUpdate, Math.abs(delta) / lengthScale);
      }
    }
    copyBoundary(nt);
    validate(next);
    return { nodes: next, maxUpdate };
  };
  // Fixed-boundary harmonic initialization can reach coordinate precision
  // before its derivative residual reaches an absolute tolerance. Estimate
  // that residual's sensitivity to rounding the stored coordinates, using
  // the frozen nine-point operator. This is a backward-error diagnostic,
  // not a bound on the global grid error or an Euler convergence test.
  const coordinateRoundoff = (grid, tolerance) => {
    if (angleFeedback || mapped || controlled || stretched || Object.values(boundaries).some(mode => mode !== 'fixed')) return null;
    const coefficients = metrics(grid, false), state = residuals(grid, coefficients);
    let maximumRatio = 0, maximumBound = 0;
    for (const row of state.rows) {
      const { i, j } = row, { alpha, beta, gamma } = coefficients[i][j];
      const scale = (alpha + gamma) * lengthScale;
      for (const key of ['x', 'y']) {
        let bound = 0;
        for (let di = -1; di <= 1; di++) for (let dj = -1; dj <= 1; dj++) {
          const a = alpha * xiDerivatives[i].second[di + 1] * Number(dj === 0)
            - 2 * beta * xiDerivatives[i].first[di + 1] * derivatives[j].first[dj + 1]
            + gamma * derivatives[j].second[dj + 1] * Number(di === 0);
          bound += Math.abs(a * grid[i + di][j + dj][key]) * Number.EPSILON / scale;
        }
        if (!Number.isFinite(bound)) return null;
        maximumBound = Math.max(maximumBound, bound);
        maximumRatio = Math.max(maximumRatio, Math.abs(row[key]) / Math.max(tolerance, bound));
      }
    }
    return { limited: maximumRatio <= 1, maximumRatio, maximumBound, residual: state.residual,
      method: 'coordinate-rounding sensitivity of frozen harmonic operator' };
  };
  const system = { initial, nx, nt, xi, eta, lengthScale, metrics, residuals, quality, sweep, linearizeRows, coordinateRoundoff,
    coordinateEquations: { xi: physicalAngleSource ? 'Poisson with implicit boundary-angle physical source' : angleFeedback ? 'Poisson with implicit boundary-angle streamwise stretch' : mapped ? 'Poisson from prescribed harmonic-map coefficients' : stretched ? 'Poisson with prescribed metric-scaled streamwise stretch' : controlled ? 'Poisson with prescribed fixed P' : 'Laplace', eta: 'Laplace',
      sourceUnits: stretched || angleFeedback && !physicalAngleSource ? 'dimensionless F multiplying squared physical gradient of normalized xi'
        : mapped ? 'metric contraction of computational map coefficients; physical Laplacian has inverse length squared units' : 'inverse physical length squared', exactMsetSpacingLaw: false, discretization,
      stationSpacing: streamwiseCoordinates ? 'prescribed nonuniform xi' : 'uniform xi', streamwiseSourceDiscretization,
      boundaryConditions: boundaries,
      ...(angleFeedback ? { lineSolve: lineGrouping === 'boundary-pairs' ? 'block Thomas with coupled boundary-adjacent rows' : '2x2 block Thomas with analytic boundary-control feedback',
        lineLinearization, lineGrouping, lineSearch, rowGroups: rowGroups.map(rows => rows.slice()), boundaryAngle: Math.PI / 2,
        controlledSides: (orthogonalBoundaryControl.sides ?? ['lower', 'upper']).slice(), controlExtension: orthogonalBoundaryControl.decay
          ? { method: 'matched exponential tails in normalized mass', decay: { ...orthogonalBoundaryControl.decay } } : 'linear in normalized mass',
        ...(orthogonalBoundaryControl.activeStations ? { activeBoundaryStations: Object.fromEntries(
          (orthogonalBoundaryControl.sides ?? ['lower', 'upper']).map(side => [side,
            Array.from({ length: nx + 1 }, (_, i) => Boolean(i && i < nx && (orthogonalBoundaryControl.activeStations[side]?.[i] ?? true))) ])) } : {}),
        ...(orthogonalBoundaryControl.corners ? { cornerControls: { method: `arithmetic mean of two neighboring ${physicalAngleSource ? 'P' : 'F'} controls`,
          stations: Object.fromEntries(Object.entries(orthogonalBoundaryControl.corners).map(([side, row]) => [side, row.slice()])) } } : {}),
        ...(orthogonalBoundaryControl.previous ? { controlUpdate: { method: 'implicit damped and limited update from fixed previous controls',
          relaxation: orthogonalBoundaryControl.relaxation ?? .3, changeLimit: orthogonalBoundaryControl.changeLimit ?? 1 } } : {}) } : {}),
      ...(Object.keys(curves).length ? { boundaryCurves: Object.fromEntries(Object.entries(curves)
        .map(([side, curve]) => [side, { points: curve.points, slopes: curve.slopes }])) } : {}),
      metricUpdate: discretization === 'giles-1985' ? 'each line' : 'each sweep' } };
  const pairedKernel = lineGrouping === 'boundary-pairs' && lineSearch === 'armijo' && !orthogonalBoundaryControl.previous;
  const scalarKernel = discretization === 'giles-1985' && !angleFeedback && !mapped && Object.values(boundaries).every(mode => mode === 'fixed');
  if (pairedKernel || scalarKernel)
    Object.defineProperty(system, 'wasmSmoothingConfig', { value: {
      nx, nt, xi: xi.slice(), eta: eta.slice(), lengthScale, rowGroups: rowGroups.map(row => row.slice()),
      xiDerivatives, derivatives, streamwiseSourceDiscretization,
      scalar: scalarKernel, mode: pairedKernel ? (physicalAngleSource ? 2 : 1) : stretched ? 3 : poissonP ? 4 : 0,
      boundary: pairedKernel ? structuredClone(orthogonalBoundaryControl) : { background: stretch ?? poissonP ?? initial.map(row => row.map(() => 0)), sides: [] },
    } });
  return system;
}

export function smoothEllipticStreamtubeGrid(system, { maxSweeps = 200, tolerance = 1e-9, omega = 1.3, requireConvex = false, detectCoordinateRoundoff = false, onSweep, backend = 'wasm' } = {}) {
  if (!Number.isInteger(maxSweeps) || maxSweeps < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isFinite(omega) || !(omega > 0 && omega < 2) || typeof requireConvex !== 'boolean'
    || typeof detectCoordinateRoundoff !== 'boolean') throw new Error('Invalid elliptic SLOR controls.');
  if (!['wasm', 'javascript'].includes(backend)) throw new Error('Unknown smoothing backend.');
  const kernel = backend === 'wasm' && system.wasmSmoothingConfig?.scalar ? createWasmSmoother(system) : null;
  let nodes = copy(system.initial), reason = 'sweep limit', converged = false, coordinatePrecision;
  let termination = { origin: 'solver', termination: 'sweep-limit' }; const history = [];
  for (let iteration = 0; iteration <= maxSweeps; iteration++) {
    try {
      const quality = system.quality(nodes);
      if (requireConvex && !quality.valid) { reason = 'Invalid starting state: folded grid'; termination = { origin: 'solver', termination: 'invalid-grid' }; break; }
      const state = kernel ? kernel.evaluate(nodes) : system.residuals(nodes);
      const { residual } = state, interiorResidual = kernel ? residual : state.interiorResidual, boundaryResidual = kernel ? 0 : state.boundaryResidual;
      if (!Number.isFinite(residual)) throw new Error('Nonfinite elliptic SLOR residual.');
      history.push({ iteration, residual, interiorResidual, boundaryResidual,
        minCornerSine: quality.minCornerSine, invalidCells: quality.invalidCells.length });
      onSweep?.(structuredClone(history.at(-1)), copy(nodes));
      if (residual <= tolerance) { converged = quality.valid; reason = converged ? 'converged' : 'folded converged grid'; break; }
      if (iteration === maxSweeps) break;
      let next = (kernel ?? system).sweep(nodes, omega);
      // Inspect the UNDAMPED correction, before convexity backtracking could
      // make an unconverged update artificially small. Retain the current
      // grid and its measured residual; never label this a tolerance root.
      if (detectCoordinateRoundoff && omega === 1 && quality.valid && next.maxUpdate <= Number.EPSILON) {
        const estimate = system.coordinateRoundoff?.(nodes, tolerance);
        if (estimate?.limited) {
          coordinatePrecision = { ...estimate, undampedUpdate: next.maxUpdate };
          reason = 'coordinate precision limit';
          termination = { origin: 'solver', termination: 'coordinate-precision-limit' };
          break;
        }
      }
      if (requireConvex) {
        const full = next, trials = []; let accepted = false;
        for (let halving = 0; halving <= 30; halving++) {
          const fraction = 2 ** -halving;
          const candidate = halving ? nodes.map((row, i) => row.map((p, j) => ({
            x: p.x + fraction * (full.nodes[i][j].x - p.x),
            y: p.y + fraction * (full.nodes[i][j].y - p.y),
          }))) : full.nodes;
          const q = system.quality(candidate);
          trials.push({ fraction, accepted: q.valid, invalidCells: q.invalidCells.length, minCornerSine: q.minCornerSine });
          if (!q.valid) continue;
          let maxUpdate = 0;
          candidate.forEach((row, i) => row.forEach((p, j) => {
            maxUpdate = Math.max(maxUpdate, Math.abs(p.x - nodes[i][j].x) / system.lengthScale,
              Math.abs(p.y - nodes[i][j].y) / system.lengthScale);
          }));
          next = { ...full, nodes: candidate, maxUpdate };
          Object.assign(history.at(-1), { fraction, convexTrials: trials }); accepted = true; break;
        }
        if (!accepted) throw slorNumericalError('No admissible convex SLOR correction.', 'no-admissible-convex-correction');
      }
      nodes = next.nodes; history.at(-1).maxUpdate = next.maxUpdate;
      if (next.maxUpdate === 0) { reason = 'roundoff stagnation'; termination = { origin: 'solver', termination: 'roundoff-stagnation' }; break; }
    } catch (error) {
      if (error?.code === 'slor-observer-failed') throw error;
      reason = error.message; termination = slorErrorTermination(error); break;
    }
  }
  const boundaryDescription = Object.values(system.coordinateEquations?.boundaryConditions ?? {}).includes('normal-curve')
    ? 'fixed walls/cuts/ends; selected farfield stations move on prescribed Hermite curves with normal adjacent segments'
    : Object.values(system.coordinateEquations?.boundaryConditions ?? {}).includes('giles-indexed-y')
    ? 'fixed inlet/outlet and each original indexed boundary y; selected boundary x copied from adjacent interior rows, without a fixed physical curve constraint'
    : Object.values(system.coordinateEquations?.boundaryConditions ?? {}).includes('giles-vertical')
    ? 'fixed inlet/outlet and boundary y; selected horizontal boundary x copied from adjacent interior rows' : 'fixed boundaries';
  const result = tagSlorTermination({ nodes, converged, reason, history, quality: system.quality(nodes), omega, tolerance,
    ...(coordinatePrecision ? { coordinatePrecision } : {}), ...(requireConvex ? { requireConvex } : {}),
    coordinateEquations: system.coordinateEquations,
    formulation: `${system.coordinateEquations?.xi === 'Poisson from prescribed harmonic-map coefficients'
      ? 'Inverse Poisson xi from prescribed harmonic-map coefficients and harmonic mass eta'
      : system.coordinateEquations?.xi === 'Poisson with implicit boundary-angle streamwise stretch' ? 'Inverse Poisson xi with implicit boundary-angle F and harmonic mass eta'
      : system.coordinateEquations?.xi === 'Poisson with implicit boundary-angle physical source' ? 'Inverse Poisson xi with implicit boundary-angle P and harmonic mass eta'
      : system.coordinateEquations?.xi === 'Poisson with prescribed metric-scaled streamwise stretch' ? 'Inverse Poisson xi with alpha-scaled prescribed F and harmonic mass eta'
          : system.coordinateEquations?.xi.startsWith('Poisson') ? 'Inverse Poisson xi with prescribed P and harmonic mass eta' : 'Inverse Laplace grid equations'}; streamwise-line SLOR with metrics frozen for ${system.coordinateEquations?.metricUpdate}; ${boundaryDescription}, ${system.coordinateEquations?.stationSpacing} and mass-fraction eta. Grid initialization only.` }, termination);
  Object.defineProperty(result, 'backend', { value: kernel ? 'wasm' : 'javascript' });
  Object.defineProperty(result, 'referenceReplays', { value: kernel?.replays ?? 0 });
  return result;
}
