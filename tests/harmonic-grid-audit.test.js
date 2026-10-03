import test from 'node:test';
import assert from 'node:assert/strict';
import { quadLaplaceMatrix, solveHarmonicGridReference, auditHarmonicGrid } from '../src/geometry/harmonic-grid-audit.js';

const rectangular = (nx, nt, perturb = 0) => ({ massFlows: Array(nt).fill(1 / nt),
  nodes: Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
    const x = i / nx, y = j / nt, d = !i || !j || i === nx || j === nt ? 0 : perturb * Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
    return { x: x + d, y: y - .6 * d };
  })) });

test('nonuniform streamwise labels preserve the exact affine physical field under reference refinement', () => {
  const xi = [0, .05, .2, .6, 1], eta = [0, .2, .5, 1];
  const region = { streamwiseCoordinates: xi, massFlows: [.2, .3, .5],
    nodes: xi.map(x => eta.map(y => ({ x, y }))) }, saved = structuredClone(region);
  const result = solveHarmonicGridReference(region, { refinement: 2, includeRefinedField: true });
  assert.deepEqual(result.streamwiseCoordinates, xi);
  assert.notEqual(result.streamwiseCoordinates, xi);
  for (const row of result.refinedField) for (const p of row) {
    assert.ok(Math.abs(p.xi - p.x) < 3e-14);
    assert.ok(Math.abs(p.eta - p.y) < 3e-14);
  }
  assert.equal(result.refinedField[1][0].xi, .025);
  assert.equal(result.refinedField[3][0].xi, .125);
  assert.ok(result.maximum.crosslineIntervals < 1e-12);
  assert.ok(result.maximum.tubeIntervals < 1e-12);
  assert.deepEqual(region, saved);
});

test('nonuniform crossline errors use the smaller adjacent supplied Xi interval', () => {
  const xi = [0, .05, .2, .6, 1], eta = [0, .2, .5, 1];
  const region = { streamwiseCoordinates: xi, massFlows: [.2, .3, .5], nodes: xi.map((x, i) => eta.map((y, j) => {
    const interior = i > 0 && i < xi.length - 1 && j > 0 && j < eta.length - 1;
    return { x: x + (interior ? .04 * Math.sin(Math.PI * x) * Math.sin(Math.PI * y) : 0), y };
  })) };
  for (const refinement of [1, 2]) {
    const result = solveHarmonicGridReference(region, { refinement });
    region.nodes.forEach((row, i) => row.forEach((p, j) => {
      // Boundary labels are physical x/y, so their exact harmonic extension
      // stays x/y even though the interior nodes carry different labels.
      const width = Math.min(i ? xi[i] - xi[i - 1] : Infinity, i < xi.length - 1 ? xi[i + 1] - xi[i] : Infinity);
      assert.ok(Math.abs(result.values[i][j].xi - p.x) < 3e-14);
      assert.ok(Math.abs(result.errors[i][j].xi - (p.x - xi[i])) < 3e-14);
      assert.ok(Math.abs(result.errors[i][j].crosslineIntervals - (p.x - xi[i]) / width) < 5e-13);
    }));
  }
});

test('nonuniform streamwise labels reject unresolved or mismatched coordinates', () => {
  const region = rectangular(4, 3);
  for (const xi of [[0, .2, .7, 1], [.01, .1, .3, .7, 1], [0, .1, .3, .7, .99],
    [0, .1, .1, .7, 1], [0, .3, .2, .7, 1], [0, .1, NaN, .7, 1], [0, .1, , .7, 1], null])
    assert.throws(() => solveHarmonicGridReference({ ...region, streamwiseCoordinates: xi }), /streamwise coordinates/);
  const result = solveHarmonicGridReference(region);
  assert.deepEqual(result.streamwiseCoordinates, [0, .25, .5, .75, 1]);
  assert.ok(result.maximum.crosslineIntervals < 1e-12);
});

test('quadrilateral Laplace matrix matches exact unit-square integrals and preserves rigid coordinates', () => {
  const p = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 1 }, { x: 0, y: 1 }];
  const k = quadLaplaceMatrix(p), exact = [4, -1, -2, -1, -1, 4, -1, -2, -2, -1, 4, -1, -1, -2, -1, 4];
  k.forEach((v, i) => assert.ok(Math.abs(v - exact[i] / 6) < 5e-16));
  const a = .7, mapped = p.map(q => ({ x: 4 + 2.7 * (q.x * Math.cos(a) - q.y * Math.sin(a)), y: -3 + 2.7 * (q.x * Math.sin(a) + q.y * Math.cos(a)) }));
  const transformed = quadLaplaceMatrix(mapped); k.forEach((v, i) => assert.ok(Math.abs(v - transformed[i]) < 2e-15));
  assert.throws(() => quadLaplaceMatrix([p[0], p[2], p[1], p[3]]), /Invalid/);
});

test('independent harmonic reference detects known coordinate errors in positive distorted grids', () => {
  const region = rectangular(12, 8, .06), saved = structuredClone(region);
  for (const refinement of [1, 2, 4]) {
    const r = solveHarmonicGridReference(region, { refinement });
    region.nodes.forEach((row, i) => row.forEach((p, j) => {
      // Exact continuum harmonic coordinates are physical x and y, whatever
      // the interior node placement. Q1 reproduces these linear fields.
      assert.ok(Math.abs(r.values[i][j].xi - p.x) < 2e-12);
      assert.ok(Math.abs(r.values[i][j].eta - p.y) < 2e-12);
    }));
  }
  const audit = auditHarmonicGrid(region);
  assert.equal(audit.passed, false); assert.equal(audit.referenceResolved, true);
  assert.ok(audit.levels.at(-1).maximum.crosslineIntervals > .7);
  assert.ok(audit.levels.at(-1).maximum.tubeIntervals > .28);
  assert.deepEqual(region, saved);
  const exact = auditHarmonicGrid(rectangular(8, 6));
  assert.equal(exact.passed, true);
  assert.throws(() => solveHarmonicGridReference(region, { refinement: 4, maxUnknowns: 10 }), /budget/);
});

test('independent quadrilateral reference approaches analytic vortex harmonic coordinates under grid refinement', t => {
  let previous = Infinity; const results = [];
  // Refine the polygonal boundary as well as the interior. The exact log(r)
  // is not linear along its chordal edges; the first very coarse grid is
  // outside the observed asymptotic range at the inlet-adjacent nodes.
  for (const [nx, nt] of [[16, 8], [32, 16], [64, 32]]) {
    const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
      const angle = -.4 + .8 * i / nx, radius = 2 * Math.exp(-.5 * j / nt);
      return { x: radius * Math.cos(angle), y: radius * Math.sin(angle) };
    }));
    const r = solveHarmonicGridReference({ nodes, massFlows: Array(nt).fill(1 / nt) }, { refinement: 2 });
    let error = 0;
    r.values.forEach((row, i) => row.forEach((p, j) => { error = Math.max(error, Math.abs(p.xi - i / nx), Math.abs(p.eta - j / nt)); }));
    assert.ok(error < .32 * previous); previous = error; results.push({ nx, nt, error });
  }
  assert.ok(previous < 4e-5); t.diagnostic(JSON.stringify(results));
});

test('mixed harmonic reference reproduces affine fields and solves genuinely free boundary xi values', () => {
  const xi = [0, .2, .5, .8, 1], eta = [0, .18, .56, 1];
  for (const boundaryConditions of [{ lower: 'giles-vertical' }, { upper: 'giles-vertical' },
    { lower: 'giles-vertical', upper: 'giles-vertical' }]) {
    const lower = boundaryConditions.lower === 'giles-vertical', upper = boundaryConditions.upper === 'giles-vertical';
    for (const perturb of [0, .05]) {
      const region = { streamwiseCoordinates: xi, massFlows: [.18, .38, .44], boundaryConditions,
        nodes: xi.map(u => eta.map(v => {
          const factor = lower && upper ? .5 + .5 * v : lower ? 1 - v : v;
          return { x: 3 + 2 * (u + perturb * Math.sin(Math.PI * u) * factor), y: -.7 + .4 * v };
        })) };
      const before = structuredClone(region);
      const result = solveHarmonicGridReference(region, { refinement: 2, includeRefinedField: true });
      // The domain is always the same rectangle. Its exact solution is
      // xi=(x-3)/2, eta=(y+.7)/.4, including the natural xi boundaries.
      for (const row of result.refinedField) for (const p of row) {
        assert.ok(Math.abs(p.xi - (p.x - 3) / 2) < 2e-13);
        assert.ok(Math.abs(p.eta - (p.y + .7) / .4) < 2e-13);
      }
      result.errors.forEach((row, i) => row.forEach((error, j) => {
        if (i === 0 || i === xi.length - 1) assert.equal(error.xi, 0);
        if (i === 0 || i === xi.length - 1 || j === 0 || j === eta.length - 1) assert.equal(error.eta, 0);
        if ((j === 0 && !lower) || (j === eta.length - 1 && !upper)) assert.equal(error.xi, 0);
      }));
      const freeJ = upper ? eta.length - 1 : 0;
      if (perturb) assert.ok(result.errors[2][freeJ].xi > .049, 'the natural boundary value must actually be solved');
      else assert.ok(result.maximum.crosslineIntervals < 1e-12);
      assert.deepEqual(result.unknownsByCoordinate, { xi: 7 * (5 + Number(lower) + Number(upper)), eta: 35 });
      assert.equal(result.unknowns, result.unknownsByCoordinate.xi);
      assert.deepEqual(result.boundaryConditions, { lower: 'fixed', upper: 'fixed', ...boundaryConditions });
      assert.deepEqual(region, before);
    }
  }
});

test('mixed reference preserves the existing fixed-boundary path and enforces the larger xi budget', () => {
  const region = rectangular(4, 3, .025);
  const defaults = solveHarmonicGridReference(region, { refinement: 2, includeRefinedField: true });
  const explicit = solveHarmonicGridReference({ ...region, boundaryConditions: { lower: 'fixed', upper: 'fixed' } },
    { refinement: 2, includeRefinedField: true });
  assert.deepEqual(defaults, explicit);
  assert.deepEqual(defaults.unknownsByCoordinate, { xi: 35, eta: 35 });
  const mixed = { ...region, boundaryConditions: { lower: 'giles-vertical', upper: 'giles-vertical' } };
  assert.throws(() => solveHarmonicGridReference(mixed, { refinement: 2, maxUnknowns: 35 }), /xi.*49.*budget/i);
  const allowed = solveHarmonicGridReference(mixed, { refinement: 2, maxUnknowns: 49 });
  assert.equal(allowed.unknowns, 49);
  assert.deepEqual(allowed.unknownsByCoordinate, { xi: 49, eta: 35 });
  // Eta has identical unknowns, geometry, quadrature and boundary data.
  assert.deepEqual(allowed.values.map(row => row.map(p => p.eta)), defaults.values.map(row => row.map(p => p.eta)));
});

test('mixed harmonic audit retains free-boundary errors in refinement and acceptance reports', () => {
  const region = rectangular(4, 3);
  region.boundaryConditions = { lower: 'giles-vertical', upper: 'giles-vertical' };
  region.nodes.forEach((row, i) => {
    if (i > 0 && i < 4) for (const j of [0, 3]) row[j].x += .04 * Math.sin(Math.PI * i / 4);
  });
  // The exact Xi=x differs from the assigned coordinate only on the free
  // boundaries. Dropping those corrections would falsely pass this audit.
  const result = auditHarmonicGrid(region);
  assert.equal(result.referenceResolved, true);
  assert.equal(result.passed, false);
  assert.ok(Math.abs(result.errors[2][3].xi - .04) < 1e-13);
  for (let i = 1; i < 4; i++) for (let j = 1; j < 3; j++) assert.ok(Math.abs(result.errors[i][j].xi) < 1e-13);
  assert.ok(result.levels.at(-1).maximum.crosslineIntervals > .15);
  assert.match(result.method, /normal-Neumann xi/);
  assert.match(result.scope, /mixed coordinate boundary conditions/);
  assert.ok(result.levels.every(level => level.unknownsByCoordinate.xi > level.unknownsByCoordinate.eta));
});

test('mixed harmonic reference approaches a known nonlinear harmonic field with a zero normal farfield derivative', t => {
  const amplitude = .08, results = [];
  const exactXi = (x, y) => x + amplitude * Math.sin(Math.PI * x) * Math.cosh(Math.PI * (1 - y)) / Math.cosh(Math.PI);
  // xi(x,0) increases strictly: its derivative is at least 1-amplitude*pi.
  // Bisection fixes the physical bottom hits for the prescribed Xi labels.
  const bottomX = u => {
    if (u === 0 || u === 1) return u;
    let lo = 0, hi = 1;
    for (let k = 0; k < 60; k++) {
      const mid = (lo + hi) / 2;
      if (exactXi(mid, 0) < u) lo = mid; else hi = mid;
    }
    return (lo + hi) / 2;
  };
  for (const [nx, nt] of [[8, 4], [16, 8], [32, 16]]) {
    const region = { massFlows: Array(nt).fill(1 / nt), boundaryConditions: { upper: 'giles-vertical' },
      nodes: Array.from({ length: nx + 1 }, (_, i) => {
        const u = i / nx, hit = bottomX(u);
        return Array.from({ length: nt + 1 }, (_, j) => {
          const y = j / nt;
          return { x: hit + y * (u - hit), y };
        });
      }) };
    // The rectangle is exact at every resolution. Its bottom Dirichlet
    // interpolant is refined with the mesh; this is continuum convergence,
    // not a claim that nested refinement removes fixed boundary-data error.
    const result = solveHarmonicGridReference(region, { maxUnknowns: 600 });
    let xiError = 0, etaError = 0;
    result.values.forEach((row, i) => row.forEach((value, j) => {
      const p = region.nodes[i][j];
      xiError = Math.max(xiError, Math.abs(value.xi - exactXi(p.x, p.y)));
      etaError = Math.max(etaError, Math.abs(value.eta - p.y));
    }));
    assert.ok(etaError < 1e-13);
    assert.ok(result.errors[nx / 2][nt].xi > .006, 'nonconstant free-boundary correction remains measurable');
    assert.ok(result.linear.xi.relativeResidual < 1e-12);
    results.push({ nx, nt, xiError, etaError, unknowns: result.unknowns });
  }
  for (let k = 1; k < results.length; k++) {
    const ratio = results[k].xiError / results[k - 1].xiError;
    assert.ok(ratio > .15 && ratio < .4, `expected second-order refinement, error ratio ${ratio}`);
  }
  assert.ok(results.at(-1).xiError < 1e-4);
  t.diagnostic(JSON.stringify({ results, scope: 'Independent Q1 mixed-boundary harmonic reference; no SLOR or airfoil solution.' }));
});

test('mixed harmonic reference rejects invalid conditions and nonhorizontal selected farfields', () => {
  const unspecified = rectangular(4, 4);
  assert.deepEqual(solveHarmonicGridReference({ ...unspecified, boundaryConditions: { lower: undefined } }),
    solveHarmonicGridReference(unspecified));
  const region = rectangular(4, 3);
  for (const boundaryConditions of [null, [], 'giles-vertical', { lower: 'unknown' }, { upper: true }, { left: 'fixed' }])
    assert.throws(() => solveHarmonicGridReference({ ...region, boundaryConditions }), /boundary conditions/i);
  for (const side of ['lower', 'upper']) {
    const j = side === 'lower' ? 0 : 3;
    const bowed = structuredClone(region); bowed.nodes[2][j].y += 1e-15;
    assert.throws(() => solveHarmonicGridReference({ ...bowed, boundaryConditions: { [side]: 'giles-vertical' } }), /horizontal/i);
    const sloped = structuredClone(region); sloped.nodes.forEach(row => { row[j].y += .1 * row[j].x; });
    assert.throws(() => solveHarmonicGridReference({ ...sloped, boundaryConditions: { [side]: 'giles-vertical' } }), /horizontal/i);
    assert.doesNotThrow(() => solveHarmonicGridReference({ ...sloped, boundaryConditions: { [side]: 'fixed' } }));
  }
});

test('normal-curve natural boundary conditions reproduce a rotated affine harmonic field', () => {
  const slope = .3, nx = 4, nt = 3;
  for (const sides of [['lower'], ['upper'], ['lower', 'upper']]) {
    const region = { massFlows: Array(nt).fill(1 / nt), boundaryConditions: {}, boundaryCurves: {},
      nodes: Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
        const v = j / nt, free = sides.includes(j === 0 ? 'lower' : j === nt ? 'upper' : 'interior');
        const u = i / nx + (i > 0 && i < nx && free ? .035 * Math.sin(Math.PI * i / nx) : 0);
        return { x: u - slope * v, y: slope * u + v };
      })) };
    for (const side of sides) {
      const j = side === 'lower' ? 0 : nt;
      region.boundaryConditions[side] = 'normal-curve';
      region.boundaryCurves[side] = { points: [region.nodes[0][j], region.nodes[nx][j]], slopes: [slope, slope] };
    }
    const before = structuredClone(region);
    const result = solveHarmonicGridReference(region, { refinement: 2, includeRefinedField: true });
    for (const row of result.refinedField) for (const p of row) {
      assert.ok(Math.abs(p.xi - (p.x + slope * p.y) / (1 + slope * slope)) < 2e-13);
      assert.ok(Math.abs(p.eta - (p.y - slope * p.x) / (1 + slope * slope)) < 2e-13);
    }
    for (const side of sides) {
      const j = side === 'lower' ? 0 : nt;
      assert.ok(result.errors[2][j].xi > .034, 'free boundary Xi must be solved, not left at its seed label');
      assert.equal(result.boundaryCurves[side].kind, 'piecewise-cubic-hermite-graph');
      assert.notEqual(result.boundaryCurves[side].points, region.boundaryCurves[side].points);
    }
    assert.deepEqual(region, before);
    assert.match(result.geometryRefinement, /No exact curved-cell integration/);
  }
});

test('normal-curve refinement follows the fixed graphs and preserves original nodes and fixed ends', () => {
  const nx = 4, nt = 3, lower = x => -.12 * x * (1 - x), upper = x => 1 + .24 * x * (1 - x);
  const points = f => [0, .5, 1].map(x => ({ x, y: f(x) }));
  const region = { massFlows: [.15, .25, .6], boundaryConditions: { lower: 'normal-curve', upper: 'normal-curve' },
    boundaryCurves: { lower: { points: points(lower), slopes: [-.12, 0, .12] },
      upper: { points: points(upper), slopes: [.24, 0, -.24] } },
    nodes: Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
      const x = i / nx, fraction = j / nt;
      return { x, y: (1 - fraction) * lower(x) + fraction * upper(x) };
    })) };
  const before = structuredClone(region), result = solveHarmonicGridReference(region, { refinement: 2, includeRefinedField: true });
  result.refinedField.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.abs(p.y - ((1 - j / (2 * nt)) * lower(p.x) + j / (2 * nt) * upper(p.x))) < 3e-15);
    if (i % 2 === 0 && j % 2 === 0) {
      assert.equal(p.x, region.nodes[i / 2][j / 2].x);
      assert.equal(p.y, region.nodes[i / 2][j / 2].y);
    }
  }));
  for (const i of [0, 2 * nx]) for (let j = 0; j <= 2 * nt; j++) {
    assert.equal(result.refinedField[i][j].x, i / (2 * nx));
    assert.ok(Math.abs(result.refinedField[i][j].y - j / (2 * nt)) < 2e-15);
  }
  assert.deepEqual(region, before);
  assert.deepEqual(solveHarmonicGridReference({ ...region, boundaryCurves: result.boundaryCurves }, { refinement: 2 }),
    solveHarmonicGridReference(region, { refinement: 2 }), 'normalized descriptors must round-trip');
  const audit = auditHarmonicGrid(region, { refinements: [1, 2] });
  assert.match(audit.method, /Q1 polygon approximations of fixed cubic Hermite graph/);
  assert.match(audit.method, /no exact curved-cell integration/);
  assert.match(audit.scope, /not the graph approximation or fixed end-data error/);
  assert.ok(audit.levels.every(level => level.boundaryCurves.upper.kind === 'piecewise-cubic-hermite-graph'));
});

test('natural-Xi refinement on fixed Hermite farfields approaches an analytic conformal field', t => {
  const nx = 8, nt = 4, height = .6, amplitude = .025;
  // z=w+a*sin(pi*w), w=u+i*v. Lines v=constant are streamlines;
  // their normals have zero derivative of harmonic Xi=u. Re(z) increases
  // along each boundary since 1-a*pi*cosh(pi*height)>0.
  const at = (u, v) => ({ x: u + amplitude * Math.sin(Math.PI * u) * Math.cosh(Math.PI * v),
    y: v + amplitude * Math.cos(Math.PI * u) * Math.sinh(Math.PI * v) });
  const parameter = (i, j) => i / nx + (i && i < nx ? .04 * Math.sin(Math.PI * i / nx) * j / nt : 0);
  const region = { massFlows: Array(nt).fill(1 / nt), boundaryConditions: { lower: 'normal-curve', upper: 'normal-curve' },
    boundaryCurves: {}, nodes: Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => at(parameter(i, j), height * j / nt))) };
  for (const [side, j] of [['lower', 0], ['upper', nt]]) {
    const v = height * j / nt;
    region.boundaryCurves[side] = { points: region.nodes.map(row => row[j]), slopes: region.nodes.map((row, i) => {
      const u = parameter(i, j);
      return -amplitude * Math.PI * Math.sin(Math.PI * u) * Math.sinh(Math.PI * v)
        / (1 + amplitude * Math.PI * Math.cos(Math.PI * u) * Math.cosh(Math.PI * v));
    }) };
  }
  const results = [];
  for (const refinement of [1, 2, 4]) {
    const result = solveHarmonicGridReference(region, { refinement, maxUnknowns: 600 });
    let xiError = 0;
    result.values.forEach((row, i) => row.forEach((p, j) => { xiError = Math.max(xiError, Math.abs(p.xi - parameter(i, j))); }));
    assert.ok(result.errors[nx / 2][nt].xi > .03, 'free farfield labels must change toward the analytic values');
    assert.ok(result.linear.xi.relativeResidual < 1e-12);
    results.push({ refinement, xiError, unknowns: result.unknowns });
  }
  for (let k = 1; k < results.length; k++) assert.ok(results[k].xiError < .6 * results[k - 1].xiError,
    `Expected decreasing geometric/field discretization error: ${JSON.stringify(results)}`);
  assert.ok(results.at(-1).xiError < 5e-4);
  t.diagnostic(JSON.stringify({ results, scope: 'Xi only: fixed Hermite graphs approximate the analytic conformal boundaries. Nested refinement removes chord/FE error, not the fixed Hermite approximation or eta end-data error.' }));
});

test('normal-curve validation rejects mismatched data and folded refined geometry without changing old modes', () => {
  const region = rectangular(4, 3, .015), flat = { points: [{ x: 0, y: 1 }, { x: 1, y: 1 }], slopes: [0, 0] };
  const normal = { ...region, boundaryConditions: { upper: 'normal-curve' }, boundaryCurves: { upper: flat } };
  for (const boundaryCurves of [null, [], 'curve', { left: flat }, { lower: flat }, {}, { upper: null }, { upper: [] }])
    assert.throws(() => solveHarmonicGridReference({ ...normal, boundaryCurves }), /curve|graph/i);
  assert.throws(() => solveHarmonicGridReference({ ...region, boundaryCurves: { upper: flat } }), /curves/i);
  for (const i of [0, 2, 4]) {
    const invalid = structuredClone(normal); invalid.nodes[i][3].y += 1e-7;
    assert.throws(() => solveHarmonicGridReference(invalid), /node.*curve/i);
  }
  const extended = structuredClone(normal); extended.boundaryCurves.upper.points[0].x = -.1;
  assert.throws(() => solveHarmonicGridReference(extended), /endpoints/i);
  const reversed = structuredClone(normal); reversed.nodes[2][3].x = reversed.nodes[1][3].x;
  assert.throws(() => solveHarmonicGridReference(reversed), /ordered/i);
  const folding = { ...rectangular(2, 2), boundaryConditions: { upper: 'normal-curve' },
    boundaryCurves: { upper: { points: [{ x: 0, y: 1 }, { x: .5, y: 1 }, { x: 1, y: 1 }], slopes: [-20, 20, 0] } } };
  assert.throws(() => solveHarmonicGridReference(folding, { refinement: 2 }), /Invalid quadrilateral/i);
  for (const boundaryConditions of [{}, { upper: 'giles-vertical' }]) {
    const old = { ...region, boundaryConditions };
    assert.deepEqual(solveHarmonicGridReference(old, { refinement: 2 }),
      solveHarmonicGridReference({ ...old, boundaryCurves: {} }, { refinement: 2 }));
  }
});
