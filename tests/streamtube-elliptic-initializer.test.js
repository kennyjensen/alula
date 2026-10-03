import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { relaxStreamtubeInitialGrid } from '../src/euler/streamtube-elliptic-initializer.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createEllipticStreamtubeGrid, gilesStreamwiseCoordinates } from '../src/geometry/elliptic-streamtube-grid.js';
import { createNormalGraphBoundary } from '../src/geometry/normal-graph-boundary.js';

test('accepted wall-angle adapter preserves the supplied grid and material boundaries and publishes independent snapshots', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 }));
  const before = system.decode(system.initial), chart = system.geometryChart(), snapshots = [];
  const result = relaxStreamtubeInitialGrid({ system }, {
    seed: 'supplied', boundaryControl: 'wall-angle', maxSweeps: 1,
    onSweep: (history, nodes) => snapshots.push({ history, nodes }),
  });
  assert.deepEqual(result.seeds, before.nodes);
  assert.deepEqual(system.decode(system.initial), before);
  assert.deepEqual(system.geometryChart(), chart);
  assert.equal(result.boundaryControl, 'wall-angle');
  for (let g = 0; g < before.nodes.length; g++) {
    const first = snapshots.find(s => s.history.region === g);
    assert.ok(first, `Region ${g} must publish its actual starting grid`);
    assert.deepEqual(first.nodes, before.nodes[g]);
    const equations = result.regions[g].coordinateEquations;
    assert.equal(equations.lineGrouping, 'boundary-pairs');
    for (const side of ['lower', 'upper']) {
      const body = system.layout.bodies[side === 'lower' ? g - 1 : g];
      const active = equations.activeBoundaryStations[side];
      active.forEach((value, i) => assert.equal(value, Boolean(i && i < system.layout.nx
        && (!body || i > body.leadingIndex && i < body.trailingIndex))));
    }
    result.nodes[g].forEach((row, i) => row.forEach((p, j) => {
      if (!i || i === system.layout.nx || !j || j === row.length - 1) assert.deepEqual(p, before.nodes[g][i][j]);
    }));
    first.nodes[1][1].x += 100;
    assert.notEqual(first.nodes[1][1].x, result.nodes[g][1][1].x);
  }
  const fresh = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 }));
  const adopted = fresh.decode(fresh.adoptGeometry(fresh.initial, result.nodes)).nodes;
  adopted.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.hypot(p.x - result.nodes[g][i][j].x, p.y - result.nodes[g][i][j].y) < 2e-13);
  })));
  for (const overrides of [{ seed: 'linear' }, { discretization: 'quadratic' }, { farfieldBoundary: 'giles-vertical' }])
    assert.throws(() => relaxStreamtubeInitialGrid({ system }, { seed: 'supplied', boundaryControl: 'wall-angle', ...overrides }), /Wall-angle SLOR/);
});

test('multielement elliptic initialization preserves every displaced wall, wake bank, cut and endpoint without changing the Euler chart', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 });
  input.displacement = { surfaces: input.bodies.map(b => ({
    upper: Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => .0001 * (1 + .2 * i)),
    lower: Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => .0001 * (1 + .1 * i)) })),
    wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(.0005)) };
  const system = createStreamtubeBodySystem(input), before = system.decode(system.initial), chart = system.geometryChart();
  const notifications = [];
  const r = relaxStreamtubeInitialGrid({ system }, { maxSweeps: 1, onSweep: h => notifications.push(h) });
  assert.equal(r.regions.length, 3); assert.equal(r.converged, false); assert.match(r.status, /not solved/);
  assert.deepEqual(system.decode(system.initial), before); assert.deepEqual(system.geometryChart(), chart);
  assert.deepEqual([...new Set(notifications.map(h => h.region))], [0, 1, 2]);
  const xi = gilesStreamwiseCoordinates(before.nodes[system.layout.primaryBody + 1].map(row => row[0]));
  assert.deepEqual(r.stationCoordinate.xi, xi);
  assert.equal(r.stationCoordinate.sharedAcrossPassages, true);
  assert.equal(r.stationCoordinate.exactModernMsetGauge, false);
  for (let g = 0; g < 3; g++) {
    const mass = before.allocation.groups[g].map(t => t.massFlow), total = mass.reduce((s, m) => s + m, 0);
    assert.equal(r.regions[g].coordinateEquations.discretization, 'giles-1985');
    assert.equal(r.regions[g].coordinateEquations.metricUpdate, 'each line');
    const direct = createEllipticStreamtubeGrid({ nodes: r.seeds[g], massFlows: mass,
      streamwiseCoordinates: xi, discretization: 'giles-1985' }).sweep(r.seeds[g]).nodes;
    assert.deepEqual(r.nodes[g], direct);
    let eta = 0;
    for (let j = 0; j <= mass.length; j++) {
      for (let i = 0; i <= system.layout.nx; i++) {
        const p = before.nodes[g][i][j], q = r.nodes[g][i][j];
        if (!i || i === system.layout.nx || !j || j === mass.length) assert.deepEqual(q, p);
        else for (const key of ['x', 'y']) assert.ok(Math.abs(r.seeds[g][i][j][key]
          - ((1 - eta) * before.nodes[g][i][0][key] + eta * before.nodes[g][i].at(-1)[key])) < 1e-13);
      }
      eta += (mass[j] ?? 0) / total;
    }
  }
  // Fresh chart adoption checks all two-bank displacement constraints,
  // independently of the smoother's fixed-boundary logic.
  const fresh = createStreamtubeBodySystem(input), adopted = fresh.adoptGeometry(fresh.initial, r.nodes), restored = fresh.decode(adopted).nodes;
  restored.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.hypot(p.x - r.nodes[g][i][j].x, p.y - r.nodes[g][i][j].y) < 2e-13);
  })));
});

test('elliptic adapter retains an explicit legacy discretization and supports prescribed common station labels', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 }));
  const legacy = relaxStreamtubeInitialGrid({ system }, { maxSweeps: 0, discretization: 'quadratic' });
  assert.equal(legacy.stationCoordinate.source, 'uniform station index');
  assert.ok(legacy.regions.every(r => r.coordinateEquations.discretization === 'quadratic'));
  const xi = Array.from({ length: system.layout.nx + 1 }, (_, i) => (i / system.layout.nx) ** 1.2);
  const prescribed = relaxStreamtubeInitialGrid({ system }, { maxSweeps: 1, streamwiseCoordinates: xi });
  assert.deepEqual(prescribed.stationCoordinate.xi, xi);
  xi[1] = .5;
  assert.notDeepEqual(prescribed.stationCoordinate.xi, xi);
  assert.throws(() => relaxStreamtubeInitialGrid({ system }, { discretization: 'unknown' }), /discretization/);
});

test('horizontal Giles farfields apply only to the two exterior rows and preserve the multielement chart and cuts', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 });
  const system = createStreamtubeBodySystem(input), before = system.decode(system.initial), chart = system.geometryChart();
  const result = relaxStreamtubeInitialGrid({ system }, { maxSweeps: 1, farfieldBoundary: 'giles-vertical' });
  assert.equal(result.farfieldBoundary, 'giles-vertical');
  assert.match(result.status, /Euler\/BL equations not solved/);
  assert.deepEqual(system.decode(system.initial), before);
  assert.deepEqual(system.geometryChart(), chart);
  const conditions = [{ lower: 'giles-vertical', upper: 'fixed' }, { lower: 'fixed', upper: 'fixed' },
    { lower: 'fixed', upper: 'giles-vertical' }];
  let moved = 0;
  for (let g = 0; g < 3; g++) {
    assert.deepEqual(result.regions[g].coordinateEquations.boundaryConditions, conditions[g]);
    const massFlows = before.allocation.groups[g].map(t => t.massFlow);
    const direct = createEllipticStreamtubeGrid({ nodes: result.seeds[g], massFlows,
      streamwiseCoordinates: result.stationCoordinate.xi, discretization: 'giles-1985',
      boundaryConditions: conditions[g] }).sweep(result.seeds[g]).nodes;
    assert.deepEqual(result.nodes[g], direct);
    result.nodes[g].forEach((row, i) => row.forEach((p, j) => {
      const original = before.nodes[g][i][j], nt = row.length - 1;
      if (!i || i === system.layout.nx) assert.deepEqual(p, original);
      else if (!j || j === nt) {
        if ((g === 0 && j === 0) || (g === 2 && j === nt)) {
          assert.equal(p.y, original.y);
          if (Math.abs(p.x - original.x) > 1e-10) moved++;
        } else assert.deepEqual(p, original);
      }
    }));
  }
  assert.ok(moved > 0, 'This fixture must actually exercise free outer x coordinates');
  const fresh = createStreamtubeBodySystem(input), adopted = fresh.adoptGeometry(fresh.initial, result.nodes);
  fresh.decode(adopted).nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.hypot(p.x - result.nodes[g][i][j].x, p.y - result.nodes[g][i][j].y) < 2e-13);
  })));
  const curved = structuredClone(before.nodes); curved[0][2][0].y += .001;
  assert.throws(() => relaxStreamtubeInitialGrid({ system, nodes: curved }, { farfieldBoundary: 'giles-vertical' }), /horizontal/);
  assert.throws(() => relaxStreamtubeInitialGrid({ system }, { farfieldBoundary: 'unknown' }), /boundary/);
});

test('indexed-y Giles adapter permits sloped exterior rows while retaining their station y and every other boundary', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 });
  for (const path of [input.outerLower, input.outerUpper]) for (const point of path) point.y += .04 * point.x;
  const system = createStreamtubeBodySystem(input), before = system.decode(system.initial), chart = system.geometryChart();
  const result = relaxStreamtubeInitialGrid({ system }, { maxSweeps: 1, farfieldBoundary: 'giles-indexed-y' });
  assert.equal(result.farfieldBoundary, 'giles-indexed-y');
  assert.match(result.status, /Giles indexed-y farfield/);
  assert.deepEqual(system.decode(system.initial), before);
  assert.deepEqual(system.geometryChart(), chart);
  const conditions = [{ lower: 'giles-indexed-y', upper: 'fixed' }, { lower: 'fixed', upper: 'fixed' },
    { lower: 'fixed', upper: 'giles-indexed-y' }];
  let moved = 0;
  result.nodes.forEach((group, g) => {
    assert.deepEqual(result.regions[g].coordinateEquations.boundaryConditions, conditions[g]);
    const massFlows = before.allocation.groups[g].map(t => t.massFlow);
    const direct = createEllipticStreamtubeGrid({ nodes: result.seeds[g], massFlows,
      streamwiseCoordinates: result.stationCoordinate.xi, discretization: 'giles-1985',
      boundaryConditions: conditions[g] }).sweep(result.seeds[g]).nodes;
    assert.deepEqual(group, direct);
    group.forEach((row, i) => row.forEach((point, j) => {
      const original = before.nodes[g][i][j], nt = row.length - 1;
      if (!i || i === system.layout.nx) assert.deepEqual(point, original);
      else if (!j || j === nt) {
        if ((g === 0 && j === 0) || (g === 2 && j === nt)) {
          assert.equal(point.y, original.y);
          if (Math.abs(point.x - original.x) > 1e-10) moved++;
        } else assert.deepEqual(point, original);
      }
    }));
  });
  assert.ok(moved > 0, 'The adapter fixture must exercise the indexed-y outer x copies.');
  const fresh = createStreamtubeBodySystem(input), state = fresh.adoptGeometry(fresh.initial, result.nodes);
  fresh.decode(state).nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((point, j) => {
    assert.ok(Math.hypot(point.x - result.nodes[g][i][j].x, point.y - result.nodes[g][i][j].y) < 2e-13);
  })));
  assert.throws(() => relaxStreamtubeInitialGrid({ system }, { farfieldBoundary: 'giles-vertical' }), /horizontal/);
});

test('curved farfield adapter preserves the prescribed graphs, cuts, walls and adoptable geometry', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 8, tubes: 3 });
  for (const path of [input.outerLower, input.outerUpper]) for (const point of path) point.y += .04 * point.x ** 2;
  const system = createStreamtubeBodySystem(input), initial = system.decode(system.initial), chart = system.geometryChart();
  const result = relaxStreamtubeInitialGrid({ system, guideField: { velocityAt: p => ({ u: 1, v: .08 * p.x }) } },
    { farfieldBoundary: 'normal-curve', maxSweeps: 1 });
  assert.deepEqual(system.decode(system.initial), initial); assert.deepEqual(system.geometryChart(), chart);
  let moved = 0;
  result.nodes.forEach((group, g) => {
    const nt = group[0].length - 1, side = g === 0 ? 'lower' : 'upper', jOuter = g === 0 ? 0 : nt;
    const spec = result.regions[g].coordinateEquations.boundaryCurves?.[side];
    const curve = spec ? createNormalGraphBoundary(spec) : null;
    if (curve) assert.deepEqual(spec.points, initial.nodes[g].map(row => row[jOuter]));
    group.forEach((row, i) => row.forEach((point, j) => {
      const before = initial.nodes[g][i][j];
      if (!i || i === system.layout.nx || (!j || j === nt) && !(curve && j === jOuter)) assert.deepEqual(point, before);
      if (curve && j === jOuter) {
        assert.ok(Math.abs(point.y - curve.evaluate(point.x).point.y) < 1e-13);
        if (Math.hypot(point.x - before.x, point.y - before.y) > 1e-10) moved++;
      }
    }));
  });
  assert.ok(moved > 0);
  const fresh = createStreamtubeBodySystem(input), state = fresh.adoptGeometry(fresh.initial, result.nodes);
  fresh.decode(state).nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((point, j) => {
    assert.ok(Math.hypot(point.x - result.nodes[g][i][j].x, point.y - result.nodes[g][i][j].y) < 2e-13);
  })));
  assert.throws(() => relaxStreamtubeInitialGrid({ system }, { farfieldBoundary: 'normal-curve' }), /velocity field/);
});
