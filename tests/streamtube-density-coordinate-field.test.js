import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { prolongPhysicalLogDensity as prolong, prolongStreamtubeDensities } from '../src/euler/streamtube-density-prolongation.js';

const archive = new URL('../docs/surface-pchip-refinement/before/src/euler/streamtube-density-prolongation.js', import.meta.url);
const original = fs.readFileSync(archive, 'utf8'), dependency = new URL('../src/euler/streamtube-cell.js', import.meta.url).href;
const archivedSource = original.replace("'./streamtube-cell.js'", JSON.stringify(dependency));
assert.equal(archivedSource.replace(JSON.stringify(dependency), "'./streamtube-cell.js'"), original);
const { prolongPhysicalLogDensity: archived } = await import(`data:text/javascript;base64,${Buffer.from(archivedSource).toString('base64')}`);
const close = (a, b, tolerance = 4e-14) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const centers = masses => masses.map(row => { let sum = 0; return row.map(m => { const center = sum + m / 2; sum += m; return center; }); });
const affine = (s, m, g) => .03125 * s + .015625 * m + .125 * g;

function fixture({ mapped = true, value = affine, gamma = 1.4 } = {}) {
  const nx = 4, masses = [[.4, .9, 1.7], [.6, 1.4]], subdivisions = [[2, 1, 3], [1, 2]], mc = centers(masses);
  const logDensity = Array.from({ length: nx }, (_, i) => masses.map((row, g) => row.map((_, j) => value(i + .5, mc[g][j], g))));
  const source = { masses, logDensity,
    normalAreas: logDensity.map(row => row.map((group, g) => group.map((v, j) => masses[g][j] / (Math.exp(v) * .8)))) };
  const childMasses = masses.map((row, g) => row.flatMap((m, j) => Array(subdivisions[g][j]).fill(m / subdivisions[g][j])));
  const nodeCoordinates = [0, .3, .8, 1.4, 2, 2.7, 3.2, 3.7, 4];
  const target = { nodeCoordinates, subdivisions, masses: childMasses,
    normalAreas: nodeCoordinates.slice(1).map(() => childMasses.map(row => row.map(m => m / .75))) };
  if (mapped) target.nodeCoordinatesByGroup = childMasses.map((row, g) => nodeCoordinates.map((s, i) =>
    Array.from({ length: row.length + 1 }, (_, j) => i === 0 || i === nodeCoordinates.length - 1 ? s
      : s + .12 * Math.sin(Math.PI * s / nx) * (j / row.length - .35 + .2 * g))));
  return { source, target, h0: 8, gamma };
}

test('affine source log density follows each four-corner coordinate and unequal physical mass center', () => {
  const f = fixture(), before = structuredClone(f), r = prolong(f), parentMass = centers(f.source.masses), childMass = centers(f.target.masses);
  let maximumPlaneDifference = 0;
  r.logDensity.forEach((row, i) => row.forEach((group, g) => group.forEach((v, j) => {
    const field = f.target.nodeCoordinatesByGroup[g];
    const s = (field[i][j] + field[i][j + 1] + field[i + 1][j] + field[i + 1][j + 1]) / 4;
    const boundedS = Math.max(.5, Math.min(3.5, s));
    const m = Math.max(parentMass[g][0], Math.min(parentMass[g].at(-1), childMass[g][j]));
    close(r.targetSectionCoordinates[i][g][j], s);
    close(v, affine(boundedS, m, g));
    maximumPlaneDifference = Math.max(maximumPlaneDifference, Math.abs(v - affine(Math.max(.5, Math.min(3.5,
      (f.target.nodeCoordinates[i] + f.target.nodeCoordinates[i + 1]) / 2)), m, g)));
    const rho = Math.exp(v), mass = f.target.masses[g][j], area = f.target.normalAreas[i][g][j];
    const q = mass / (rho * area), h = f.h0 - q * q / 2, p = (f.gamma - 1) / f.gamma * rho * h;
    close(rho * q * area, mass); close(f.gamma / (f.gamma - 1) * p / rho + q * q / 2, f.h0);
    assert.ok(rho > 0 && h > 0 && p > 0);
  })));
  assert.ok(maximumPlaneDifference > .001, 'The coordinate field must affect sampled density.');
  assert.equal(r.diagnostics.commonIsentropeInversions, 0); assert.deepEqual(f, before);
  r.nodeCoordinatesByGroup[0][1][0] = 99; r.targetSectionCoordinates[0][0][0] = 99; r.logDensity[0][0][0] = 99;
  assert.deepEqual(f, before);
});

test('common station planes retain exact old center arithmetic and complete transfer values', () => {
  const f = fixture({ mapped: false }), baseline = archived(f);
  f.target.nodeCoordinatesByGroup = f.target.masses.map(row => f.target.nodeCoordinates.map(s => Array(row.length + 1).fill(s)));
  const actual = prolong(f);
  for (const [key, value] of Object.entries(baseline)) {
    if (key === 'targetSectionCoordinates') actual[key].forEach((row, i) => row.forEach(group => group.forEach(s => assert.equal(s, value[i]))));
    else assert.deepEqual(actual[key], value, key);
  }
});

test('omitted field preserves every archived result field and default failure exactly', () => {
  for (const value of [affine, () => .125, (s, m, g) => .11 * Math.sin(s + m) + .03 * g]) for (const gamma of [1.4, 5 / 3]) {
    const f = fixture({ mapped: false, value, gamma });
    assert.deepEqual(prolong(f), archived(f));
    assert.equal(Object.hasOwn(prolong(f), 'nodeCoordinatesByGroup'), false);
  }
  const f = fixture({ mapped: false }); f.target.masses[0][0] *= 1.01;
  const failure = fn => { try { fn(f); } catch (e) { return { message: e.message, name: e.name }; } assert.fail('Expected failure'); };
  assert.deepEqual(failure(prolong), failure(archived));
  const manifest = JSON.parse(fs.readFileSync(new URL('../docs/surface-pchip-refinement/before/density-manifest.json', import.meta.url)));
  assert.equal(createHash('sha256').update(original).digest('hex'), manifest.sha256['src/euler/streamtube-density-prolongation.js']);
});

test('field rejects malformed dimensions, holes, nonfinite/bounded values and nonincreasing rails atomically', () => {
  const mutations = [
    f => { f.target.nodeCoordinatesByGroup = null; },
    f => f.target.nodeCoordinatesByGroup.pop(),
    f => { delete f.target.nodeCoordinatesByGroup[0]; },
    f => f.target.nodeCoordinatesByGroup[0].pop(),
    f => f.target.nodeCoordinatesByGroup[0][1].pop(),
    f => { delete f.target.nodeCoordinatesByGroup[0][1][0]; },
    f => { f.target.nodeCoordinatesByGroup[0][1][0] = NaN; },
    f => { f.target.nodeCoordinatesByGroup[0][1][0] = Infinity; },
    f => { f.target.nodeCoordinatesByGroup[0][1][0] = -.1; },
    f => { f.target.nodeCoordinatesByGroup[0][1][0] = 4.1; },
    f => { f.target.nodeCoordinatesByGroup[0][0][0] = .01; },
    f => { f.target.nodeCoordinatesByGroup[0].at(-1)[0] = 3.99; },
    f => { f.target.nodeCoordinatesByGroup[1][2][1] = f.target.nodeCoordinatesByGroup[1][1][1]; },
    f => { f.target.nodeCoordinatesByGroup[1][2][1] = .1; },
  ];
  for (const mutate of mutations) {
    const f = fixture(); mutate(f); const before = structuredClone(f);
    assert.throws(() => prolong(f), /per-node density source-coordinate field/); assert.deepEqual(f, before);
  }
});

test('opt-in coordinates preserve the unchanged mass-partition and physical gas rejection gates', () => {
  for (const [mutate, message] of [
    [f => { f.target.masses[0][0] *= 1.01; }, /parent tube mass/],
    [f => { f.target.normalAreas[2][0][1] = 1e-8; }, /inadmissible/],
    [f => { f.target.normalAreas[2][0][1] = -1; }, /inadmissible/],
  ]) {
    const f = fixture(); mutate(f); const before = structuredClone(f);
    assert.throws(() => prolong(f), message); assert.deepEqual(f, before);
  }
});

test('layout adapter forwards the field, changes only density slots and performs no residual evaluation', () => {
  const f = fixture(), sourceCounts = f.source.masses.map(row => row.length), counts = f.target.masses.map(row => row.length);
  const layout = (nx, tubes) => ({ nx, tubes,
    densityIndex: (i, g, j) => i * tubes.reduce((a, b) => a + b, 0) + tubes.slice(0, g).reduce((a, b) => a + b, 0) + j });
  const conditions = { mach: .3, gamma: f.gamma, h0: f.h0, lengthScale: 1, massScale: 5 };
  const source = { conditions, layout: layout(4, sourceCounts) }, sourceState = Float64Array.from(f.source.logDensity.flat(2));
  const sourceFlow = { allocation: { groups: f.source.masses.map(row => row.map(massFlow => ({ massFlow }))) },
    cells: Array.from({ length: 3 }, (_, i) => sourceCounts.map((n, g) => Array.from({ length: n }, (_, j) =>
      ({ geometry: { normalAreas: [f.source.normalAreas[i][g][j], f.source.normalAreas[i + 1][g][j]] } })))) };
  const nx = f.target.nodeCoordinates.length - 1, targetLayout = layout(nx, counts);
  const nodes = f.target.masses.map((row, g) => {
    const ys = [g * 4]; for (const m of row) ys.push(ys.at(-1) + m);
    return f.target.nodeCoordinates.map(x => ys.map(y => ({ x, y })));
  });
  let decodes = 0;
  const target = { conditions, layout: targetLayout, decode: () => { decodes++; return { nodes,
    allocation: { groups: f.target.masses.map(row => row.map(massFlow => ({ massFlow }))) } }; },
    evaluate: () => assert.fail('No residual evaluation is authorized') };
  const targetState = new Float64Array(nx * counts.reduce((a, b) => a + b, 0) + 3).fill(.4567), before = targetState.slice();
  const result = prolongStreamtubeDensities(source, sourceState, sourceFlow, target, targetState, f.target);
  const expected = prolong({ ...f, target: { ...f.target, normalAreas: Array.from({ length: nx }, () =>
    f.target.masses.map(row => row.slice())) } });
  result.state.subarray(0, -3).forEach((v, i) => close(v, expected.logDensity.flat(2)[i]));
  assert.deepEqual(result.state.subarray(-3), before.subarray(-3)); assert.deepEqual(targetState, before);
  assert.equal(decodes, 1); assert.equal(result.diagnostics.commonIsentropeInversions, 0);
});
