// SPDX-License-Identifier: GPL-2.0-or-later
// Public input-routing checks. Panel solves are linear; the viscous adapter
// is limited to its initial state with zero global Newton/wake updates.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { parseCoordinates } from '../src/geometry/parse.js';
import { prepareAirfoilElement } from '../src/geometry/airfoil-element.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { solveInviscid } from '../src/inviscid/linear-vortex.js';
import { createDisplacementOperator } from '../src/inviscid/displacement.js';
import { solveViscousAssembly } from '../src/viscous/result.js';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';

const native = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/finite-base-panel.json', import.meta.url)));
const rows = points => points.map(p => `${p.x} ${p.y}`).join('\n');
const surface = () => structuredClone(native.cases[1].points);
const explicit = e => ({ ...e, points: [...e.points, { ...e.points[0] }],
  trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: e.points.length - 1 } });
const meshControls = { surfaceIntervals: 8, tubes: 3, inletIntervals: 4, outletIntervals: 4 };
const zeroUpdates = { maxIterations: 0, maxWakeIterations: 0 };
const blInput = elements => ({ elements, alpha: 4, reynolds: 1e6, trips: [1, 1],
  wakeLength: 1, wakeCount: 12, initialization: 'march', maxIterations: 1 });
const beforeDir = new URL('../docs/airfoil-formats/before/public-entry/', import.meta.url);
async function archived(file) {
  const url = new URL('../' + file, import.meta.url);
  const code = fs.readFileSync(new URL(file.replaceAll('/', '__') + '.txt', beforeDir), 'utf8')
    .replace(/from '(\.[^']+)'/g, (_, p) => `from '${new URL(p, url).href}'`);
  return import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
}
function sameOuter(a, b) {
  for (const key of ['bodies', 'wakes', 'panels', 'sourcePanels', 'sourceMatrix', 'q0', 'influence', 'field', 'diagnostics'])
    assert.deepEqual(a[key], b[key], key);
  const mass = Float64Array.from({ length: a.total }, (_, i) => .0001 * Math.sin(.3 * i));
  assert.deepEqual(a.evaluate(mass), b.evaluate(mass));
  assert.deepEqual(a.sources(mass), b.sources(mass));
}
function assembly() {
  return parseCoordinates('Open assembly\n-3 4 -2 2\n' + rows(surface()) + '\n999 999\n'
    + rows(transform(surface(), { chord: .3, x: 1.08, y: -.15 })));
}

test('plain/labeled open XFOIL endpoints reach the existing finite-base panel equations in either orientation', async () => {
  const old = await archived('src/inviscid/linear-vortex.js');
  for (const reverse of [false, true]) {
    const points = surface(); if (reverse) points.reverse();
    const document = parseCoordinates((reverse ? 'Finite source\n' : '') + rows(points));
    const snapshot = structuredClone(document), prepared = document.elements.map(prepareAirfoilElement);
    const result = solveInviscid({ elements: document.elements, alpha: 4 });
    const reference = old.solveInviscid({ elements: prepared, alpha: 4 });
    assert.deepEqual(result, reference);
    assert.equal(result.status, 'solved'); assert.equal(result.field.basePanels.length, 1);
    assert.deepEqual(result.elements[0].points, prepared[0].points);
    assert.deepEqual(document, snapshot);
    assert.deepEqual(prepared[0].sourcePoints, points);
  }
});

test('MSES open multielement surfaces retain their frame and match explicit finite panel/displacement inputs', () => {
  const document = assembly(), snapshot = structuredClone(document);
  const input = { elements: document.elements, alpha: 4, wakeCount: 8 };
  const canonical = { ...input, elements: document.elements.map(explicit) };
  assert.deepEqual(solveInviscid(input), solveInviscid(canonical));
  const a = createDisplacementOperator(input), b = createDisplacementOperator(canonical);
  sameOuter(a, b);
  assert.equal(a.bodies.length, 2);
  a.bodies.forEach((body, i) => {
    assert.equal(body.basePanels.length, 1);
    assert.deepEqual(body.points, document.elements[i].points);
    assert.deepEqual(body.solidPoints, canonical.elements[i].points);
  });
  assert.deepEqual(document, snapshot);
});

test('initial streamtube topology accepts raw MSES elements and retains both distinct TE corners', async () => {
  const document = assembly(), snapshot = structuredClone(document);
  const old = await archived('src/geometry/streamtube-topology.js');
  const input = { elements: document.elements, alpha: 4, mach: .2 };
  const result = createInitialStreamtubeTopology(input, meshControls);
  assert.deepEqual(result, old.createInitialStreamtubeTopology({ ...input, elements: document.elements.map(explicit) }, meshControls));
  for (const body of result.bodies) {
    const points = document.elements[body.element].points;
    assert.deepEqual(body.points.slice(0, -1), points);
    assert.deepEqual(body.trailingEdge, { kind: 'finite-base', upperIndex: 0, lowerIndex: points.length - 1 });
    assert.notDeepEqual(body.points[0], body.points[body.trailingEdge.lowerIndex]);
  }
  assert.deepEqual(document, snapshot);
});

test('viscous public normalization routes raw finite XFOIL input identically without a Newton update', async () => {
  const document = parseCoordinates('Finite BL source\n' + rows(transform(surface(), { chord: 1.7, x: .1, y: -.02 })));
  const snapshot = structuredClone(document), old = await archived('src/viscous/result.js');
  const input = { ...blInput(document.elements), referenceChord: 1.7 };
  const result = solveViscousAssembly(input, zeroUpdates);
  const expected = old.solveViscousAssembly({ ...input, elements: document.elements.map(explicit) }, zeroUpdates);
  assert.deepEqual(result, expected);
  assert.equal(result.status, 'unconverged');
  assert.equal(result.diagnostics.iterations, 0);
  assert.equal(result.wakeHistory.length, 0);
  assert.deepEqual(result.elements[0].points.slice(0, -1), document.elements[0].points);
  assert.deepEqual(document, snapshot);
});

test('all four public sharp entry results retain their complete archived data and arithmetic', async () => {
  const input = { elements: [{ name: 'Sharp', points: naca4('2412', 40) }], alpha: 4, wakeCount: 12 };
  const snapshot = structuredClone(input);
  const oldPanel = await archived('src/inviscid/linear-vortex.js');
  for (const boundaryCondition of ['normal-velocity', 'streamfunction']) {
    assert.deepEqual(solveInviscid({ ...input, boundaryCondition }), oldPanel.solveInviscid({ ...input, boundaryCondition }));
  }
  const oldOuter = await archived('src/inviscid/displacement.js');
  sameOuter(createDisplacementOperator(input), oldOuter.createDisplacementOperator(input));
  const oldMesh = await archived('src/geometry/streamtube-topology.js');
  assert.deepEqual(createInitialStreamtubeTopology(input, meshControls), oldMesh.createInitialStreamtubeTopology(input, meshControls));
  const oldBL = await archived('src/viscous/result.js');
  assert.deepEqual(solveViscousAssembly(blInput(input.elements), zeroUpdates), oldBL.solveViscousAssembly(blInput(input.elements), zeroUpdates));
  assert.deepEqual(input, snapshot);
});

test('unsupported repeated surface-corner markers are rejected explicitly at each public entry', () => {
  const points = surface(); points.splice(12, 0, { ...points[12] });
  const input = { elements: [{ points }] }, check = /Repeated successive.*surface corners/;
  assert.throws(() => solveInviscid(input), check);
  assert.throws(() => createDisplacementOperator(input), check);
  assert.throws(() => createInitialStreamtubeTopology(input, meshControls), check);
  assert.throws(() => solveViscousAssembly(blInput(input.elements), zeroUpdates), check);
});
