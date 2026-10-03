// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { prepareContour, signedArea, validateAssembly } from '../src/geometry/airfoil.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';

const hash = path => createHash('sha256').update(fs.readFileSync(path)).digest('hex');
const sourceRows = text => text.split(/\r?\n/).map(line => line.trim().split(/\s+/).map(Number))
  .filter(row => row.length === 3 && row.every(Number.isFinite));
test('RAE preset preserves NASA measured ordinates and the converter lower sign', () => {
  const preset = getBenchmarkAirfoil('rae2822');
  const root = 'third_party/airfoils/rae2822/';
  assert.equal(hash(root + 'geom.txt'), preset.provenance.sourceHash);
  assert.equal(hash(root + 'geom.f90'), preset.provenance.auxiliarySourceHash);
  assert.match(fs.readFileSync(root + 'geom.f90', 'utf8'), /yl\s*=\s*-\s*yl/);
  const rows = sourceRows(fs.readFileSync(root + 'geom.txt', 'utf8'));
  assert.equal(rows.length, 65);
  const upper = rows.map(([x, , y]) => ({ x, y }));
  const lower = rows.map(([x, y]) => ({ x, y: y === 0 ? 0 : -y }));
  assert.deepEqual(preset.elements[0].points, upper.toReversed().concat(lower.slice(1)));
  assert.equal(preset.elements[0].points.length, 129);
  assert.equal(preset.referenceChord, 1); assert.equal(preset.solverUnsupportedReason, undefined);
  assert.deepEqual(preset.elements[0].points[64], { x: 0, y: 0 });
  assert.ok(preset.elements[0].points.some(p => p.x === .952 && p.y === .00125), 'aft lower negative table ordinates must become positive y');
});
test('RAE exact sharp-TE contour passes current contour, curve and assembly preparation', () => {
  const p = getBenchmarkAirfoil('rae2822').elements[0].points;
  assert.deepEqual(prepareContour(p), p); assert.ok(signedArea(p) > 0);
  assert.deepEqual(p[0], { x: 1, y: 0 }); assert.deepEqual(p.at(-1), p[0]);
  const curve = createContourCurve(p);
  assert.deepEqual(curve.evaluate(0).point, p[0]);
  assert.deepEqual(curve.evaluate(curve.length).point, p.at(-1));
  assert.doesNotThrow(() => validateAssembly([p]));
});
test('30P30N preserves exact deployed source coordinates, section order and reference frame', () => {
  const p = getBenchmarkAirfoil('30p30n');
  assert.equal(hash(p.provenance.sourcePath), p.provenance.sourceHash);
  const text = fs.readFileSync(p.provenance.sourcePath, 'utf8'), groups = new Map(); let current;
  for (const line of text.split(/\r?\n/)) {
    if (['# Slat', '# Main Element', '# Flap'].includes(line.trim())) {
      current = line.trim().slice(2); groups.set(current, []); continue;
    }
    const row = line.trim().split(/\s+/).map(Number);
    if (current && row.length === 2 && row.every(Number.isFinite)) groups.get(current).push({ x: row[0], y: row[1] });
  }
  assert.deepEqual(p.elements.map(e => e.name), ['Slat', 'Main Element', 'Flap']);
  assert.deepEqual(p.elements.map(e => e.points.length), [201, 221, 242]);
  for (const e of p.elements) assert.deepEqual(e.points, groups.get(e.name));
  assert.equal(p.referenceChord, 1);
  assert.ok(Math.min(...p.elements[0].points.map(p => p.x)) < 0, 'deployed slat is ahead of stowed LE');
  assert.ok(Math.max(...p.elements[2].points.map(p => p.x)) > 1, 'deployed flap remains aft of stowed TE');
});
test('30P30N source endpoints remain distinct without a preset solver restriction', () => {
  const p = getBenchmarkAirfoil('30p30n'), flap = p.elements[2].points;
  assert.deepEqual(flap[0], { x: 1.128307, y: -.145799 });
  assert.deepEqual(flap.at(-1), { x: 1.130859, y: -.140497 });
  const gap = Math.hypot(flap.at(-1).x - flap[0].x, flap.at(-1).y - flap[0].y);
  assert.equal(gap, .005884208357969709);
  assert.equal(p.solverSupportedModels, undefined);
  assert.equal(p.solverUnsupportedReason, undefined);
  for (const e of p.elements.slice(0, 2)) assert.doesNotThrow(() => prepareContour(e.points));
  assert.throws(() => prepareContour(flap), /closed, sharp trailing edge/);
  assert.doesNotThrow(() => validateAssembly(p.elements.map(e => e.points)), 'supplied element segments do not touch');
});
test('Preset calls own independent geometry and provenance; unknown identifiers return null', () => {
  for (const id of ['rae2822', '30p30n']) {
    const original = getBenchmarkAirfoil(id), changed = getBenchmarkAirfoil(id);
    changed.elements[0].points[0].x = 999; changed.elements.push({ name: 'bad', points: [] });
    changed.provenance.sourceUrls[0] = 'bad'; changed.provenance.configuration = 'bad';
    assert.deepEqual(getBenchmarkAirfoil(id), original);
    assert.ok(fs.existsSync(original.provenance.sourcePath));
    assert.equal(typeof original.provenance.documentationPath, 'string');
    assert.equal(typeof original.provenance.description, 'string');
  }
  for (const id of [undefined, null, '', 'unknown', '__proto__', 'constructor']) assert.equal(getBenchmarkAirfoil(id), null);
});
