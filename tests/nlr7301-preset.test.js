// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { execFileSync } from 'node:child_process';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { prepareContour, signedArea, validateAssembly } from '../src/geometry/airfoil.js';

const root = 'third_party/airfoils/nlr7301/';
const hash = path => createHash('sha256').update(fs.readFileSync(path)).digest('hex');
const geometry = JSON.parse(fs.readFileSync(root + 'geometry.json'));

test('NLR wall extraction exactly reproduces original SAAB archive boundaries without fitting', () => {
  assert.equal(hash(root + 'UFR3-01_nlr7301grid.zip'), '28aa28b47d981681df3b6300a41a6a9e1c029514ca28cb91bda699b32b7bfc62');
  const stdout = execFileSync('python3', [root + 'extract.py', '--check'], { encoding: 'utf8' });
  assert.match(stdout, /"unpairedVertexDegrees"/);
  assert.deepEqual(geometry.edgeMultiplicity, { 1: 1920, 2: 1440 });
  assert.deepEqual(geometry.unpairedVertexDegrees, { 2: 1920 });
  assert.deepEqual(geometry.loopPointCounts, { main: 456, flap: 248, farfield: 1216 });
});

test('NLR static preset and downloadable coordinates retain all exact wall and base points', () => {
  const p = getBenchmarkAirfoil('nlr7301');
  assert.equal(p.referenceChord, 1); assert.equal(hash(p.provenance.sourcePath), p.provenance.sourceHash);
  assert.equal(p.provenance.sourceArchiveHash, geometry.sourceArchiveHash);
  assert.equal(p.provenance.sourceMemberHash, geometry.sourceMemberHash);
  const groups = [[]];
  for (const line of fs.readFileSync(p.provenance.sourcePath, 'utf8').trim().split(/\r?\n/).slice(2)) {
    const [x, y] = line.split(/\s+/).map(Number);
    if (x === 999 && y === 999) groups.push([]); else groups.at(-1).push({ x, y });
  }
  assert.deepEqual(p.elements.map(e => e.points.length), [457, 249]);
  p.elements.forEach((e, k) => {
    assert.deepEqual(e.points, geometry.elements[k].points.map(([x, y]) => ({ x, y })));
    assert.deepEqual(e.points, groups[k]); assert.deepEqual(e.points[0], e.points.at(-1));
    assert.ok(signedArea(e.points) > 0);
    const base = geometry.elements[k].sourceEdges.find(s => s.side === 'imin');
    assert.equal(base.count, 32, 'The original finite base remains fully represented.');
    assert.ok(geometry.elements[k].trailingEdgeBaseLength > 0);
  });
  assert.doesNotThrow(() => validateAssembly(p.elements.map(e => e.points)));
});

test('Closed finite-base NLR geometry has no preset-specific solver lock', () => {
  const p = getBenchmarkAirfoil('nlr7301');
  for (const key of ['solverSupportedModels','solverSupportedQuadModes','solverUnsupportedReason']) assert.equal(p[key], undefined);
  for (const e of p.elements) assert.throws(() => prepareContour(e.points), /sharp trailing edge/);
  assert.match(fs.readFileSync('index.html', 'utf8'), /value="nlr7301">NLR 7301 \+ flap<\/option>/);
  assert.ok(fs.existsSync(p.provenance.documentationPath));
});

test('NLR callers own their data and earlier geometry/provenance remain exactly unchanged', async () => {
  const a = getBenchmarkAirfoil('nlr7301'), b = getBenchmarkAirfoil('nlr7301');
  b.elements[0].points[0].x = 999; b.provenance.sourceUrls[0] = 'changed';
  assert.deepEqual(getBenchmarkAirfoil('nlr7301'), a);
  const archived = fs.readFileSync('docs/benchmark-presets/nlr-before-src-geometry-benchmark-airfoils.js.txt', 'utf8');
  const before = await import('data:text/javascript;base64,' + Buffer.from(archived).toString('base64'));
  // Only retired availability metadata differs; original geometry and
  // provenance remain exact for both earlier presets.
  const withoutAvailability = ({ solverSupportedModels, solverUnsupportedReason, ...data }) => data;
  for (const id of ['rae2822','30p30n']) assert.deepEqual(withoutAvailability(getBenchmarkAirfoil(id)), withoutAvailability(before.getBenchmarkAirfoil(id)));
});
