// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { parseCoordinates } from '../src/geometry/parse.js';
import { prepareContour, validateAssembly } from '../src/geometry/airfoil.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';

test('MSES RAE preset exactly preserves the public MIT sample, independently of NASA ordinates', () => {
  const p = getBenchmarkAirfoil('rae2822-mses'), bytes = fs.readFileSync(p.provenance.sourcePath);
  assert.equal(createHash('sha256').update(bytes).digest('hex'), p.provenance.sourceHash);
  assert.deepEqual(bytes, fs.readFileSync('docs/rae2822/mses-website-reference/blade.rae'));
  const source = parseCoordinates(bytes.toString());
  assert.deepEqual(source.domain, { xMin: -2, xMax: 3, yMin: -3, yMax: 3.5 });
  assert.equal(p.referenceChord, 1); assert.equal(p.elements[0].points.length, 129);
  assert.deepEqual(p.elements[0].points, source.elements[0].points);
  const nasa = getBenchmarkAirfoil('rae2822').elements[0].points;
  const differences = p.elements[0].points.map((q, i) => { assert.equal(q.x, nasa[i].x); return Math.abs(q.y - nasa[i].y); });
  assert.ok(Math.abs(Math.max(...differences.slice(0, 65)) - .00031) < 1e-16);
  assert.ok(Math.abs(Math.max(...differences.slice(64)) - .00027) < 1e-16);
});

test('MSES RAE sharp contour remains supported and returned geometry/provenance is detached', () => {
  const a = getBenchmarkAirfoil('rae2822-mses'), p = a.elements[0].points;
  assert.equal(a.solverUnsupportedReason, undefined); assert.deepEqual(prepareContour(p), p);
  assert.doesNotThrow(() => createContourCurve(p)); assert.doesNotThrow(() => validateAssembly([p]));
  const b = getBenchmarkAirfoil('rae2822-mses'); b.elements[0].points[0].x = 3; b.provenance.sourceUrls[0] = 'modified';
  assert.deepEqual(getBenchmarkAirfoil('rae2822-mses'), a);
  assert.ok(fs.existsSync(a.provenance.documentationPath));
  assert.match(fs.readFileSync('index.html', 'utf8'), /value="rae2822-mses">RAE 2822<\/option>/);
});
