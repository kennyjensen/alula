// SPDX-License-Identifier: GPL-2.0-or-later
// node scripts/validation/benchmark-geometry-kernels.js BEFORE_ROOT UNUSED_OUT
import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { createContourCurve } from '../../src/geometry/contour-curve.js';
import { streamtubeCellGeometry } from '../../src/euler/streamtube-cell.js';
import { naca4 } from '../../src/geometry/airfoil.js';

const [root, output] = process.argv.slice(2);
if (!root || !output || fs.existsSync(output)) throw new Error('Supply baseline source root and unused output filename.');
const oldCurve = await import(pathToFileURL(path.resolve(root, 'src/geometry/contour-curve.js')));
const oldCell = await import(pathToFileURL(path.resolve(root, 'src/euler/streamtube-cell.js')));
const points = naca4('2412', 160), beforeCurve = oldCurve.createContourCurve(points), afterCurve = createContourCurve(points);
const positions = Array.from({ length: 1000 }, (_, i) => i / 999 * beforeCurve.length);
const stencils = Array.from({ length: 128 }, (_, i) => {
  const width = 10 ** (-6 * i / 127), offset = .1 * Math.sin(i) * width;
  return { lower: [{ x: 0, y: 0 }, { x: .45, y: .025 }, { x: 1.1, y: .08 }],
    upper: [{ x: offset, y: width }, { x: .45 + offset, y: .025 + width }, { x: 1.1 + offset, y: .08 + width }] };
});
const kernels = [
  { name: 'contour', inputs: positions, repeats: 1000, before: s => beforeCurve.evaluate(s), after: s => afterCurve.evaluate(s), scalar: v => v.point.x },
  { name: 'surface-branch', inputs: positions.map(s => s / beforeCurve.length), repeats: 1000,
    before: f => beforeCurve.branch('upper', f, .48 * beforeCurve.length),
    after: f => afterCurve.branch('upper', f, .48 * afterCurve.length), scalar: v => v.point.x },
  // With a contour-only change these also provide unchanged timing controls.
  ...['convex', 'positive-simple'].map(geometryDomain => ({ name: `cell-${geometryDomain}`, inputs: stencils, repeats: 1000,
    before: p => oldCell.streamtubeCellGeometry(p.lower, p.upper, { geometryDomain }),
    after: p => streamtubeCellGeometry(p.lower, p.upper, { geometryDomain }), scalar: v => v.area })),
];
const samples = [], hashes = {};
for (const k of kernels) {
  const expected = k.inputs.map(k.before), actual = k.inputs.map(k.after);
  assert.deepEqual(actual, expected);
  hashes[k.name] = createHash('sha256').update(JSON.stringify(actual)).digest('hex');
  for (const label of ['before', 'after', 'after', 'before']) {
    for (let i = 0; i < 10; i++) for (const input of k.inputs) k[label](input);
    let checksum = 0;
    const start = performance.now();
    for (let i = 0; i < k.repeats; i++) for (const input of k.inputs) checksum += k.scalar(k[label](input));
    samples.push({ kernel: k.name, label, calls: k.inputs.length * k.repeats, milliseconds: performance.now() - start, checksum });
  }
}
fs.writeFileSync(output, JSON.stringify({ runtime: process.version, scope: 'Isolated warm ABBA kernels; not whole-solve speedups.', hashes, samples }, null, 2), { flag: 'wx' });
