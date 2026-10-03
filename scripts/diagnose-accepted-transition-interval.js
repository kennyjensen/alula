// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';

const file = 'docs/accepted-coupled-trust-region.json', saved = JSON.parse(fs.readFileSync(file)), started = performance.now();
const report = { date: new Date().toISOString(), input: { path: file, sha256: sha256(file) },
  sourceHashes: numericalSourceHashes(['scripts/diagnose-accepted-transition-interval.js']),
  scope: 'Four retained local transition blocks and the exact centered perturbations used by their Jacobian. No coupled iteration.', blocks: [] };
const seed = saved.restart, system = createCoupledStreamtubeBody(seed.input, {
  ...seed.options, initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
const value = system.evaluate(system.initial), { states, geometry } = value.layers;
for (const [key, number] of Object.entries(saved.result.families))
  assert.ok(Math.abs(value.families[key] - number) < 1e-8 * Math.max(1, Math.abs(number)), `Restore ${key}`);
system.bl.surfaces.forEach((surface, index) => {
  const upstream = surface.ids[surface.transition - 1], downstream = surface.ids[surface.transition];
  const input = { upstream: states[upstream], downstream: states[downstream], regime: 'transition', tripS: geometry.surfaceData[index].tripS };
  const base = system.bl.kernel.interval(input), block = { body: surface.body, side: surface.side,
    upstream, downstream, input, base: { residual: base.residual, transition: base.transition }, perturbations: [] };
  for (const side of ['upstream', 'downstream']) ['aux', 'theta', 'deltaStar', 'ue', 's'].forEach((key, k) => {
    const v = input[side][key], h = Math.cbrt(Number.EPSILON) * Math.max(Math.abs(v), k === 0 ? .01 : k < 3 ? 1e-7 : 1e-6);
    for (const sign of [-1, 1]) {
      const point = { ...input, [side]: { ...input[side], [key]: v + sign * h } };
      let outcome;
      try { const r = system.bl.kernel.interval(point); outcome = { residual: r.residual, transition: r.transition }; }
      catch (error) { outcome = { error: error.message }; }
      block.perturbations.push({ side, key, sign, h, ...outcome });
    }
  });
  report.blocks.push(block);
});
assert.deepEqual(changedSources(report.sourceHashes), []); assert.equal(sha256(file), report.input.sha256);
report.seconds = (performance.now() - started) / 1000;
fs.writeFileSync('docs/accepted-transition-interval-diagnosis.json', JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ seconds: report.seconds, blocks: report.blocks.map(b => ({ body: b.body, side: b.side,
  input: b.input, transition: b.base.transition, failed: b.perturbations.filter(p => p.error || !p.transition.forced) })) }));
