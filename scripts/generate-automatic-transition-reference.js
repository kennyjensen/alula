// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import { mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, relative } from 'node:path';
import { spawnSync } from 'node:child_process';
import assert from 'node:assert/strict';
import { buildBLReference } from './reference/build-bl.js';
import { sha256 } from './validation/provenance.js';

const output = process.argv[2] ?? 'tests/fixtures/fortran/automatic-transition.json';
const directory = await mkdtemp(join(tmpdir(), 'mses-auto-transition-'));
console.log(JSON.stringify({ stage: 'build', directory }));
const native = await buildBLReference(directory, 'automatic-transition');
const run = cases => {
  const lines = [cases.length];
  for (const c of cases) {
    lines.push([c.parameters.reynolds, c.parameters.mach, c.parameters.ncrit, c.tripS ?? 1e30].join(' '));
    for (const s of [c.upstream, c.downstream]) lines.push([s.s, s.aux, s.theta, s.deltaStar, s.ue].join(' '));
  }
  const r = spawnSync(native.executable, [], { input: lines.join('\n') + '\n', encoding: 'utf8', timeout: 10000 });
  if (r.error || r.status !== 0) throw new Error(`Native automatic-transition oracle failed: ${r.error?.message ?? r.stderr}`);
  const results = cases.map(() => ({}));
  for (const line of r.stdout.trim().split('\n')) {
    const [tag, id, ...values] = line.trim().split(/\s+/), c = results[Number(id) - 1];
    assert.ok(c, line);
    if (tag === 'CHECK') Object.assign(c, { transition: values[0] === 'T', forced: values[1] === 'T', s: Number(values[2]), amplification: Number(values[3]) });
    else if (tag === 'R') c.residual = values.map(Number);
    else if (tag === 'SHEAR') c.shear = Number(values[0]);
    else throw new Error(`Unexpected Fortran output: ${line}`);
  }
  results.forEach(c => {
    assert.ok(Number.isFinite(c.amplification) && Number.isFinite(c.s) && c.shear > 0);
    if (c.transition) assert.ok(c.residual?.length === 3 && c.residual.every(Number.isFinite));
  });
  return results;
};
const natural = { upstream: { s: .2, aux: 8.8, theta: .0007, deltaStar: .0024, ue: 1.2 },
  downstream: { s: .24, aux: .035, theta: .0009, deltaStar: .002, ue: 1.15 } };
const cases = [];
for (const mach of [0, .2, .3]) for (const [name, tripS] of [['natural', undefined], ['earlier-trip', .205], ['later-trip', .235]])
  cases.push({ name: `${name}-M${mach}`, parameters: { reynolds: 1e6, mach, ncrit: 9 }, ...structuredClone(natural), ...(tripS ? { tripS } : {}) });
cases.push({ name: 'near-node-forced', parameters: { reynolds: 1e6, mach: .2, ncrit: 9 },
  ...structuredClone(natural), upstream: { ...natural.upstream, aux: 0 }, tripS: .24 - 1e-10 });
const base = run(cases), samples = [], sampleMetadata = [];
// Preserve the native TRFORC flag. Independently run each case without its
// trip: an earlier trip that pins XT is physically forced even when the
// original strict floating-point comparison labels it free.
const untripped = run(cases.map(c => ({ ...c, tripS: undefined })));
for (let i = 0; i < cases.length; i++) {
  const c = cases[i], r = base[i];
  r.untripped = untripped[i];
  r.selectedForced = c.tripS !== undefined && (!untripped[i].transition || c.tripS < untripped[i].s);
  if (r.selectedForced) assert.ok(Math.abs(r.s - c.tripS) < 1e-15);
}
for (const [i, c] of cases.entries()) {
  assert.equal(base[i].transition, true, c.name);
  const fields = ['upstream', 'downstream'].flatMap(side => ['aux', 'theta', 'deltaStar', 'ue', 's'].map(key => ({ side, key })));
  fields.push({ key: 'tripS' });
  for (const field of fields) {
    if (field.key === 'tripS' && !base[i].selectedForced) continue;
    const v = field.side ? c[field.side][field.key] : c.tripS;
    const h = 2e-5 * Math.max(Math.abs(v), field.key === 'aux' ? .01 : 1e-7);
    // Independent fourth-order native stencils. Near the retained endpoint
    // use five one-sided points; production uses three second-order points.
    const sign = c.name === 'near-node-forced' && field.key === 's' && field.side === 'downstream' ? 1
      : c.name === 'near-node-forced' && field.key === 'tripS' ? -1 : 0;
    const offsets = sign ? [0, 1, 2, 3, 4].map(k => sign * k) : [-2, -1, 1, 2];
    const indices = offsets.map(k => {
      const q = structuredClone(c);
      if (field.side) q[field.side][field.key] = v + k * h; else q.tripS = v + k * h;
      samples.push(q); return samples.length - 1;
    });
    sampleMetadata.push({ case: i, ...field, step: h, sign, offsets, indices });
  }
}
const sampled = run(samples);
for (const m of sampleMetadata) {
  const values = m.indices.map(i => sampled[i]);
  values.forEach(v => {
    assert.equal(v.transition, true, JSON.stringify(m));
  });
  const coefficients = m.sign ? [-25, 48, -36, 16, -3] : [1, -8, 8, -1];
  const packed = values.map(v => [...v.residual, v.s]);
  m.expected = Array.from({ length: 4 }, (_, k) => coefficients.reduce((s, a, j) => s + a * packed[j][k], 0) / (12 * m.step * (m.sign || 1)));
  m.samples = values;
  delete m.indices;
}

// The oracle propagates its OWN N values, without importing JS transition
// decisions. These are local prescribed-state tests, not airfoil predictions.
const states = [.01, .05, .1, .15, .2, .25, .3, .4, .6, .8, 1].map(s => ({ s, aux: .03,
  theta: .0007 * Math.sqrt(s / .2), deltaStar: 3.4 * .0007 * Math.sqrt(s / .2), ue: 1.25 - .25 * s }));
const profiles = [];
for (const ncrit of [3, 6, 9, 12]) for (const tripS of [undefined, .18]) {
  const parameters = { reynolds: 1e6, mach: .2, ncrit }, amplification = [0], checks = [];
  let expected = { index: null, kind: 'laminar', s: null };
  for (let j = 1; j < states.length; j++) {
    const upstream = { ...states[j - 1], aux: amplification[j - 1] };
    const r = run([{ parameters, upstream, downstream: states[j], tripS }])[0];
    checks.push(r); amplification.push(r.amplification);
    if (r.transition) { expected = { index: j, kind: r.forced ? 'forced' : 'natural', s: r.s }; break; }
  }
  profiles.push({ parameters, states, ...(tripS ? { tripS } : {}), expected: { ...expected, amplification }, checks });
}
const sources = [...native.files, 'third_party/Xfoil/src/XBL.INC', 'third_party/Xfoil/src/BLPAR.INC',
  'scripts/reference/build-bl.js', 'scripts/generate-automatic-transition-reference.js'];
const report = { provenance: { date: new Date().toISOString(), compiler: native.compiler, flags: native.flags,
  scope: 'Executed unmodified original Fortran TRCHEK/BLSYS/BLVAR. Local transition roots, residuals, fourth-order finite differences and sequential prescribed-state interval selection; not a coupled airfoil solve.',
  sha256: Object.fromEntries(sources.map(p => [relative(process.cwd(), p), sha256(p)])) },
  cases: cases.map((input, i) => ({ input, expected: base[i] })), derivatives: sampleMetadata, profiles };
fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ output, cases: cases.length, nativePerturbations: samples.length, derivativeColumns: sampleMetadata.length, profiles: profiles.length, compiler: native.compiler }));
