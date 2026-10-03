import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';
const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/coupled-inverse-bl.json', import.meta.url)));
const sha = value => createHash('sha256').update(value).digest('hex');

test('every prepared inverse BL block and selected inverse derivatives agree with executed Fortran', t => {
  for (const [file, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(sha(fs.readFileSync(new URL('../' + file, import.meta.url))), hash, file);
  let native = fs.readFileSync(new URL('../third_party/Xfoil/src/xblsys.f', import.meta.url), 'utf8');
  for (const p of fixture.provenance.variant.patches) { assert.equal(native.split(p.from).length, 2); native = native.replace(p.from, p.to); }
  assert.equal(sha(native), fixture.provenance.variant.sourceSHA256);
  const differentiated = new Set(fixture.derivatives.map(d => d.case)); let maxDifference = 0, maxResidual = 0, maxDerivative = 0;
  const evaluated = fixture.cases.map((c, i) => {
    const k = createIntegralKernel(c.parameters), q = c.input;
    const v = q.regime === 'te' ? k.trailingEdge(q.upper, q.lower, q.downstream)
      : q.regime === 'transition' ? evaluateTransitionInterval(k, q, { jacobian: differentiated.has(i) }) : k.interval(q);
    v.residual.forEach((a, row) => {
      maxDifference = Math.max(maxDifference, Math.abs(a - c.expected[row]));
      maxResidual = Math.max(maxResidual, Math.abs(a * c.rowScale[row]));
    });
    return v.partials?.downstream ?? v.downstream;
  });
  for (const d of fixture.derivatives) d.expected.forEach((v, row) => {
    const a = evaluated[d.case][row][d.col]; maxDerivative = Math.max(maxDerivative, Math.abs(a - v) / Math.max(1, Math.abs(a), Math.abs(v)));
  });
  assert.ok(maxDifference < 3e-10); assert.ok(maxResidual < 3e-11); assert.ok(maxDerivative < 2e-6);
  t.diagnostic(JSON.stringify({ blocks: fixture.cases.length, columns: fixture.derivatives.length, maxDifference, maxResidual, maxDerivative }));
});
