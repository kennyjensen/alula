import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/coupled-auxiliary-roots.json', import.meta.url)));
const sha = value => createHash('sha256').update(value).digest('hex');

test('conditional auxiliary roots and remaining equations match the independently executed native blocks', t => {
  for (const [file, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(sha(fs.readFileSync(new URL('../' + file, import.meta.url))), hash, file);
  let source = fs.readFileSync(new URL('../third_party/Xfoil/src/xblsys.f', import.meta.url), 'utf8');
  for (const p of fixture.provenance.variant.patches) { assert.equal(source.split(p.from).length, 2); source = source.replace(p.from, p.to); }
  assert.equal(sha(source), fixture.provenance.variant.sourceSHA256);
  let maxDifference = 0, maxAuxiliary = 0;
  for (const c of fixture.cases) {
    const k = createIntegralKernel(c.parameters), input = c.input;
    const value = input.regime === 'te' ? k.trailingEdge(input.upper, input.lower, input.downstream)
      : input.regime === 'transition' ? evaluateTransitionInterval(k, input) : k.interval(input);
    value.residual.forEach((v, i) => { const error = Math.abs(v - c.expected.residual[i]); assert.ok(error < 3e-10); maxDifference = Math.max(maxDifference, error); });
    maxAuxiliary = Math.max(maxAuxiliary, Math.abs(value.residual[0] * c.rowScale[0]));
  }
  assert.ok(maxAuxiliary < 2e-12);
  for (const c of fixture.brackets) {
    assert.ok(c.lower.residual[0] > 0 && c.upper.residual[0] < 0);
    assert.ok(c.lower.downstreamShearDerivative < 0 && c.upper.downstreamShearDerivative < 0);
  }
  t.diagnostic(JSON.stringify({ intervals: fixture.cases.length, brackets: fixture.brackets.length, maxDifference, maxAuxiliary }));
});
