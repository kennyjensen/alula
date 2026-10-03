import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { axset } from '../src/viscous/xfoil/xblsys.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';

const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/transition-root-separated-interval.json', import.meta.url)));
const pack = r => [...r.residual, r.transition.s];

test('separated-interval transition failures resolve the unchanged natural amplification equation', () => {
  assert.equal(fixture.cases.length, 10);
  for (const [i, { parameters, input }] of fixture.cases.entries()) {
    const before = structuredClone(input), kernel = createIntegralKernel(parameters);
    const root = kernel.transitionCheck(input);
    assert.equal(root.transition, true);
    assert.equal(root.forced, false);
    // Reconstruct f from amplification. Subtracting two absolute locations
    // loses bits for transition near the upstream end of this long interval.
    const fraction = (parameters.ncrit - input.upstream.aux) / (root.amplification - input.upstream.aux);
    assert.ok(Math.abs(fraction - fixture.independentFractionRoots[i].fraction) < 3e-14);
    assert.ok(Math.abs(root.s - (input.upstream.s * (1 - fraction) + input.downstream.s * fraction)) < 1e-15);
    const at = { s: root.s, aux: parameters.ncrit };
    for (const key of ['theta', 'deltaStar', 'ue'])
      at[key] = input.upstream[key] * (1 - fraction) + input.downstream[key] * fraction;
    const up = kernel.station(input.upstream), transition = kernel.station(at);
    const ax = axset(up.hk, input.upstream.theta, up.reTheta, input.upstream.aux,
      transition.rawHk, at.theta, transition.reTheta, parameters.ncrit, parameters.ncrit, 0).ax;
    const residual = root.amplification - input.upstream.aux - ax * (input.downstream.s - input.upstream.s);
    assert.ok(Math.abs(residual) < 1e-11, `case ${i}: amplification residual ${residual}`);
    const interval = kernel.interval({ ...input, regime: 'transition' });
    assert.equal(interval.transition.s, root.s);
    assert.equal(interval.transition.amplification, root.amplification);
    assert.ok([...interval.residual, ...interval.upstream.flat(), ...interval.downstream.flat()].every(Number.isFinite));
    // Repeated checks after interpolated station and mixed-interval work must
    // replay the same root: the bracket's scratch context cannot leak out.
    assert.deepEqual(kernel.transitionCheck(input), root);
    assert.deepEqual(input, before);
  }
});

test('a bracket-recovered mixed interval retains resolved-root Jacobian agreement', t => {
  const { parameters, input } = fixture.cases[0], kernel = createIntegralKernel(parameters);
  const result = evaluateTransitionInterval(kernel, input, { jacobian: true });
  let maximum = 0;
  for (const side of ['upstream', 'downstream']) for (const [column, key] of ['aux', 'theta', 'deltaStar', 'ue', 's'].entries()) {
    const h = 1e-4 * (key === 's' ? input.downstream.s - input.upstream.s : Math.max(Math.abs(input[side][key]), key === 'aux' ? .01 : 1e-7));
    const samples = [-2, -1, 1, 2].map(offset => pack(evaluateTransitionInterval(kernel,
      { ...input, [side]: { ...input[side], [key]: input[side][key] + offset * h } })));
    const reference = samples[0].map((_, row) => (samples[0][row] - 8 * samples[1][row] + 8 * samples[2][row] - samples[3][row]) / (12 * h));
    const actual = [...result.partials[side].map(row => row[column]), result.partials.location[side][column]];
    const norm = Math.max(1, ...reference.map(Math.abs));
    const error = Math.max(...actual.map((value, row) => Math.abs(value - reference[row]))) / norm;
    maximum = Math.max(maximum, error);
    assert.ok(error < 2e-5, `${side}.${key}: relative column error ${error}`);
  }
  t.diagnostic(JSON.stringify({ maximumRelativeColumnError: maximum }));
});
