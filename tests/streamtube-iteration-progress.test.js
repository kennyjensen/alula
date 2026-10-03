import test from 'node:test';
import assert from 'node:assert/strict';
import { createIterationProgress, residualMerit, requireCoupledResidualDecrease } from '../src/euler/streamtube-iteration-progress.js';
import { quadResidualProgressRows } from '../src/ui/quad-residual-progress.js';
const sample = (vector, overrides = {}) => ({ vector, before: 1, after: 1, residual: 1,
  tolerance: 1e-10, phase: 'a', mcrit: .99, ...overrides });
test('dynamic-pressure merit cannot hide a worsening pressure balance at low Mach', () => {
  const controls = { step: .5, tolerance: 1e-10 };
  for (const mach of [.01, .2, .74]) {
    const pressureScale = 1 / (1.4 * mach ** 2), weights = [pressureScale, 1];
    // Row zero is a pressure equation divided by pInf; row one is geometry.
    // Doubling the physical pressure error outweighs the small geometry gain.
    assert.throws(() => requireCoupledResidualDecrease([1 / pressureScale, 1],
      [2 / pressureScale, .9], { ...controls, weights }), { code: 'COUPLED_RESIDUAL_DECREASE' });
    const accepted = requireCoupledResidualDecrease([1 / pressureScale, 1],
      [.5 / pressureScale, .5], { ...controls, weights });
    assert(Math.abs(accepted.beforeSquaredNorm - 2) < 1e-14);
    assert(Math.abs(accepted.afterSquaredNorm - .5) < 1e-14);
  }
});
test('merit weights preserve the governing convergence tolerance and reject invalid scales', () => {
  const controls = { step: 1, tolerance: 1e-10 };
  const root = requireCoupledResidualDecrease([1e-13, 1e-13], [5e-11, 5e-11],
    { ...controls, weights: [1e4, 1] });
  assert.equal(root.converged, true);
  assert.equal(root.maximumResidual, 5e-11);
  for (const weights of [[1], [0, 1], [-1, 1], [NaN, 1], [Infinity, 1]])
    assert.throws(() => requireCoupledResidualDecrease([1, 1], [.5, .5], { ...controls, weights }));
});
test('stagnation detects actual negligible motion, never declares convergence', () => {
  const p = createIterationProgress(); p.seed([1, 2], 'a', .99);
  assert.equal(p.observe(sample([1, 2])).cause, null);
  assert.equal(p.observe(sample([1, 2])).cause, null);
  assert.equal(p.observe(sample([1, 2])).cause, 'negligible-state-change');
  assert.equal(p.observe(sample([1, 2], { residual: 1e-12 })).cause, null);
  const improving = createIterationProgress(); improving.seed([1, 2], 'a', .99);
  for (let i = 0; i < 4; i++) assert.equal(improving.observe(sample([1, 2], { after: .1 })).cause, null,
    'small physical motion alone must not stop a rapidly improving residual');
});
test('two-state cycle requires physical recurrence, not alternating residual numbers', () => {
  const p = createIterationProgress(); p.seed([0], 'a', .99);
  let r;
  for (let i = 1; i <= 4; i++) r = p.observe(sample([i % 2], { residual: i % 2 ? 2 : 1,
    phase: i % 2 ? 'b' : 'a', phaseChanged: true, mcrit: i % 2 ? .8 : .99 }));
  assert.equal(r.cause, 'two-state-cycle'); assert.equal(r.transitionChanged, true); assert.equal(r.mcritChanged, true);
  const q = createIterationProgress(); q.seed([0], 'a', .99);
  for (let i = 1; i < 6; i++) assert.notEqual(q.observe(sample([i], { residual: i % 2 ? 2 : 1 })).cause, 'two-state-cycle');
  const changing = createIterationProgress(); changing.seed([0], '0', .99);
  for (let i = 1; i <= 4; i++) assert.equal(changing.observe(sample([i % 2], {
    phase: String(i), phaseChanged: true })).cause, null, 'auxiliary coordinates with different phase meanings are not the same physical state');
});
test('plateau uses paired fixed-equation merits and ignores transition relabeling', () => {
  for (const changed of [false, true]) {
    const p = createIterationProgress(); p.seed([0], 'a', .99); let r;
    for (let i = 1; i <= 6; i++) r = p.observe(sample([i], { after: .995, phaseChanged: changed }));
    assert.equal(r.cause, changed ? null : 'residual-plateau');
  }
});
test('extra budget requires sustained rapid improvement, with bounded stable norm calculation', () => {
  const p = createIterationProgress(); p.seed([0], 'a', .99);
  for (let i = 1; i <= 3; i++) p.observe(sample([i], { after: .2 }));
  assert.equal(p.canExtend(), true);
  p.observe(sample([4], { after: 1 })); assert.equal(p.canExtend(), false);
  assert.equal(residualMerit(new Float64Array(200000).fill(1e200)), 1e200);
});
test('UI distinguishes active and prescribed residuals, including unsampled and inadmissible checks', () => {
  const h = { euler: .2, boundaryLayer: .3, edgeMatching: .1, residualContext: { mcrit: .8 } };
  assert(quadResidualProgressRows(h).some(r => r.value === 'Not sampled this iteration'));
  const rows = quadResidualProgressRows({ ...h, prescribedResidual: { euler: .4, boundaryLayer: .3, edgeMatching: .1 } });
  assert(rows.some(r => r.label === 'Prescribed Euler equation error'.replace('Euler', 'euler') && r.value === '4.00e-1'));
  assert(!rows.some(r => r.value === 'Not sampled this iteration'));
  assert(quadResidualProgressRows({ ...h, prescribedResidual: { unavailable: 'invalid' } }).some(r => r.value.includes('inadmissible')));
});
