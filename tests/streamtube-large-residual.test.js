import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import vm from 'node:vm';

// Read the actual internal reducers, including restart and certification.
// No global flow solve is needed to exercise engine argument-count limits.
for (const file of ['streamtube-coupled-ises.js', 'streamtube-coupled.js',
  'streamtube-coupled-flow-restart.js', 'streamtube-coupled-hybrid-certification.js',
  'streamtube-independent-wake-restart.js']) {
  test(`${file}: million-entry residual preserves extrema and nonfinite semantics`, () => {
    const testOnly = file !== 'streamtube-coupled-ises.js' && file !== 'streamtube-coupled.js';
    const source = fs.readFileSync(new URL(`../src/euler/${testOnly ? 'tests/' : ''}${file}`, import.meta.url), 'utf8');
    const declaration = source.match(/const maximum = [^;]+;/)[0];
    const maximum = vm.runInNewContext(`${declaration}\nmaximum`);
    const seed = file === 'streamtube-independent-wake-restart.js' ? 0 : -Infinity;
    for (const a of [[], [0, -0], [-3, 2, -4], [1e-300, -1e300], [Infinity], [NaN, 1], [1, NaN]])
      for (const values of [a, Float64Array.from(a)])
        assert.equal(maximum(values), Math.max(seed, ...a.map(Math.abs)));
    const large = new Float64Array(1_000_000); large[600_000] = -7; large[999_999] = 6;
    assert.equal(maximum(large), 7);
    large[999_999] = NaN; assert.ok(Number.isNaN(maximum(large)));
  });
}

test('real ISES initial report accepts a million-entry residual without a Newton solve', async () => {
  const location = new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url);
  let source = fs.readFileSync(location, 'utf8')
    .replace('import { createCoupledStreamtubeBody,', 'import { createCoupledStreamtubeBody as originalCreate,')
    .replace(/from '([^']+)'/g, (_, relative) => `from '${new URL(relative, location).href}'`);
  source += `\nfunction createCoupledStreamtubeBody(...args) {
    const s = originalCreate(...args), original = s.admissibleValue;
    s.admissibleValue = (...a) => {
      const value = original(...a);
      if (!value) return value;
      const residual = new Float64Array(1_000_000); residual.set(value.residual);
      return { ...value, residual };
    };
    return s;
  }\n`;
  const { solveCoupledStreamtubeIses } = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
  const cp = JSON.parse(fs.readFileSync(new URL('../docs/solver-reliability/rae-wake-contact/cold-checkpoint.json', import.meta.url)));
  const reports = [], before = structuredClone(cp);
  const result = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0,
    onIteration: h => reports.push(h) });
  assert.equal(reports.length, 1);
  assert.equal(reports[0].residual, Math.max(...Object.values(cp.families)));
  assert.deepEqual(result.families, cp.families);
  assert.deepEqual(cp, before);
});
