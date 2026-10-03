import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';

test('browser transition components match the executed native interval and derivative fixtures', async ({ page }) => {
  const fixturePath = 'tests/fixtures/fortran/automatic-transition.json';
  const sourceHashes = numericalSourceHashes([fixturePath, 'tests/browser/automatic-transition-kernel.spec.js']);
  await page.route('**/transition-kernel-harness', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Transition kernel test</title>' }));
  await page.goto('/transition-kernel-harness');
  const result = await page.evaluate(async () => {
    const [{ createIntegralKernel }, { evaluateTransitionInterval }, { selectSurfaceTransition }, fixture] = await Promise.all([
      import('/src/viscous/integral.js'), import('/src/viscous/transition-interval.js'), import('/src/viscous/transition-selection.js'),
      fetch('/tests/fixtures/fortran/automatic-transition.json').then(r => r.json()),
    ]);
    const values = fixture.cases.map(({ input }) => evaluateTransitionInterval(createIntegralKernel({ ...input.parameters, exactJacobian: true }), input, { jacobian: true }));
    const error = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
    const residualError = Math.max(...values.flatMap((v, i) => v.residual.map((r, j) => error(r, fixture.cases[i].expected.residual[j]))));
    const locationError = Math.max(...values.map((v, i) => error(v.transition.s, fixture.cases[i].expected.s)));
    const derivativeError = Math.max(...fixture.derivatives.flatMap(c => {
      const p = values[c.case].partials, col = ['aux', 'theta', 'deltaStar', 'ue', 's'].indexOf(c.key);
      const actual = c.side ? [...p[c.side].map(row => row[col]), p.location[c.side][col]] : [...p.trip, p.location.trip];
      return actual.map((v, i) => error(v, c.expected[i]));
    }));
    const profileErrors = fixture.profiles.map(c => {
      const r = selectSurfaceTransition(createIntegralKernel(c.parameters), c.states, c);
      return { sameIndex: r.index === c.expected.index, sameKind: r.kind === c.expected.kind, locationError: error(r.s, c.expected.s) };
    });
    return { cases: values.length, derivativeColumns: fixture.derivatives.length, residualError, locationError, derivativeError,
      flagsAgreeWithSelectedRoots: values.every((v, i) => v.transition.forced === fixture.cases[i].expected.selectedForced), profileErrors };
  });
  expect(result.cases).toBe(10); expect(result.derivativeColumns).toBe(104);
  expect(result.residualError).toBeLessThan(3e-10); expect(result.locationError).toBeLessThan(3e-10);
  expect(result.derivativeError).toBeLessThan(2e-6); expect(result.flagsAgreeWithSelectedRoots).toBe(true);
  for (const p of result.profileErrors) {
    expect(p.sameIndex).toBe(true); expect(p.sameKind).toBe(true); expect(p.locationError).toBeLessThan(5e-5);
  }
  expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-automatic-transition-browser.json', JSON.stringify({ date: new Date().toISOString(), sourceHashes,
    fixture: { path: fixturePath, sha256: sha256(fixturePath) }, result, physicalAcceptance: false,
    scope: 'Chromium execution of browser-native local transition components and original Fortran comparisons. This is not a coupled Euler/BL solve and does not enable automatic transition in the GUI.' }, null, 2) + '\n');
  console.log(JSON.stringify(result));
});
