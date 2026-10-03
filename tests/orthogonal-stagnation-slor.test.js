import test from 'node:test';
import assert from 'node:assert/strict';
import { runOrthogonalStagnationControl } from '../scripts/validation/orthogonal-stagnation-slor.js';

for (const sourceForm of ['metric-stretch', 'poisson']) test(`${sourceForm} corner control refines physical stagnation mass and regular-wall angle`, t => {
  const evidence = [8, 16].map(nt => runOrthogonalStagnationControl({ nt, controlUpdate: 'damped', sourceForm }));
  for (const r of evidence) {
    assert.equal(r.converged, true, r.reason); assert.equal(r.quality.valid, true);
    assert.ok(r.history.every(row => row.invalidCells === 0));
    assert.equal(r.boundariesExactlyFixed, true);
    assert.equal(r.sourceForm, sourceForm);
    assert.ok(r.residual < 1e-10 && r.history.at(-1).controlResidual < 1e-9);
  }
  assert.ok(evidence[1].massError < (sourceForm === 'poisson' ? .3 : .4) * evidence[0].massError);
  assert.ok(evidence[1].maximumWallDerivativeShear < .25 * evidence[0].maximumWallDerivativeShear);
  // At the marked corner, no 90-degree angle is asserted. The entire
  // physical mass field, including that neighborhood, is still measured.
  assert.ok(evidence[1].massError < .0007);
  t.diagnostic(JSON.stringify(evidence.map(r => ({ nt: r.nt, nx: r.nx, sweeps: r.sweeps, massError: r.massError,
    derivativeShear: r.maximumWallDerivativeShear, firstRowShear: r.firstRowShear, maximumAdjacentRatio: r.maximumAdjacentRatio,
    minimumCornerSine: r.quality.minCornerSine, seconds: r.seconds }))));
});
