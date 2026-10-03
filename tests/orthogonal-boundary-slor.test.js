import test from 'node:test';
import assert from 'node:assert/strict';
import { smoothOrthogonalBoundaryGrid } from '../src/geometry/tests/orthogonal-boundary-slor.js';
import { runOrthogonalStagnationControl } from '../scripts/validation/orthogonal-stagnation-slor.js';

test('damped-control convergence requires the undamped PDE and the boundary control state', () => {
  const nodes = Array.from({ length: 6 }, (_, i) => Array.from({ length: 5 }, (_, j) => ({ x: i / 5, y: j / 4 })));
  const make = background => ({ nodes, massFlows: [1, 1, 1, 1], discretization: 'giles-1985', orthogonalBoundaryControl: { background } });
  const uniform = smoothOrthogonalBoundaryGrid(make(nodes.map(row => row.map(() => 0))), { maxSweeps: 0 });
  assert.equal(uniform.converged, true); assert.equal(uniform.sweeps, 0); assert.deepEqual(uniform.nodes, nodes);
  // The boundary corrections cancel this constant background in the
  // undamped equations. The initial, nonzero control state is still wrong.
  const unfinished = smoothOrthogonalBoundaryGrid(make(nodes.map(row => row.map(() => 2))), { maxSweeps: 0 });
  assert.ok(unfinished.history[0].residual < 1e-12);
  assert.equal(unfinished.history[0].controlResidual, 2);
  assert.equal(unfinished.converged, false);
  assert.equal(unfinished.reason, 'sweep limit');
});

test('damped corner controls converge to the same stagnation-grid root as direct controls', t => {
  const direct = runOrthogonalStagnationControl({ nt: 8 }), damped = runOrthogonalStagnationControl({ nt: 8, controlUpdate: 'damped' });
  assert.equal(direct.converged, true, direct.reason); assert.equal(damped.converged, true, damped.reason);
  let difference = 0;
  direct.nodes.forEach((row, i) => row.forEach((p, j) => { difference = Math.max(difference,
    Math.hypot(p.x - damped.nodes[i][j].x, p.y - damped.nodes[i][j].y)); }));
  assert.ok(difference < 1e-9, String(difference));
  assert.ok(damped.history.at(-1).controlResidual < 1e-9);
  assert.ok(damped.history.every(row => row.invalidCells === 0));
  t.diagnostic(JSON.stringify({ directSweeps: direct.sweeps, dampedSweeps: damped.sweeps, difference,
    massError: damped.massError, derivativeShear: damped.maximumWallDerivativeShear }));
});
