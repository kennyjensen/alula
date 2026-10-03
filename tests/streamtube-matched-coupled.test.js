import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { directChannelConservation } from './oracles/streamtube.js';

const maximum = a => Math.max(...Array.from(a, Math.abs));

for (const [name, file, unknowns, stations] of [
  ['coarse', 'streamtube-matched-coupled-root.json', 3185, 79],
  ['refined', 'streamtube-matched-refined-root.json', 7919, 137],
]) test(`saved ${name} attached coupled root retains conservation and both Jacobian coupling directions`, t => {
  const fixture = JSON.parse(readFileSync(new URL(`./fixtures/${file}`, import.meta.url)));
  // Reuse the attained root: no panel trace, BL marching or nonlinear solve.
  // Native profile/refinement acceptance is a separate, currently failing gate.
  const system = createCoupledStreamtubeBody(fixture.input, { ...fixture.options,
    initialEuler: fixture.initialEuler, initialBL: fixture.initialBL });
  const x = system.initial, value = system.evaluate(x), j = system.jacobian(x);
  assert.equal(system.n, unknowns); assert.equal(system.bl.stations.length, stations);
  assert.ok(maximum(value.residual) < 1e-10);
  const errors = [];
  for (const mode of ['euler', 'bl', 'both']) {
    const d = x.map((v, i) => i < system.ne
      ? (mode === 'bl' ? 0 : .001 * Math.sin(i * .73 + .2))
      : (mode === 'euler' ? 0 : Math.max(.01, Math.abs(v)) * Math.sin(i * .43 + .4)));
    const h = 2e-6, exact = sparseProduct(j, d);
    const plus = system.residual(x.map((v, i) => v + h * d[i]));
    const minus = system.residual(x.map((v, i) => v - h * d[i]));
    let worst = { error: 0 };
    for (let row = 0; row < system.n; row++) {
      const fd = (plus[row] - minus[row]) / (2 * h);
      const error = Math.abs(fd - exact[row]) / Math.max(1, Math.abs(fd), Math.abs(exact[row]));
      if (error > worst.error) worst = { error, row, fd, exact: exact[row] };
    }
    assert.ok(worst.error < 5e-6, JSON.stringify({ mode, ...worst }));
    // Check actual cross-coupling, so a block-diagonal derivative cannot pass.
    if (mode === 'bl') assert.ok(maximum(exact.slice(0, system.ne)) > 1e-6);
    if (mode === 'euler') assert.ok(maximum(exact.slice(system.ne)) > 1e-6);
    errors.push({ mode, ...worst });
  }
  for (const [g, nodes] of value.outer.nodes.entries()) {
    const c = directChannelConservation({ nodes, sections: value.outer.sections.map(row => row[g]),
      cells: value.outer.cells.map(row => row[g]) }, system.euler.conditions.gamma);
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(maximum(c[key]) < 2e-9, `${g}: ${key}`);
  }
  const rebased = system.rebase(x), after = system.evaluate(rebased);
  assert.ok(maximum(after.residual) < 1e-10);
  value.outer.nodes.forEach((region, g) => region.forEach((row, i) => row.forEach((p, k) => {
    const q = after.outer.nodes[g][i][k]; assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 5e-13);
  })));
  t.diagnostic(JSON.stringify({ residual: maximum(value.residual), errors }));
});
