import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { redistributeStreamtubeTangentially, assembleTangentialCoordinate } from '../src/geometry/streamtube-tangential-redistribution.js';
import { solveAlternatingScalarLines } from '../src/numerics/alternating-scalar-lines.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

// Inject a controlled candidate-domain rejection after the real tiny Euler
// evaluation. This tests initial retry/commit control, without a Newton solve
// or tuning a manufactured flow to a particular gas rejection threshold.
async function driver({ rejectedCandidates = 0, failMovement = false } = {}) {
  const url = new URL('../src/euler/streamtube-ises-update.js', import.meta.url), calls = [], systems = [];
  const key = `ises-initial-redistribution-test-${Math.random()}`;
  globalThis[key] = {
    createStreamtubeBodySystem(input) {
      const system = createStreamtubeBodySystem(input), index = systems.length, evaluate = system.evaluate;
      systems.push(system);
      system.evaluate = x => {
        const result = evaluate(x);
        if (index > 0 && index <= rejectedCandidates) throw new Error('Controlled candidate interface pressure rejection.');
        return result;
      };
      return system;
    },
    solveStreamtubeBody,
    redistributeStreamtubeTangentially(nodes, options) {
      calls.push({ nodes: structuredClone(nodes), options: structuredClone(options) });
      if (failMovement && options.correctionScale > .5) throw new Error('Controlled nonpositive coordinate span.');
      return redistributeStreamtubeTangentially(nodes, options);
    },
  };
  let text = fs.readFileSync(url, 'utf8');
  text = text.replace("import { createStreamtubeBodySystem, solveStreamtubeBody } from './streamtube-body.js';",
    `const { createStreamtubeBodySystem, solveStreamtubeBody } = globalThis[${JSON.stringify(key)}];`);
  text = text.replace("import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';",
    `const { redistributeStreamtubeTangentially } = globalThis[${JSON.stringify(key)}];`);
  text = text.replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try {
    const { solveStreamtubeIses } = await import(`data:text/javascript;base64,${Buffer.from(text).toString('base64')}`);
    return { solve: solveStreamtubeIses, calls, systems };
  } finally { delete globalThis[key]; }
}

test('coordinate scaling changes the full SMOVE coordinate field, including its span, while retaining the original five-pair solution', () => {
  const nx = 7, nt = 4;
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => ({
    x: i / nx + .075 * Math.sin(Math.PI * i / nx) * Math.sin(Math.PI * j / nt),
    y: j / nt + .025 * Math.sin(2 * Math.PI * i / nx) * Math.sin(Math.PI * j / nt),
  })));
  const before = structuredClone(nodes), scale = .25;
  const full = redistributeStreamtubeTangentially(nodes), explicitFull = redistributeStreamtubeTangentially(nodes, { correctionScale: 1 });
  assert.deepEqual(explicitFull, full);
  const a = assembleTangentialCoordinate(nodes), solved = solveAlternatingScalarLines(a.matrix, a.rhs, { ...a.lineDimensions, pairs: 5 });
  const actual = redistributeStreamtubeTangentially(nodes, { correctionScale: scale });
  assert.deepEqual(actual.solution, full.solution);
  const c = (i, j) => i > 0 && i < nx && j > 0 && j < nt ? solved.x[(i - 1) * (nt - 1) + j - 1] : 0;
  let differsFromBlendedFull = false;
  for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) {
    if (i === 0 || i === nx || j === 0 || j === nt) assert.deepEqual(actual.nodes[i][j], nodes[i][j]);
    else {
      const span = a.increments[i] + a.increments[i - 1] + scale * c(i + 1, j) - scale * c(i - 1, j);
      const fraction = -scale * c(i, j) / span;
      for (const key of ['x', 'y']) {
        const expected = nodes[i][j][key] + fraction * (nodes[i + 1][j][key] - nodes[i - 1][j][key]);
        assert.equal(actual.nodes[i][j][key], expected);
        if (Math.abs(expected - (nodes[i][j][key] + scale * (full.nodes[i][j][key] - nodes[i][j][key]))) > 1e-6)
          differsFromBlendedFull = true;
      }
    }
  }
  assert.equal(differsFromBlendedFull, true);
  assert.deepEqual(nodes, before);
  for (const correctionScale of [0, -1, 1.01, Infinity, NaN, null])
    assert.throws(() => redistributeStreamtubeTangentially(nodes, { correctionScale }), /correction scale/);
});

test('admissible initial SMOVE retries all three passages from the original physical state and commits only a fully evaluated candidate', async () => {
  const { solve, calls, systems } = await driver({ rejectedCandidates: 2 });
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), before = structuredClone(input), publications = [];
  const result = solve(input, { maxIterations: 0, stepAcceptance: 'admissible', maxBacktracks: 2,
    onCheckpoint: (cp, details) => publications.push({ cp, details }) });
  assert.equal(result.initialRedistribution.accepted, true, result.reason);
  assert.equal(result.initialRedistribution.correctionScale, .25);
  assert.equal(result.initialRedistribution.backtracks, 2);
  assert.deepEqual(result.initialRedistribution.rejections.map(r => [r.correctionScale, r.stage]), [[1, 'admissibility'], [.5, 'admissibility']]);
  assert.deepEqual(calls.map(c => c.options.correctionScale), [1, 1, 1, .5, .5, .5, .25, .25, .25]);
  for (let k = 3; k < calls.length; k++) assert.deepEqual(calls[k].nodes, calls[k % 3].nodes);
  assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(result.finalQuality.valid, true);
  assert.equal(publications.length, 1);
  assert.deepEqual(publications[0].cp.initialEuler.nodes, result.nodes);
  assert.deepEqual(publications[0].cp.residual, Array.from(result.residual));
  assert.deepEqual(publications[0].cp.initialEuler.x.slice(0, systems[0].layout.densityCount), Array.from(systems[0].initial.slice(0, systems[0].layout.densityCount)));
  assert.deepEqual(input, before);
});

test('span rejection uses the same bounded initial policy; listing and exhausted retries retain the original state without publishing a checkpoint', async () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), before = structuredClone(input);
  const span = await driver({ failMovement: true });
  const recovered = span.solve(input, { maxIterations: 0, stepAcceptance: 'admissible', maxBacktracks: 1 });
  assert.equal(recovered.initialRedistribution.accepted, true);
  assert.equal(recovered.initialRedistribution.correctionScale, .5);
  assert.equal(recovered.initialRedistribution.rejections[0].stage, 'SMOVE');
  for (const [policy, expectedAttempts] of [['listing', 1], ['admissible', 2]]) {
    const h = await driver({ rejectedCandidates: 10 });
    const result = h.solve(input, { maxIterations: 0, stepAcceptance: policy, maxBacktracks: 1,
      onCheckpoint: () => assert.fail('Rejected initial movement must not publish a restart.') });
    assert.equal(result.initialRedistribution.accepted, false);
    assert.match(result.reason, /initial redistribution rejected/);
    assert.equal(h.systems.length, 1 + expectedAttempts);
    assert.equal(result.linearDiagnostics.solves, 0);
    assert.deepEqual(result.nodes, h.systems[0].evaluate(h.systems[0].initial).nodes);
    assert.deepEqual(result.residual, h.systems[0].evaluate(h.systems[0].initial).residual);
  }
  assert.deepEqual(input, before);
});

test('already valid full SMOVE and zero-update resume are unchanged between listing and admissible policies', async () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), a = await driver(), b = await driver();
  const listing = a.solve(input, { maxIterations: 0, stepAcceptance: 'listing', retainCheckpoint: true });
  const admissible = b.solve(input, { maxIterations: 0, stepAcceptance: 'admissible', retainCheckpoint: true });
  assert.deepEqual(listing.initialRedistribution, admissible.initialRedistribution);
  for (const key of ['nodes', 'residual', 'history', 'linearDiagnostics']) assert.deepEqual(admissible[key], listing[key]);
  assert.equal('correctionScale' in admissible.initialRedistribution, false);
  const c = await driver({ failMovement: true });
  const resumed = c.solve(undefined, { resume: admissible.checkpoint, maxIterations: 0, stepAcceptance: 'admissible' });
  assert.equal(c.calls.length, 0);
  assert.equal(resumed.initialRedistribution.resumed, true);
  assert.deepEqual(resumed.nodes, admissible.nodes);
  assert.deepEqual(resumed.residual, admissible.residual);
});
