import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { streamtubeFlowSnapshot, compareStreamtubeFlow } from '../src/euler/streamtube-flow-preview.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('displayed centerline velocities are the Euler conservation velocities in every tube', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  const flow = system.evaluate(system.initial), display = streamtubeFlowSnapshot(flow);
  assert.equal(display.lines.length, 6); assert.equal(display.iteration, 0);
  assert.equal(display.residual, flow.diagnostics.residual);
  assert.equal(display.speedUnit, 'U_infinity');
  let count = 0;
  for (const line of display.lines) {
    assert.equal(line.points.length, line.speedRatios.length + 1);
    line.speedRatios.forEach((speed, i) => {
      const { group, tube } = line, cell = flow.cells[Math.max(0, i - 1)][group][tube];
      const side = i === 0 ? 0 : 1, direction = cell.geometry.directions[side];
      const a = line.points[i], b = line.points[i + 1], length = Math.hypot(b.x - a.x, b.y - a.y);
      assert.equal(speed, cell.states[side].q);
      assert.equal(line.machNumbers[i], Math.sqrt(cell.states[side].machSquared));
      assert.ok(Math.abs((b.x - a.x) / length - direction.x) < 1e-14);
      assert.ok(Math.abs((b.y - a.y) / length - direction.y) < 1e-14);
      // Independently recover speed from the conserved tube mass flux.
      const mass = flow.allocation.groups[group][tube].massFlow;
      assert.equal(speed, mass / (cell.states[side].rho * cell.geometry.normalAreas[side]));
      count++;
    });
  }
  assert.equal(count, system.layout.nx * 6);
  assert.deepEqual(structuredClone(display), display);
  const retained = structuredClone(display);
  flow.nodes[0][0][0].x += 100; flow.sections[0][0][0].q = 123;
  assert.deepEqual(display, retained, 'snapshot must own its coordinates and speeds');
});

test('live flow snapshots match accepted iterates and do not change the numerical trajectory', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const baseline = solveStreamtubeBody(createStreamtubeBodySystem(input), { maxIterations: 3, stepMethod: 'dogleg' });
  const frames = [], result = solveStreamtubeBody(createStreamtubeBodySystem(input), {
    maxIterations: 3, stepMethod: 'dogleg', onMesh: state => {
      const snapshot = streamtubeMeshSnapshot(state);
      assert.equal(snapshot.flow.residual, state.iteration.residual);
      assert.equal(snapshot.flow.iteration, state.iteration.iteration);
      frames.push(snapshot.flow);
    },
  });
  assert.deepEqual(result.x, baseline.x); assert.deepEqual(result.history, baseline.history);
  assert.deepEqual(result.nodes, baseline.nodes); assert.deepEqual(result.sections, baseline.sections);
  assert.deepEqual(frames.map(f => f.iteration), [1, 2, 3]);
  for (let k = 1; k < frames.length; k++) {
    assert.notDeepEqual(frames[k].lines.map(l => l.speedRatios), frames[k - 1].lines.map(l => l.speedRatios));
    assert.notDeepEqual(frames[k].lines.map(l => l.points), frames[k - 1].lines.map(l => l.points));
  }
  assert.deepEqual(frames.at(-1), streamtubeFlowSnapshot(result, 3));
  const retained = structuredClone(frames), compared = compareStreamtubeFlow(frames[2], frames[0], frames[1]);
  const differences = frames[2].lines.flatMap((line, k) => line.speedRatios.map((q, i) => q - frames[0].lines[k].speedRatios[i]));
  assert.deepEqual(compared.lines.flatMap(line => line.speedChanges), differences);
  assert.equal(compared.minimumSpeedChange, Math.min(...differences));
  assert.equal(compared.maximumSpeedChange, Math.max(...differences));
  assert.equal(compared.maximumSpeedChangeFromPrevious, Math.max(...frames[2].lines.flatMap((line, k) =>
    line.speedRatios.map((q, i) => Math.abs(q - frames[1].lines[k].speedRatios[i])))));
  assert.deepEqual(frames, retained, 'comparing must not overwrite the original flow fields');
  assert.ok(compareStreamtubeFlow(frames[0], frames[0]).lines.every(line => line.speedChanges.every(q => q === 0)));
  assert.throws(() => compareStreamtubeFlow(frames[2], { ...frames[0], lines: frames[0].lines.slice(1) }), /corresponding/);
});
