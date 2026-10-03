// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

const read = relative => JSON.parse(fs.readFileSync(new URL(relative, import.meta.url)));
const serial = value => JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v);

test('MSES RAE warm MRCHDU profile moves the interior before it overtakes the first tube', () => {
  const parent = read('../docs/rae2822/mses-fine047-stage-tail/accepted-parent.json');
  const unextended = read('../docs/rae2822/mses047-mrchdu-fold/target-before-evaluate.json');
  const checkpoint = parent.checkpoint, before = serial(checkpoint);
  const prepared = initializeCoupledStreamtubeFromFlow(.74, checkpoint, { blPredictor: 'xfoil-mrchdu' });
  const { system, initial, value } = prepared;
  assert.equal(streamtubeMeshSnapshot({ system: system.euler, nodes: unextended.initialEuler.nodes }).quality.valid, false,
    'The saved original target must still expose the regression.');
  assert.equal(streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes }).quality.valid, true);
  const independent = directStreamtubeVolumeGeometry(value.outer.nodes);
  assert.equal(independent.valid, true, 'All primal, half and dual volumes and global embedding must pass.');
  assert.equal(serial(initial.subarray(system.ne)), serial(unextended.initialBL), 'Retain the native predicted profile.');
  assert.equal(serial(system.bl.snapshotActive()), serial(unextended.options.transitionState));
  assert.equal(serial(checkpoint), before, 'Do not modify the accepted source checkpoint.');
  assert.equal(serial(prepared.checkpoint.continuation), serial(checkpoint.continuation));
  for (const key of ['captured', 'stagnation', 'strengths'])
    assert.equal(serial(value.outer[key]), serial(parent.flow[key]));
  for (let i = 0; i < value.outer.sections.length; i++) for (let g = 0; g < system.euler.layout.tubes.length; g++)
    for (let j = 0; j < system.euler.layout.tubes[g]; j++)
      assert.equal(value.outer.sections[i][g][j].rho, parent.flow.sections[i][g][j].rho);
  for (let g = 0; g < system.euler.layout.tubes.length; g++) for (let j = 0; j < system.euler.layout.tubes[g]; j++)
    assert.equal(value.outer.allocation.groups[g][j].massFlow, parent.flow.allocation.groups[g][j].massFlow);
  for (let b = 0; b < system.euler.layout.elements; b++) {
    const body = system.euler.layout.bodies[b];
    for (const [g, j] of [[b, system.euler.layout.tubes[b]], [b + 1, 0]]) {
      for (let i = body.leadingIndex; i <= body.trailingIndex; i++)
        assert.deepEqual(value.outer.undisplacedNodes[g][i][j], parent.flow.undisplacedNodes[g][i][j]);
      for (let i = body.trailingIndex + 1; i <= system.euler.layout.nx; i++)
        assert.deepEqual(value.outer.nodes[g][i][j], parent.flow.nodes[g][i][j]);
    }
  }
  assert.ok(prepared.diagnostics.targetBLPrediction.interiorGridExtension.interiorGeometryDeparture > 0);
  assert.equal(prepared.diagnostics.operations.newtonUpdates, 0);
  assert.equal(prepared.diagnostics.targetConverged, false, 'An admissible initial guess is not a solved flow.');
  const f = JSON.parse(serial(prepared.checkpoint)).restart;
  const replay = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const replayed = replay.evaluate(replay.initial);
  assert.deepEqual(replay.initial, initial);
  assert.deepEqual(replayed.residual, value.residual);
  assert.deepEqual(replayed.outer.nodes, value.outer.nodes);
});
