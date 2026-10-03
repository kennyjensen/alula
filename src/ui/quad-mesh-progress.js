// SPDX-License-Identifier: GPL-2.0-or-later
// Compare accepted grids within the same solve phase/attempt. Starting a
// new coupled guess resets the reference instead of comparing different grids.
import { compareStreamtubeFlow } from '../euler/streamtube-flow-preview.js';
import { quadCoupledNcrit } from './quad-coupled-ncrit.js';
export function createQuadMeshProgress(referenceChord) {
  let phase = 'euler', startupAttempt = 0, initial, previous, initialization, ncrit = {};
  return {
    stage(s) { phase = s.stage; startupAttempt = s.startupAttempt ?? 0; ncrit = quadCoupledNcrit(s); initial = previous = null; },
    mesh(mesh) {
      if (phase === 'euler') initialization = mesh.initialization;
      const actual = { ...ncrit, ...quadCoupledNcrit(mesh), ...quadCoupledNcrit(mesh.iteration) };
      const output = { ...mesh, ...actual, initialization: { ...initialization, ...mesh.initialization,
        ...(initialization?.gridSmoothing ? { gridSmoothing: initialization.gridSmoothing } : {}) } };
      if (!mesh.flow) return output;
      if (!previous || previous.vertices.length !== mesh.vertices.length) initial = previous = mesh;
      const movement = reference => mesh.vertices.reduce((m, p, i) => Math.max(m, Math.hypot(p.x - reference.vertices[i].x, p.y - reference.vertices[i].y)), 0);
      output.flow = compareStreamtubeFlow(mesh.flow, initial.flow, previous.flow);
      output.iteration = { ...mesh.iteration, ...actual, stage: phase, startupAttempt, movementReferenceChord: referenceChord,
        maximumNodeMovement: movement(previous), maximumNodeMovementFromInitial: movement(initial) };
      previous = output; return output;
    },
  };
}
