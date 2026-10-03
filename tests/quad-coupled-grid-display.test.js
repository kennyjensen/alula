// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledGrid, quadCoupledGridLabel, quadCoupledGridResult } from '../src/ui/quad-coupled-grid.js';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';

const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json', import.meta.url)));
const input = { ...saved.input, gridIntervals: 128 };
const retained = () => ({ ...structuredClone(saved.result), boundaryLayer: structuredClone(saved.result.numericalBoundaryLayer),
  gridSequence: { kind: 'coarse-to-fine', reachedTarget: false, actualGridIntervals: 16, requestedGridIntervals: 128 },
  reason: 'Controlled grid-level transfer failure' });

test('a retained grid cannot be promoted to the requested solve by a stale converged flag', () => {
  const raw = retained(), before = structuredClone(raw), next = quadCoupledGridResult(raw, input);
  assert.equal(raw.converged, true);
  assert.equal(next.converged, false); assert.equal(next.stateConverged, true);
  assert.equal(next.coefficientStatus, 'unconverged');
  assert.equal(next.sourceCase.gridIntervals, 16); assert.equal(next.requestedCase.gridIntervals, 128);
  assert.equal(next.checkpoint, raw.checkpoint); assert.equal(next.flow, raw.flow);
  assert.equal(quadCoupledGridLabel(next), ' · Grid 16 → target 128 intervals/side');
  assert.deepEqual(raw, before);
  const legacy = { converged: true };
  assert.equal(quadCoupledGridResult(legacy, input), legacy);
});

test('an unresolved target-grid solve is unsuccessful without inventing a coarser grid', () => {
  const raw = retained(); delete raw.gridSequence.actualGridIntervals;
  const next = quadCoupledGridResult(raw, input);
  assert.equal(next.converged, false);
  assert.equal(quadCoupledGrid(next, input).differentGrid, false);
  assert.equal(quadCoupledGridLabel(next), '');
  assert.equal(next.sourceCase, undefined);
  const complete = retained();
  complete.gridSequence = { ...complete.gridSequence, reachedTarget: true, actualGridIntervals: 128 };
  assert.equal(quadCoupledGridResult(complete, input), complete);
  const legacySequence = retained();
  legacySequence.automaticRefinement = legacySequence.gridSequence; delete legacySequence.gridSequence;
  assert.equal(quadCoupledGridResult(legacySequence, input).converged, false);
});

test('retained display preserves its Cp and rejects a failed grid overlay even at the same Mach and Ncrit', () => {
  const raw = retained(), before = structuredClone(raw);
  const unsequenced = { ...raw }; delete unsequenced.gridSequence;
  const expected = quadCoupledResultForDisplay(unsequenced, saved.input);
  // Deliberately no Mach/Ncrit continuation: grid provenance alone guards this.
  assert.equal(raw.machContinuation, undefined); assert.equal(raw.actualNcrit, undefined);
  for (const variant of ['changed-vertex', 'changed-topology', 'missing-mach']) {
    const latest = { vertices: structuredClone(raw.mesh.vertices), mach: raw.mach,
      flow: { targetMarker: true, lines: [] }, iteration: { iteration: 999, targetMarker: true },
      initialization: { targetMarker: true } };
    if (variant === 'changed-vertex') latest.vertices[0].x += .01;
    if (variant === 'changed-topology') latest.vertices.pop();
    if (variant === 'missing-mach') delete latest.mach;
    const view = quadCoupledResultForDisplay(raw, input, latest);
    assert.equal(view.converged, false); assert.equal(view.stateConverged, true);
    assert.equal(view.sourceCase.gridIntervals, 16); assert.equal(view.requestedCase.gridIntervals, 128);
    assert.equal(view.coefficientStatus, 'unconverged');
    assert.deepEqual(view.mesh.vertices, raw.mesh.vertices);
    assert.deepEqual(view.mesh.flow, raw.mesh.flow); assert.deepEqual(view.mesh.iteration, raw.mesh.iteration);
    assert.equal(view.mesh.initialization.targetMarker, undefined);
    assert.equal(view.mesh.initialization.flowSolved, true);
    assert.deepEqual(view.elements.map(e => e.cp), expected.elements.map(e => e.cp));
    assert.deepEqual(view.coefficients, saved.result.coefficients);
    assert.ok(view.warnings.some(w => /128 surface intervals.*retained 16-interval/.test(w)));
  }
  assert.deepEqual(raw, before);
});

test('intermediate-grid results display without Euler startup metadata and preserve actual mesh smoothing', () => {
  for (const initialization of [undefined, { method: 'Complete-state nested grid sequencing' }]) {
    const raw = retained();
    raw.initialization = initialization;
    const smoothing = { enabled: true, converged: true };
    raw.mesh.initialization.gridSmoothing = smoothing;
    const before = structuredClone(raw), view = quadCoupledResultForDisplay(raw, input);
    assert.deepEqual(view.mesh.initialization.gridSmoothing, smoothing);
    assert.equal(view.stateConverged, true);
    assert.equal(view.converged, false);
    assert.ok(['cl', 'cd', 'cm'].every(k => Number.isFinite(view[k])));
    assert.deepEqual(view.coefficients, saved.result.coefficients);
    assert.deepEqual(raw, before);
    delete raw.mesh.initialization.gridSmoothing;
    assert.equal(quadCoupledResultForDisplay(raw, input).mesh.initialization.gridSmoothing, undefined);
  }
});
