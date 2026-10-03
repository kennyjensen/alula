// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4Standard } from '../src/geometry/airfoil.js';
import { prepareAirfoilElement } from '../src/geometry/airfoil-element.js';
import { prepareStreamtubeAssembly } from '../src/euler/streamtube-result.js';
import { initializeStreamtubeStartup } from '../src/euler/streamtube-startup.js';

test('finite NACA0012 128/11 builds a convex smoothed mesh and physical gas state', () => {
  const element = prepareAirfoilElement({ name: 'Main element', points: naca4Standard('0012', 160) });
  const before = structuredClone(element);
  const prepared = prepareStreamtubeAssembly({ elements: [element], alpha: 4, mach: .2,
    referenceChord: 1, gridIntervals: 128, gridTubes: 11, eulerIsmom: 4,
    gridCrosslinePlacement: 'potential', gridChordExponent: 0, gridSurfaceSpacing: 'supplied',
    gridEllipticSmoothing: true, gridSmoothingMethod: 'elliptic' });
  assert.deepEqual(element, before);
  const connector = prepared.diagnostics.geometricStagnationConnectors[0];
  assert.equal(connector.trigger, 'STAGNATION_GUIDE_REFERENCE_MISMATCH');
  assert.equal(connector.sourcePanel, 83);
  assert.ok(Math.abs(connector.joinPotentialError) < 8 * connector.potentialTolerance);
  assert.match(prepared.diagnostics.finiteBaseInteriorTreatment, /Retain panel-traced/);
  assert.equal(prepared.mesh.quality.valid, true);
  assert.deepEqual(prepared.mesh.quality.invalidCells, []);
  assert.ok(prepared.system.layout.tubes.every(n => n === 14));
  const startup = initializeStreamtubeStartup(prepared.system, prepared.initial);
  assert.ok(startup.initial.every(Number.isFinite));
  assert.ok(startup.flow);
  console.log(JSON.stringify({ cells: prepared.mesh.cells.length, quality: prepared.mesh.quality,
    smoothing: prepared.diagnostics.gridSmoothing.converged, gas: startup.diagnostics }));
});
