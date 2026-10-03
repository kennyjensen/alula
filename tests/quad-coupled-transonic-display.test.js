// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';
import { quadCoupledTransonicCoefficients } from '../src/ui/quad-coupled-transonic-coefficients.js';

// Execute the exact archived pre-change adapter with its original imports.
// This compares the complete default output, not a selected coefficient list.
const archiveURL = new URL('../docs/transonic-gui/quad-coupled-result.before.js.txt', import.meta.url);
const originalURL = new URL('../src/ui/quad-coupled-result.js', import.meta.url);
const archive = fs.readFileSync(archiveURL, 'utf8');
assert.equal(createHash('sha256').update(archive).digest('hex'), '0b801823ca1704dd3a3c7e71de136d2d10c7ef0e422ec333912474645d5ec2ef');
const rewritten = archive.replace(/from\s+(['"])(\.{1,2}\/[^'"]+)\1/g,
  (_, quote, specifier) => `from ${quote}${new URL(specifier, originalURL).href}${quote}`);
const old = (await import(`data:text/javascript;base64,${Buffer.from(rewritten).toString('base64')}`)).quadCoupledResultForDisplay;
const read = path => JSON.parse(fs.readFileSync(new URL(path, import.meta.url)));
const fixed = read('../docs/current-coupled-startup-default-browser.json');
const automatic = read('../docs/current-multielement-automatic-16x9-slor-browser.json');
const automaticRaw = () => ({ ...structuredClone(automatic.result), boundaryLayer: structuredClone(automatic.result.numericalBoundaryLayer) });
const historicalRaw = () => {
  const raw = automaticRaw();
  raw.checkpoint.restart.input.streamwiseMode = 'hybrid';
  raw.checkpoint.restart.input.hybrid = { epsilonP: 1e-5 };
  raw.checkpoint.restart.input.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
  raw.checkpoint.restart.options.blThermodynamics = 'historical-common-isentrope';
  raw.solverInput = structuredClone(raw.checkpoint.restart.input);
  raw.conditions.blThermodynamics = 'historical-common-isentrope';
  return raw;
};

test('omitted-mode fixed-trip display is entirely identical to the archived adapter', () => {
  const before = JSON.stringify(fixed.result);
  assert.deepEqual(quadCoupledResultForDisplay(fixed.result, fixed.input), old(fixed.result, fixed.input));
  assert.equal(JSON.stringify(fixed.result), before);
});

test('omitted-mode automatic display and current captured public coefficients retain exact output', () => {
  const raw = automaticRaw(), before = JSON.stringify(raw), mesh = automatic.result.mesh;
  const view = quadCoupledResultForDisplay(raw, automatic.input, mesh);
  assert.deepEqual(view, old(raw, automatic.input, mesh));
  assert.deepEqual(view.coefficients, automatic.result.coefficients);
  assert.equal(JSON.stringify(raw), before);
});

test('explicit historical display uses physical Cp and the identical live coefficient helper at actual Mach', () => {
  const raw = automaticRaw(), checkpoint = raw.checkpoint, input = structuredClone(automatic.input);
  checkpoint.restart.input.streamwiseMode = 'hybrid';
  checkpoint.restart.input.hybrid = { epsilonP: 1e-5 };
  checkpoint.restart.input.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
  checkpoint.restart.options.blThermodynamics = 'historical-common-isentrope';
  raw.solverInput = structuredClone(checkpoint.restart.input);
  raw.conditions.blThermodynamics = 'historical-common-isentrope';
  raw.mach = input.mach = .4; // Deliberately retain a different requested Mach.
  raw.converged = false; raw.continuation = { currentMach: .2, targetMach: .4, reachedTarget: false };
  const before = JSON.stringify(raw), view = quadCoupledResultForDisplay(raw, input);
  assert.equal(view.mach, .2); assert.equal(view.coefficients.conditions.mach, .2);
  assert.match(view.pressureKind, /physical Euler interface/);
  assert.ok(view.warnings.some(w => /Requested Mach 0.4.*Mach 0.2/.test(w)));
  const pInf = 1 / (1.4 * .2 * .2);
  for (const surface of view.boundaryLayer.surfaces) for (const p of surface.stations) {
    const pressure = surface.side === 'upper' ? raw.flow.cells[p.i - 1][p.body + 1][0].interfacePressure.lower
      : raw.flow.cells[p.i - 1][p.body].at(-1).interfacePressure.upper;
    assert.equal(p.cp, 2 * (pressure - pInf));
  }
  const live = quadCoupledTransonicCoefficients({ checkpoint, flow: raw.flow,
    bl: { ...raw.boundaryLayer, scale: 1 / Math.sqrt(raw.kernelReynolds) }, bodies: raw.solverInput.bodies,
    solverLength: raw.solverLength, referenceChord: raw.referenceChord, momentReference: input.momentReference,
    mach: .2, alpha: raw.alpha, gamma: 1.4 });
  assert.deepEqual(view.coefficients, live);
  assert.equal(view.cd, live.viscousDragCoefficient + live.eulerWaveDragCoefficient);
  assert.equal(view.coefficients.eulerExit.sections.length, 39);
  assert.equal(view.coefficients.wakes.length, 2);
  assert.equal(JSON.stringify(raw), before);
});

test('failed target overlay is rejected unless both actual Mach and every physical vertex match the retained historical result', () => {
  const raw = historicalRaw(); raw.converged = false; raw.stateConverged = true;
  raw.continuation = { method: 'freestream-mach', currentMach: .2, targetMach: .4, reachedTarget: false };
  raw.mesh.flow.sourceMarker = 'retained-source';
  const before = JSON.stringify(raw);
  for (const variant of ['wrong-mach-and-grid', 'wrong-mach', 'wrong-grid', 'missing-mach', 'conflicting-mach']) {
    const latest = { vertices: structuredClone(raw.mesh.vertices), mach: .2,
      flow: { targetMarker: 'failed-target', lines: [] },
      iteration: { iteration: 999, maximumNodeMovement: 999, targetMarker: 'failed-target' },
      initialization: { targetMarker: 'failed-target' } };
    if (variant.includes('wrong-mach')) latest.mach = .4;
    if (variant.includes('grid')) latest.vertices[0].x += .01;
    if (variant === 'missing-mach') delete latest.mach;
    if (variant === 'conflicting-mach') latest.actualMach = .4;
    const latestBefore = structuredClone(latest), view = quadCoupledResultForDisplay(raw, automatic.input, latest);
    assert.deepEqual(view.mesh.vertices, raw.mesh.vertices, variant);
    assert.deepEqual(view.mesh.flow, raw.mesh.flow, variant);
    assert.deepEqual(view.mesh.iteration, raw.mesh.iteration, variant);
    assert.equal(view.mesh.initialization.targetMarker, undefined, variant);
    assert.equal(view.mesh.initialization.flowSolved, true, 'The retained state itself converged.');
    assert.equal(view.mach, .2); assert.deepEqual(latest, latestBefore);
  }
  assert.equal(JSON.stringify(raw), before);
});

test('generic Mach-continuation fallback creates only its own actual flow snapshot when no stored overlay exists', () => {
  const raw = automaticRaw(); raw.converged = false; raw.stateConverged = true; raw.mach = .4;
  raw.machContinuation = { actualMach: .2, targetMach: .4, reachedTarget: false };
  delete raw.mesh.flow; delete raw.mesh.iteration;
  const latest = { mach: .4, vertices: raw.mesh.vertices.map(p => ({ x: p.x + .1, y: p.y })),
    flow: { targetMarker: 'failed-target', changeReference: 'failed target', lines: [] },
    iteration: { iteration: 999, maximumNodeMovement: 999 }, initialization: { targetMarker: 'failed-target' } };
  const view = quadCoupledResultForDisplay(raw, automatic.input, latest);
  assert.equal(view.mach, .2); assert.deepEqual(view.mesh.vertices, raw.mesh.vertices);
  assert.equal(view.mesh.flow.targetMarker, undefined); assert.equal(view.mesh.flow.changeReference, undefined);
  assert.equal(view.mesh.initialization.targetMarker, undefined);
  assert.equal(Object.hasOwn(view.mesh, 'iteration'), false);
  assert.equal(view.mesh.flow.lines.length, 39);
  assert.equal(view.mesh.flow.lines[0].speedRatios[0], raw.flow.sections[0][0][0].q);
  const a = raw.flow.nodes[0][0][0], b = raw.flow.nodes[0][0][1];
  assert.deepEqual(view.mesh.flow.lines[0].points[0], { x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
  assert.deepEqual(view.coefficients, automatic.result.coefficients, 'Generic fallback coefficients use actual source Mach.');
  assert.ok(view.warnings.some(w => /Requested Mach 0.4.*Mach 0.2/.test(w)));
});

test('a matching actual-Mach and physical-grid progress snapshot retains its comparison overlay', () => {
  const nonhistoricalCache = automaticRaw();
  nonhistoricalCache.machContinuation = { actualMach: .2, targetMach: .2, reachedTarget: true };
  for (const raw of [historicalRaw(), nonhistoricalCache]) {
    const latest = { mach: .2, actualMach: .2, vertices: structuredClone(raw.mesh.vertices),
      flow: { ...raw.mesh.flow, acceptedMarker: true }, iteration: { iteration: 7, mach: .2, maximumNodeMovement: .001 },
      initialization: { acceptedMarker: true } };
    const view = quadCoupledResultForDisplay(raw, automatic.input, latest);
    assert.equal(view.mesh.flow, latest.flow); assert.equal(view.mesh.iteration, latest.iteration);
    assert.equal(view.mesh.initialization.acceptedMarker, true);
    assert.equal(view.mesh.flow.lines.length, 39);
    if (raw === nonhistoricalCache) assert.deepEqual(view.coefficients, automatic.result.coefficients);
  }
});

test('Ncrit continuation retains actual pressure and rejects an overlay from a different Ncrit', () => {
  const raw = automaticRaw(), originalView = quadCoupledResultForDisplay(raw, automatic.input), request = { ...automatic.input, ncrit: 10 };
  raw.actualNcrit = 9; raw.targetNcrit = 10;
  raw.ncritContinuation = { actualNcrit: 9, targetNcrit: 10, reachedTarget: false };
  const original = JSON.stringify(raw), latest = { mach: raw.mach, actualNcrit: 10,
    vertices: structuredClone(raw.mesh.vertices), flow: { wrongNcrit: true, lines: [] },
    iteration: { iteration: 999, actualNcrit: 10 }, initialization: { wrongNcrit: true } };
  const view = quadCoupledResultForDisplay(raw, request, latest);
  assert.equal(view.converged, false); assert.equal(view.stateConverged, true);
  assert.equal(view.actualNcrit, 9); assert.equal(view.targetNcrit, 10);
  assert.equal(view.sourceCase.ncrit, 9); assert.equal(view.requestedCase.ncrit, 10);
  assert.deepEqual(view.coefficients, automatic.result.coefficients);
  assert.deepEqual(view.elements.map(e => e.cp), originalView.elements.map(e => e.cp));
  assert.equal(view.mesh.flow.wrongNcrit, undefined);
  assert.equal(view.mesh.initialization.wrongNcrit, undefined);
  assert.equal(view.coefficientStatus, 'unconverged');
  assert(view.warnings.some(w => /Requested Ncrit 10.*Ncrit 9/.test(w)));
  assert.equal(JSON.stringify(raw), original);
});
