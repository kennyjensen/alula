// SPDX-License-Identifier: GPL-2.0-or-later
// Read-only presentation: dimensional BL profiles in original element order.
// This never integrates displaced-wall traction or changes a solved state.
import { createContourCurve } from '../geometry/contour-curve.js';
import { prepareContour } from '../geometry/airfoil.js';
import { createIntegralKernel } from '../viscous/integral.js';
import { isentropicState } from '../potential/isentropic.js';
import { quadCoupledCoefficients } from './quad-coupled-coefficients.js';
import { createContourArc } from '../geometry/contour-arc.js';
import { quadCoupledTransonicCoefficientsFromResult } from './quad-coupled-transonic-coefficients.js';
import { streamtubeFlowSnapshot } from '../euler/streamtube-flow-preview.js';
import { createSurfaceContourCurve } from '../geometry/contour-topology.js';
import { quadCoupledGridFailure } from './quad-coupled-failure.js';
import { quadCoupledNcrit, quadCoupledNcritResult } from './quad-coupled-ncrit.js';
import { quadCoupledGrid, quadCoupledGridResult } from './quad-coupled-grid.js';

// The transition root uses physical distance from the solved stagnation
// point, normalized by solverLength. Convert it back onto the solid contour.
function automaticTransitionLocation(curve, stagnation, surface, transition, length) {
  const arc = createContourArc(curve), origin = arc.at(stagnation), upper = surface.side === 'upper';
  const total = upper ? origin : arc.length - origin;
  const target = origin + (upper ? -1 : 1) * transition.s * length;
  const tolerance = 32 * Number.EPSILON * arc.length;
  if (target < -tolerance || target > arc.length + tolerance) throw new Error('Solved transition location is outside the solid contour.');
  let lo = 0, hi = curve.length;
  for (let i = 0; i < 52; i++) { const mid = .5 * (lo + hi); if (arc.at(mid) < target) lo = mid; else hi = mid; }
  return { transitionPoint: curve.evaluate(.5 * (lo + hi)).point, transitionFraction: transition.s * length / total,
    transitionDistance: transition.s * length, transitionKind: transition.kind === 'trailing-edge' ? 'laminar-to-te' : transition.kind,
    forced: transition.kind === 'forced' };
}

// A failed Mach stage may leave the last progress frame on another physical
// grid. Progress overlays are usable only with evidence that they describe
// this returned state; a requested-Mach label is not that evidence.
function compatibleMachMesh(latest, actual, mach, actualNcrit) {
  const tags = [latest?.actualMach, latest?.mach, latest?.iteration?.actualMach, latest?.iteration?.mach,
    latest?.flow?.actualMach, latest?.flow?.mach].filter(value => value !== undefined);
  const ncritTags = [latest?.actualNcrit, latest?.iteration?.actualNcrit, latest?.flow?.actualNcrit]
    .filter(value => value !== undefined);
  return tags.length > 0 && tags.every(value => Number.isFinite(value) && value === mach)
    && (actualNcrit === undefined || ncritTags.length > 0 && ncritTags.every(value => value === actualNcrit))
    && Array.isArray(latest?.vertices) && Array.isArray(actual?.vertices)
    && latest.vertices.length === actual.vertices.length
    && latest.vertices.every((p, i) => p && Number.isFinite(p.x) && Number.isFinite(p.y)
      && p.x === actual.vertices[i]?.x && p.y === actual.vertices[i]?.y);
}

export function quadCoupledResultForDisplay(raw, input, latestMesh) {
  raw = quadCoupledGridResult(quadCoupledNcritResult(raw, input), input);
  const grid = quadCoupledGrid(raw, input);
  const ncrit = quadCoupledNcrit(raw), ncritNotReached = ncrit.actualNcrit !== ncrit.targetNcrit
    || raw.ncritContinuation?.reachedTarget === false;
  const length = raw.solverLength, referenceChord = raw.referenceChord;
  if (raw.model !== 'research-streamtube-euler-bl' || !(length > 0) || !(referenceChord > 0))
    throw new Error('Expected a normalized coupled quad result.');
  const historical = raw.conditions.blThermodynamics === 'historical-common-isentrope';
  const failureDiagnostic = quadCoupledGridFailure(raw, input);
  const machContinuation = raw.machContinuation ?? (raw.continuation?.method === 'freestream-mach' ? raw.continuation : null);
  const guardedMesh = historical || machContinuation !== null || Object.keys(ncrit).length > 0 || Object.keys(grid).length > 0;
  const gas = { mach: guardedMesh ? raw.conditions.mach : raw.mach, gamma: raw.solverInput.gamma ?? 1.4 };
  const kernel = createIntegralKernel({ ...gas, reynolds: raw.kernelReynolds, ncrit: raw.conditions.ncrit });
  const bodies = raw.solverInput.bodies, original = raw.boundaryLayer;
  if (bodies.length !== input.elements.length || new Set(bodies.map(b => b.element)).size !== bodies.length)
    throw new Error('Coupled display requires the original element mapping.');
  const points = new Map(), surfaces = original.surfaces.map(surface => {
    const transition = original.transitions?.find(t => t.body === surface.body && t.side === surface.side);
    const body = bodies[surface.body], element = body.element;
    const group = surface.side === 'upper' ? surface.body + 1 : surface.body;
    const bank = row => surface.side === 'upper' ? row[0] : row.at(-1);
    const stations = surface.ids.map(id => {
      const state = original.stations[id], solid = bank(raw.flow.undisplacedNodes[group][state.i]);
      const displacement = bank(raw.flow.nodes[group][state.i]);
      const turbulent = ['leading-transition', 'transition', 'turbulent'].includes(state.regime);
      const terminal = transition?.kind === 'trailing-edge' && id === surface.ids.at(-1);
      const properties = kernel.station(state, turbulent && !terminal ? 'turbulent' : 'laminar');
      const p = { ...state, ...solid, index: id, s: state.s * length, theta: state.theta * length,
        deltaStar: state.deltaStar * length, h: state.deltaStar / state.theta, hk: properties.hk,
        cf: properties.cf * properties.rho * state.ue ** 2,
        cp: historical ? 2 * ((surface.side === 'upper'
          ? raw.flow.cells[state.i - 1][surface.body + 1][0].interfacePressure.lower
          : raw.flow.cells[state.i - 1][surface.body].at(-1).interfacePressure.upper) - 1 / (gas.gamma * gas.mach * gas.mach))
          : isentropicState(state.ue, 0, gas).cp, displacement: { ...displacement },
        amplification: terminal ? transition.amplification : turbulent ? null : state.aux, ctau: turbulent ? state.aux : null,
        ...(terminal ? { regime: 'laminar' } : {}) };
      points.set(id, p); return p;
    });
    const curve = body.trailingEdge?.kind === 'finite-base'
      ? createSurfaceContourCurve(body.points, body) : createContourCurve(prepareContour(body.points));
    return { element, name: input.elements[element].name ?? `Element ${element + 1}`, side: surface.side,
      ...(transition ? automaticTransitionLocation(curve, raw.flow.stagnation[surface.body], surface, transition, length)
        : { forced: true, transitionKind: 'forced', transitionPoint: curve.evaluate(surface.tripParameter).point }), stations };
  });
  surfaces.sort((a, b) => a.element - b.element || (a.side === 'upper' ? -1 : 1));
  const wakes = original.wakes.map(wake => {
    const element = bodies[wake.body].element;
    return { element, name: input.elements[element].name ?? `Element ${element + 1}`, stations: wake.ids.map(id => {
      const state = original.stations[id], a = raw.flow.nodes[wake.body][state.i].at(-1), b = raw.flow.nodes[wake.body + 1][state.i][0];
      const p = { ...state, index: id, x: .5 * (a.x + b.x), y: .5 * (a.y + b.y), s: state.s * length,
        theta: state.theta * length, deltaStar: state.deltaStar * length, h: state.deltaStar / state.theta, cf: 0,
        ...(state.wakeGap === undefined ? {} : { wakeGap: state.wakeGap * length,
          fluidDeltaStar: (state.deltaStar - state.wakeGap) * length, fluidH: (state.deltaStar - state.wakeGap) / state.theta }) };
      points.set(id, p); return p;
    }) };
  }).sort((a, b) => a.element - b.element);
  const elements = input.elements.map((e, element) => ({ ...e, cp: [
    ...surfaces.find(s => s.element === element && s.side === 'upper').stations.slice().reverse(),
    ...surfaces.find(s => s.element === element && s.side === 'lower').stations,
  ].map(({ x, y, cp }) => ({ x, y, cp })) }));
  const progressMesh = !guardedMesh || compatibleMachMesh(latestMesh, raw.mesh, gas.mach, ncrit.actualNcrit) ? latestMesh : undefined;
  const mesh = { ...raw.mesh, initialization: { ...progressMesh?.initialization, ...raw.mesh.initialization,
    gridSmoothing: raw.mesh.initialization?.gridSmoothing ?? raw.initialization?.euler?.gridSmoothing,
    flowSolved: guardedMesh ? raw.stateConverged ?? raw.converged : raw.converged },
    ...(progressMesh?.iteration ? { iteration: progressMesh.iteration } : {}), ...(progressMesh?.flow ? { flow: progressMesh.flow } : {}) };
  if (guardedMesh && !mesh.flow) mesh.flow = streamtubeFlowSnapshot(raw.flow, raw.history.at(-1)?.iteration ?? 0);
  const momentReference=input.momentReference??{x:referenceChord/4,y:0},stagnation=[];
  bodies.forEach((body,b)=>{stagnation[body.element]=raw.flow.undisplacedNodes[b+1][body.leadingIndex][0];});
  const coefficients=historical ? quadCoupledTransonicCoefficientsFromResult(raw,{momentReference})
    : quadCoupledCoefficients({surfaces,wakes,stagnation,bodies,alpha:raw.alpha,mach:gas.mach,
      gamma:gas.gamma,referenceChord,momentReference});
  return { ...raw, elements, mesh, numericalBoundaryLayer: original,
    ...(failureDiagnostic ? { solverReason: raw.reason, reason: failureDiagnostic.message, failureDiagnostic } : {}),
    ...(guardedMesh ? { mach: gas.mach } : {}),
    cl:coefficients.cl,cm:coefficients.cm,cd:coefficients.cd,coefficients,momentReference,
    coefficientStatus:raw.converged?'research-unvalidated':'unconverged',
    forceStatus:historical ? 'Unvalidated physical Euler pressure lift/moment and summed viscous-wake plus Euler-entropy exit drag; shear lift/moment omitted.'
      : 'Unvalidated estimates: solid-contour BL-edge pressure lift/moment and separate-wake momentum drag.',
    boundaryLayer: { lengthUnit: 'case-coordinate', surfaces, wakes, stations: original.stations.map(s => points.get(s.id)) },
    pressureKind: historical ? coefficients.pressureKind : 'isentropic BL-edge pressure from solved edge speed',
    diagnostics: { ...raw.flow.diagnostics, equationResidual: Math.max(...Object.values(raw.families)),
      ...(raw.convergence ? { convergence: raw.convergence, residualConverged: raw.residualConverged } : {}),
      eulerResidual: raw.families.euler, boundaryLayerResidual: raw.families.boundaryLayer, edgeMatchingResidual: raw.families.edgeMatching,
      iterations: raw.history.length - 1, unknowns: raw.x.length, cells: mesh.cells.length,
      ...(failureDiagnostic ? { failure: failureDiagnostic } : {}) },
    warnings: [...(raw.alphaContinuation?.reachedTarget === false
      ? [`Requested alpha ${raw.targetAlpha}° was not reached. Displayed coefficients belong to alpha ${raw.actualAlpha}°.`] : []),
      ...(grid.differentGrid ? [`Requested grid (${grid.requestedGridIntervals} surface intervals per side) was not reached. Displayed coefficients and profiles belong to the retained ${grid.actualGridIntervals}-interval grid.`] : []),
      ...(ncritNotReached ? [`Requested Ncrit ${ncrit.targetNcrit} was not reached. Displayed coefficients and profiles belong to Ncrit ${ncrit.actualNcrit ?? 'unknown'}.`] : []),
      ...(!raw.converged ? [(machContinuation ?? (historical ? raw.continuation : null))?.reachedTarget === false
      ? `Requested Mach ${(machContinuation ?? raw.continuation).targetMach} was not reached. Displayed profiles are at Mach ${gas.mach}.`
      : `Euler/BL equations are unconverged: ${failureDiagnostic?.message ?? raw.reason}. Displayed profiles are provisional.`] : []),
      'Research Euler/BL model: physical accuracy and aerodynamic forces remain unvalidated.',
      ...coefficients.warnings,
      ...(original.transitions ? ['Automatic eᴺ transition is solved with earlier prescribed trips taking precedence. Transition/profile accuracy remains unvalidated; wake merging is not implemented.']
        : ['Fixed material trips are used. Select Automatic (eᴺ) to solve for natural transition. Wake merging is not implemented.'])] };
}
