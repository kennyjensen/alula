// SPDX-License-Identifier: GPL-2.0-or-later
// Observables from an existing hybrid Euler/BL state; no flow or BL residual
// evaluation. Physical Euler pressure/entropy and historical BL gas are kept
// separate. MSES total drag = viscous wake defect + inviscid entropy defect.
import { streamtubeSolidPressureForces } from '../euler/streamtube-forces.js';
import { streamtubeExitDefect } from '../euler/streamtube-exit-defect.js';
import { streamtubeViscousExitDefect } from '../euler/streamtube-viscous-exit-defect.js';
import { quadCoupledWakeGap } from './quad-coupled-wake-gap.js';

const positive = value => Number.isFinite(value) && value > 0;
const vector = value => Array.isArray(value) || ArrayBuffer.isView(value) && !(value instanceof DataView);
const require = (condition, message) => { if (!condition) throw new Error(message); };

// Identical physical-speed gas convention to integral.js / native BLKIN.
// Gas is independent of s, theta and Re; no dummy BL arc coordinate is used.
export function historicalCommonIsentropeWakeGas(ue, { mach, gamma = 1.4 }) {
  require(positive(ue) && positive(mach) && mach < 1 && Number.isFinite(gamma) && gamma > 1,
    'Historical wake gas requires positive edge speed and valid compressible reference conditions.');
  const gm1 = gamma - 1, totalRatio = 1 + .5 * gm1 * mach * mach;
  const hstinv = gm1 * mach * mach / totalRatio, rhoTotal = totalRatio ** (1 / gm1);
  const h0 = 1 / hstinv, enthalpyRatio = 1 - .5 * ue * ue * hstinv;
  const enthalpy = h0 - .5 * ue * ue;
  require(positive(enthalpyRatio) && positive(enthalpy), 'Historical wake edge has nonpositive static enthalpy.');
  const density = rhoTotal * enthalpyRatio ** (1 / gm1), pressure = gm1 / gamma * density * enthalpy;
  const machSquared = ue * ue / (gm1 * enthalpy);
  require([density, pressure, machSquared].every(positive), 'Historical wake gas is nonfinite or nonpositive.');
  return { density, pressure, speed: ue, enthalpy, stagnationEnthalpy: h0, machSquared,
    gasSource: 'historical-common-isentrope' };
}

export function quadCoupledTransonicCoefficients({ checkpoint, flow, bl, bodies,
  solverLength, referenceChord, momentReference, alpha, mach, gamma } = {}) {
  const restart = checkpoint?.restart, input = restart?.input, options = restart?.options;
  require(checkpoint?.version === 1 && input && options && vector(restart.initialBL)
    && input.streamwiseMode === 'hybrid' && input.upwind
    && (input.flowModel === undefined || input.flowModel === 'compressible')
    && options.blThermodynamics === 'historical-common-isentrope' && options.edgeMatching === 'section-velocity',
  'Transonic coefficients require the explicit hybrid/historical-BL checkpoint.');
  const actualMach = input.mach, actualAlpha = input.alpha ?? 0, actualGamma = input.gamma ?? 1.4;
  require((mach === undefined || mach === actualMach) && (alpha === undefined || alpha === actualAlpha)
    && (gamma === undefined || gamma === actualGamma), 'Coefficient conditions disagree with the actual checkpoint state.');
  mach = actualMach; alpha = actualAlpha; gamma = actualGamma;
  require(positive(mach) && mach < 1 && Number.isFinite(alpha) && Number.isFinite(gamma) && gamma > 1
    && positive(solverLength) && positive(referenceChord), 'Invalid transonic coefficient reference conditions.');
  momentReference ??= { x: referenceChord / 4, y: 0 };
  require([momentReference?.x, momentReference?.y].every(Number.isFinite), 'Invalid transonic moment reference.');
  require(Array.isArray(bodies) && Array.isArray(input.bodies) && bodies.length > 0 && bodies.length === input.bodies.length
    && bodies.every((b, i) => b.leadingIndex === input.bodies[i].leadingIndex && b.trailingIndex === input.bodies[i].trailingIndex
      && b.element === input.bodies[i].element && Number.isInteger(b.element) && b.element >= 0 && b.element < bodies.length)
    && new Set(bodies.map(b => b.element)).size === bodies.length, 'Coefficient body mapping differs from the checkpoint.');
  require(Array.isArray(bl?.stations) && Array.isArray(bl?.wakes) && positive(options.reynolds)
    && bl.scale === 1 / Math.sqrt(options.reynolds) && restart.initialBL.length === 4 * bl.stations.length,
  'Missing complete packed BL metadata or inconsistent thickness normalization.');
  require(Array.isArray(input.outerLower) && Array.isArray(input.weights), 'Missing coefficient grid dimensions.');
  const layout = { nx: input.outerLower.length - 1, tubes: input.weights.map(w => w.length), bodies };
  const pInf = 1 / (gamma * mach * mach), h0 = 1 / ((gamma - 1) * mach * mach) + .5;
  const conditions = { flowModel: 'compressible', mach, alpha, gamma, pInf, h0 };
  const pressure = streamtubeSolidPressureForces({ flow, layout, conditions, referenceChord, momentReference });
  const angle = alpha * Math.PI / 180;
  const freestream = { density: 1, speed: 1, pressure: pInf, stagnationEnthalpy: h0,
    direction: { x: Math.cos(angle), y: Math.sin(angle) } };
  const warnings = [
    'Research transonic Euler/BL coefficients: force accuracy and grid refinement remain unvalidated.',
    'CL and Cm use physical Euler pressure on the solid contour; shear lift and shear moment are omitted.',
    'Viscous wake loss uses the historical common-isentrope BL gas. One merged wake per element approximates unequal-entropy Euler wake banks.',
  ];
  let eulerExit = null, viscousExit = null, eulerExitError = null, viscousExitError = null;
  const wakeInputs = [];
  try {
    require(flow.sections?.length === layout.nx && Array.isArray(flow.sections.at(-1))
      && flow.sections.at(-1).length === layout.tubes.length && Array.isArray(flow.allocation?.groups)
      && flow.allocation.groups.length === layout.tubes.length, 'Incomplete physical Euler exit.');
    const sections = [];
    flow.sections.at(-1).forEach((group, g) => {
      require(group.length === layout.tubes[g] && flow.allocation.groups[g]?.length === layout.tubes[g],
        'Incomplete physical Euler exit tube allocation.');
      group.forEach((s, j) => {
        const row0 = flow.nodes[g][layout.nx - 1], row1 = flow.nodes[g][layout.nx];
        const dx = .5 * ((row1[j].x - row0[j].x) + (row1[j + 1].x - row0[j + 1].x));
        const dy = .5 * ((row1[j].y - row0[j].y) + (row1[j + 1].y - row0[j + 1].y));
        // Saved massFlow is already rhoInfinity*UInfinity*physical length.
        // Never multiply it by solverLength a second time.
        sections.push({ density: s.rho, speed: s.q, pressure: s.p,
          stagnationEnthalpy: s.enthalpy + .5 * s.q * s.q,
          massFlow: flow.allocation.groups[g][j].massFlow, direction: { x: dx, y: dy } });
      });
    });
    eulerExit = streamtubeExitDefect({ sections, freestream, gamma, referenceChord });
  } catch (error) { eulerExitError = error.message; warnings.push(`Euler exit contribution unavailable: ${error.message}`); }
  try {
    require(bl.wakes.length === bodies.length && new Set(bl.wakes.map(w => w.body)).size === bodies.length,
      'Each element must have exactly one complete viscous wake.');
    for (const wake of bl.wakes) {
      require(Number.isInteger(wake.body) && wake.body >= 0 && wake.body < bodies.length
        && Array.isArray(wake.ids) && wake.ids.length > 0, 'Invalid viscous wake mapping.');
      const id = wake.ids.at(-1), station = bl.stations[id], packed = restart.initialBL;
      require(Number.isInteger(id) && id >= 0 && station?.body === wake.body && station.kind === 'wake'
        && station.i === layout.nx, 'Viscous wake must reach its final outlet station.');
      const gap = quadCoupledWakeGap({ body: bodies[wake.body], bodyIndex: wake.body, flow, station, solverLength });
      require(gap === 0, 'Wake exit has a nonzero prescribed dead-air gap; extend the wake before drag extrapolation.');
      const theta = packed[4 * id + 1] * bl.scale * solverLength;
      const deltaStar = packed[4 * id + 2] * bl.scale * solverLength;
      const gas = historicalCommonIsentropeWakeGas(packed[4 * id + 3], { mach, gamma });
      const rawHk = (deltaStar / theta - .29 * gas.machSquared) / (1 + .113 * gas.machSquared);
      require(positive(theta) && positive(deltaStar) && Number.isFinite(rawHk) && rawHk > 1,
        `Wake ${bodies[wake.body].element + 1} has inadmissible historical BL thickness or shape.`);
      wakeInputs.push({ ...gas, theta, deltaStar, rawHk, body: wake.body,
        element: bodies[wake.body].element, station: id });
    }
    viscousExit = streamtubeViscousExitDefect({ wakes: wakeInputs, freestream, gamma, referenceChord,
      gasSource: 'historical-common-isentrope' });
  } catch (error) { viscousExitError = error.message; warnings.push(`Viscous exit contribution unavailable: ${error.message}`); }
  const eulerWaveDragCoefficient = eulerExit?.eulerWaveDragCoefficient ?? null;
  const viscousDragCoefficient = viscousExit?.viscousDragCoefficient ?? null;
  const sum = eulerWaveDragCoefficient !== null && viscousDragCoefficient !== null
    ? eulerWaveDragCoefficient + viscousDragCoefficient : null;
  const cd = Number.isFinite(sum) ? sum : null;
  if (sum !== null && cd === null) warnings.push('Total exit drag is nonfinite; no partial sum is shown.');
  const thinGrid = { minNormalArea: null, maxDisplacementOverAdjacentTubeWidth: null };
  // Descriptive resolution diagnostic only; it is not a thin-layer or mesh
  // acceptance gate and does not change the coefficients.
  for (const row of flow.cells ?? []) for (const group of row) for (const cell of group)
    for (const width of cell.geometry?.normalAreas ?? []) if (positive(width))
      thinGrid.minNormalArea = Math.min(thinGrid.minNormalArea ?? width, width);
  for (const branch of bl.surfaces ?? []) for (const id of branch.ids ?? []) {
    const s = bl.stations[id], g = branch.side === 'upper' ? branch.body + 1 : branch.body;
    const row = flow.nodes[g]?.[s?.i]; if (!row || row.length < 2) continue;
    const j = branch.side === 'upper' ? 0 : row.length - 1, k = branch.side === 'upper' ? 1 : row.length - 2;
    const width = Math.hypot(row[j].x - row[k].x, row[j].y - row[k].y);
    const delta = restart.initialBL[4 * id + 2] * bl.scale * solverLength;
    if (positive(width) && positive(delta)) thinGrid.maxDisplacementOverAdjacentTubeWidth = Math.max(
      thinGrid.maxDisplacementOverAdjacentTubeWidth ?? 0, delta / width);
  }
  warnings.push('Pressure-integral drag is a separate diagnostic; it is not added to total exit drag. Thin-tube resolution and force accuracy remain unvalidated.');
  return { cl: pressure.cl, cm: pressure.cm, cx: pressure.cx, cy: pressure.cy, cd,
    pressureIntegralDrag: pressure.pressureIntegralDrag, eulerWaveDragCoefficient, viscousDragCoefficient,
    elements: pressure.perBody.map(p => ({ ...p })).sort((a, b) => a.element - b.element),
    wakes: viscousExit ? viscousExit.wakes.map((w, i) => ({ ...w, element: wakeInputs[i].element,
      station: wakeInputs[i].station })).sort((a, b) => a.element - b.element) : [],
    decomposition: { total: cd, viscousWake: viscousDragCoefficient, eulerEntropy: eulerWaveDragCoefficient,
      pressureIntegralSeparate: pressure.pressureIntegralDrag, complete: cd !== null,
      includesPressureIntegral: false, includesAdditionalSkinFriction: false },
    physicalPressure: pressure, eulerExit, viscousExit, wakeGas: wakeInputs,
    errors: { eulerExit: eulerExitError, viscousExit: viscousExitError }, thinGrid,
    conditions, referenceChord, solverLength, momentReference: { ...momentReference }, warnings,
    pressureKind: 'physical Euler interface pressure on the undisplaced solid contour',
    method: { lift: 'Physical Euler interface pressure on the solid contour; shear lift omitted',
      moment: 'Physical Euler interface pressure on the solid contour; positive nose-up; shear moment omitted',
      drag: 'MSES manual 2.7.2: summed historical-BL viscous wake defects plus all physical Euler exit entropy defects',
      pressure: 'Physical Euler interface pressure; stagnation uses the mean of its two physical banks' },
    physicalAcceptance: false };
}

export function quadCoupledTransonicCoefficientsFromResult(raw, { momentReference } = {}) {
  require(raw?.conditions?.blThermodynamics === 'historical-common-isentrope' && raw.checkpoint,
    'Missing explicit historical-mode result checkpoint.');
  const reynolds = raw.checkpoint.restart.options.reynolds;
  require(raw.kernelReynolds === reynolds, 'Result and checkpoint BL Reynolds normalizations differ.');
  return quadCoupledTransonicCoefficients({ checkpoint: raw.checkpoint, flow: raw.flow,
    bl: { ...raw.boundaryLayer, scale: 1 / Math.sqrt(reynolds) }, bodies: raw.solverInput.bodies,
    solverLength: raw.solverLength, referenceChord: raw.referenceChord, momentReference,
    alpha: raw.solverInput.alpha ?? 0, mach: raw.conditions.mach, gamma: raw.solverInput.gamma ?? 1.4 });
}
