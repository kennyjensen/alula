// SPDX-License-Identifier: GPL-2.0-or-later
// A shared panel field supplies passage order and inlet cuts for reentrant
// lower surfaces. Physical contour points and panel strengths are retained.
import { createInitialStreamtubeTopology } from '../geometry/streamtube-topology.js';
import { solveInviscid, velocityAt } from '../inviscid/linear-vortex.js';
import { streamfunctionAt } from '../inviscid/streamfunction.js';
import { streamtubeCellGeometry } from './streamtube-cell.js';

const seeds = new WeakMap();
const descriptor = input => JSON.stringify({ alpha: input.alpha ?? 0, bodies: input.bodies.map(b => ({
  points: b.points.map(p => [p.x, p.y]), trailingEdge: b.trailingEdge ?? null, element: b.element })) });

export function validatePanelTopologySeed(panel, input) {
  const record = seeds.get(panel);
  if (!record || record.descriptor !== descriptor(input) || record.panel !== JSON.stringify(panel))
    throw new Error('Prepared panel topology seed differs from its original geometry, incidence or field.');
  return panel;
}

function inletCut(field, x, lowY, highY, level, tolerance) {
  const sample = y => ({ y, psi: streamfunctionAt({ x, y }, field), u: velocityAt({ x, y }, field).u });
  let low = sample(lowY), high = sample(highY);
  if (!(low.u > 0 && high.u > 0 && low.psi < level && level < high.psi))
    throw new Error('Panel body streamfunction is outside the ordered physical inlet.');
  let point = sample(lowY + (highY - lowY) * (level - low.psi) / (high.psi - low.psi));
  for (let iteration = 0; iteration < 60; iteration++) {
    if (!(Number.isFinite(point.u) && point.u > 0)) throw new Error('Panel inlet cut has nonpositive streamwise velocity.');
    if (Math.abs(point.psi - level) <= tolerance) return { x, y: point.y, residual: point.psi - level, iterations: iteration };
    if (point.psi < level) low = point; else high = point;
    const next = point.y - (point.psi - level) / point.u;
    point = sample(next > low.y && next < high.y ? next : .5 * (low.y + high.y));
  }
  throw new Error('Panel inlet cut did not close its original streamfunction level.');
}

export function preparePanelStreamtubeTopology(input, controls, { crosslinePlacement = 'potential' } = {}) {
  let originalFailure;
  try { return { topology: createInitialStreamtubeTopology(input, controls) }; }
  catch (error) {
    if (crosslinePlacement !== 'potential' || error.code !== 'STREAMTUBE_NONMONOTONE_BODY' || error.diagnostics?.side !== 'lower') throw error;
    originalFailure = { code: error.code, message: error.message, ...error.diagnostics };
  }
  const topology = createInitialStreamtubeTopology(input, { ...controls, parametricLowerBranches: true });
  const oldBodies = topology.bodies, panel = solveInviscid({ elements: oldBodies.map(b => ({ points: b.points,
    ...(b.trailingEdge ? { trailingEdge: b.trailingEdge } : {}) })), alpha: topology.alpha, boundaryCondition: 'streamfunction' });
  if (panel.status !== 'solved' || panel.boundaryCondition !== 'streamfunction')
    throw new Error('Potential topology requires a solved common streamfunction panel field.');
  const levels = panel.diagnostics.surfaceStreamfunctions;
  if (!Array.isArray(levels) || levels.length !== oldBodies.length || !levels.every(Number.isFinite))
    throw new Error('Potential topology is missing finite body streamfunction levels.');
  const charts = oldBodies.map((body, element) => {
    const surface = panel.elements[element].cp.filter(p => !p.base), incoming = [];
    for (let i = 1; i < surface.length; i++) if (surface[i - 1].qt < 0 && surface[i].qt >= 0) incoming.push(i);
    if (incoming.length !== 1) throw new Error('Panel topology requires one resolved incoming stagnation point per body.');
    // This locates the incoming branch only. The initializer must still
    // certify its complete contour-potential coordinate before accepting a grid.
    return { element: body.element, incomingCollocationInterval: incoming[0], potentialChartValidated: false };
  });
  const order = oldBodies.map((_, i) => i).sort((a, b) => levels[a] - levels[b]);
  if (order.some((j, k) => k && !(levels[j] > levels[order[k - 1]]))) throw new Error('Coincident panel dividing-streamline levels need a different block topology.');
  const oldToNew = order.map(() => 0); order.forEach((old, current) => { oldToNew[old] = current; });
  const x = topology.outerLower[0].x, low = topology.outerLower[0].y, high = topology.outerUpper[0].y;
  const tolerance = 1e-12 * (high - low), cuts = order.map(old => inletCut(panel.field, x, low, high, levels[old], tolerance));
  if (cuts.some((p, i) => i && !(p.y > cuts[i - 1].y))) throw new Error('The panel field does not define ordered inlet cuts.');
  const lowerLevel = streamfunctionAt(topology.outerLower[0], panel.field), upperLevel = streamfunctionAt(topology.outerUpper[0], panel.field);
  topology.bodies = order.map(old => oldBodies[old]);
  topology.primaryBody = oldToNew[topology.primaryBody];
  // These paths supply provisional coordinates only. Tracing replaces them
  // before any gas inversion or physical mesh is accepted.
  topology.cutPaths = cuts.map(cut => topology.outerLower.map(p => ({ x: p.x, y: cut.y })));
  topology.captureLevels = [lowerLevel, ...order.map(old => levels[old]), upperLevel];
  if (topology.captureLevels.some((p, i, a) => i && !(p > a[i - 1]))) throw new Error('Panel topology has nonpositive captured inlet mass.');
  for (const key of ['surfaceIntervalsByElement', 'requestedSurfaceIntervalsByElement'])
    topology.gridSpacing[key] = topology.bodies.map(body => topology.gridSpacing[key].find(row => row.element === body.element));
  const remapPanel = p => ({ ...p, element: oldToNew[p.element] });
  const prepared = { ...panel, elements: order.map(old => panel.elements[old]),
    diagnostics: { ...panel.diagnostics, surfaceStreamfunctions: order.map(old => levels[old]),
      teProbes: panel.diagnostics.teProbes?.map(remapPanel) },
    field: { ...panel.field, panels: panel.field.panels.map(remapPanel),
      ...(panel.field.basePanels ? { basePanels: panel.field.basePanels.map(remapPanel) } : {}) } };
  const diagnostics = { method: 'Shared panel streamfunction passage order and physical inlet cuts', originalFailure,
    originalElementOrder: oldBodies.map(b => b.element), elementOrder: topology.bodies.map(b => b.element),
    bodyStreamfunctions: order.map(old => levels[old]), inletCuts: cuts, surfaceCharts: order.map(old => charts[old]),
    panelSolves: 1, panelStrengthsChanged: false, solidGeometryChanged: false };
  topology.gridSpacing.potentialTopology = diagnostics;
  seeds.set(prepared, { descriptor: descriptor(topology), panel: JSON.stringify(prepared) });
  return { topology, panelSolution: prepared, diagnostics };
}

// Uncorrected panel velocity seed on the final Euler geometry. This satisfies
// the eliminated mass and energy equations, not the Euler momentum/entropy rows.
function sections(system, state, visit) {
  const { nodes, allocation } = system.decode(state), { nx, tubes } = system.layout;
  for (let g = 0; g < tubes.length; g++) for (let j = 0; j < tubes[g]; j++) {
    for (let i = 1; i < nx; i++) {
      const geometry = streamtubeCellGeometry(
        [nodes[g][i-1][j], nodes[g][i][j], nodes[g][i+1][j]],
        [nodes[g][i-1][j+1], nodes[g][i][j+1], nodes[g][i+1][j+1]]);
      for (const k of i === 1 ? [0, 1] : [1]) {
        const row = i - 1 + k;
        const corners = [nodes[g][row][j], nodes[g][row+1][j], nodes[g][row][j+1], nodes[g][row+1][j+1]];
        const point = { x: corners.reduce((s,p)=>s+p.x,0)/4, y: corners.reduce((s,p)=>s+p.y,0)/4 };
        visit({ i: row, g, j, point, direction: geometry.directions[k], area: geometry.normalAreas[k],
          mass: allocation.groups[g][j].massFlow, column: system.layout.densityIndex(row,g,j) });
      }
    }
  }
}

export function sampleStreamtubePanelVelocity(system, state, velocityAt) {
  if (typeof velocityAt !== 'function') throw new Error('Panel velocity sampler is required.');
  const speeds = new Array(system.layout.densityCount);
  let maximumDirectionMismatch = 0, worstSection;
  sections(system,state,({ i,g,j,point,direction,column })=>{
    const v = velocityAt(point), q = v.u*direction.x + v.v*direction.y;
    if (![v.u,v.v,q].every(Number.isFinite) || q <= 0)
      throw new Error(`Panel velocity is not forward and finite at Euler section ${column}.`);
    speeds[column] = q;
    const mismatch = Math.abs(v.u*direction.y-v.v*direction.x)/Math.hypot(v.u,v.v);
    if(mismatch>maximumDirectionMismatch){maximumDirectionMismatch=mismatch;worstSection={i,group:g,tube:j,point,direction,velocity:v};}
  });
  return { method: 'uncorrected-panel-velocity', speeds, maximumDirectionMismatch, worstSection };
}

export function initializeStreamtubePanelVelocity(system, state, sample) {
  const { densityCount } = system.layout;
  if (system.conditions.flowModel !== 'compressible' || sample?.method !== 'uncorrected-panel-velocity'
    || !Array.isArray(sample.speeds) || sample.speeds.length !== densityCount
    || !sample.speeds.every(q=>Number.isFinite(q)&&q>0)) throw new Error('Invalid panel velocity seed.');
  const initial = Float64Array.from(state);
  sections(system,state,({ column,area,mass })=>{
    const rho = mass/(area*sample.speeds[column]);
    if (!(rho>0) || !Number.isFinite(rho)) throw new Error('Invalid panel-seeded Euler density.');
    initial[column] = Math.log(rho);
  });
  // No pressure clipping, branch substitution or geometry/mass alteration.
  const flow = system.evaluate(initial);
  let maximumSpeedError = 0, maximumEntropyDeparture = 0;
  const { gamma,pInf } = system.conditions;
  for (let i=0;i<system.layout.nx;i++) for(let g=0;g<system.layout.tubes.length;g++)
    for(let j=0;j<system.layout.tubes[g];j++) {
      const s=flow.sections[i][g][j], q=sample.speeds[system.layout.densityIndex(i,g,j)];
      maximumSpeedError=Math.max(maximumSpeedError,Math.abs(s.q-q));
      maximumEntropyDeparture=Math.max(maximumEntropyDeparture,Math.abs(Math.log(s.p/pInf)-gamma*Math.log(s.rho)));
    }
  return { initial,flow,diagnostics:{method:sample.method,initialGuessOnly:true,targetEquationsUnchanged:true,
    pgCorrection:false,maximumSpeedError,maximumDirectionMismatch:sample.maximumDirectionMismatch,
    maximumEntropyDeparture,residual:flow.diagnostics.residual,residualByFamily:flow.diagnostics.residualByFamily} };
}
