// SPDX-License-Identifier: GPL-2.0-or-later
import { createInitialStreamtubeTopology } from '../geometry/streamtube-topology.js';
import { createPanelStreamtubeGrid } from './streamtube-body-initializer.js';
import { streamtubeMeshConnectivity, validateStreamtubeMeshConnectivity } from '../geometry/streamtube-mesh-connectivity.js';
import { repairSubsonicStreamtubeGrid, repairStreamtubeGrid } from './streamtube-grid-repair.js';
import { relaxStreamtubeInitialGrid } from './streamtube-elliptic-initializer.js';
import { streamtubeFlowSnapshot } from './streamtube-flow-preview.js';
import { requireStreamtubeInitialGridDomain } from './streamtube-grid-repair.js';

export const STREAMTUBE_MINIMUM_CORNER_SINE = 1e-12;

// Scan in published cell order with the same arithmetic and corner threshold.
// Reuse four point references instead of allocating arrays for every cell.
function logicalMeshQuality(layout, nodes) {
  let minArea = Infinity, minCornerSine = Infinity; const invalidCells = [];
  const p = new Array(4); let index = 0;
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < layout.nx; i++) for (let tube = 0; tube < layout.tubes[g]; tube++, index++) {
    p[0] = nodes[g][i][tube]; p[1] = nodes[g][i + 1][tube];
    p[2] = nodes[g][i + 1][tube + 1]; p[3] = nodes[g][i][tube + 1];
    let twiceArea = 0, valid = true;
    for (let j = 0; j < 4; j++) {
      const a = p[j], b = p[(j + 1) % 4], c = p[(j + 2) % 4];
      twiceArea += (a.x - p[0].x) * (b.y - p[0].y) - (a.y - p[0].y) * (b.x - p[0].x);
      const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
      const sine = (ux * vy - uy * vx) / (Math.hypot(ux, uy) * Math.hypot(vx, vy));
      minCornerSine = Math.min(minCornerSine, sine); if (!(sine > STREAMTUBE_MINIMUM_CORNER_SINE)) valid = false;
    }
    minArea = Math.min(minArea, .5 * twiceArea); if (!valid || !(twiceArea > 0)) invalidCells.push(index);
  }
  return { valid: invalidCells.length === 0, invalidCells, minArea, minCornerSine };
}

export function streamtubeMeshQuality({ system, nodes }) {
  validateStreamtubeMeshConnectivity(system.layout, nodes);
  return logicalMeshQuality(system.layout, nodes);
}

// Actual quadrilateral connectivity from the intrinsic solver's logical
// grid. Shared cuts and the two LE/TE endpoints use shared vertex indices.
export function streamtubeMeshSnapshot({ system, nodes, diagnostics, initial, iteration, flow, flowSolved = false }) {
  const { layout } = system, { vertices, cells } = streamtubeMeshConnectivity(layout, nodes);
  return { topology: 'intrinsic-quadrilateral-streamtubes', vertices, cells,
    ...(flow ? { flow: streamtubeFlowSnapshot(flow, iteration?.iteration ?? 0) } : {}),
    ...(iteration ? { iteration: { ...iteration } } : {}),
    quality: logicalMeshQuality(layout, nodes),
    initialization: { ...diagnostics, groups: layout.elements + 1, streamwiseSegments: layout.nx, tubes: layout.tubes,
      surfaceIntervals: layout.bodies.map(body => ({ element: body.element, intervals: body.trailingIndex - body.leadingIndex })),
      ...(initial ? { massFlows: system.decode(initial).allocation.groups.map(group => group.map(t => t.massFlow)) } : {}),
      elementOrder: layout.bodies.map(b => b.element), flowSolved } };
}

export function buildStreamtubeMeshPreview(input, controls = {}) {
  const topology = createInitialStreamtubeTopology(input, controls);
  return prepareStreamtubeMesh(createPanelStreamtubeGrid(topology, controls), controls).mesh;
}

// Optional elliptic SLOR, otherwise a bounded fallback for invalid traced
// grids. Publishing precedes relaxation and gas inversion, even on failure.
export function prepareStreamtubeMesh(prepared, { onMesh, gridRepair = {}, ellipticSmoothing = false } = {}) {
  let mesh = streamtubeMeshSnapshot(prepared);
  onMesh?.(mesh, 'initial');
  if (ellipticSmoothing) {
    const { system, guideField } = prepared, preview = onMesh ? structuredClone(prepared.nodes) : null;
    let nodes = prepared.nodes, initial = prepared.initial, diagnostics = prepared.diagnostics, report;
    try {
      const r = relaxStreamtubeInitialGrid(prepared, { omega: 1.3, maxSweeps: 600,
        farfieldBoundary: 'fixed', allowPartialInitialGuess: true,
        ...(typeof ellipticSmoothing === 'object' ? ellipticSmoothing : {}),
        seed: ellipticSmoothing.boundaryControl === 'wall-angle' ? 'supplied' : 'linear',
        onSweep: onMesh ? (h, region) => {
          preview[h.region] = region;
          if (onMesh && h.iteration % 20 === 0) onMesh(streamtubeMeshSnapshot({ system, nodes: preview,
            diagnostics: { ...diagnostics, gridSmoothing: { attempted: true, converged: false, current: h } } }), 'smoothing');
        } : undefined });
      nodes = r.partialInitialGuess?.nodes ?? r.nodes;
      report = { attempted: true, method: 'elliptic SLOR', converged: false, regions: r.regions,
        boundaryControl: r.boundaryControl, stationCoordinate: r.stationCoordinate, farfieldBoundary: r.farfieldBoundary,
        ...(r.regions.some(region => region.spacingFallback?.attempted) ? {
          fixedSpacingPassages: r.regions.flatMap((region, g) => region.spacingFallback?.converged ? [g] : []),
        } : {}),
        ...(r.regions.some(region => region.harmonicFallback?.attempted) ? {
          harmonicPassages: r.regions.flatMap((region, g) => region.harmonicFallback?.converged ? [g] : []),
        } : {}),
        ...(r.farfieldBoundary === 'normal-curve' ? {
          farfieldGeometry: 'Fixed C1 Hermite graph through the original traced stations with panel-velocity tangents; a between-node approximation to the panel streamline.',
          exactModernMsetBoundaryRule: false,
        } : {}) };
      if (!r.converged) {
        report.reason = r.regions.filter(q => !q.converged).map(q => q.reason).join('; ');
        if (!r.partialInitialGuess) throw new Error(report.reason);
        report.partialPassages = r.partialInitialGuess.passages;
        report.partialInitialGuess = r.partialInitialGuess.reports;
        report.geometryAcceptance = requireStreamtubeInitialGridDomain(nodes);
      }
      if (!streamtubeMeshSnapshot({ system, nodes }).quality.valid) throw new Error('Smoothed grid failed the independent corner check.');
      if (!guideField?.admissibleNode) throw new Error('Body-clearance checks are required for elliptic initialization.');
      for (let g = 0; g < nodes.length; g++) for (let i = 0; i <= system.layout.nx; i++) for (let j = 1; j < nodes[g][i].length - 1; j++)
        if (!guideField.admissibleNode(nodes[g][i][j])) throw new Error(`Smoothed interior node enters a body at region ${g + 1}, station ${i}, streamline ${j}.`);
      diagnostics = guideField.diagnosticsForNodes(nodes);
      if (r.farfieldBoundary === 'normal-curve' && typeof guideField.streamfunctionAt === 'function') {
        report.farfieldStreamfunctionDrift = [0, nodes.length - 1].map((g, side) => {
          const level = prepared.input.captureLevels[side ? prepared.input.captureLevels.length - 1 : 0];
          return Math.max(...nodes[g].map(row => Math.abs(guideField.streamfunctionAt(row[side ? row.length - 1 : 0]) - level)));
        });
      }
      initial = system.adoptGeometry(prepared.initial, nodes); report.converged = r.converged;
      if (r.partialInitialGuess) report.initialGuessAccepted = true;
    } catch (error) {
      if (error?.code === 'slor-observer-failed') throw error;
      // A rejected smoother must not replace a valid user-visible base grid.
      // Keep its candidate separately in the exported diagnostics.
      const rejectedMesh = streamtubeMeshSnapshot({ system, nodes });
      report = { ...(report ?? { attempted: true, method: 'elliptic SLOR', converged: false }),
        converged: false, reason: error.message, retainedOriginal: true, rejectedMesh };
      delete report.initialGuessAccepted;
      nodes = prepared.nodes; initial = prepared.initial; diagnostics = prepared.diagnostics;
    }
    const updated = { ...prepared, nodes, initial, diagnostics: { ...diagnostics, gridSmoothing: report } };
    mesh = streamtubeMeshSnapshot(updated); onMesh?.(mesh, 'initial');
    return { ...updated, mesh };
  }
  if (mesh.quality.valid || gridRepair === false) return { ...prepared, mesh };
  const initialQuality = structuredClone(mesh.quality), { system, guideField } = prepared;
  let report, updated = prepared;
  try {
    if (!guideField) throw new Error('Panel guide geometry checks are required for initial-grid repair.');
    const repair = system.conditions.flowModel === 'compressible' ? repairSubsonicStreamtubeGrid : repairStreamtubeGrid;
    const r = repair(prepared, { minimumCorner: .2, maxSweeps: 60, ...gridRepair });
    report = { attempted: true, converged: false, reason: r.reason, initialQuality, quality: r.quality,
      maxDisplacement: r.maxDisplacement, sweeps: r.history.length - 1,
      mach: system.conditions.mach, massFluxFraction: r.massFluxFraction, massFluxLimit: r.massFluxLimit };
    if (r.converged) {
      for (let g = 0; g < r.nodes.length; g++) for (let i = 0; i < r.nodes[g].length; i++) for (let j = 0; j < r.nodes[g][i].length; j++) {
        const a = prepared.nodes[g][i][j], b = r.nodes[g][i][j];
        if ((a.x !== b.x || a.y !== b.y) && !guideField.admissibleMove(a, b))
          throw new Error(`Grid repair would cross a body at region ${g + 1}, station ${i}, streamline ${j}.`);
      }
      const candidate = streamtubeMeshSnapshot({ system, nodes: r.nodes });
      if (!candidate.quality.valid) throw new Error('Repaired initial grid failed the independent mesh snapshot check.');
      const diagnostics = guideField.diagnosticsForNodes(r.nodes);
      // All rejection checks precede adoption. No partial repair changes
      // the input chart, capture masses, walls or inlet/outlet constraints.
      const initial = system.adoptGeometry(prepared.initial, r.nodes);
      report.converged = true;
      updated = { ...prepared, initial, nodes: r.nodes, diagnostics,
        status: 'repaired initial geometry; flow equations not solved' };
    }
  } catch (error) {
    report = { ...(report ?? { attempted: true, initialQuality }), converged: false, reason: error.message };
  }
  updated = { ...updated, diagnostics: { ...updated.diagnostics, gridRepair: report } };
  mesh = streamtubeMeshSnapshot(updated);
  onMesh?.(mesh, 'initial');
  return { ...updated, mesh };
}
