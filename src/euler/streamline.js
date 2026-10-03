// SPDX-License-Identifier: GPL-2.0-or-later
import { freestream } from './gas.js';
import { buildMesh } from './mesh.js';
import { eulerResidual } from './residual.js';
import { encodeStates, decodeStates, summarizeEuler } from './solve.js';
import { solveNewton } from '../numerics/newton.js';

// A square coupled reference system: 4C conservative flow unknowns/equations
// plus one y displacement and one zero shared-face mass-flux equation for
// every internal streamline segment. Inlet streamline positions are fixed.
// This channel experiment verifies flow/grid off-diagonal Newton blocks. Body
// dividing-streamline capture and stagnation position are NOT closed here.
export function createStreamlineSystem(mesh, { mach = 0.3, alpha = 0, gamma = 1.4, initial } = {}) {
  if (!mesh.structured || !mesh.streamlineNodes?.length) throw new Error('Moving-streamline reference requires an unmapped channel with internal rows.');
  if (mesh.cells.length > 200) throw new Error('Moving-streamline dense reference is limited to 200 cells.');
  const reference = freestream({ mach, alpha, gamma }); const count = mesh.cells.length;
  const nodes = mesh.streamlineNodes;
  const connectivity = mesh.cells.map(c => c.vertices);
  const boundaries = mesh.faces.filter(f => f.boundary).map(f => ({ ...f.boundary, a: f.a, b: f.b }));
  const faceIds = nodes.map(node => mesh.faces.findIndex(f => (f.a === node.previous && f.b === node.vertex) || (f.b === node.previous && f.a === node.vertex)));
  if (faceIds.some(i => i < 0 || mesh.faces[i].neighbor === null)) throw new Error('Missing shared streamline face.');
  const deform = x => {
    const vertices = mesh.vertices.map(v => ({ ...v }));
    nodes.forEach((node, i) => { vertices[node.vertex].y += node.scale * x[4 * count + i]; });
    return buildMesh(vertices, connectivity, boundaries);
  };
  const admissible = x => {
    if (!decodeStates(x, count, reference).every(s => Number.isFinite(s.rho) && s.rho > 0 && s.p > 0)) return false;
    try { deform(x); return true; } catch { return false; }
  };
  const residual = x => {
    const current = deform(x); const evaluated = eulerResidual(current, decodeStates(x, count, reference), reference);
    return Float64Array.from([...evaluated.residual, ...faceIds.map(id => evaluated.faceFluxes[id][0] / (current.faces[id].length * reference.rho * Math.hypot(reference.u, reference.v)))]);
  };
  const guess = new Float64Array(4 * count + nodes.length);
  if (initial && initial.length !== count) throw new Error('Invalid Euler initial-state size.');
  guess.set(encodeStates(initial ?? mesh.cells.map(() => reference), reference));
  return { initial: guess, admissible, residual, meshAt: deform, reference, faceIds,
    unknowns: { flow: 4 * count, grid: nodes.length }, statesAt: x => decodeStates(x, count, reference) };
}

export function solveStreamlineChannel(mesh, options = {}) {
  const { tolerance = 1e-9, maxIterations = 30, onIteration, onMesh } = options;
  const system = createStreamlineSystem(mesh, options);
  const { reference, faceIds } = system;
  let previous = mesh;
  const solve = solveNewton({ initial: system.initial, admissible: system.admissible, residual: system.residual, tolerance, maxIterations, onIteration,
    onState: onMesh ? ({ x, iteration }) => {
      const current = system.meshAt(x);
      const maximumNodeMovement = current.vertices.reduce((maximum, p, i) => Math.max(maximum,
        Math.hypot(p.x - previous.vertices[i].x, p.y - previous.vertices[i].y)), 0);
      previous = current;
      onMesh({ mesh: current, iteration: { ...iteration, maximumNodeMovement } });
    } : undefined });
  const current = system.meshAt(solve.x); const states = system.statesAt(solve.x);
  const result = summarizeEuler(current, states, reference, solve);
  result.model = 'euler-moving-streamline-channel';
  result.displacements = [...solve.x.slice(system.unknowns.flow)].map((v, i) => v * mesh.streamlineNodes[i].scale);
  result.streamlineFaces = faceIds;
  result.diagnostics.streamlineMassResidual = Math.max(...faceIds.map(id => Math.abs(result.faceFluxes[id][0]) / current.faces[id].length));
  return result;
}
