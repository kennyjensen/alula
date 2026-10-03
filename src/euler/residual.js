// SPDX-License-Identifier: GPL-2.0-or-later
import { hllc, totalConditions } from './gas.js';
import { normInf } from '../numerics/linear.js';

// Subsonic reservoir inlet: specify entropy, total enthalpy and direction;
// retain the outgoing acoustic invariant from the interior. The inlet normal
// need not coincide with the prescribed velocity direction.
export function boundaryState(interior, face, reference, gamma = 1.4) {
  const { nx, ny, boundary } = face;
  const un = interior.u * nx + interior.v * ny;
  if (boundary.type === 'wall') return { ...interior, u: interior.u - 2 * un * nx, v: interior.v - 2 * un * ny };
  if (boundary.type === 'farfield') return boundary.state ?? reference;
  const a = Math.sqrt(gamma * interior.p / interior.rho);
  if (Math.abs(un) >= a) throw new Error('Subsonic inlet/outlet boundary became sonic.');
  if (boundary.type === 'outlet') {
    const p = boundary.pressure ?? reference.p;
    const rho = interior.rho * (p / interior.p) ** (1 / gamma);
    const ab = Math.sqrt(gamma * p / rho);
    const ub = un + 2 * (a - ab) / (gamma - 1);
    if (ub <= 0 || ub >= ab) throw new Error('Outlet reversal or choking is outside the reference boundary model.');
    return { rho, p, u: interior.u + (ub - un) * nx, v: interior.v + (ub - un) * ny };
  }
  if (boundary.type === 'inlet') {
    const total = totalConditions(reference, gamma);
    const speed = Math.hypot(reference.u, reference.v);
    const dx = reference.u / speed; const dy = reference.v / speed;
    const d = dx * nx + dy * ny;
    if (d >= 0) throw new Error('Reservoir inlet direction must enter the domain.');
    const invariant = un + 2 * a / (gamma - 1);
    const sound = q => Math.sqrt((gamma - 1) * (total.h0 - q * q / 2));
    const equation = q => d * q + 2 * sound(q) / (gamma - 1) - invariant;
    // Monotonic equation on the subsonic branch; no quadratic root ambiguity.
    let lo = 0; let hi = Math.sqrt(2 * (gamma - 1) * total.h0 / (gamma + 1));
    if (equation(lo) < 0 || equation(hi) > 0) throw new Error('Inlet reversal or choking is outside the reference boundary model.');
    for (let k = 0; k < 52; k++) { const mid = (lo + hi) / 2; if (equation(mid) > 0) lo = mid; else hi = mid; }
    const q = (lo + hi) / 2; const ab = sound(q);
    const rho = (ab * ab / (gamma * Math.exp(total.entropy))) ** (1 / (gamma - 1));
    return { rho, p: rho * ab * ab / gamma, u: q * dx, v: q * dy };
  }
  throw new Error(`Unknown Euler boundary type: ${boundary.type}.`);
}

export function integratedFaceFlux(left, right, face, reference, gamma = reference.gamma ?? 1.4) {
  return hllc(left, right ?? boundaryState(left, face, reference, gamma), face.nx, face.ny, gamma).map(value => value * face.length);
}

export function eulerResidual(mesh, states, reference, { gamma = reference.gamma ?? 1.4,reconstruction } = {}) {
  if (states.length !== mesh.cells.length) throw new Error('Euler state/mesh size mismatch.');
  const integrated = new Float64Array(4 * states.length);
  const boundaryFlux = new Float64Array(4); const faceFluxes = [];
  const wallForces = new Map(); let boundaryMassMagnitude = 0; let wallLeakage = 0;
  for (const [index,face] of mesh.faces.entries()) {
    const left = reconstruction?reconstruction.sample(reconstruction.faces[index].left,states):states[face.owner];
    const right = reconstruction?reconstruction.sample(reconstruction.faces[index].right,states):face.neighbor === null ? null : states[face.neighbor];
    const flux = integratedFaceFlux(left, right, face, reference, gamma);
    faceFluxes.push(flux);
    for (let k = 0; k < 4; k++) {
      integrated[4 * face.owner + k] += flux[k];
      if (face.neighbor !== null) integrated[4 * face.neighbor + k] -= flux[k];
      else boundaryFlux[k] += flux[k];
    }
    if (face.neighbor === null) {
      boundaryMassMagnitude += Math.abs(flux[0]);
      if (face.boundary.type === 'wall') {
        wallLeakage = Math.max(wallLeakage, Math.abs(flux[0]), Math.abs(flux[3]));
        // Force on solid = flux leaving the fluid; remove uniform reference
        // pressure locally to avoid cancellation in integrated closed contours.
        const key = face.boundary.element ?? face.boundary.side ?? 'wall';
        const force = wallForces.get(key) ?? { x: 0, y: 0 };
        force.x += flux[1] - reference.p * face.nx * face.length;
        force.y += flux[2] - reference.p * face.ny * face.length;
        wallForces.set(key, force);
      }
    }
  }
  const speed = Math.hypot(reference.u, reference.v); const mass = reference.rho * speed;
  const scales = [mass, mass * speed, mass * speed, mass * totalConditions(reference, gamma).h0];
  const residual = integrated.map((v, i) => v / (mesh.cells[Math.floor(i / 4)].perimeter * scales[i % 4]));
  const summed = Array.from({ length: 4 }, (_, k) => mesh.cells.reduce((sum, _, i) => sum + integrated[4 * i + k], 0));
  return { residual, integrated, faceFluxes, boundaryFlux: [...boundaryFlux], wallForces: Object.fromEntries(wallForces),
    diagnostics: { residual: normInf(residual), wallLeakage,
      relativeMassImbalance: Math.abs(boundaryFlux[0]) / Math.max(boundaryMassMagnitude / 2, Number.MIN_VALUE),
      sharedFluxCancellation: Math.max(...summed.map((v, k) => Math.abs(v - boundaryFlux[k]) / scales[k])) } };
}
