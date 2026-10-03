// SPDX-License-Identifier: GPL-2.0-or-later
// Explicit local XFOIL MRCHDU preparation; not a simultaneous Euler/BL solve.
// Original xbl.f:880–1204 uses current Ue/Hk, changes all four BL variables,
// and checks transition during its mixed-mode station iterations. SETBL
// invokes it before assembling its global system (xbl.f:93).
import { mrchdu, blpini } from './xfoil/xbl.js';
import { ensureCtx, blprv, blkin, blvar, blmid, trchek, tesys, blsys, hkin, syncComToVars } from './xfoil/xblsys.js';
import { createIntegralKernel } from './integral.js';

export function relaxXfoilProfiles({ surfaces, wake = [], phases, tripS, normalGap = 0 }, parameters = {}) {
  const kernel = createIntegralKernel(parameters);
  // Native SETBL calls MRCHDU on the old physical profile after changing gas
  // parameters. BLVAR and the local DSLIM update can repair raw Hk <= 1 in
  // that starting guess. This exception is INPUT-only: retain every strict
  // output check below, as native extrapolation need not return usable states.
  const inputShapeRecoveries = [];
  const inputStation = (station, regime, location) => {
    try { return kernel.station(station, regime); }
    catch (error) {
      const d = error?.diagnostics;
      if (error?.code !== 'BL_EDGE_STATE_DOMAIN' || d?.condition !== 'raw-hk'
        || !Array.isArray(d.failedChecks) || d.failedChecks.length !== 1 || d.failedChecks[0] !== 'raw-hk'
        || !Number.isFinite(d.rawHk) || d.rawHk > 1
        || ![d.enthalpyRatio, d.density, d.viscosity, d.reTheta, d.machSquared, d.physicalUe].every(Number.isFinite)
        || !(d.enthalpyRatio > 0 && d.density > 0 && d.viscosity > 0 && d.reTheta > 0
          && d.machSquared >= 0 && d.physicalUe > 0)) throw error;
      inputShapeRecoveries.push({ ...location, rawHk: d.rawHk, fluidShapeFactor: d.fluidShapeFactor,
        machSquared: d.machSquared, enthalpyRatio: d.enthalpyRatio });
    }
  };
  const { reynolds, mach, gamma, ncrit, velocityConvention } = kernel.parameters;
  if (velocityConvention !== 'physical') throw new Error('MRCHDU requires physical edge velocities.');
  const ordered = p => Array.isArray(p) && p.length >= 2 && p.every((a, i) =>
    ['s', 'ue', 'aux', 'theta', 'deltaStar'].every(k => Number.isFinite(a[k])) && a.s > 0 && a.ue > 0
    && a.theta > 0 && a.deltaStar - (a.wakeGap ?? 0) > a.theta && (!i || a.s > p[i - 1].s));
  if (!Array.isArray(surfaces) || surfaces.length !== 2 || !surfaces.every(ordered)
    || !Array.isArray(wake) || (wake.length && !ordered(wake)))
    throw new Error('Supply complete ordered physical surface and wake states.');
  if (!Number.isFinite(normalGap) || normalGap < 0 || !Array.isArray(phases) || phases.length !== 2
    || phases.some((p, side) => !Number.isInteger(p) || p < 1 || p >= surfaces[side].length))
    throw new Error('Invalid MRCHDU base gap or previous transition intervals.');
  // XIFSET's geometry-free branch places a terminal trip at the TE. An
  // earlier physical trip requires a separately verified geometry adapter.
  const terminalTrips = surfaces.map(p => p.at(-1).s);
  if (tripS !== undefined && (!Array.isArray(tripS) || tripS.length !== 2 || tripS.some((s, i) => s !== terminalTrips[i])))
    throw new Error('Geometry-free MRCHDU currently supports terminal trips only.');
  surfaces.forEach((profile, side) => profile.forEach((p, i) => {
    if ((p.wakeGap ?? 0) !== 0) throw new Error('A solid surface cannot contain a wake gap.');
    // UPDATE can leave finite laminar N below zero or above Ncrit: its
    // additive N correction is not clamped (xbl.f:1433,1514). MRCHDU retains
    // ITROLD and reintegrates N through TRCHEK (xbl.f:918,941,972). Admit
    // that incoming iterate, retaining positive turbulent shear and the
    // strict returned N/phase checks below. ordered() requires finite N.
    if (i >= phases[side] && p.aux <= 0) throw new Error('MRCHDU auxiliary state does not match its previous phase.');
    inputStation(p, i < phases[side] ? 'laminar' : 'turbulent', { part: 'surface', side: side === 0 ? 'upper' : 'lower', index: i });
  }));
  if (wake.length) {
    if (wake[0].s !== terminalTrips[1]) throw new Error('The first wake station must share the lower TE arclength.');
    wake.forEach((p, i) => {
      if (!Number.isFinite(p.wakeGap ?? 0) || (p.wakeGap ?? 0) < 0 || !(p.aux > 0)) throw new Error('Invalid MRCHDU wake state.');
      inputStation(p, 'wake', { part: 'wake', index: i });
    });
    if ((wake[0].wakeGap ?? 0) !== normalGap) throw new Error('The first wake gap must equal the physical base gap.');
  } else if (normalGap !== 0) throw new Error('A finite base requires the merged wake states.');

  const gm1 = gamma - 1, totalRatio = 1 + .5 * gm1 * mach ** 2;
  const hstinv = gm1 * mach ** 2 / totalRatio, hstinvMs = gm1 / totalRatio ** 2;
  const herat = 1 - .5 * hstinv, heratMs = -.5 * hstinvMs, hvrat = .35;
  const viscousFactor = herat ** 1.5 * (1 + hvrat) / (herat + hvrat), messages = [];
  const ctx = { QINFBL: 1, TKBL: 0, TKBL_MS: 0,
    RSTBL: totalRatio ** (1 / gm1), RSTBL_MS: .5 * totalRatio ** (1 / gm1 - 1),
    HSTINV: hstinv, HSTINV_MS: hstinvMs, REYBL: reynolds * viscousFactor, REYBL_RE: viscousFactor,
    REYBL_MS: reynolds * viscousFactor * (1.5 / herat - 1 / (herat + hvrat)) * heratMs,
    GAMBL: gamma, GM1BL: gm1, HVRAT: hvrat, AMCRIT: ncrit, BULE: 1, IDAMPV: 0, ANTE: normalGap,
    blprv, blkin, blvar, blmid, trchek, tesys, blsys, hkin, syncComToVars, log: m => messages.push(m) };
  const profiles = [surfaces[0], [...surfaces[1], ...wake]], size = Math.max(...profiles.map(p => p.length)) + 2;
  for (const name of ['XSSI', 'UEDG', 'CTAU', 'THET', 'DSTR', 'MASS', 'TAU', 'DIS', 'CTQ', 'DELT', 'TSTR'])
    ctx[name] = Array.from({ length: size }, () => new Float64Array(3));
  for (const name of ['NBL', 'IBLTE', 'ITRAN', 'XSSITR', 'TFORCE']) ctx[name] = new Float64Array(3);
  ctx.ACRIT = [0, ncrit, ncrit]; ctx.XSTRIP = [0, 1, 1]; ctx.WGAP = Float64Array.from([0, ...wake.map(p => p.wakeGap ?? 0)]);
  profiles.forEach((p, side) => {
    const is = side + 1; ctx.NBL[is] = p.length + 1; ctx.IBLTE[is] = surfaces[side].length + 1; ctx.ITRAN[is] = phases[side] + 2;
    p.forEach((a, i) => {
      ctx.XSSI[i + 2][is] = a.s; ctx.UEDG[i + 2][is] = a.ue; ctx.CTAU[i + 2][is] = a.aux;
      ctx.THET[i + 2][is] = a.theta; ctx.DSTR[i + 2][is] = a.deltaStar; ctx.MASS[i + 2][is] = a.deltaStar * a.ue;
    });
  });
  ensureCtx(ctx); blpini(ctx); mrchdu(ctx);
  const read = (p, side) => p.map((a, i) => ({ s: a.s, ue: ctx.UEDG[i + 2][side + 1], aux: ctx.CTAU[i + 2][side + 1],
    theta: ctx.THET[i + 2][side + 1], deltaStar: ctx.DSTR[i + 2][side + 1], ...(a.wakeGap === undefined ? {} : { wakeGap: a.wakeGap }) }));
  const states = profiles.map(read);
  const result = { surfaces: surfaces.map((p, side) => ({ states: states[side].slice(0, p.length),
    transition: ctx.ITRAN[side + 1] - 2, s: ctx.XSSITR[side + 1], forced: Boolean(ctx.TFORCE[side + 1]) })),
    wake: states[1].slice(surfaces[1].length), normalGap, tripS: terminalTrips, parameters: { ...kernel.parameters },
    messages, localConvergenceWarnings: messages.filter(m => m.includes('Convergence failed')),
    method: 'XFOIL MRCHDU current-profile preparation', flowSolved: false,
    ...(inputShapeRecoveries.length ? { inputShapeRecovery: {
      method: 'native-mrchdu-local-shape', count: inputShapeRecoveries.length,
      minimumRawHk: Math.min(...inputShapeRecoveries.map(s => s.rawHk)), stations: inputShapeRecoveries,
      interpretation: 'Starting-profile shape only; unchanged native MRCHDU and strict returned-state checks.' } } : {}),
    changes: { surfaces: [], wake: [] } };
  const differences = (before, after) => after.map((a, i) => ({ s: a.s,
    ...Object.fromEntries(['ue', 'aux', 'theta', 'deltaStar'].map(k => [k, a[k] - before[i][k]])) }));
  result.changes.surfaces = result.surfaces.map((s, i) => differences(surfaces[i], s.states));
  result.changes.wake = differences(wake, result.wake);
  try {
    result.surfaces.forEach(s => {
      if (!ordered(s.states) || !Number.isInteger(s.transition) || s.transition < 1 || s.transition >= s.states.length)
        throw new Error('MRCHDU returned an unusable surface or transition.');
      s.states.forEach((p, i) => {
        if (i < s.transition ? p.aux < 0 || p.aux >= ncrit : !(p.aux > 0))
          throw new Error('MRCHDU returned auxiliary state inconsistent with its new phase.');
        kernel.station(p, i < s.transition ? 'laminar' : 'turbulent');
      });
    });
    if (wake.length && !ordered(result.wake)) throw new Error('MRCHDU returned an unusable wake.');
    result.wake.forEach(p => {
      if (!(p.aux > 0)) throw new Error('MRCHDU returned nonpositive wake shear.');
      kernel.station(p, 'wake');
    });
  } catch (error) { error.diagnostics = result; throw error; }
  return result;
}
