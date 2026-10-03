// SPDX-License-Identifier: GPL-2.0-or-later
// Research adapter for original MRCHDU transition re-establishment.
// Takes two complete surface states and their old transition indices.
// No wake or outer-flow correction; not used by the production quad solver.
import { mrchdu, blpini } from '../../src/viscous/xfoil/xbl.js';
import { ensureCtx, blprv, blkin, blvar, blmid, trchek, tesys, blsys, hkin, syncComToVars } from '../../src/viscous/xfoil/xblsys.js';
import { createIntegralKernel } from '../../src/viscous/integral.js';

export function relaxXfoilSurfaces(profiles, parameters = {}, phases) {
  const { reynolds, mach, gamma, ncrit, velocityConvention } = createIntegralKernel(parameters).parameters;
  if (velocityConvention !== 'physical') throw new Error('MRCHDU surface initialization requires physical edge velocities.');
  if (!Array.isArray(profiles) || profiles.length !== 2 || profiles.some(profile =>
    !Array.isArray(profile) || profile.length < 2 || profile.some((p, i) =>
      !Number.isFinite(p.s) || !Number.isFinite(p.ue) || !(p.s > 0) || !(p.ue > 0) || (i && p.s <= profile[i - 1].s))))
    throw new Error('Supply two ordered surface profiles with positive physical edge speed.');
  const gm1 = gamma - 1, totalRatio = 1 + .5 * gm1 * mach ** 2;
  const hstinv = gm1 * mach ** 2 / totalRatio, hstinvMs = gm1 / totalRatio ** 2;
  const herat = 1 - .5 * hstinv, heratMs = -.5 * hstinvMs, hvrat = .35;
  const viscousFactor = herat ** 1.5 * (1 + hvrat) / (herat + hvrat);
  const messages = [], ctx = { QINFBL: 1, TKBL: 0, TKBL_MS: 0,
    RSTBL: totalRatio ** (1 / gm1), RSTBL_MS: .5 * totalRatio ** (1 / gm1 - 1),
    HSTINV: hstinv, HSTINV_MS: hstinvMs, REYBL: reynolds * viscousFactor, REYBL_RE: viscousFactor,
    REYBL_MS: reynolds * viscousFactor * (1.5 / herat - 1 / (herat + hvrat)) * heratMs,
    GAMBL: gamma, GM1BL: gm1, HVRAT: hvrat, AMCRIT: ncrit, BULE: 1, IDAMPV: 0, ANTE: 0,
    blprv, blkin, blvar, blmid, trchek, tesys, blsys, hkin, syncComToVars, log: message => messages.push(message) };
  const size = Math.max(...profiles.map(p => p.length)) + 2;
  for (const name of ['XSSI', 'UEDG', 'CTAU', 'THET', 'DSTR', 'MASS', 'TAU', 'DIS', 'CTQ', 'DELT', 'TSTR'])
    ctx[name] = Array.from({ length: size }, () => new Float64Array(3));
  for (const name of ['NBL', 'IBLTE', 'ITRAN', 'XSSITR', 'TFORCE']) ctx[name] = new Float64Array(3);
  ctx.ACRIT = [0, ncrit, ncrit]; ctx.XSTRIP = [0, 1, 1];
  profiles.forEach((profile, side) => {
    const is = side + 1; ctx.NBL[is] = ctx.IBLTE[is] = profile.length + 1;
    profile.forEach((p, i) => { ctx.XSSI[i + 2][is] = p.s; ctx.UEDG[i + 2][is] = p.ue; });
  });
  profiles.forEach((profile, side) => {
    const is = side + 1; ctx.ITRAN[is] = phases[side] + 2;
    profile.forEach((p, i) => {
      ctx.CTAU[i + 2][is] = p.aux; ctx.THET[i + 2][is] = p.theta; ctx.DSTR[i + 2][is] = p.deltaStar;
      ctx.MASS[i + 2][is] = p.deltaStar * p.ue;
    });
  });
  ensureCtx(ctx); blpini(ctx); mrchdu(ctx);
  const surfaces = profiles.map((profile, side) => {
    const is = side + 1, transition = ctx.ITRAN[is] - 2;
    const states = profile.map((p, i) => ({ s: p.s, ue: ctx.UEDG[i + 2][is], aux: ctx.CTAU[i + 2][is],
      theta: ctx.THET[i + 2][is], deltaStar: ctx.DSTR[i + 2][is] }));
    if (!Number.isInteger(transition) || transition < 1 || transition >= profile.length
      || states.some(p => !Object.values(p).every(Number.isFinite) || p.ue <= 0 || p.theta <= 0 || p.deltaStar <= p.theta))
      throw new Error('MRCHDU returned an unusable surface guess.');
    return { states, transition, s: ctx.XSSITR[is], forced: Boolean(ctx.TFORCE[is]),
      targetHK: null };
  });
  return { surfaces, messages, method: 'XFOIL MRCHDU surface relaxation', flowSolved: false,
    localConvergenceWarnings: messages.filter(message => message.includes('Convergence failed')) };
}
