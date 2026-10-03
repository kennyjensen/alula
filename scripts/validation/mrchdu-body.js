// SPDX-License-Identifier: GPL-2.0-or-later
// Research preparation through both surfaces and their merged sharp-TE wake.
// Native IBLPAN solves the wake only in side 2; XICALC puts its first point
// at the same arclength as the lower TE. Not a coupled flow solution.
import { mrchdu, blpini } from '../../src/viscous/xfoil/xbl.js';
import { ensureCtx, blprv, blkin, blvar, blmid, trchek, tesys, blsys, hkin, syncComToVars } from '../../src/viscous/xfoil/xblsys.js';
import { createIntegralKernel } from '../../src/viscous/integral.js';

export function relaxXfoilBody({ surfaces, wake, phases }, parameters = {}) {
  const { reynolds, mach, gamma, ncrit, velocityConvention } = createIntegralKernel(parameters).parameters;
  if (velocityConvention !== 'physical') throw new Error('MRCHDU requires physical edge velocities.');
  const ordered = p => Array.isArray(p) && p.length >= 2 && p.every((v, i) =>
    ['s', 'ue', 'aux', 'theta', 'deltaStar'].every(k => Number.isFinite(v[k]))
    && v.s > 0 && v.ue > 0 && v.theta > 0 && v.deltaStar > v.theta && (!i || v.s > p[i - 1].s));
  if (!Array.isArray(surfaces) || surfaces.length !== 2 || !surfaces.every(ordered) || !ordered(wake))
    throw new Error('Supply ordered upper, lower and wake states with positive speed and thickness.');
  if (!Array.isArray(phases) || phases.length !== 2 || phases.some((p, i) => !Number.isInteger(p) || p < 1 || p >= surfaces[i].length))
    throw new Error('Supply the old transition interval on each surface.');
  if (wake[0].s !== surfaces[1].at(-1).s)
    throw new Error('The first wake state must share the lower trailing-edge arclength.');
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
  const profiles = [surfaces[0], [...surfaces[1], ...wake]], size = Math.max(...profiles.map(p => p.length)) + 2;
  for (const name of ['XSSI', 'UEDG', 'CTAU', 'THET', 'DSTR', 'MASS', 'TAU', 'DIS', 'CTQ', 'DELT', 'TSTR'])
    ctx[name] = Array.from({ length: size }, () => new Float64Array(3));
  for (const name of ['NBL', 'IBLTE', 'ITRAN', 'XSSITR', 'TFORCE']) ctx[name] = new Float64Array(3);
  ctx.ACRIT = [0, ncrit, ncrit]; ctx.XSTRIP = [0, 1, 1]; ctx.WGAP = new Float64Array(wake.length + 1);
  profiles.forEach((profile, side) => {
    const is = side + 1; ctx.NBL[is] = profile.length + 1;
    ctx.IBLTE[is] = surfaces[side].length + 1; ctx.ITRAN[is] = phases[side] + 2;
    profile.forEach((p, i) => {
      ctx.XSSI[i + 2][is] = p.s; ctx.UEDG[i + 2][is] = p.ue; ctx.CTAU[i + 2][is] = p.aux;
      ctx.THET[i + 2][is] = p.theta; ctx.DSTR[i + 2][is] = p.deltaStar; ctx.MASS[i + 2][is] = p.deltaStar * p.ue;
    });
  });
  ensureCtx(ctx); blpini(ctx); mrchdu(ctx);
  const states = profiles.map((profile, side) => profile.map((p, i) => ({ s: p.s, ue: ctx.UEDG[i + 2][side + 1],
    aux: ctx.CTAU[i + 2][side + 1], theta: ctx.THET[i + 2][side + 1], deltaStar: ctx.DSTR[i + 2][side + 1] })));
  const result = { surfaces: surfaces.map((profile, side) => ({ states: states[side].slice(0, profile.length),
    transition: ctx.ITRAN[side + 1] - 2, s: ctx.XSSITR[side + 1], forced: Boolean(ctx.TFORCE[side + 1]) })),
  wake: states[1].slice(surfaces[1].length), messages, method: 'XFOIL MRCHDU surfaces and sharp-TE wake', flowSolved: false,
  localConvergenceWarnings: messages.filter(message => message.includes('Convergence failed')) };
  if (!result.surfaces.every(s => ordered(s.states) && Number.isInteger(s.transition) && s.transition >= 1 && s.transition < s.states.length)
    || !ordered(result.wake)) throw new Error('MRCHDU returned an unusable body/wake guess.');
  return result;
}
