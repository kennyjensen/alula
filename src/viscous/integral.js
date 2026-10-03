// SPDX-License-Identifier: GPL-2.0-or-later
// Topology-independent Drela integral-BL blocks. The native kernels retain
// their provenance in xfoil/; this adapter has no panel or single-body arrays.
import { blpini } from './xfoil/xbl.js';
import { ensureCtx, blprv, blkin, blvar, blsys, trchek, tesys, copyCom, hkin, XFOIL_TRANSITION_ROOT_ACCURACY } from './xfoil/xblsys.js';

// Copy only public primitive inputs, only after a root failure. Native context
// arrays, aliases and caller metadata do not belong in a replay fixture.
const transitionRootStation = station => station && Object.fromEntries(
  ['s', 'aux', 'theta', 'deltaStar', 'ue', 'wakeGap', 'amplification']
    .filter(key => station[key] !== undefined && (station[key] === null || Number.isFinite(station[key])))
    .map(key => [key, station[key]]));

export function createIntegralKernel({ reynolds = 1e6, mach = 0, gamma = 1.4, ncrit = 9,
  velocityConvention = 'physical',exactJacobian=false, transitionTolerance = 5e-5,
  hkFloorLinearization = 'exact' } = {}) {
  if (![reynolds, mach, gamma, ncrit].every(Number.isFinite) || reynolds <= 0 || mach < 0 || mach >= 1 || gamma <= 1
    || ncrit < 0 || !['physical', 'xfoil'].includes(velocityConvention)) throw new Error('Invalid integral boundary-layer parameters.');
  if (!Number.isFinite(transitionTolerance) || transitionTolerance <= 0 || transitionTolerance > 5e-5)
    throw new Error('Invalid transition amplification root tolerance.');
  if (!['exact', 'native'].includes(hkFloorLinearization)) throw new Error('Invalid Hk-floor linearization policy.');
  const gm1 = gamma - 1; const totalRatio = 1 + .5 * gm1 * mach ** 2;
  const hstinv = gm1 * mach ** 2 / totalRatio;
  const hstinvMs = gm1 / totalRatio ** 2;
  const herat = 1 - .5 * hstinv; const heratMs = -.5 * hstinvMs;
  const hvrat = .35;
  const viscousFactor = herat ** 1.5 * (1 + hvrat) / (herat + hvrat);
  const beta = Math.sqrt(1 - mach * mach);
  const tk = velocityConvention === 'xfoil' ? (1 - beta) / (1 + beta) : 0;
  const tkMs = velocityConvention === 'xfoil' ? 1 / (beta * (1 + beta) ** 2) : 0;
  const ctx = { QINFBL: 1, TKBL: tk, TKBL_MS: tkMs,
    RSTBL: totalRatio ** (1 / gm1), RSTBL_MS: .5 * totalRatio ** (1 / gm1 - 1),
    HSTINV: hstinv, HSTINV_MS: hstinvMs,
    REYBL: reynolds * viscousFactor, REYBL_RE: viscousFactor,
    REYBL_MS: reynolds * viscousFactor * (1.5 / herat - 1 / (herat + hvrat)) * heratMs,
    GAMBL: gamma, GM1BL: gm1, HVRAT: hvrat, AMCRIT: ncrit, BULE: 1, IDAMPV: 0, EXACT_BL_JACOBIAN: exactJacobian,
    TRANSITION_TOLERANCE: transitionTolerance, NATIVE_HK_FLOOR_LINEARIZATION: hkFloorLinearization === 'native' };
  ensureCtx(ctx); blpini(ctx);
  const captureTransitionRootFailure = (error, method, upstream, downstream, tripS, regime, similarityExponent) => {
    if (error?.code === XFOIL_TRANSITION_ROOT_ACCURACY) error.transitionRootFailure = {
      version: 1, code: XFOIL_TRANSITION_ROOT_ACCURACY,
      parameters: { reynolds, mach, gamma, ncrit, velocityConvention, exactJacobian, transitionTolerance,
        ...(hkFloorLinearization !== 'exact' ? { hkFloorLinearization } : {}) },
      method, input: { upstream: transitionRootStation(upstream), downstream: transitionRootStation(downstream), tripS,
        ...(method === 'interval' ? { regime, similarityExponent } : {}) },
    };
    throw error;
  };
  const fill = station => {
    const { s, theta, deltaStar, ue, aux = 0, wakeGap = 0 } = station;
    if (![s, theta, deltaStar, ue, aux, wakeGap].every(Number.isFinite) || s <= 0 || theta <= 0 || ue <= 0 || deltaStar - wakeGap <= theta || wakeGap < 0) throw Object.assign(new Error('Inadmissible integral BL station.'),{station:{...station}});
    blprv(s, station.amplification ?? aux, aux, theta, deltaStar, wakeGap, ue, ctx); blkin(ctx);
    if (!(ctx.HK2 > 1) || ![ctx.R2, ctx.V2, ctx.RT2, ctx.M2].every(Number.isFinite) || ctx.R2 <= 0 || ctx.V2 <= 0) {
      // Diagnose only after the unchanged BLKIN domain check has failed.
      // In particular, Hk <= 1 is a shape-domain failure, not a sonic gate.
      const enthalpyRatio = 1 - .5 * ctx.U2 * ctx.U2 * ctx.HSTINV;
      const failedChecks = [];
      if (!(ctx.HK2 > 1)) failedChecks.push('raw-hk');
      if (!Number.isFinite(ctx.R2) || ctx.R2 <= 0) failedChecks.push('density');
      if (!Number.isFinite(ctx.V2) || ctx.V2 <= 0) failedChecks.push('viscosity');
      if (!Number.isFinite(ctx.RT2)) failedChecks.push('re-theta');
      if (!Number.isFinite(ctx.M2)) failedChecks.push('mach-squared');
      const condition = !Number.isFinite(enthalpyRatio) || enthalpyRatio <= 0 ? 'thermal-energy' : failedChecks[0];
      const detail = {
        'thermal-energy': `edge thermal-energy ratio must be finite and positive (h/h0=${enthalpyRatio}, Ue=${ctx.U2}).`,
        'raw-hk': `raw kinematic shape factor Hk must exceed 1 (Hk=${ctx.HK2}, fluid H=${ctx.H2}, edge Mach squared=${ctx.M2}).`,
        density: `edge density must be finite and positive (rho=${ctx.R2}).`,
        viscosity: `edge viscosity must be finite and positive (nu=${ctx.V2}).`,
        're-theta': `momentum-thickness Reynolds number must be finite (ReTheta=${ctx.RT2}).`,
        'mach-squared': `edge Mach squared must be finite (MachSquared=${ctx.M2}).`
      }[condition];
      throw Object.assign(new Error(`Inadmissible compressible BL edge state. ${detail}`), {
        code: 'BL_EDGE_STATE_DOMAIN', station: { ...station },
        diagnostics: { condition, failedChecks, rawHk: ctx.HK2, fluidShapeFactor: ctx.H2,
          enthalpyRatio, density: ctx.R2, viscosity: ctx.V2, reTheta: ctx.RT2,
          machSquared: ctx.M2, physicalUe: ctx.U2 }
      });
    }
  };
  const properties = () => ({ ue: ctx.U2, rho: ctx.R2, viscosity: ctx.V2, machSquared: ctx.M2,
    // BLVAR bounds Hk for closure evaluation. Inverse kinematic constraints
    // must instead use BLKIN/HKIN's unbounded shape relation; its derivative
    // remains nonzero inside the closure's constant branch.
    h: ctx.H2, rawHk: hkin(ctx.H2, ctx.M2).hk, hk: ctx.HK2, hStar: ctx.HS2, reTheta: ctx.RT2, cf: ctx.CF2,
    dissipation: ctx.DI2, equilibriumShear: ctx.CQ2,
    transitionShear: ctx.CTRCON*Math.exp(-ctx.CTRCEX/(ctx.HK2-1))*ctx.CQ2 });
  const output = () => {
    // BLDIF stores the Newton RHS (-R), but VS1/VS2 differentiate R.
    const residual = Array.from(ctx.VSREZ).slice(1, 4).map(v => -v);
    const upstream = ctx.VS1.slice(1, 4).map(r => [...r.slice(1, 6)]);
    const downstream = ctx.VS2.slice(1, 4).map(r => [...r.slice(1, 6)]);
    if (![...residual, ...upstream.flat(), ...downstream.flat()].every(Number.isFinite)) throw new Error('Nonfinite integral BL residual/Jacobian.');
    return { residual, upstream, downstream, reynoldsDerivative: [...ctx.VSR.slice(1, 4)],
      machSquaredDerivative: [...ctx.VSM.slice(1, 4)], tripDerivative: [...ctx.VSX.slice(1, 4)], properties: properties() };
  };
  const interval = ({ upstream, downstream, regime = 'laminar', tripS = Number.MAX_VALUE, similarityExponent = 1 },
    { hkFloorCorrection = false } = {}) => {
    if (!['similarity', 'laminar', 'turbulent', 'wake', 'transition'].includes(regime)) throw new Error('Unknown BL interval regime.');
    if (!Number.isFinite(tripS) || !Number.isFinite(similarityExponent)) throw new Error('Invalid BL transition/similarity controls.');
    ctx.SIMI = regime === 'similarity'; ctx.WAKE = regime === 'wake';
    ctx.HK_FLOOR_LINEARIZATION_ACTIVE = false;
    ctx.TRAN = regime === 'transition'; ctx.TURB = regime === 'turbulent' || regime === 'wake';
    ctx.BULE = similarityExponent; ctx.XIFORC = tripS;
    // Set every primitive and derived upstream value; interval evaluation is
    // independent of which other element/surface was evaluated previously.
    fill(upstream ?? downstream);
    blvar(regime === 'wake' ? 3 : regime === 'turbulent' ? 2 : 1, ctx); copyCom(ctx, 2, 1);
    fill(regime === 'transition' ? { ...downstream, amplification: upstream.aux } : downstream);
    if (!ctx.SIMI && !(ctx.X2 > ctx.X1)) throw new Error('BL interval stations must increase downstream.');
    let transition = null;
    if (regime === 'transition') {
      try { trchek(ctx); }
      catch (error) { captureTransitionRootFailure(error, 'interval', upstream, downstream, tripS, regime, similarityExponent); }
      if (!ctx.TRAN) throw new Error('Transition is outside this active interval.');
      transition = { s: ctx.XT, forced: ctx.TRFORC, amplification: ctx.AMPL2 };
    }
    blsys(ctx);
    const value = { ...output(), transition };
    if (!ctx.HK_FLOOR_LINEARIZATION_ACTIVE) return value;
    value.hkFloorLinearizationUsed = true;
    if (hkFloorCorrection) {
      // Isolate the native Hk-floor extension from every other analytic
      // partial. Re-evaluate the identical interval and local root with
      // only that policy disabled; no primitive or residual is changed.
      let baseline;
      ctx.NATIVE_HK_FLOOR_LINEARIZATION = false;
      try { baseline = interval({ upstream, downstream, regime, tripS, similarityExponent }); }
      finally { ctx.NATIVE_HK_FLOOR_LINEARIZATION = true; }
      value.hkFloorJacobianCorrection = {
        upstream: value.upstream.map((row, r) => row.map((v, k) => v - baseline.upstream[r][k])),
        downstream: value.downstream.map((row, r) => row.map((v, k) => v - baseline.downstream[r][k])),
        trip: value.tripDerivative.map((v, r) => v - baseline.tripDerivative[r]),
      };
    }
    return value;
  };
  const station = (state, regime = 'laminar') => {
    ctx.SIMI = false; fill(state); blvar(regime === 'wake' ? 3 : regime === 'turbulent' ? 2 : 1, ctx);
    return properties();
  };
  const transitionCheck = ({upstream,downstream,tripS=Number.MAX_VALUE}) => {
    ctx.SIMI=false;ctx.WAKE=false;ctx.TURB=false;ctx.TRAN=false;ctx.XIFORC=tripS;
    fill(upstream);blvar(1,ctx);copyCom(ctx,2,1);
    fill({...downstream,amplification:upstream.aux});
    try { trchek(ctx); }
    catch (error) { captureTransitionRootFailure(error, 'transitionCheck', upstream, downstream, tripS); }
    const result={transition:ctx.TRAN,amplification:ctx.AMPL2,s:ctx.XT,forced:ctx.TRFORC};
    if(!result.transition){
      // TRCHEK stops its internal N iteration at a looser tolerance than our
      // global residual target. Refine the same native amplification equation
      // before projecting N during an active-set update; otherwise that
      // projection can impose a ~1e-7 floor on the full Newton residual.
      for(let iteration=0;iteration<8;iteration++){
        const block=interval({upstream,downstream:{...downstream,aux:result.amplification},regime:'laminar'});
        if(Math.abs(block.residual[0])<1e-12)break;
        result.amplification-=block.residual[0]/block.downstream[0][0];
      }
    }
    return result;
  };
  const trailingEdge = (upper, lower, firstWake, gap = 0) => {
    if (!(gap >= 0) || !Number.isFinite(gap)) throw new Error('Invalid trailing-edge gap.');
    const theta = upper.theta + lower.theta;
    const deltaStar = upper.deltaStar + lower.deltaStar + gap;
    const aux = (upper.aux * upper.theta + lower.aux * lower.theta) / theta;
    ctx.SIMI = false; fill(firstWake); tesys(aux, theta, deltaStar, ctx);
    return { ...output(), matched: { theta, deltaStar, aux } };
  };
  return { interval, station, transitionCheck, trailingEdge, parameters: { reynolds, mach, gamma, ncrit, velocityConvention,exactJacobian,
    ...(transitionTolerance !== 5e-5 ? { transitionTolerance } : {}),
    ...(hkFloorLinearization !== 'exact' ? { hkFloorLinearization } : {}) } };
}
