// SPDX-License-Identifier: GPL-2.0-or-later
// Isolated station-pair assembly for native parity and derivative verification.
import { blpini } from '../xfoil/xbl.js';
import { ensureCtx, blprv, blkin, blvar, blmid, bldif, copyCom } from '../xfoil/xblsys.js';

export function assembleInterval({ type, reynolds, stations }) {
  const ctx = { QINFBL: 1, TKBL: 0, TKBL_MS: 0.25, RSTBL: 1, RSTBL_MS: 0.5,
    HSTINV: 0, HSTINV_MS: 0.4, REYBL: reynolds, REYBL_MS: 0, REYBL_RE: 1,
    GAMBL: 1.4, GM1BL: 0.4, HVRAT: 0.35, AMCRIT: 9, BULE: 1, IDAMPV: 0 };
  ensureCtx(ctx); blpini(ctx);
  for (let side = 0; side < 2; side++) {
    blprv(...stations[side], ctx); blkin(ctx); blvar(Math.max(1, type), ctx);
    if (side === 0) copyCom(ctx, 2, 1);
  }
  blmid(Math.max(1, type), ctx); bldif(type, ctx);
  return { residual: Array.from(ctx.VSREZ).slice(1), reynoldsDerivative: Array.from(ctx.VSR).slice(1),
    machSquaredDerivative: Array.from(ctx.VSM).slice(1), tripDerivative: Array.from(ctx.VSX).slice(1),
    upstream: ctx.VS1.slice(1).map(r => Array.from(r).slice(1)),
    downstream: ctx.VS2.slice(1).map(r => Array.from(r).slice(1)) };
}
