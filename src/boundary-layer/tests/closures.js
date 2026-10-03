// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from XFOIL src/xblsys.f, Copyright (C) 2000 Mark Drela.
// Modified 2026-09-05: isolated JS laminar correlations and exact derivatives.
// GPL v2 or later; distributed without warranty. See ../../../LICENSE.
function check(hk, rt = 1) {
  if (!Number.isFinite(hk) || hk <= 1 || !Number.isFinite(rt) || rt <= 0) {
    throw new Error('Boundary-layer closures require Hk > 1 and Rtheta > 0.');
  }
}

export function hkin(h, msq) {
  if (!Number.isFinite(h) || !Number.isFinite(msq) || msq < 0) throw new Error('Invalid shape factor or Mach squared.');
  const hk = (h - 0.29 * msq) / (1 + 0.113 * msq);
  return { hk, hkH: 1 / (1 + 0.113 * msq), hkMsq: (-0.29 - 0.113 * hk) / (1 + 0.113 * msq) };
}

// H* = kinetic-energy thickness / momentum thickness. HSL is not shear stress.
export function hsl(hk) {
  check(hk);
  let hs; let hsHk;
  if (hk < 4.35) {
    const t = hk - 4.35;
    hs = 0.0111 * t ** 2 / (hk + 1) - 0.0278 * t ** 3 / (hk + 1) + 1.528 - 0.0002 * (t * hk) ** 2;
    hsHk = 0.0111 * (2 * t - t ** 2 / (hk + 1)) / (hk + 1)
      - 0.0278 * (3 * t ** 2 - t ** 3 / (hk + 1)) / (hk + 1) - 0.0004 * t * hk * (t + hk);
  } else {
    hs = 0.015 * (hk - 4.35) ** 2 / hk + 1.528;
    hsHk = 0.030 * (hk - 4.35) / hk - 0.015 * (hk - 4.35) ** 2 / hk ** 2;
  }
  return { hs, hsHk, hsRt: 0, hsMsq: 0 };
}

// Cf = tau_wall / (0.5 rho_edge U_edge^2), not tau/(rho U^2).
export function cfl(hk, rt) {
  check(hk, rt);
  let cf; let cfHk;
  if (hk < 5.5) {
    const t = (5.5 - hk) ** 3 / (hk + 1);
    cf = (0.0727 * t - 0.07) / rt;
    cfHk = (-0.2181 * (5.5 - hk) ** 2 / (hk + 1) - 0.0727 * t / (hk + 1)) / rt;
  } else {
    const t = 1 - 1 / (hk - 4.5);
    cf = (0.015 * t ** 2 - 0.07) / rt;
    cfHk = 0.030 * t / (hk - 4.5) ** 2 / rt;
  }
  return { cf, cfHk, cfRt: -cf / rt, cfMsq: 0 };
}

// DI = 2 C_dissipation / H*. This is NOT an airfoil drag coefficient.
export function dil(hk, rt) {
  check(hk, rt);
  let di; let diHk;
  if (hk < 4) {
    di = (0.00205 * (4 - hk) ** 5.5 + 0.207) / rt;
    diHk = -0.00205 * 5.5 * (4 - hk) ** 4.5 / rt;
  } else {
    const t = hk - 4; const den = 1 + 0.02 * t * t;
    di = (-0.0016 * t * t / den + 0.207) / rt;
    diHk = -0.0032 * t * (1 / den - 0.02 * t * t / den ** 2) / rt;
  }
  return { di, diHk, diRt: -di / rt };
}
