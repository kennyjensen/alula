// SPDX-License-Identifier: GPL-2.0-or-later
import { readFileSync } from 'node:fs';

export function nativeSurfaceLimits() {
  const include = readFileSync('third_party/Xfoil/src/XFOIL.INC', 'utf8');
  const solve = readFileSync('third_party/Xfoil/src/xsolve.f', 'utf8');
  const iqx = Number(include.match(/\bIQX\s*=\s*(\d+)/)?.[1]);
  const nvx = Number(solve.match(/\bNVX\s*=\s*(\d+)/)?.[1]);
  if (!(iqx > 6 && nvx > 1)) throw new Error('Cannot read original native surface-array limits.');
  // XFOIL.INC reserves six IQX slots. xpanel.f factors N+1 rows through
  // LUDCMP, whose local VV array is separately dimensioned to NVX.
  return { iqx, nvx, maximumSurfaceNodes: Math.min(iqx - 6, nvx - 1) };
}

export function assertNativeSurfaceCount(count, limits = nativeSurfaceLimits()) {
  if (!Number.isInteger(count) || count < 4 || count > limits.maximumSurfaceNodes)
    throw new Error(`Native surface node count ${count} exceeds its supported range 4–${limits.maximumSurfaceNodes} (IQX=${limits.iqx}, LUDCMP NVX=${limits.nvx}, matrix N+1).`);
}
