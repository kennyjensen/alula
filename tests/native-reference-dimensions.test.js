import test from 'node:test';
import assert from 'node:assert/strict';
import { nativeSurfaceLimits, assertNativeSurfaceCount } from '../scripts/reference/dimensions.js';

test('native reference preflight respects the separate N+1 LUDCMP work array', () => {
  const limits = nativeSurfaceLimits();
  assert.deepEqual(limits, { iqx: 601, nvx: 515, maximumSurfaceNodes: 514 });
  for (const count of [321, 351, 411, 514]) assertNativeSurfaceCount(count, limits);
  for (const count of [515, 531, 595]) assert.throws(() => assertNativeSurfaceCount(count, limits), /LUDCMP NVX=515/);
  assert.throws(() => assertNativeSurfaceCount(3, limits), /supported range/);
});
