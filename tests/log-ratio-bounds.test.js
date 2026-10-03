import test from 'node:test';
import assert from 'node:assert/strict';
import { fitLogRatios } from '../src/numerics/log-ratio-fit.js';
import { passageStationTargets } from '../src/geometry/passage-density-stations.js';
const close = (a, b) => assert.ok(Math.abs(a - b) < 2e-13, `${a} versus ${b}`);

test('one-sided endpoint bounds give the analytic constrained gauge without changing ratio residuals', () => {
  const connections = [{ left: 0, right: 1, ratio: 2 }];
  const a = fitLogRatios({ natural: [2, 8], connections });
  const b = fitLogRatios({ natural: [2, 8], connections, maximumValues: [1, Infinity] });
  close(a.values[0], Math.sqrt(8)); close(b.values[0], 1); close(b.values[1], 2);
  close(a.residuals[0], b.residuals[0]); assert.deepEqual(b.limitedGroups, [0]);
});

test('incompatible parallel ratios retain their optimum and endpoint units scale together', () => {
  const connections = [{ left: 0, right: 1, ratio: 2 }, { left: 0, right: 1, ratio: 8 }];
  const a = fitLogRatios({ natural: [3, 3, 9], connections, maximumValues: [Infinity, 4, 2] });
  [1, 4, 2].forEach((v, i) => close(a.values[i], v));
  close(a.residuals[0], Math.log(2)); close(a.residuals[1], -Math.log(2));
  const b = fitLogRatios({ natural: [15, 15, 45], connections, maximumValues: [Infinity, 20, 10] });
  a.values.forEach((v, i) => close(b.values[i], 5 * v));
  assert.throws(() => fitLogRatios({ natural: [1], connections: [], maximumValues: [0] }), /upper bounds/);
});

test('passage endpoint limits preserve the original density request at a coarse cut to fine wall join', () => {
  const uniform = n => Array.from({ length: n + 1 }, (_, i) => i / n);
  const components = [
    { start: 0, end: 1, intervals: 16, boundaries: [{ kind: 'cut', body: 0, end: 'upstream', length: 4, progress: uniform(16) }] },
    ...['upper', 'lower'].map(side => ({ start: 1, end: 2, intervals: 4,
      boundaries: [{ kind: 'wall', body: 0, side, length: 1, range: [0, 1], requested: uniform(8) }] })),
    { start: 2, end: 3, intervals: 16, boundaries: [{ kind: 'cut', body: 0, end: 'wake', length: 4, progress: uniform(16) }] },
  ];
  const a = passageStationTargets({ components }), b = passageStationTargets({ components, retainEndpointResolution: true });
  close(a.endpointFit.values[2], .25); close(b.endpointFit.values[2], .125);
  close(b.endpointFit.values[1] * 4, .125);
  a.endpointFit.residuals.forEach((v, i) => close(v, b.endpointFit.residuals[i]));
  for (const k of [2, 3, 4, 5]) assert.ok(b.endpointFit.values[k] <= .125 + 1e-14);
});
