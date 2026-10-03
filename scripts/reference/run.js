// SPDX-License-Identifier: GPL-2.0-or-later
import { spawnSync } from 'node:child_process';

export function runReference(executable, points, options = {}) {
  const { alpha = 0, reynolds = 1e6, ncrit = 9, trips = [1, 1], maxIterations = 100 } = options;
  // Only the separate compressible driver reads this final field. Existing
  // callers and the original fixture input remain byte-for-byte unchanged.
  const extra=options.mach===undefined?[]:[options.mach];
  const input = [points.length, alpha, reynolds, ncrit, ncrit, ...trips, maxIterations,...extra].join(' ')
    + '\n' + points.map(p => `${p.x} ${p.y}`).join('\n') + '\n';
  const run = spawnSync(executable, [], { input, encoding: 'utf8', timeout: 60_000, maxBuffer: 8 * 1024 ** 2 });
  if (run.error || run.status !== 0) throw new Error(`Native XFOIL failed: ${run.error?.message ?? ''}\n${run.stderr}\n${run.stdout.slice(-3000)}`);
  const result = { cp: [], bl: [], log: run.stdout };
  for (const line of run.stdout.split('\n')) {
    const [tag, ...words] = line.trim().split(/\s+/);
    const v = words.map(Number);
    if (tag === 'RESULT') {
      result.converged = words[0] === 'T';
      [result.cl, result.cm, result.cd, result.cdf, result.pressureIntegralDrag,
        result.rmsUpdate, ...result.transition] = v.slice(1);
    }
    if (tag === 'CP') {
      const [index, x, y, cpInviscid, cp, ue] = v;
      result.cp.push({ index, x, y, cpInviscid, cp, ue });
    }
    if (tag === 'BL') {
      const [side, station, index, s, ue, theta, deltaStar, ctau, mass, tau, dissipation] = v;
      result.bl.push({ side, station, index, s, ue, theta, deltaStar, ctau, mass, tau, dissipation });
    }
  }
  if (typeof result.converged !== 'boolean') throw new Error('Missing native RESULT record.');
  return result;
}
