// SPDX-License-Identifier: GPL-2.0-or-later
import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { resolve } from 'node:path';
import { spawnSync } from 'node:child_process';

export async function buildReference(directory, driver = 'driver') {
  await mkdir(directory, { recursive: true });
  const source = resolve('third_party/Xfoil/src');
  const xfoil = await readFile(`${source}/xfoil.f`, 'utf8');
  // Remove the interactive PROGRAM, retaining every subroutine unchanged.
  await writeFile(`${directory}/xfoil-core.f`, xfoil.slice(xfoil.indexOf('      SUBROUTINE INIT')));
  const files = ['xpanel', 'xoper', 'xsolve', 'xbl', 'xblsys', 'xgeom', 'xutils', 'spline', 'xtcam', 'xgdes', 'userio'];
  const executable = `${directory}/xfoil-${driver}`;
  const args = ['-O2', '-std=legacy', '-fdefault-real-8', '-fallow-argument-mismatch',
    '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', `-I${source}`,
    resolve(`scripts/reference/${driver}.f`), `${directory}/xfoil-core.f`,
    ...files.map(f => `${source}/${f}.f`), '-o', executable];
  const fc = process.env.FC || 'gfortran';
  const built = spawnSync(fc, args, { encoding: 'utf8', timeout: 120_000 });
  if (built.error || built.status !== 0) throw new Error(`Fortran build failed (${fc}): ${built.error?.message ?? ''}\n${built.stderr}`);
  return { executable, compiler: spawnSync(fc, ['--version'], { encoding: 'utf8' }).stdout.split('\n')[0], flags: args.slice(0, 7), files: ['xfoil', ...files].map(f => `${source}/${f}.f`) };
}
