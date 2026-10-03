// SPDX-License-Identifier: GPL-2.0-or-later
// Small native BL oracle build: link only the unchanged BL source files.
import { mkdir } from 'node:fs/promises';
import { resolve, join } from 'node:path';
import { spawnSync } from 'node:child_process';

export async function buildBLReference(directory, driver) {
  await mkdir(directory, { recursive: true });
  const source = resolve('third_party/Xfoil/src'), compiler = process.env.FC ?? 'gfortran';
  const flags = ['-O2', '-std=legacy', '-fdefault-real-8', '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections'];
  const files = [`scripts/reference/${driver}.f`, 'third_party/Xfoil/src/xblsys.f', 'third_party/Xfoil/src/xbl.f'];
  const executable = join(directory, driver), built = spawnSync(compiler, [...flags, `-I${source}`, ...files, '-o', executable], { encoding: 'utf8', timeout: 60000 });
  if (built.error || built.status !== 0) throw new Error(`Native BL build failed: ${built.error?.message ?? built.stderr}`);
  return { executable, flags, files, compiler: spawnSync(compiler, ['--version'], { encoding: 'utf8' }).stdout.split('\n')[0] };
}
