// SPDX-License-Identifier: GPL-2.0-or-later
import { spawnSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
const directory = fileURLToPath(new URL('.', import.meta.url));
const compiler = fileURLToPath(new URL('../../../node_modules/assemblyscript/bin/asc.js', import.meta.url));
const options = ['--optimize', '--runtime', 'incremental'];
const result = spawnSync(process.execPath, [compiler, directory + 'smoothing.ts', '-o', directory + 'smoothing.wasm', ...options], { stdio: 'inherit' });
if (result.error) throw result.error;
if (result.status) process.exit(result.status);
const hash = name => createHash('sha256').update(readFileSync(directory + name)).digest('hex');
writeFileSync(directory + 'BUILD.json', JSON.stringify({ compiler: 'assemblyscript', version: '0.28.9', options,
  sourceSha256: hash('smoothing.ts'), wasmSha256: hash('smoothing.wasm'), floatingPoint: 'Float64; no fast-math' }, null, 2) + '\n');
