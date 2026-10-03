// SPDX-License-Identifier: GPL-2.0-or-later
import { cp, mkdir, rm } from 'node:fs/promises';
import { sep } from 'node:path';
await rm('dist', { recursive: true, force: true });
await mkdir('dist', { recursive: true });
for (const path of ['index.html', 'euler.html', 'src', 'public', 'CNAME', 'LICENSE', 'NOTICE.md', 'third_party/ises', 'third_party/grape', 'third_party/klu', 'third_party/cdt2d']) {
  await cp(path, `dist/${path}`, {
    recursive: true,
    ...(path === 'src' ? { filter: source => !source.split(sep).includes('tests') } : {}),
  });
}
// Include notices and the coordinate sources embedded in the benchmark presets.
await cp('third_party/README.md', 'dist/third_party/README.md');
for (const dataset of ['rae2822-mses', 'rae2822', 'nlr7301', '30p30n']) {
  const directory = `third_party/airfoils/${dataset}`;
  await mkdir(`dist/${directory}`, { recursive: true });
  await cp(`${directory}/README.md`, `dist/${directory}/README.md`);
}
for (const file of [
  'rae2822-mses/blade.rae', 'rae2822/geom.txt',
  'nlr7301/nlr7301-gap26.dat', 'nlr7301/provenance.json',
  '30p30n/30p-30n.dat', '30p30n/provenance.json',
]) await cp(`third_party/airfoils/${file}`, `dist/third_party/airfoils/${file}`);
console.log('Static browser app built in dist/. No runtime dependencies or server-side solver.');
