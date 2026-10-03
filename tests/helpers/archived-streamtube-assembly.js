// SPDX-License-Identifier: GPL-2.0-or-later
// Keep all four changed assembly modules on their recorded implementation.
// Other physical kernels remain shared; full-size source-bound matrix receipts
// independently compare the complete assembly before and after optimization.
import fs from 'node:fs';

const entries = [
  ['numerics/sparse.js', 'sparse.js.txt'],
  ['euler/streamtube-body-jacobian.js', 'streamtube-body-jacobian.js.txt'],
  ['euler/streamtube-body.js', 'streamtube-body.js.txt'],
  ['euler/streamtube-coupled.js', 'streamtube-coupled.js.txt'],
];
const urls = new Map();
for (const [relative, archive] of entries) {
  const original = new URL(`../../src/${relative}`, import.meta.url);
  const text = fs.readFileSync(new URL(`../../docs/performance-assembly/before/${archive}`, import.meta.url), 'utf8');
  const source = text.replace(/from\s+(['"])(\.[^'"]+)\1/g, (_, quote, path) => {
    const absolute = new URL(path, original).href;
    return `from ${quote}${urls.get(absolute) ?? absolute}${quote}`;
  });
  urls.set(original.href, `data:text/javascript;base64,${Buffer.from(`${source}\n//# sourceURL=mses-assembly-before/${relative}\n`).toString('base64')}`);
}

export const archivedBody = await import(urls.get(new URL('../../src/euler/streamtube-body.js', import.meta.url).href));
export const archivedCoupled = await import(urls.get(new URL('../../src/euler/streamtube-coupled.js', import.meta.url).href));
