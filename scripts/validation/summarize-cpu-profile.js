// SPDX-License-Identifier: GPL-2.0-or-later
// Inclusive entries overlap; self entries and the phase categories do not.
import fs from 'node:fs';
import { pathToFileURL } from 'node:url';

export function summarizeCpuProfile(profile) {
  const nodes = new Map(profile.nodes.map(n => [n.id, n])), parents = new Map();
  for (const n of profile.nodes) for (const id of n.children ?? []) parents.set(id, n.id);
  const self = new Map(), inclusive = new Map(), categories = new Map(); let elapsed = 0;
  const add = (map, key, ms) => map.set(key, (map.get(key) ?? 0) + ms);
  const key = n => `${n.callFrame.functionName || '(anonymous)'} — ${n.callFrame.url}:${n.callFrame.lineNumber + 1}`;
  for (let i = 0; i < (profile.samples?.length ?? 0); i++) {
    const milliseconds = profile.timeDeltas[i] / 1000, stack = []; elapsed += milliseconds;
    let id = profile.samples[i]; add(self, key(nodes.get(id)), milliseconds);
    while (id !== undefined) { stack.push(nodes.get(id)); id = parents.get(id); }
    for (const frame of new Set(stack.map(key))) add(inclusive, frame, milliseconds);
    const has = (name, suffix = '') => stack.some(n => n.callFrame.functionName === name && n.callFrame.url.endsWith(suffix));
    const category = has('(garbage collector)') ? 'garbage collection'
      : has('adopt', '/streamtube-coupled-ises.js') ? 'trial reconstruction and admission'
        : has('jacobian', '/streamtube-coupled.js') ? 'coupled Jacobian'
          : has('solveSparseDirect') || has('solveSparseDirectAligned') || has('solveCoupledLinearSystem') ? 'sparse linear solves'
            : 'other';
    add(categories, category, milliseconds);
  }
  const table = map => [...map].sort((a, b) => b[1] - a[1]).map(([name, milliseconds]) =>
    ({ name, milliseconds, percent: 100 * milliseconds / elapsed }));
  return { sampledMilliseconds: elapsed, categories: table(categories), self: table(self).slice(0, 30),
    inclusive: table(inclusive).slice(0, 30), note: 'Inclusive function times overlap; categories are disjoint. Sampling is approximate.' };
}
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href)
  console.log(JSON.stringify(summarizeCpuProfile(JSON.parse(fs.readFileSync(process.argv[2]))), null, 2));
