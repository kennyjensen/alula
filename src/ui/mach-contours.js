// SPDX-License-Identifier: GPL-2.0-or-later
// Display-only linear interpolation on the dual lattice of section/tube
// centers. Keep passages separate: never interpolate through a solid or wake.
export function machContours(flow, step = .05) {
  if (!(step > 0) || !Number.isFinite(step)) throw new Error('Positive contour spacing required.');
  const groups = new Map(), contours = new Map();
  for (const line of flow?.lines ?? []) {
    if (!line.machNumbers) continue;
    if (!groups.has(line.group)) groups.set(line.group, []);
    groups.get(line.group).push(line);
  }
  const sample = (line, i) => ({ x: (line.points[i].x + line.points[i + 1].x) / 2,
    y: (line.points[i].y + line.points[i + 1].y) / 2, m: line.machNumbers[i] });
  function triangle(v) {
    if (!v.every(p => Number.isFinite(p.m) && Number.isFinite(p.x) && Number.isFinite(p.y))) return;
    const low = Math.min(...v.map(p => p.m)), high = Math.max(...v.map(p => p.m));
    // Bound rendering work for pathological unconverged states.
    for (let k = Math.max(1, Math.ceil(low / step)); k <= Math.min(200, Math.floor(high / step)); k++) {
      const level = k * step, hits = [];
      for (let j = 0; j < 3; j++) {
        const a = v[j], b = v[(j + 1) % 3];
        if ((a.m < level) === (b.m < level)) continue;
        const t = (level - a.m) / (b.m - a.m);
        hits.push({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
      }
      if (hits.length === 2) {
        if (!contours.has(k)) contours.set(k, { level, segments: [] });
        contours.get(k).segments.push(hits);
      }
    }
  }
  for (const lines of groups.values()) {
    lines.sort((a, b) => a.tube - b.tube);
    for (let j = 1; j < lines.length; j++) {
      const a = lines[j - 1], b = lines[j];
      if (b.tube !== a.tube + 1) continue;
      for (let i = 1; i < Math.min(a.machNumbers.length, b.machNumbers.length); i++) {
        const p = [sample(a, i - 1), sample(a, i), sample(b, i), sample(b, i - 1)];
        triangle([p[0], p[1], p[2]]); triangle([p[0], p[2], p[3]]);
      }
    }
  }
  return [...contours.values()].sort((a, b) => a.level - b.level);
}
