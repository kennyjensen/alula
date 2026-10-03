// Independent physical polygon audit. No Euler kernel or mesh-quality imports.
// The adaptive orientation predicate evaluates signs for the stored doubles;
// no absolute-area repair or intersection tolerance hides crossed edges.
import orientation from 'robust-orientation';

const turn = (a, b, c) => -orientation([a.x, a.y], [b.x, b.y], [c.x, c.y]);
const midpoint = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
const onBox = (a, b, p) => p.x >= Math.min(a.x, b.x) && p.x <= Math.max(a.x, b.x)
  && p.y >= Math.min(a.y, b.y) && p.y <= Math.max(a.y, b.y);
const opposite = (a, b) => a < 0 && b > 0 || a > 0 && b < 0;
const intersects = (a, b, c, d) => {
  const u = turn(a, b, c), v = turn(a, b, d), w = turn(c, d, a), z = turn(c, d, b);
  if (opposite(u, v) && opposite(w, z)) return 'proper';
  if (u === 0 && onBox(a, b, c) || v === 0 && onBox(a, b, d)
    || w === 0 && onBox(c, d, a) || z === 0 && onBox(c, d, b)) return 'touch-or-overlap';
  return null;
};

export function directPolygonGeometry(p) {
  if (!Array.isArray(p) || p.length < 3 || p.some(a => !Number.isFinite(a?.x) || !Number.isFinite(a?.y)))
    throw new Error('Finite polygon vertices required.');
  let twiceArea = 0, compensation = 0, minCornerSine = Infinity;
  const intersections = [], degenerateEdges = [];
  for (let i = 1; i < p.length - 1; i++) {
    const term = turn(p[0], p[i], p[i + 1]) - compensation, sum = twiceArea + term;
    compensation = (sum - twiceArea) - term; twiceArea = sum;
  }
  for (let i = 0; i < p.length; i++) {
    const a = p[i], b = p[(i + 1) % p.length], c = p[(i + 2) % p.length], t = turn(a, b, c);
    const ab = Math.hypot(b.x - a.x, b.y - a.y), bc = Math.hypot(c.x - b.x, c.y - b.y);
    if (ab === 0 || t === 0 && (b.x - a.x) * (c.x - b.x) + (b.y - a.y) * (c.y - b.y) < 0) degenerateEdges.push(i);
    if (ab * bc > 0) minCornerSine = Math.min(minCornerSine, t / (ab * bc));
    for (let j = i + 2; j < p.length; j++) {
      if (i === 0 && j === p.length - 1) continue;
      const type = intersects(a, b, p[j], p[(j + 1) % p.length]);
      if (type) intersections.push({ edges: [i, j], type });
    }
  }
  return { valid: twiceArea > 0 && degenerateEdges.length === 0 && intersections.length === 0,
    area: twiceArea / 2, minCornerSine, degenerateEdges, intersections };
}

export function directStreamtubeVolumeGeometry(nodes) {
  const create = () => ({ count: 0, minArea: Infinity, failures: [] });
  const primal = create(), halves = create(), dual = create(), concavePrimal = [];
  const edges = new Map();
  const collect = (summary, vertices, location) => {
    const q = directPolygonGeometry(vertices); summary.count++; summary.minArea = Math.min(summary.minArea, q.area);
    if (!q.valid) summary.failures.push({ ...location, ...q, vertices });
    return q;
  };
  const id = p => `${p.x},${p.y}`;
  const edge = (a, b, cell) => {
    const ids = [id(a), id(b)].sort(), key = ids.join('|'), sign = id(a) === ids[0] ? 1 : -1;
    if (edges.has(key)) { edges.get(key).cells.push(cell); edges.get(key).directions.push(sign); return; }
    edges.set(key, { a, b, ids, cells: [cell], directions: [sign], minX: Math.min(a.x, b.x), maxX: Math.max(a.x, b.x),
      minY: Math.min(a.y, b.y), maxY: Math.max(a.y, b.y) });
  };
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length - 1; i++) for (let j = 0; j < nodes[g][i].length - 1; j++) {
    const location = { group: g, interval: i, tube: j };
    const p = [nodes[g][i][j], nodes[g][i + 1][j], nodes[g][i + 1][j + 1], nodes[g][i][j + 1]];
    const q = collect(primal, p, location);
    if (q.minCornerSine < 0) concavePrimal.push({ ...location, ...q, vertices: p });
    for (let k = 0; k < 4; k++) edge(p[k], p[(k + 1) % 4], location);
    const lower = midpoint(p[0], p[1]), upper = midpoint(p[2], p[3]);
    collect(halves, [p[0], lower, upper, p[3]], { ...location, half: 'upstream' });
    collect(halves, [lower, p[1], p[2], upper], { ...location, half: 'downstream' });
    if (i > 0) {
      const beforeLower = midpoint(nodes[g][i - 1][j], p[0]), beforeUpper = midpoint(nodes[g][i - 1][j + 1], p[3]);
      collect(dual, [beforeLower, p[0], lower, upper, p[3], beforeUpper], { group: g, station: i, tube: j });
    }
  }
  // Sweep unique displayed edges. Shared topological endpoints are allowed;
  // repeated edges represent adjacent cells and are collected only once.
  let active = [], candidates = 0; const crossings = [];
  for (const current of [...edges.values()].sort((a, b) => a.minX - b.minX)) {
    active = active.filter(e => e.maxX >= current.minX);
    for (const e of active) {
      if (e.maxY < current.minY || current.maxY < e.minY) continue;
      candidates++;
      const shared = e.ids.find(v => current.ids.includes(v));
      if (shared !== undefined) {
        const a = id(e.a) === shared ? e.a : e.b, b = id(e.a) === shared ? e.b : e.a;
        const c = id(current.a) === shared ? current.b : current.a;
        if (turn(a, b, c) === 0 && (b.x - a.x) * (c.x - a.x) + (b.y - a.y) * (c.y - a.y) > 0)
          crossings.push({ type: 'shared-endpoint-overlap', first: e, second: current });
        continue;
      }
      const type = intersects(e.a, e.b, current.a, current.b);
      if (type) crossings.push({ type, first: e, second: current });
    }
    active.push(current);
  }
  const nonmanifoldEdges = [...edges.values()].filter(e => e.cells.length > 2
    || e.cells.length === 2 && e.directions[0] === e.directions[1]);
  const embedding = { edges: edges.size, testedPairs: candidates, crossings, nonmanifoldEdges };
  const valid = [primal, halves, dual].every(s => s.failures.length === 0) && crossings.length === 0 && nonmanifoldEdges.length === 0;
  return { valid, primal, halves, dual, concavePrimal, embedding };
}
