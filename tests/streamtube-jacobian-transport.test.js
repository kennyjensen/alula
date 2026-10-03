// SPDX-License-Identifier: GPL-2.0-or-later
// Local real transport chains only; no mesh generation, global Jacobian or LU.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeJacobianTransport } from '../src/euler/streamtube-jacobian-transport.js';
import { prepareStreamtubeTransportChain } from '../src/euler/streamtube-transport-chain.js';
import { createStreamtubeBodyLayout, allocateStreamtubeMasses } from '../src/euler/streamtube-body-layout.js';

const zero = Object.freeze({ x: 0, y: 0 });
function fixture(regime, mucon = .75) {
  const nx = 12, tubes = [2, 1, 2], layout = createStreamtubeBodyLayout({ segments: nx, tubes,
    bodies: [{ leadingIndex: 3, trailingIndex: 9 }, { leadingIndex: 2, trailingIndex: 5 }], independentWakeBanks: true });
  const massScale = .73, allocation = allocateStreamtubeMasses([-1, 0, .6, 2], [[2, 1], [1], [1, 3]]);
  const massDerivatives = allocation.groups.map(group => group.map(m => new Map([...m.derivatives]
    .map(([body, d]) => [layout.globals.capture[body], d * massScale]))));
  const geometry = layout.nodes.map(group => group.map((row, k) => row.map(node => {
    const col = node.column ?? layout.n + 2 * node.body + (node.side === 'upper' ? 1 : 0);
    return new Map([[col, { x: .02 * Math.sin(k + col), y: .03 * Math.cos(k + col) }],
      [layout.globals.stagnation[0], { x: .003 * Math.cos(k), y: .002 * Math.sin(k) }]]);
  })));
  const x = Array.from({ length: nx + 1 }, (_, k) => .1 * k + .005 * k * k), width = .42;
  const speeds = Array.from({ length: nx }, (_, k) => regime === 'subsonic' ? .4 + .03 * Math.cos(k)
    : regime === 'supersonic' ? 1.3 + .04 * Math.cos(k) : .98 + .3 * Math.cos(k));
  const chains = tubes.map((count, g) => Array.from({ length: count }, (_, j) => {
    const massFlow = allocation.groups[g][j].massFlow * massScale;
    return prepareStreamtubeTransportChain({ lower: x.map(x => ({ x, y: g + j * width })),
      upper: x.map(x => ({ x, y: g + (j + 1) * width })), densities: speeds.map(q => massFlow / (width * q)),
      massFlow, stagnationEnthalpy: 2.7, gamma: 1.4,
      upwind: { mucon, mcrit: .9, boundary: { kind: 'unfiltered-first-two' } } }, { linearize: true });
  }));
  return { chains, layout, geometry, massDerivatives };
}
// Literal pre-cache per-cell callback, kept independently of the new helper.
function original(f, i, g, j, col, counts) {
  const chain = f.chains[g][j], sections = new Map(), dm = f.massDerivatives[g][j].get(col) ?? 0;
  const at = k => {
    if (!sections.has(k)) {
      counts.sections++;
      sections.set(k, chain.sectionTangent(k, {
        lower: [f.geometry[g][k][j].get(col) ?? zero, f.geometry[g][k + 1][j].get(col) ?? zero],
        upper: [f.geometry[g][k][j + 1].get(col) ?? zero, f.geometry[g][k + 1][j + 1].get(col) ?? zero],
        density: f.layout.densityIndex(k, g, j) === col ? chain.sections[k].rho : 0, massFlow: dm,
      }));
    }
    return sections.get(k);
  };
  counts.speeds += 2;
  return [chain.transportTangent(i - 1, at), chain.transportTangent(i, at)];
}
function columns(f, i, g, j) {
  const result = new Set([...f.massDerivatives[g][j].keys(), f.layout.n + 500]);
  for (let k = Math.max(0, i - 3); k <= i; k++) result.add(f.layout.densityIndex(k, g, j));
  for (let k = Math.max(0, i - 3); k <= i + 1; k++) for (const side of [j, j + 1])
    for (const col of f.geometry[g][k][side].keys()) result.add(col);
  return [...result];
}

test('rolling section/speed reuse exactly preserves subsonic, transonic and supersonic derivatives in every passage', t => {
  const counts = { sections: 0, speeds: 0 }, cached = { sections: 0, speeds: 0 }; let comparisons = 0;
  for (const regime of ['subsonic', 'transonic', 'supersonic']) for (const mucon of [-.75, .75]) {
    const f = fixture(regime, mucon), before = structuredClone({ geometry: f.geometry, massDerivatives: f.massDerivatives });
    const machSquared = f.chains.flat().flatMap(c => c.sections.map(s => s.machSquared));
    assert.equal(machSquared.some(m => m > 1), regime !== 'subsonic');
    assert.equal(machSquared.some(m => m < 1), regime !== 'supersonic');
    let currentCol; const seenSections = new Map(), seenSpeeds = new Map();
    const record = (seen, g, j, k) => { const key = `${g}/${j}/${k}/${currentCol}`; seen.set(key, (seen.get(key) ?? 0) + 1); };
    const chains = f.chains.map((group, g) => group.map((chain, j) => ({ ...chain,
      sectionTangent(k, d) { cached.sections++; record(seenSections, g, j, k); return chain.sectionTangent(k, d); },
      transportTangent(k, at) { cached.speeds++; record(seenSpeeds, g, j, k); return chain.transportTangent(k, at); },
    })));
    const cache = createStreamtubeJacobianTransport({ ...f, chains });
    for (let i = 1; i < f.layout.nx; i++) for (const g of [2, 0, 1]) for (let j = f.layout.tubes[g] - 1; j >= 0; j--)
      for (const col of columns(f, i, g, j)) {
        currentCol = col;
        const expected = original(f, i, g, j, col, counts), result = cache.forCell(i, g, j, col);
        assert.deepEqual(result, expected, `${regime}/${mucon}: ${i},${g},${j},${col}`); comparisons += 2;
        result[0] = Infinity; // Returned pairs cannot mutate the scalar cache.
        assert.deepEqual(cache.forCell(i, g, j, col), expected);
      }
    assert.ok([...seenSections.values(), ...seenSpeeds.values()].every(n => n === 1));
    assert.deepEqual({ geometry: f.geometry, massDerivatives: f.massDerivatives }, before);
  }
  assert.ok(cached.sections < counts.sections); assert.ok(cached.speeds < counts.speeds);
  t.diagnostic(JSON.stringify({ exactScalarComparisons: comparisons, originalKernelCalls: counts, cachedKernelCalls: cached }));
});

test('explicit zero and signed zero survive cache hits and physical tube mass derivatives remain separate', () => {
  let calls = 0;
  const nx = 6, layout = { nx, tubes: [2], densityIndex: (k, g, j) => 2 * k + j };
  const geometry = [Array.from({ length: nx + 1 }, () => [new Map(), new Map(), new Map()])];
  const makeChain = () => ({ sections: Array.from({ length: nx }, () => ({ rho: 1 })),
    sectionTangent(k, { massFlow }) { calls++; return { q: massFlow }; }, transportTangent: (k, at) => at(k).q });
  const f = { layout, geometry, chains: [[makeChain(), makeChain()]],
    massDerivatives: [[new Map([[20, -0], [21, 0], [22, -.7]]), new Map([[20, 0], [21, -0], [22, .3]])]] };
  const cache = createStreamtubeJacobianTransport(f);
  for (let i = 1; i < nx; i++) for (const j of [1, 0]) for (const col of [20, 21, 22, 99]) {
    const expected = f.massDerivatives[0][j].get(col) ?? 0, pair = cache.forCell(i, 0, j, col), count = calls;
    assert.ok(pair.every(v => Object.is(v, expected)));
    assert.ok(cache.forCell(i, 0, j, col).every(v => Object.is(v, expected))); assert.equal(calls, count);
  }
  const count = calls, next = createStreamtubeJacobianTransport(f);
  assert.deepEqual(next.forCell(1, 0, 0, 20), [-0, -0]); assert.equal(calls, count + 2);
});

test('backward and invalid queries reject before corrupting a current window; tubes advance independently', () => {
  const f = fixture('transonic'), cache = createStreamtubeJacobianTransport(f), col = f.layout.globals.capture[1];
  const expected = cache.forCell(5, 0, 1, col);
  assert.throws(() => cache.forCell(4, 0, 1, col), /advance monotonically/);
  for (const query of [[0, 0, 1, col], [12, 0, 1, col], [5.5, 0, 1, col], [5, -1, 1, col],
    [5, 3, 1, col], [5, 0, 2, col], [5, 0, 1, -1], [5, 0, 1, Infinity]])
    assert.throws(() => cache.forCell(...query), /topology or query/);
  assert.deepEqual(cache.forCell(5, 0, 1, col), expected);
  assert.deepEqual(cache.forCell(1, 0, 0, col), original(f, 1, 0, 0, col, { sections: 0, speeds: 0 }));
  assert.deepEqual(cache.forCell(11, 0, 1, col), original(f, 11, 0, 1, col, { sections: 0, speeds: 0 }));
  assert.throws(() => cache.forCell(5, 0, 1, col), /advance monotonically/);
  for (const bad of [{}, { ...f, geometry: [] }, { ...f, massDerivatives: [] },
    { ...f, layout: { ...f.layout, densityIndex: undefined } }])
    assert.throws(() => createStreamtubeJacobianTransport(bad), /topology or query/);
});

test('kernel failures propagate unchanged and a failed scalar is never cached', () => {
  const f = fixture('subsonic'), chain = f.chains[0][0], error = Object.assign(new Error('Controlled derivative rejection'), { code: 'TEST_DOMAIN' });
  let failures = 1;
  f.chains[0][0] = { ...chain, sectionTangent(k, d) { if (failures-- > 0) throw error; return chain.sectionTangent(k, d); } };
  const cache = createStreamtubeJacobianTransport(f), col = f.layout.globals.stagnation[0];
  assert.throws(() => cache.forCell(1, 0, 0, col), e => e === error);
  assert.deepEqual(cache.forCell(1, 0, 0, col), original(f, 1, 0, 0, col, { sections: 0, speeds: 0 }));
  const invalidStencil = fixture('subsonic');
  invalidStencil.chains[0][0].transportTangent = (k, at) => at(k + 2).q;
  assert.throws(() => createStreamtubeJacobianTransport(invalidStencil).forCell(1, 0, 0, col), /topology or query/);
});
