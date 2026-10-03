import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodyLayout, allocateStreamtubeMasses, streamtubeEquationOrder } from '../src/euler/streamtube-body-layout.js';

test('intrinsic body topology has a square multielement count including capture, stagnation and farfield constraints', () => {
  for (let elements = 1; elements <= 4; elements++) for (const rounded of [true, false]) {
    const nx = 20, tubes = Array.from({ length: elements + 1 }, (_, g) => g + 2);
    const bodies = Array.from({ length: elements }, (_, b) => ({ leadingIndex: 2 + b, trailingIndex: 14 + b, roundedLeadingEdge: rounded || b % 2 === 0 }));
    const l = createStreamtubeBodyLayout({ segments: nx, tubes, bodies, primaryBody: elements - 1 });
    const nt = tubes.reduce((a, b) => a + b), roundCount = bodies.filter(b => b.roundedLeadingEdge).length;
    assert.equal(l.densityCount, nx * nt);
    assert.equal(l.positionCount, (nx + 1) * (nt + 1) - bodies.reduce((s, b) => s + b.trailingIndex - b.leadingIndex + 1, 0));
    assert.equal(l.globalCount, elements + roundCount + 3);
    assert.equal(l.rows.length, l.n);
    assert.equal(l.rowCounts.trailingKutta, elements);
    assert.equal(l.rowCounts.leadingKutta, roundCount);
    assert.equal(l.rowCounts.farfieldMatch, 3);
    assert.equal(l.rowCounts.endTangency, 2 * (nt + 1));
    assert.equal(l.globals.capture.filter(c => c !== null).length, elements - 1);
    assert.equal(new Set(l.positions.map(p => p.column)).size, l.positionCount);
    assert.deepEqual(l.rows.map(r => r.index), Array.from({ length: l.n }, (_, i) => i));
  }
});

test('cut geometry is identified once while physical walls and Kutta pressure rows remain distinct', () => {
  const l = createStreamtubeBodyLayout({ segments: 12, tubes: [2, 3, 2], bodies: [{ leadingIndex: 2, trailingIndex: 6 }, { leadingIndex: 5, trailingIndex: 9 }] });
  for (let b = 0; b < l.elements; b++) for (let i = 0; i <= l.nx; i++) {
    const lower = l.nodes[b][i].at(-1), upper = l.nodes[b + 1][i][0];
    if (l.active(b, i)) {
      assert.equal(lower.kind, 'wall'); assert.equal(upper.kind, 'wall');
      assert.equal(lower.column, null); assert.equal(upper.column, null);
      assert.equal(lower.side, 'lower'); assert.equal(upper.side, 'upper');
      assert.equal(l.rows.filter(r => r.kind === 'cutPressure' && r.body === b && r.i === i).length, 0);
    } else {
      assert.equal(lower, upper); assert.equal(lower.kind, 'cut');
      if (i > 0 && i < l.nx) assert.equal(l.rows.filter(r => r.kind === 'cutPressure' && r.body === b && r.i === i).length, 1);
    }
  }
  assert.throws(() => createStreamtubeBodyLayout({ segments: 12, tubes: [2, 2], bodies: [{ leadingIndex: 0, trailingIndex: 6 }] }), /strictly inside/);
  assert.throws(() => createStreamtubeBodyLayout({ segments: 12, tubes: [2, 0], bodies: [{ leadingIndex: 2, trailingIndex: 6 }] }), /dimensions/);
});

test('captured-flow changes redistribute neighboring groups conservatively with exact derivatives', () => {
  const levels = [-2, -.9, .2, .8, 2], weights = [[1, 2], [1, 3, 2], [2, 1], [1, 2, 4]], primary = 1;
  const allocation = allocateStreamtubeMasses(levels, weights, primary), states = allocation.groups.flat();
  assert.ok(Math.abs(states.reduce((s, v) => s + v.massFlow, 0) - 4) < 1e-14);
  for (const body of [0, 2]) {
    const h = 1e-6, plus = levels.slice(), minus = levels.slice(); plus[body + 1] += h; minus[body + 1] -= h;
    const p = allocateStreamtubeMasses(plus, weights, primary).groups.flat(), m = allocateStreamtubeMasses(minus, weights, primary).groups.flat();
    let sum = 0;
    states.forEach((s, i) => {
      const derivative = s.derivatives.get(body) ?? 0; sum += derivative;
      assert.ok(Math.abs(derivative - (p[i].massFlow - m[i].massFlow) / (2 * h)) < 2e-10);
    });
    assert.ok(Math.abs(sum) < 1e-14);
  }
  assert.ok(states.every(s => !s.derivatives.has(primary)));
  const shifted = allocateStreamtubeMasses(levels.map(v => v + 3), weights, primary).groups.flat();
  states.forEach((s, i) => assert.ok(Math.abs(s.massFlow - shifted[i].massFlow) < 1e-14));
  assert.throws(() => allocateStreamtubeMasses([0, 1, .9], [[1], [1]]), /levels/);
  assert.throws(() => allocateStreamtubeMasses([0, 1, 2], [[1], [0]]), /weights/);
});

// Include sharp noses, nonprimary capture variables, and independent wake banks.
test('Euler equation alignment is a bijection across supported body topologies', () => {
  for (const densityUnknowns of [false, true]) for (const independentWakeBanks of [false, true])
    for (const elements of [1, 2, 3]) {
      const layout = createStreamtubeBodyLayout({ segments: 20, tubes: Array(elements + 1).fill(3),
        bodies: Array.from({ length: elements }, (_, i) => ({ leadingIndex: 2+i, trailingIndex: 14+i, roundedLeadingEdge: i%2 === 0 })),
        primaryBody: elements-1, densityUnknowns, independentWakeBanks });
      const order = streamtubeEquationOrder(layout);
      assert.deepEqual(Array.from(order).sort((a,b) => a-b), Array.from({length:layout.n}, (_,i) => i));
      for (const row of layout.rows.filter(r => r.kind === 'internalPressure'))
        assert.equal(order[layout.nodes[row.group][row.i][row.j].column], row.index);
    }
});
