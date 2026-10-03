// SPDX-License-Identifier: GPL-2.0-or-later
// Prescribe physical x stations BEFORE elliptic relaxation, as ISET's SINL
// and SOUT inputs do. The geometric distribution and multiblock matching here
// are our construction, not a transcription of OUTLIN or a recovered MSET law.
// Each block keeps its anchors/count; the first discrete interval is prescribed.
// Require at least half of double precision in every achieved interval.
const spacingTolerance = Math.sqrt(Number.EPSILON);
export function createPhysicalCutStations({ anchors, edgeSpacing, bodyEnd = 'first' } = {}) {
  if (!Array.isArray(anchors) || anchors.length < 2
    || !anchors.every((a, i) => Number.isSafeInteger(a.index) && a.index >= 0 && Number.isFinite(a.x)
      && (!i || a.index > anchors[i - 1].index && a.x > anchors[i - 1].x))
    || !Number.isFinite(edgeSpacing) || edgeSpacing <= 0 || !['first', 'last'].includes(bodyEnd))
    throw new Error('Physical cut stations require ordered anchors and a positive edge spacing.');
  const firstIndex = anchors[0].index, count = anchors.at(-1).index - firstIndex;
  if (count > 10000) throw new Error('Physical cut station count exceeds 10000 intervals.');
  const ordered = bodyEnd === 'first' ? anchors : anchors.toReversed();
  const direction = bodyEnd === 'first' ? 1 : -1, x = Array(count + 1), blocks = [];
  x[ordered[0].index - firstIndex] = ordered[0].x;
  let first = edgeSpacing;
  for (let b = 1; b < ordered.length; b++) {
    const a = ordered[b - 1], c = ordered[b], n = Math.abs(c.index - a.index), length = Math.abs(c.x - a.x);
    let logRatio = 0;
    if (!Number.isFinite(length) || !(first > 0)) throw new Error('Unresolved physical cut block.');
    if (n === 1) {
      if (Math.abs(length - first) > 64 * Number.EPSILON * Math.max(length, first))
        throw new Error('One-interval cut block cannot match its fixed anchors and adjoining spacing.');
    } else {
      if (first >= length) throw new Error('Cut edge spacing must be smaller than its fixed block extent.');
      // Fit the remaining length: adding a tiny tail to first/length can
      // round to one when first is almost the entire block.
      const tail = length - first, relative = first / tail;
      const logFirst = relative > 0 && Number.isFinite(relative) ? Math.log(relative) : Math.log(first) - Math.log(tail);
      const sum = z => {
        let s = 0;
        for (let k = 1; k < n; k++) { s += Math.exp(logFirst + k * z); if (s > 1) break; }
        return s;
      };
      let lo = -1, hi = Math.max(1, -logFirst / (n - 1));
      while (sum(lo) > 1 && lo > -1024) lo *= 2;
      if (sum(lo) > 1) throw new Error('Physical cut spacing has no representable geometric ratio.');
      for (let k = 0; k < 80; k++) {
        const mid = .5 * (lo + hi);
        if (mid === lo || mid === hi) break;
        if (sum(mid) < 1) lo = mid; else hi = mid;
      }
      logRatio = .5 * (lo + hi);
    }
    let previous = a.x, last, achievedFirstSpacing;
    for (let k = 0; k < n; k++) {
      const requested = !k ? first : Math.exp(Math.log(first) + k * logRatio);
      const next = k === n - 1 ? c.x : previous + direction * requested;
      last = direction * (next - previous);
      if (!(last > 0) || !Number.isFinite(next)) throw new Error('Physical cut stations collapse at floating-point resolution.');
      if (Math.abs(last / requested - 1) > spacingTolerance)
        throw new Error('Physical cut interval accuracy is unresolved at floating-point resolution.');
      if (!k) achievedFirstSpacing = last;
      x[a.index + direction * (k + 1) - firstIndex] = next;
      previous = next;
    }
    const growth = Math.exp(logRatio);
    if (!(growth > 0) || !Number.isFinite(growth)) throw new Error('Physical cut growth is not representable.');
    blocks.push({ fromIndex: a.index, toIndex: c.index, intervals: n, length,
      firstSpacing: first, achievedFirstSpacing, lastSpacing: last, growth });
    // Use the achieved discrete interval at the anchor, not a derivative of
    // a continuous rank map. Independent block counts therefore cannot add
    // another scaling jump here.
    first = last;
  }
  return { x, firstIndex, blocks, metric: 'physical x', bodyEnd, edgeSpacing, relativeSpacingTolerance: spacingTolerance };
}
