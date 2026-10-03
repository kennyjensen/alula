// SPDX-License-Identifier: GPL-2.0-or-later
// One instance per relaxation: observers may run nested solves safely.
const url = new URL('./wasm/smoothing.wasm', import.meta.url);
const bytes = url.protocol === 'file:' ? await (await import('node:fs/promises')).readFile(url)
  : await fetch(url).then(response => { if (!response.ok) throw new Error(`Cannot load grid smoother: HTTP ${response.status}`); return response.arrayBuffer(); });
const module = await WebAssembly.compile(bytes);
const rejection = code => code === 3 ? 'folded grid' : code === 1
  ? 'No resolved positive harmonic-mass normal-speed branch at this boundary.' : 'Nonfinite elliptic metric or control.';
export function createWasmSmoother(system) {
  const config = system.wasmSmoothingConfig;
  if (!config) return null;
  const { nx, nt, rowGroups, boundary, xi, eta, lengthScale } = config;
  const wasm = new WebAssembly.Instance(module, { env: { abort() { throw new Error('WASM smoothing memory or bounds failure.'); } } }).exports;
  wasm.allocate(nx, nt, rowGroups.length, Math.max(...rowGroups.map(row => row.length)), config.mode, config.streamwiseSourceDiscretization === 'grape-1980', lengthScale);
  const view = (id, size) => new Float64Array(wasm.memory.buffer, wasm.pointer(id), size);
  const put = (id, array) => view(id, array.length).set(array);
  put(1, xi); put(2, eta); put(3, boundary.background.flat());
  put(4, eta.flatMap(v => boundary.decay ? [
    Math.exp(-boundary.decay.lower * v) * Math.expm1(-(boundary.decay.lower + boundary.decay.upper) * (1 - v)) / Math.expm1(-boundary.decay.lower - boundary.decay.upper),
    Math.exp(-boundary.decay.upper * (1 - v)) * Math.expm1(-(boundary.decay.lower + boundary.decay.upper) * v) / Math.expm1(-boundary.decay.lower - boundary.decay.upper),
  ] : [1 - v, v]));
  put(5, ['lower', 'upper'].flatMap(side => Array.from({ length: nx + 1 }, (_, i) =>
    Number((boundary.sides ?? ['lower', 'upper']).includes(side) && i > 0 && i < nx && (boundary.activeStations?.[side]?.[i] ?? true)))));
  put(6, ['lower', 'upper'].flatMap(side => Array.from({ length: nx + 1 }, (_, i) => Number(boundary.corners?.[side]?.includes(i) ?? false))));
  put(7, config.xiDerivatives.flatMap(d => d ? [...d.first, ...d.second] : Array(6).fill(0)));
  put(8, config.derivatives.flatMap(d => d ? [...d.first, ...d.second] : Array(6).fill(0)));
  put(9, rowGroups.flatMap(row => [row[0], row.length]));
  let replays = 0;
  const load = nodes => {
    // Always copy: public observer/caller objects are mutable.
    const target = view(0, 2 * (nx + 1) * (nt + 1)); let k = 0;
    for (const row of nodes) for (const p of row) { target[k++] = p.x; target[k++] = p.y; }
  };
  const read = () => {
    const data = view(0, 2 * (nx + 1) * (nt + 1)); let k = 0;
    return Array.from({ length: nx + 1 }, () => Array.from({ length: nt + 1 }, () => ({ x: data[k++], y: data[k++] })));
  };
  return {
    get replays() { return replays; },
    evaluate(nodes) {
      load(nodes);
      if (!wasm.evaluate()) {
        replays++;
        // Re-evaluate in the reference to preserve detailed domain diagnostics.
        const quality = system.quality(nodes);
        if (!quality.valid && !config.scalar) return { quality, residual: Infinity, merit: Infinity };
        const state = system.residuals(nodes);
        return { quality, residual: state.residual, merit: .5 * state.rows.reduce((s, p) => s + p.x ** 2 + p.y ** 2, 0) };
      }
      return { quality: { valid: true, invalidCells: [], minCornerSine: wasm.minSine.value }, residual: wasm.residual.value, merit: wasm.merit.value };
    },
    sweep(nodes, omega, { lineTolerance = 0 } = {}) {
      load(nodes);
      if (!(config.scalar ? wasm.scalarSweep(omega) : wasm.sweep(omega, lineTolerance))) {
        // A failure must retain the reference's typed termination and detailed
        // branch/linear diagnostics. Never commit a partially completed sweep.
        replays++;
        return system.sweep(nodes, omega, { lineTolerance });
      }
      if (config.scalar) { return { nodes: read(), maxUpdate: wasm.maxUpdate.value }; }
      const events = view(10, 5 * wasm.eventCount.value), summaries = view(11, 6 * rowGroups.length);
      const rowSteps = rowGroups.map((rows, g) => {
        if (summaries[6 * g + 4]) return { rows: rows.slice(), skipped: true, residual: summaries[6 * g + 5], tolerance: lineTolerance };
        const trials = [];
        for (let e = 0; e < events.length; e += 5) if (events[e] === g)
          trials.push({ fraction: events[e + 1], merit: events[e + 2], rejected: events[e + 3] ? rejection(events[e + 3]) : null, accepted: Boolean(events[e + 4]) });
        return { rows: rows.slice(), baseMerit: summaries[6 * g], merit: summaries[6 * g + 1], fraction: summaries[6 * g + 2], minCornerSine: summaries[6 * g + 3], trials };
      });
      return { nodes: read(), maxUpdate: wasm.maxUpdate.value, rowSteps };
    },
  };
}
