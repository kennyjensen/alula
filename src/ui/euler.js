// SPDX-License-Identifier: GPL-2.0-or-later
import { buildEulerLabCase } from '../euler/cases.js';

const $ = id => document.getElementById(id);
let problem; let result; let liveGrid; let worker; let job = 0; let history = []; let lastControls;
const field = $('field'); const graph = $('history');
const sci = value => Number.isFinite(value) ? value.toExponential(2) : '—';
const colors = ['#89e5d0', '#d6ba7d', '#9ab9ee'];
const controls = () => ({ configuration: $('configuration').value, panels: +$('panels').value, rows: +$('rows').value,
  padding: +$('padding').value, mach: +$('mach').value, alpha: +$('alpha').value, maxIterations: +$('iterations').value });
const status = (text, state = '') => { $('status').textContent = text; $('status').className = `status-pill ${state}`; };
const error = message => { $('error').textContent = message ?? ''; $('error').hidden = !message; };
function stop() { worker?.terminate(); worker = null; job++; $('cancel').hidden = true; }
function resetMetrics() {
  for (const id of ['residual', 'mass', 'max-mach', 'cancellation', 'leakage', 'entropy', 'enthalpy']) $(id).textContent = '—';
  $('streamline').textContent = problem?.moving ? '—' : 'Channel case only';
}
function prepare() {
  stop(); result = null; liveGrid = null; history = []; $('export').disabled = true; error(null);
  const moving = $('configuration').value === 'channel';
  $('alpha').disabled = moving; $('padding').disabled = moving;
  if (moving) $('alpha').value = '0';
  try {
    problem = buildEulerLabCase(controls()); $('solve').disabled = false;
    $('title').textContent = problem.title; $('cells').textContent = problem.mesh.cells.length;
    const extra = moving ? problem.mesh.streamlineNodes.length : 0;
    $('unknowns').textContent = `${4 * problem.mesh.cells.length + extra} simultaneous unknowns${moving ? ' · flow + grid' : ''}`;
    $('geometry-note').textContent = moving ? 'Curved channel, 8% contraction. Internal streamline locations move during Newton iterations.'
      : 'NACA 0012 elements. Flap: 0.4 chord, LE (1.05, −0.20), 10° down. Slat: 0.25 chord, LE (−0.40, 0.22).';
    $('run-meta').textContent = moving ? 'Conservative flow + moving internal streamlines.' : 'All elements share one conservative flow system.';
    status('Ready', 'stale');
  } catch (e) { problem = null; $('cells').textContent = '—'; $('solve').disabled = true; error(e.message); status('Invalid case', 'failed'); }
  resetMetrics(); draw();
}

function solve() {
  if (!problem) return;
  stop(); result = null; liveGrid = null; history = []; resetMetrics(); error(null); $('export').disabled = true;
  lastControls = structuredClone(problem.controls); $('solve').disabled = true; $('cancel').hidden = false;
  status('Solving'); $('run-meta').textContent = 'Assembling the coupled Newton system…';
  const id = job; worker = new Worker(new URL('../worker/euler.js', import.meta.url), { type: 'module' });
  worker.onmessage = ({ data }) => {
    if (data.id !== job) return;
    if (data.type === 'mesh') {
      liveGrid = { mesh: data.mesh, iteration: data.iteration };
      drawField(); return;
    }
    if (data.type === 'iteration') {
      history.push(data.iteration); $('residual').textContent = sci(data.iteration.residual);
      $('run-meta').textContent = `Newton iteration ${data.iteration.iteration} · residual ${sci(data.iteration.residual)}`;
      drawHistory(); return;
    }
    stop(); $('solve').disabled = false;
    if (data.type === 'error') { status('Failed', 'failed'); error(data.message); draw(); return; }
    history = data.result.history;
    if (!data.result.converged) {
      status('Not converged', 'failed'); error(`Euler solve stopped: ${data.result.reason}. ${liveGrid ? 'The last accepted grid is shown. ' : ''}No converged field is available.`);
      $('run-meta').textContent = `${history.length - 1} Newton iterations · ${(data.elapsed / 1000).toFixed(2)} s`;
      draw(); return;
    }
    result = data.result; const d = result.diagnostics;
    status('Equations converged'); $('run-meta').textContent = `${history.length - 1} Newton iterations · ${(data.elapsed / 1000).toFixed(2)} s · accuracy: verification only`;
    $('residual').textContent = sci(d.residual); $('mass').textContent = sci(d.relativeMassImbalance);
    $('max-mach').textContent = d.maxMach.toFixed(3); $('cancellation').textContent = sci(d.sharedFluxCancellation);
    $('leakage').textContent = sci(d.wallLeakage); $('entropy').textContent = sci(d.maxEntropyError);
    $('enthalpy').textContent = sci(d.maxTotalEnthalpyError);
    $('streamline').textContent = d.streamlineMassResidual === undefined ? 'Channel case only' : sci(d.streamlineMassResidual);
    $('export').disabled = false; draw();
  };
  worker.onerror = e => { if (id !== job) return; stop(); $('solve').disabled = false; status('Failed', 'failed'); error(e.message || 'Euler worker failed.'); draw(); };
  worker.postMessage({ id, controls: lastControls }); draw();
}

function context(canvas) {
  const box = canvas.getBoundingClientRect(); const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.round(box.width * ratio); canvas.height = Math.round(box.height * ratio);
  const ctx = canvas.getContext('2d'); ctx.scale(ratio, ratio); ctx.clearRect(0, 0, box.width, box.height);
  return { ctx, width: box.width, height: box.height };
}
const valueAt = (s, quantity) => quantity === 'mach' ? Math.hypot(s.u, s.v) / Math.sqrt(1.4 * s.p / s.rho)
  : quantity === 'cp' ? 2 * (s.p - result.reference.p) : s.rho;
function color(value) {
  const stops = [[24, 46, 81], [42, 123, 139], [137, 229, 208], [236, 196, 111]];
  const t = Math.max(0, Math.min(1, value)) * 3; const i = Math.min(2, Math.floor(t)); const f = t - i;
  return `rgb(${stops[i].map((v, k) => Math.round(v + f * (stops[i + 1][k] - v))).join(',')})`;
}
function drawField() {
  const { ctx, width, height } = context(field);
  if (!problem) return;
  const mesh = result?.mesh ?? liveGrid?.mesh ?? problem.mesh;
  const focus = $('near-body').checked && problem.elements.length;
  const points = focus ? problem.elements.flatMap(e => e.points) : mesh.vertices;
  let xmin = Math.min(...points.map(p => p.x)); let xmax = Math.max(...points.map(p => p.x));
  let ymin = Math.min(...points.map(p => p.y)); let ymax = Math.max(...points.map(p => p.y));
  const pad = focus ? .16 * (xmax - xmin) : .03 * (xmax - xmin);
  xmin -= pad; xmax += pad; ymin -= pad; ymax += pad;
  const scale = Math.min((width - 36) / (xmax - xmin), (height - 44) / (ymax - ymin));
  const xy = p => [width / 2 + (p.x - (xmin + xmax) / 2) * scale, height / 2 - (p.y - (ymin + ymax) / 2) * scale];
  const path = points => { ctx.beginPath(); points.forEach((p, i) => { const [x, y] = xy(p); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); };
  const quantity = $('quantity').value;
  const values = result ? result.states.map(s => valueAt(s, quantity)) : [];
  const low = Math.min(...values); const high = Math.max(...values);
  for (let i = 0; i < mesh.cells.length; i++) {
    path(mesh.cells[i].vertices.map(id => mesh.vertices[id])); ctx.closePath();
    ctx.fillStyle = result ? color((values[i] - low) / Math.max(high - low, 1e-12)) : '#13252e'; ctx.fill();
  }
  if ($('show-grid').checked) {
    ctx.strokeStyle = result ? '#07172066' : '#41616e88'; ctx.lineWidth = .6;
    for (const f of mesh.faces) { path([mesh.vertices[f.a], mesh.vertices[f.b]]); ctx.stroke(); }
  }
  if ($('show-cuts').checked) {
    ctx.strokeStyle = '#f2d594'; ctx.lineWidth = 1.2; ctx.setLineDash([5, 5]);
    for (const cut of problem.mesh.cuts ?? []) { path([mesh.vertices[cut.a], mesh.vertices[cut.b]]); ctx.stroke(); }
    ctx.setLineDash([]);
  }
  problem.elements.forEach((e, i) => { path(e.points); ctx.closePath(); ctx.fillStyle = '#0c141b'; ctx.fill(); ctx.strokeStyle = colors[i]; ctx.lineWidth = 1.7; ctx.stroke(); });
  if (problem.moving) {
    ctx.strokeStyle = '#d6ba7d'; ctx.lineWidth = 2;
    for (const f of mesh.faces.filter(f => f.boundary?.type === 'wall')) { path([mesh.vertices[f.a], mesh.vertices[f.b]]); ctx.stroke(); }
  }
  $('field-label').textContent = result ? `${quantity === 'cp' ? 'Cp' : quantity === 'rho' ? 'ρ / ρ∞' : 'Mach'} · cell averages${problem.moving ? ' · solved streamlines' : ''}`
    : liveGrid ? `Moving grid · Newton iteration ${liveGrid.iteration.iteration} · unconverged` : 'Connected grid · preview';
  $('field-range').textContent = result ? `${low.toFixed(4)} → ${high.toFixed(4)} · blue → gold` : 'No solved field';
  $('grid-meta').textContent = `Minimum area ${sci(mesh.diagnostics.minArea)}${problem.mesh.cuts ? ` · ${problem.mesh.cuts.length} shared cuts` : ''}`
    + (liveGrid && !result ? ` · max node move ${sci(liveGrid.iteration.maximumNodeMovement)}` : '');
}
function drawHistory() {
  const { ctx, width, height } = context(graph);
  const left = 46; const right = width - 22; const top = 25; const bottom = height - 31;
  const min = Math.min(-10, ...history.map(h => Math.floor(Math.log10(Math.max(h.residual, 1e-16)))));
  const max = Math.max(0, ...history.map(h => Math.ceil(Math.log10(Math.max(h.residual, 1e-16)))));
  const yy = value => top + (max - value) / (max - min) * (bottom - top);
  ctx.font = '10px Consolas, monospace'; ctx.lineWidth = 1;
  for (let log = max; log >= min; log -= 2) {
    ctx.fillStyle = '#8299a6'; ctx.fillText(`${log}`, 14, yy(log) + 3);
    ctx.strokeStyle = '#26353e'; ctx.beginPath(); ctx.moveTo(left, yy(log)); ctx.lineTo(right, yy(log)); ctx.stroke();
  }
  ctx.strokeStyle = '#d6ba7d88'; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(left, yy(-9)); ctx.lineTo(right, yy(-9)); ctx.stroke(); ctx.setLineDash([]);
  ctx.fillStyle = '#8299a6'; ctx.fillText('iteration', right - 55, height - 8);
  if (!history.length) { ctx.fillText('Run a case to inspect convergence.', left + 12, height / 2); return; }
  const xx = i => left + i / Math.max(1, history.length - 1) * (right - left);
  ctx.beginPath(); history.forEach((h, i) => { const y = yy(Math.log10(Math.max(h.residual, 1e-16))); i ? ctx.lineTo(xx(i), y) : ctx.moveTo(xx(i), y); });
  ctx.strokeStyle = '#89e5d0'; ctx.lineWidth = 1.8; ctx.stroke();
  history.forEach((h, i) => { ctx.beginPath(); ctx.arc(xx(i), yy(Math.log10(Math.max(h.residual, 1e-16))), 2.7, 0, Math.PI * 2); ctx.fillStyle = '#89e5d0'; ctx.fill(); ctx.fillStyle = '#8299a6'; ctx.fillText(String(h.iteration), xx(i) - 3, bottom + 14); });
}
function draw() { drawField(); drawHistory(); }
for (const id of ['configuration', 'panels', 'rows', 'padding', 'mach', 'alpha', 'iterations']) $(id).addEventListener('input', prepare);
for (const id of ['show-grid', 'show-cuts', 'near-body', 'quantity']) $(id).addEventListener('change', drawField);
$('solve').addEventListener('click', solve);
$('cancel').addEventListener('click', () => { stop(); $('solve').disabled = false; status('Cancelled', 'stale'); $('run-meta').textContent = 'Calculation stopped. Run again to solve.'; });
$('export').addEventListener('click', () => {
  if (!result?.converged) return;
  const data = { schemaVersion: 1, controls: lastControls, warning: 'Euler verification state; no BL coupling or released airfoil force predictions.', result };
  const blob = new Blob([JSON.stringify(data, (_, value) => ArrayBuffer.isView(value) ? [...value] : value, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob); const link = document.createElement('a'); link.href = url; link.download = 'alula-euler.json'; link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
});
new ResizeObserver(draw).observe(field.parentElement);
prepare();
