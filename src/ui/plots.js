// SPDX-License-Identifier: GPL-2.0-or-later
import { flowSpeedColor, flowSpeedMaximum, flowChangeColor, flowChangeMaximum } from './flow-colors.js';
import { machContours } from './mach-contours.js';
const machPaths = new WeakMap();
export const colors = ['#89e5d0', '#d6ba7d', '#9ab9ee', '#d5a6ca', '#c1d892', '#edaa8c'];
const mono = '10px Consolas, "Liberation Mono", monospace';
const meshPaths = new WeakMap();
const flowPaths = new WeakMap();
function coloredFlowPaths(flow, mode) {
  if (!flowPaths.has(flow)) flowPaths.set(flow, new Map());
  const cached = flowPaths.get(flow), changes = mode === 'change';
  if (!cached.has(mode)) {
    // Batch segments by color: only 64 canvas strokes per snapshot, also
    // cached for pan/zoom. Never interpolate speeds between different tubes.
    const count = changes ? 129 : 64;
    const paths = Array.from({ length: count }, (_, i) => ({ path: new Path2D(), color: changes
      ? flowChangeColor((i / (count - 1) * 2 - 1) * flowChangeMaximum) : flowSpeedColor(i / (count - 1) * flowSpeedMaximum) }));
    for (const line of flow.lines) line.speedRatios.forEach((speed, i) => {
      const fraction = changes ? .5 + (line.speedChanges?.[i] ?? 0) / (2 * flowChangeMaximum) : speed / flowSpeedMaximum;
      const { path } = paths[Math.round(Math.max(0, Math.min(1, fraction)) * (count - 1))];
      path.moveTo(line.points[i].x, line.points[i].y); path.lineTo(line.points[i + 1].x, line.points[i + 1].y);
    });
    cached.set(mode, paths);
  }
  return cached.get(mode);
}
function meshPath(mesh) {
  if (!meshPaths.has(mesh)) {
    const path = new Path2D(), edges = new Set();
    if (mesh.edgeCurves) for (const edge of mesh.edgeCurves) {
      const a = mesh.vertices[edge.a], b = mesh.vertices[edge.b];
      path.moveTo(a.x, a.y); path.bezierCurveTo(edge.c1.x, edge.c1.y, edge.c2.x, edge.c2.y, b.x, b.y);
    }
    else for (const cell of mesh.cells) for (let j = 0; j < cell.length; j++) {
      const a = cell[j], b = cell[(j + 1) % cell.length], key = a < b ? `${a}:${b}` : `${b}:${a}`;
      if (edges.has(key)) continue; edges.add(key);
      path.moveTo(mesh.vertices[a].x, mesh.vertices[a].y); path.lineTo(mesh.vertices[b].x, mesh.vertices[b].y);
    }
    meshPaths.set(mesh, path);
  }
  return meshPaths.get(mesh);
}
function surface(canvas) {
  const { width, height } = canvas.getBoundingClientRect();
  const ratio = window.devicePixelRatio || 1;
  if (canvas.width !== Math.round(width * ratio) || canvas.height !== Math.round(height * ratio)) {
    canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
  }
  const ctx = canvas.getContext('2d'); ctx.setTransform(ratio, 0, 0, ratio, 0, 0); ctx.clearRect(0, 0, width, height);
  return { ctx, width, height };
}
export function assemblyBounds(elements) {
  const points = elements.flatMap(e => e.points);
  const xMin = Math.min(...points.map(p => p.x)); const xMax = Math.max(...points.map(p => p.x));
  const yMin = Math.min(...points.map(p => p.y)); const yMax = Math.max(...points.map(p => p.y));
  const span = Math.max(xMax - xMin, yMax - yMin, 0.01);
  const centerY = (yMin + yMax) / 2;
  return { xMin: xMin - span * 0.19, xMax: xMax + span * 0.19,
    yMin: centerY - span * 0.30, yMax: centerY + span * 0.30 };
}
export function drawGeometry(canvas, { elements, lines = [], bounds, zoom = 1, pan = { x: 0, y: 0 },
  showMachContours = false, showStreamlines = true, showPanels = false, showMesh = false, mesh = null, referenceMesh = null, referenceChord = 1, boundaryLayer = null, gridAudit = null, gridAuditOverlay = 'none', flowColorMode = 'speed' }) {
  const { ctx, width, height } = surface(canvas);
  if (!bounds) return;
  const scale = Math.min((width - 30) / (bounds.xMax - bounds.xMin), (height - 30) / (bounds.yMax - bounds.yMin)) * zoom;
  const cx = (bounds.xMin + bounds.xMax) / 2; const cy = (bounds.yMin + bounds.yMax) / 2;
  const screen = p => ({ x: width / 2 + (p.x - cx) * scale + pan.x, y: height / 2 - (p.y - cy) * scale + pan.y });
  const origin = screen({ x: 0, y: 0 });
  const gridStep = niceStep(65 / scale / referenceChord) * referenceChord;
  ctx.lineWidth = 1; ctx.font = '9px Consolas, monospace';
  for (let x = Math.floor(-origin.x / scale / gridStep) * gridStep; screen({ x, y: 0 }).x < width; x += gridStep) {
    const px = screen({ x, y: 0 }).x;
    ctx.strokeStyle = Math.abs(x) < gridStep * 0.01 ? '#37505c66' : '#24394555';
    ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, height); ctx.stroke();
    if (px > 15 && px < width - 30) { ctx.fillStyle = '#3f5b6a'; ctx.fillText((x / referenceChord).toFixed(1), px + 4, height - 28); }
  }
  const yBottom = (origin.y - height) / scale;
  for (let y = Math.floor(yBottom / gridStep) * gridStep; screen({ x: 0, y }).y > 0; y += gridStep) {
    const py = screen({ x: 0, y }).y;
    ctx.strokeStyle = Math.abs(y) < gridStep * 0.01 ? '#37505c66' : '#24394555';
    ctx.beginPath(); ctx.moveTo(0, py); ctx.lineTo(width, py); ctx.stroke();
  }
  if (showMesh && mesh) {
    if (referenceMesh) {
      ctx.save(); ctx.translate(origin.x, origin.y); ctx.scale(scale, -scale);
      ctx.strokeStyle = '#d6ba7dc0'; ctx.lineWidth = .9 / scale; ctx.setLineDash([4 / scale, 3 / scale]);
      ctx.stroke(meshPath(referenceMesh)); ctx.restore();
    }
    if (gridAudit && gridAuditOverlay !== 'none') {
      const flagged = gridAuditOverlay === 'spacing' ? new Set(gridAudit.spacing.flaggedCells.map(c => c.cell)) : null;
      mesh.cells.forEach((cell, id) => {
        const value = flagged ? flagged.has(id) ? 1 : 0 : gridAudit.fieldErrors[id]?.[gridAuditOverlay === 'eta' ? 'tubeIntervals' : 'crosslineIntervals'] ?? 0;
        if (!(value > (flagged ? 0 : .05))) return;
        ctx.beginPath(); cell.forEach((vertex, k) => { const p = screen(mesh.vertices[vertex]); if (k) ctx.lineTo(p.x, p.y); else ctx.moveTo(p.x, p.y); }); ctx.closePath();
        ctx.fillStyle = flagged ? '#edab4c55' : `rgba(229,91,112,${Math.min(.65, .15 + .5 * value)})`; ctx.fill();
      });
    }
    ctx.save(); ctx.translate(origin.x, origin.y); ctx.scale(scale, -scale);
    ctx.strokeStyle = '#6a9dab88'; ctx.lineWidth = .65 / scale;
    ctx.stroke(meshPath(mesh)); ctx.restore();
  }
  if (showStreamlines) {
    if (mesh?.flow) {
      ctx.save(); ctx.translate(origin.x, origin.y); ctx.scale(scale, -scale);
      ctx.lineWidth = 1.5 / scale; ctx.lineCap = 'round';
      for (const { path, color } of coloredFlowPaths(mesh.flow, flowColorMode)) { ctx.strokeStyle = color; ctx.stroke(path); }
      ctx.restore();
      // Sparse arrowheads show the local velocity direction, without a
      // fictitious time animation during a steady Newton iteration.
      ctx.lineWidth = 1.2;
      for (const line of mesh.flow.lines) {
        let distance = 0;
        for (let i = 0; i < line.speedRatios.length; i++) {
          const a = screen(line.points[i]), b = screen(line.points[i + 1]);
          const length = Math.hypot(b.x - a.x, b.y - a.y); distance += length;
          if (distance < 110 || !(length > 0)) continue;
          const x = (a.x + b.x) / 2, y = (a.y + b.y) / 2;
          if (x < 8 || x > width - 8 || y < 8 || y > height - 8) continue;
          distance = 0;
          const dx = (b.x - a.x) / length, dy = (b.y - a.y) / length;
          ctx.strokeStyle = flowColorMode === 'change' ? flowChangeColor(line.speedChanges?.[i] ?? 0) : flowSpeedColor(line.speedRatios[i]); ctx.beginPath();
          ctx.moveTo(x - 4 * dx - 3 * dy, y - 4 * dy + 3 * dx); ctx.lineTo(x, y);
          ctx.lineTo(x - 4 * dx + 3 * dy, y - 4 * dy - 3 * dx); ctx.stroke();
        }
      }
    }
    ctx.lineWidth = 0.9;
    lines.forEach((line, i) => {
      ctx.strokeStyle = i % 3 === 0 ? '#608b9966' : '#476d7e55';
      ctx.beginPath();
      line.forEach((p, j) => { const q = screen(p); if (j) ctx.lineTo(q.x, q.y); else ctx.moveTo(q.x, q.y); });
      ctx.stroke();
    });
  }
  if (showMachContours && mesh?.flow) {
    if (!machPaths.has(mesh.flow)) machPaths.set(mesh.flow, machContours(mesh.flow).map(contour => {
      const path = new Path2D();
      for (const [a, b] of contour.segments) { path.moveTo(a.x, a.y); path.lineTo(b.x, b.y); }
      return { ...contour, path };
    }));
    ctx.save(); ctx.translate(origin.x, origin.y); ctx.scale(scale, -scale);
    for (const contour of machPaths.get(mesh.flow)) {
      const sonic = Math.abs(contour.level - 1) < 1e-10;
      ctx.strokeStyle = sonic ? '#ffd16c' : '#97caff'; ctx.lineWidth = (sonic ? 2 : .85) / scale;
      ctx.stroke(contour.path);
    }
    ctx.restore();
    ctx.font = mono;
    for (const contour of machPaths.get(mesh.flow)) {
      if (Math.round(contour.level * 20) % 2) continue;
      const segment = contour.segments.find(([a, b]) => {
        const p = screen(a); return p.x > 35 && p.x < width - 45 && p.y > 25 && p.y < height - 35 && Math.hypot(a.x-b.x,a.y-b.y)*scale > 8;
      });
      if (!segment) continue;
      const p = screen(segment[0]); ctx.fillStyle = '#111f27'; ctx.fillRect(p.x - 2, p.y - 10, 30, 12);
      ctx.fillStyle = Math.abs(contour.level - 1) < 1e-10 ? '#ffd16c' : '#97caff'; ctx.fillText(contour.level.toFixed(2), p.x, p.y);
    }
  }
  elements.forEach((e, i) => {
    const color = colors[i % colors.length];
    ctx.beginPath();
    e.points.forEach((p, j) => { const q = screen(p); if (j) ctx.lineTo(q.x, q.y); else ctx.moveTo(q.x, q.y); });
    ctx.closePath(); ctx.fillStyle = '#111f27'; ctx.fill(); ctx.fillStyle = `${color}16`; ctx.fill();
    ctx.strokeStyle = color; ctx.lineWidth = 1.7; ctx.stroke();
    if (showPanels) {
      ctx.fillStyle = color;
      e.points.slice(0, -1).forEach(p => { const q = screen(p); ctx.beginPath(); ctx.arc(q.x, q.y, 1.8, 0, Math.PI * 2); ctx.fill(); });
    }
    let center = e.points.reduce((sum, p) => ({ x: sum.x + p.x, y: sum.y + p.y }), { x: 0, y: 0 });
    center = screen({ x: center.x / e.points.length, y: center.y / e.points.length });
    ctx.font = '9px Consolas, monospace'; ctx.fillStyle = `${color}b0`; ctx.textAlign = 'center';
    ctx.fillText(String(i + 1).padStart(2, '0'), center.x, center.y + 39); ctx.textAlign = 'start';
  });
  if (showStreamlines && boundaryLayer) {
    ctx.strokeStyle = '#d6ba7dcc'; ctx.lineWidth = 1.3; ctx.setLineDash([4, 3]);
    for (const s of boundaryLayer.surfaces) {
      ctx.strokeStyle = `${colors[s.element??0]}cc`;
      ctx.beginPath(); s.stations.forEach((v, i) => {
        const q = screen(v.displacement); if (i) ctx.lineTo(q.x, q.y); else ctx.moveTo(q.x, q.y);
      }); ctx.stroke();
    }
    ctx.setLineDash([]);
    for(const wake of boundaryLayer.wakes){
      ctx.strokeStyle = `${colors[wake.element]}cc`;ctx.beginPath();
      wake.stations.forEach((v, i) => { const q = screen(v); if (i) ctx.lineTo(q.x, q.y); else ctx.moveTo(q.x, q.y); }); ctx.stroke();
    }
  }
}

function niceStep(value) {
  const magnitude = 10 ** Math.floor(Math.log10(value));
  const fraction = value / magnitude;
  return (fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10) * magnitude;
}
export function drawPressure(canvas, result, hover = null, emptyMessage = 'Run analysis to calculate Cp') {
  const { ctx, width, height } = surface(canvas);
  const box = { left: 47, top: 20, right: width - 19, bottom: height - 28 };
  const points = result?.elements.flatMap(e => e.cp) ?? [];
  const reference = result?.referenceChord ?? 1;
  const xMin = Math.min(0, ...points.map(p => p.x / reference));
  const xMax = Math.max(1, ...points.map(p => p.x / reference));
  const cpValues = points.flatMap(p => Number.isFinite(p.cpInviscid) ? [p.cp, p.cpInviscid] : [p.cp]);
  const cpMin = Math.min(-0.5, ...cpValues);
  const cpMax = Math.max(1, ...cpValues);
  const yStep = niceStep((cpMax - cpMin) / 5);
  const yMin = Math.floor(cpMin / yStep) * yStep; const yMax = Math.ceil(cpMax / yStep) * yStep;
  const sx = x => box.left + (x / reference - xMin) / (xMax - xMin) * (box.right - box.left);
  const sy = cp => box.top + (cp - yMin) / (yMax - yMin) * (box.bottom - box.top);
  ctx.font = mono; ctx.lineWidth = 1;
  for (let y = yMin; y <= yMax + yStep / 10; y += yStep) {
    const py = sy(y);
    ctx.strokeStyle = Math.abs(y) < 1e-10 ? '#3b5260' : '#26394688';
    ctx.beginPath(); ctx.moveTo(box.left, py); ctx.lineTo(box.right, py); ctx.stroke();
    ctx.fillStyle = '#748e9e'; ctx.textAlign = 'right'; ctx.fillText(y.toFixed(yStep < 1 ? 1 : 0), box.left - 10, py + 3);
  }
  const xStep = niceStep((xMax - xMin) / 5);
  for (let x = Math.ceil(xMin / xStep) * xStep; x <= xMax + 1e-10; x += xStep) {
    const px = sx(x * reference);
    ctx.strokeStyle = '#26394688'; ctx.beginPath(); ctx.moveTo(px, box.top); ctx.lineTo(px, box.bottom); ctx.stroke();
    ctx.fillStyle = '#748e9e'; ctx.textAlign = 'center'; ctx.fillText(x.toFixed(1), px, box.bottom + 17);
  }
  ctx.textAlign = 'start';
  if (!result) { ctx.fillStyle = '#5f7a89'; ctx.textAlign = 'center'; ctx.fillText(emptyMessage, width / 2, height / 2); ctx.textAlign = 'start'; return []; }
  const plotted = [];
  result.elements.forEach((e, k) => {
    if (e.cp.every(p => Number.isFinite(p.cpInviscid))) {
      ctx.strokeStyle = `${colors[k]}66`; ctx.setLineDash([4, 4]); ctx.lineWidth = 1;
      ctx.beginPath(); e.cp.forEach((p, i) => { if (i) ctx.lineTo(sx(p.x), sy(p.cpInviscid)); else ctx.moveTo(sx(p.x), sy(p.cpInviscid)); });
      ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.strokeStyle = colors[k]; ctx.lineWidth = 1.6; ctx.beginPath();
    e.cp.forEach((p, i) => {
      const x = sx(p.x); const y = sy(p.cp);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      plotted.push({ ...p, screenX: x, screenY: y, name: e.name, color: colors[k] });
    });
    ctx.stroke();
  });
  if (hover) {
    ctx.fillStyle = hover.color; ctx.beginPath(); ctx.arc(hover.screenX, hover.screenY, 3.5, 0, 2 * Math.PI); ctx.fill();
    ctx.strokeStyle = '#a3bcc366'; ctx.setLineDash([3, 4]); ctx.beginPath(); ctx.moveTo(hover.screenX, box.top); ctx.lineTo(hover.screenX, box.bottom); ctx.stroke(); ctx.setLineDash([]);
  }
  return plotted;
}

export function drawBoundaryLayer(canvas, result, quantity = 'cf') {
  const { ctx, width, height } = surface(canvas);
  if (!result?.boundaryLayer || width === 0) return;
  const reference = result.referenceChord;
  const dimensional = quantity === 'theta' || quantity === 'deltaStar';
  const value = p => p[quantity] / (dimensional ? reference : 1);
  const stations = result.boundaryLayer.surfaces.flatMap(s => s.stations);
  const xmin = Math.min(...stations.map(s => s.x / reference)); const xmax = Math.max(...stations.map(s => s.x / reference));
  const ymin = Math.min(0, ...stations.map(value)); const ymax = Math.max(1e-5, ...stations.map(value));
  const step = niceStep((ymax - ymin) / 4); const lo = Math.floor(ymin / step) * step; const hi = Math.ceil(ymax / step) * step;
  const box = { left: 62, right: width - 25, top: 20, bottom: height - 30 };
  const sx = x => box.left + (x - xmin) / (xmax - xmin) * (box.right - box.left);
  const sy = y => box.bottom - (y - lo) / (hi - lo) * (box.bottom - box.top);
  ctx.font = mono; ctx.lineWidth = 1;
  for (let y = lo; y <= hi + step * 0.01; y += step) {
    ctx.strokeStyle = Math.abs(y) < 1e-12 ? '#506773' : '#26394688';
    ctx.beginPath(); ctx.moveTo(box.left, sy(y)); ctx.lineTo(box.right, sy(y)); ctx.stroke();
    ctx.fillStyle = '#748e9e'; ctx.textAlign = 'right'; ctx.fillText(y.toFixed(quantity === 'h' ? 1 : 4), box.left - 10, sy(y) + 3);
  }
  for (let k = 0; k <= 5; k++) {
    const x = xmin + k * (xmax - xmin) / 5;
    ctx.fillStyle = '#748e9e'; ctx.textAlign = 'center'; ctx.fillText(x.toFixed(2), sx(x), box.bottom + 18);
  }
  result.boundaryLayer.surfaces.forEach(s => {
    ctx.strokeStyle = colors[s.element]; ctx.lineWidth = 1.7; ctx.setLineDash(s.side==='lower'?[6,3]:[]);ctx.beginPath();
    s.stations.forEach((v, i) => { if (i) ctx.lineTo(sx(v.x / reference), sy(value(v))); else ctx.moveTo(sx(v.x / reference), sy(value(v))); }); ctx.stroke();
    const transition = s.transitionKind === 'laminar-to-te' ? null : s.stations.find(v => ['leading-transition', 'transition', 'turbulent'].includes(v.regime));
    if (transition) {
      // The transition interval supplies its interpolated physical coordinate.
      const tx = s.transitionPoint.x / reference;
      ctx.setLineDash([3, 4]); ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(sx(tx), box.top); ctx.lineTo(sx(tx), box.bottom); ctx.stroke(); ctx.setLineDash([]);
    }
  });
  ctx.setLineDash([]);
  ctx.textAlign = 'start';
}
