// SPDX-License-Identifier: GPL-2.0-or-later
// Pointer coordinates and pan offsets use CSS pixels, independent of DPR.
// Pinching keeps the plot point under the finger midpoint anchored there.
export function bindPlotNavigation(canvas, { getView, setView, minZoom = .5, maxZoom = 192, onSwipe, swipeEnabled = () => false }) {
  const pointers = new Map(); let gesture = null; let swipe = null;
  const point = event => {
    const box = canvas.getBoundingClientRect(); return { x: event.clientX - box.left, y: event.clientY - box.top };
  };
  const sample = () => {
    const [a, b] = pointers.values();
    if (!a) return null;
    return b ? { x: .5 * (a.x + b.x), y: .5 * (a.y + b.y), distance: Math.hypot(b.x - a.x, b.y - a.y) }
      : { ...a, distance: 0 };
  };
  const rebase = () => {
    const p = sample(), view = getView();
    gesture = p ? { point: p, zoom: view.zoom, pan: { ...view.pan } } : null;
  };
  const change = (view, from, to, requestedZoom) => {
    const zoom = Math.max(minZoom, Math.min(maxZoom, requestedZoom)), ratio = zoom / view.zoom;
    const box = canvas.getBoundingClientRect(), cx = box.width / 2, cy = box.height / 2;
    setView({ zoom, pan: { x: to.x - cx - ratio * (from.x - cx - view.pan.x),
      y: to.y - cy - ratio * (from.y - cy - view.pan.y) } });
  };
  const down = event => {
    if (event.pointerType === 'mouse' && event.button !== 0) return;
    event.preventDefault(); pointers.set(event.pointerId, point(event));
    canvas.setPointerCapture(event.pointerId); rebase();
    if (pointers.size > 1) swipe = null;
    else if (event.pointerType === 'touch' && swipeEnabled())
      swipe = { start: point(event), end: point(event), view: getView() };
  };
  const move = event => {
    if (!pointers.has(event.pointerId) || !gesture) return;
    event.preventDefault(); pointers.set(event.pointerId, point(event));
    if (swipe) swipe.end = point(event);
    const p = sample(), factor = pointers.size > 1 && gesture.point.distance > 0 ? p.distance / gesture.point.distance : 1;
    change(gesture, gesture.point, p, gesture.zoom * factor);
  };
  const up = event => {
    if (!pointers.delete(event.pointerId)) return;
    if (canvas.hasPointerCapture(event.pointerId)) canvas.releasePointerCapture(event.pointerId);
    if (swipe && pointers.size === 0) {
      const dx = swipe.end.x - swipe.start.x, dy = swipe.end.y - swipe.start.y;
      if (event.type === 'pointerup' && Math.abs(dx) > 70 && Math.abs(dx) > 1.5 * Math.abs(dy)) {
        setView(swipe.view); onSwipe?.(dx < 0 ? 1 : -1);
      }
      swipe = null;
    }
    rebase(); // Adding/removing a finger starts from the current view.
  };
  const clear = () => {
    const ids = [...pointers.keys()]; pointers.clear(); gesture = null; swipe = null;
    for (const id of ids) if (canvas.hasPointerCapture(id)) canvas.releasePointerCapture(id);
  };
  const wheel = event => {
    event.preventDefault(); const p = point(event);
    change(getView(), p, p, getView().zoom * Math.exp(-event.deltaY * .001)); rebase();
  };
  const listeners = { pointerdown: down, pointermove: move, pointerup: up, pointercancel: up, lostpointercapture: up, wheel };
  for (const [type, listener] of Object.entries(listeners)) canvas.addEventListener(type, listener, { passive: false });
  window.addEventListener('blur', clear);
  return () => {
    clear(); for (const [type, listener] of Object.entries(listeners)) canvas.removeEventListener(type, listener);
    window.removeEventListener('blur', clear);
  };
}
