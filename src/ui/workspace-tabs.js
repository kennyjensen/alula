// SPDX-License-Identifier: GPL-2.0-or-later
// The same panes form desktop columns and a mobile scroll-snap workspace.
export function bindWorkspaceTabs({ shell, tabs, settings, visualization, actions, actionSlot }) {
  const panes = [settings, visualization];
  const buttons = [...tabs.querySelectorAll('button')];
  const mobile = matchMedia('(max-width: 700px)');
  let selected = 1, gesture = null;
  const update = () => {
    const actionParent = mobile.matches ? document.body : actionSlot;
    if (actions.parentElement !== actionParent) actionParent.append(actions);
    buttons.forEach((button, i) => {
      button.setAttribute('aria-selected', String(i === selected));
      button.tabIndex = i === selected ? 0 : -1;
      panes[i].inert = mobile.matches && i !== selected;
      if (mobile.matches) panes[i].setAttribute('role', 'tabpanel');
      else panes[i].removeAttribute('role');
    });
  };
  const select = (index, smooth = true) => {
    selected = Math.max(0, Math.min(1, index));
    update();
    if (mobile.matches) shell.scrollTo({ left: selected * shell.clientWidth,
      behavior: smooth && !matchMedia('(prefers-reduced-motion: reduce)').matches ? 'smooth' : 'instant' });
  };
  buttons.forEach((button, index) => button.addEventListener('click', () => select(index, false)));
  tabs.addEventListener('keydown', event => {
    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    select(event.key === 'Home' ? 0 : event.key === 'End' ? 1
      : event.key === 'ArrowLeft' ? selected - 1 : selected + 1);
    buttons[selected].focus();
  });
  // Horizontal touch gestures switch panes; vertical gestures keep scrolling
  // the active pane. Controls and the zoomable canvas own their gestures.
  shell.addEventListener('pointerdown', event => {
    if (!mobile.matches || event.pointerType !== 'touch') return;
    if (gesture) { gesture = null; shell.style.scrollSnapType = ''; select(selected); return; }
    if (event.target.closest('#geometry-canvas, input, select, textarea, button, a, summary')) return;
    gesture = { id: event.pointerId, x: event.clientX, y: event.clientY,
      left: selected * shell.clientWidth, dragging: false };
  });
  shell.addEventListener('pointermove', event => {
    if (!gesture || event.pointerId !== gesture.id) return;
    const dx = event.clientX - gesture.x, dy = event.clientY - gesture.y;
    if (!gesture.dragging && Math.abs(dx) > 20 && Math.abs(dx) > 1.5 * Math.abs(dy)) {
      gesture.dragging = true; shell.setPointerCapture(event.pointerId);
      shell.style.scrollSnapType = 'none';
    }
    if (gesture.dragging) { event.preventDefault(); shell.scrollLeft = gesture.left - dx; }
  }, { passive: false });
  const finishGesture = event => {
    if (!gesture || event.pointerId !== gesture.id) return;
    const dx = event.clientX - gesture.x;
    const switchPane = event.type === 'pointerup' && gesture.dragging && Math.abs(dx) > 70;
    gesture = null; shell.style.scrollSnapType = '';
    if (shell.hasPointerCapture(event.pointerId)) shell.releasePointerCapture(event.pointerId);
    select(selected + (switchPane ? dx < 0 ? 1 : -1 : 0));
  };
  shell.addEventListener('pointerup', finishGesture);
  shell.addEventListener('pointercancel', finishGesture);
  window.addEventListener('resize', () => select(selected, false));
  mobile.addEventListener('change', () => {
    select(selected, false);
    if (!mobile.matches) shell.scrollLeft = 0;
  });
  select(selected, false);
  return { swipe: direction => { if (mobile.matches) select(selected + direction); }, mobile };
}
