/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import { modifiersOf } from '@molstar/seqx-viewer/bindings';
import { CanvasBackend, type CanvasLike, FontCache } from '@molstar/seqx-viewer/canvas';
import { createDragSession } from './drag-session.js';
import type { NavigatorModel } from './navigator.js';
import type { SeqViewer } from './viewer.js';

// The navigator part (_spec/README.md "Composition", Navigators): a part the host places, like a legend. `mountNavigator`
// draws one navigator (a declared node's id, or a host's `viewer.navigator(def)`) into an element the host owns: its
// own base and overlay canvases sized by the element, its listeners on it alone, one tab stop whose windows are sliders
// (roving: a press on a window focuses its slider, PageDown and PageUp move the focus; `aria-valuetext` set when the
// view settles after this part's own gesture or key, never a frame, never for a move made elsewhere). `navigatorsMount`, the
// default `ui.navigators`, places the view's declared navigators by their `place` around the canvas through the same
// call; `railToggleMount`, `ui.railToggle`, a button in the canvas's corner that hides and shows the ruler rails (the
// reader's session). Everything each adds is taken away by the function it returns.

/** Mounts navigator `id` of the viewer into `element`; the function it returns takes it away. */
export function mountNavigator(element: HTMLElement, viewer: SeqViewer, id: string): () => void {
  const model = viewer.navigatorModel(id);
  if (!model) return () => {};
  return mountModel(element, viewer, model);
}

function mountModel(element: HTMLElement, viewer: SeqViewer, model: NavigatorModel): () => void {
  const doc = element.ownerDocument;
  const container = doc.createElement('div');
  container.tabIndex = 0;
  container.setAttribute('role', 'slider');
  container.dataset.navigator = model.id;
  Object.assign(container.style, { position: 'relative', outline: 'none', touchAction: 'none', overflow: 'hidden' });
  const base = doc.createElement('canvas');
  const overlay = doc.createElement('canvas');
  const labels = doc.createElement('canvas');
  Object.assign(base.style, { display: 'block' });
  for (const c of [overlay, labels])
    Object.assign(c.style, { position: 'absolute', left: '0', top: '0', pointerEvents: 'none' });
  // the readout is the chrome's surface (the hover card's inset and backdrop, so it reads over data), one line cut with an
  // ellipsis before it would wrap, and shown only while it says something
  const readout = doc.createElement('div');
  readout.className = 'msui-surface seqx-navigator-readout';
  readout.setAttribute('aria-live', 'off');
  Object.assign(readout.style, {
    position: 'absolute',
    right: 'var(--msui-space-1, 4px)',
    top: 'var(--msui-space-1, 4px)',
    display: 'none',
    maxWidth: 'calc(100% - 2 * var(--msui-space-1, 4px))',
    boxSizing: 'border-box',
    overflow: 'hidden',
    whiteSpace: 'nowrap',
    textOverflow: 'ellipsis',
    font: 'var(--msui-text-xs, 10px) / var(--msui-leading-xs, 14px) var(--msui-font-sans, sans-serif)',
    color: 'var(--msui-text-muted, gray)',
    pointerEvents: 'none',
  });
  container.append(base, overlay, labels, readout);
  element.append(container);
  const backend = new CanvasBackend(
    base as unknown as CanvasLike,
    overlay as unknown as CanvasLike,
    labels as unknown as CanvasLike,
    { fonts: new FontCache(viewer.tokens.fonts), fade: viewer.context.renderConfig.fade },
  );

  let frame: number | undefined;
  const draw = () => {
    frame = undefined;
    const f = model.draw();
    box();
    if (!f) return;
    const dirty = backend.resize(model.width, model.height, viewer.dpr);
    backend.draw(f.lists, f.rebuilt || f.base || dirty);
    const said = model.readout.value ?? '';
    readout.textContent = said;
    readout.style.display = said ? 'block' : 'none';
    container.style.cursor = model.cursor;
  };
  const schedule = () => {
    if (frame === undefined) frame = requestFrameCallback(draw);
  };
  // the focused slider's value, as the model says it (on settle, never a frame)
  container.setAttribute('aria-label', `Navigator ${model.id}`);
  const say = () => {
    const v = model.values.value[model.focused.value];
    container.setAttribute('aria-valuetext', v?.text ?? '');
    if (!v || v.disabled) container.setAttribute('aria-disabled', 'true');
    else container.removeAttribute('aria-disabled');
  };
  say();
  const unsubscribe = [
    model.version.listen(schedule),
    viewer.drawn.listen(schedule),
    model.values.listen(say),
    model.focused.listen(say),
  ];

  const local = (e: MouseEvent) => {
    const r = container.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top] as const;
  };
  // a drag ends wherever the button is released: the window's own listeners while it runs (drag-session.ts)
  const drag = createDragSession(container, () => model.up());
  const touch = (e: PointerEvent) => ({ touch: e.pointerType === 'touch' });
  const listeners: [string, (e: any) => void, AddEventListenerOptions | undefined][] = [
    ['pointermove', (e: PointerEvent) => model.move(...local(e), touch(e)), undefined],
    ['pointerleave', () => model.leave(), undefined],
    [
      'pointerdown',
      (e: PointerEvent) => {
        if (e.button !== 0) return;
        // a release that never arrived: the earlier drag ends before this press begins
        drag.end();
        const [x, y] = local(e);
        if (model.down(x, y, touch(e))) {
          container.setPointerCapture?.(e.pointerId);
          drag.start();
        }
        container.focus({ preventScroll: true });
      },
      undefined,
    ],
    [
      'pointerup',
      (e: PointerEvent) => {
        drag.end();
        model.up();
        if (container.hasPointerCapture?.(e.pointerId)) container.releasePointerCapture(e.pointerId);
      },
      undefined,
    ],
    ['pointercancel', () => (drag.end(), model.cancel()), undefined],
    ['lostpointercapture', () => (drag.end(), model.up()), undefined],
    ['click', (e: MouseEvent) => e.button === 0 && model.click(...local(e), e.detail), undefined],
    ['dblclick', (e: MouseEvent) => model.doubleClick(...local(e)), undefined],
    ['focus', () => model.setFocus(true), undefined],
    ['blur', () => model.setFocus(false), undefined],
    [
      'wheel',
      (e: WheelEvent) => {
        const [x, y] = local(e);
        if (model.wheel(x, y, e.deltaX, e.deltaY, modifiersOf(e, viewer.input.mac))) e.preventDefault();
      },
      { passive: false },
    ],
    [
      'keydown',
      (e: KeyboardEvent) => {
        if (model.key(e.key, e)) e.preventDefault();
      },
      undefined,
    ],
  ];
  for (const [type, fn, options] of listeners) container.addEventListener(type, fn, options);

  // its box (_spec/README.md "Composition", Navigators): a strip across the element, as wide as it; a side part and a
  // miniature at the form's box (an overview the view's height), never the element's width (an empty cell in a row
  // measures its content: a canvas's default 300 px would hold itself); the height the form's unless the element sets
  // one. The container's size is set to the model's, so the cell beside the canvas is the part's width and no more.
  const across = () => {
    const nav = model.navigator;
    return !!nav && (nav.place === 'top' || nav.place === 'bottom') && nav.axes !== 'both';
  };
  const box = () => {
    container.style.width = across() ? '' : `${model.width}px`;
    container.style.height = `${model.height || 0}px`;
  };
  let observer: ResizeObserver | undefined;
  const measure = () => {
    if (!model.navigator) return;
    model.resize(across() ? element.clientWidth : 0, element.style.height ? element.clientHeight : 0);
    box();
    schedule();
  };
  if (typeof ResizeObserver !== 'undefined') {
    observer = new ResizeObserver(measure);
    observer.observe(element);
  }
  measure();
  schedule();

  return () => {
    for (const [type, fn, options] of listeners) container.removeEventListener(type, fn, options);
    for (const u of unsubscribe) u();
    drag.dispose();
    observer?.disconnect();
    if (frame !== undefined) cancelFrameCallback(frame);
    container.remove();
  };
}

/**
 * The default `ui.navigators` (_spec/README.md "Composition", Navigators): the view's shown navigators placed by their
 * `place` around the canvas (top before it, bottom after it, left and right beside it in a row), through
 * `mountNavigator`, as a host would; placed again when the view's navigators change.
 */
export function navigatorsMount(viewer: SeqViewer, target: HTMLElement): () => void {
  const doc = target.ownerDocument;
  const canvas = target.querySelector<HTMLElement>('[role="application"]');
  const row = doc.createElement('div');
  Object.assign(row.style, { display: 'flex', alignItems: 'stretch', width: '100%' });
  const middle = doc.createElement('div');
  Object.assign(middle.style, { flex: '1 1 auto', minWidth: '0' });
  // a narrow frame (core's "Render", A narrow frame): each side part's room is the row less the other side parts; the
  // canvas keeps the layout's least width and the middle scrolls past it (only then: a card is never clipped)
  let sides: { cell: HTMLElement; id: string }[] = [];
  const fit = () => {
    const width = row.clientWidth;
    for (const s of sides) {
      let others = 0;
      for (const o of sides) if (o !== s) others += o.cell.offsetWidth;
      viewer.navigatorModel(s.id)?.setRoom(width - others);
    }
    const least = viewer.frame?.layout.narrow?.minWidth ?? 0;
    middle.style.overflowX = least > middle.clientWidth + 0.5 ? 'auto' : '';
  };
  const top = doc.createElement('div');
  const bottom = doc.createElement('div');
  if (canvas) {
    canvas.replaceWith(row);
    middle.append(canvas);
  } else target.append(row);
  row.append(middle);
  middle.before(doc.createComment('left'));
  row.before(top);
  row.after(bottom);
  let mounted: { key: string; dispose: () => void }[] = [];
  const sync = () => {
    const navs = (viewer.view.resolved.value?.navigators ?? []).filter((n) => !n.hidden);
    const key = navs.map((n) => `${n.id}:${n.place}`).join(',');
    if (mounted.map((m) => m.key).join(',') === key && mounted.length === navs.length) return;
    for (const m of mounted) m.dispose();
    sides = [];
    mounted = navs.map((n) => {
      const cell = doc.createElement('div');
      if (n.place === 'left' || n.place === 'right') {
        cell.style.flex = '0 0 auto';
        sides.push({ cell, id: n.id });
      }
      if (n.place === 'top') top.append(cell);
      else if (n.place === 'bottom') bottom.append(cell);
      else if (n.place === 'left') middle.before(cell);
      else row.append(cell);
      const dispose = mountNavigator(cell, viewer, n.id);
      return {
        key: `${n.id}:${n.place}`,
        dispose: () => {
          dispose();
          cell.remove();
        },
      };
    });
  };
  const unsubscribe = viewer.view.resolved.listen(sync);
  sync();
  const unfit = viewer.drawn.listen(fit);
  let observer: ResizeObserver | undefined;
  if (typeof ResizeObserver !== 'undefined') {
    observer = new ResizeObserver(fit);
    observer.observe(row);
  }
  fit();
  return () => {
    unsubscribe();
    unfit();
    observer?.disconnect();
    for (const m of mounted) m.dispose();
    mounted = [];
    if (canvas) row.replaceWith(canvas);
    else row.remove();
    top.remove();
    bottom.remove();
  };
}

/**
 * The default `ui.railToggle` (_spec/README.md "Composition", Navigators): a button in the canvas's corner that hides
 * and shows the ruler rails, the reader's session, never the spec.
 */
export function railToggleMount(viewer: SeqViewer, target: HTMLElement): () => void {
  const doc = target.ownerDocument;
  const holder = target.querySelector<HTMLElement>('[role="application"]') ?? target;
  const button = doc.createElement('button');
  button.type = 'button';
  Object.assign(button.style, {
    position: 'absolute',
    right: '2px',
    top: '2px',
    font: '10px var(--msui-font-sans, sans-serif)',
    padding: '0 4px',
    color: 'var(--msui-text-muted, gray)',
    background: 'var(--msui-bg, transparent)',
    border: '1px solid var(--msui-line, currentColor)',
    borderRadius: '3px',
    cursor: 'pointer',
  });
  const say = () => {
    const shown = viewer.railsShown.value;
    button.textContent = shown ? 'Hide rails' : 'Show rails';
    button.setAttribute('aria-pressed', String(shown));
  };
  const onClick = (e: MouseEvent) => {
    e.stopPropagation();
    viewer.setRailsShown(!viewer.railsShown.value);
  };
  button.addEventListener('click', onClick);
  button.addEventListener('pointerdown', (e) => e.stopPropagation());
  holder.append(button);
  const unsubscribe = viewer.railsShown.listen(say);
  say();
  return () => {
    unsubscribe();
    button.removeEventListener('click', onClick);
    button.remove();
  };
}
