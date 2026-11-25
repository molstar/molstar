/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import { modifiersOf } from '@molstar/seqx-viewer/bindings';
import { CanvasBackend, type CanvasLike, FontCache } from '@molstar/seqx-viewer/canvas';
import { createDragSession } from './drag-session.js';
import { isSquareNavigator } from '@molstar/seqx-core/render/navigator-square-kind';
import { type NavigatorModel, ReadoutLinePx } from './navigator.js';
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
export function mountNavigator(
  element: HTMLElement,
  viewer: SeqViewer,
  id: string,
  options: MountNavigatorOptions = {},
): () => void {
  const model = viewer.navigatorModel(id);
  if (!model) return () => {};
  return mountModel(element, viewer, model, options);
}

/** `mountNavigator`'s options: `hint: false` leaves the square's line of gestures out (default `navigator.squareHint`). */
export interface MountNavigatorOptions {
  readonly hint?: boolean;
}

/** What the square says of its gestures (_spec/README.md "Composition", Navigators, The square's hint). */
export const SquareHintBase = 'Drag to move · drag an edge to resize';
/** With a switch (`include`) a double-click is the way home. */
export const SquareHint = `${SquareHintBase} · double-click for the whole genome`;

function mountModel(
  element: HTMLElement,
  viewer: SeqViewer,
  model: NavigatorModel,
  options: MountNavigatorOptions,
): () => void {
  const doc = element.ownerDocument;
  const container = doc.createElement('div');
  container.tabIndex = 0;
  container.setAttribute('role', 'slider');
  container.dataset.navigator = model.id;
  // the part's box is the readout's line above the canvases: padding, so a point's y is the picture's less the line
  Object.assign(container.style, {
    position: 'relative',
    outline: 'none',
    touchAction: 'none',
    boxSizing: 'content-box',
    paddingTop: `${ReadoutLinePx}px`,
  });
  const base = doc.createElement('canvas');
  const overlay = doc.createElement('canvas');
  const labels = doc.createElement('canvas');
  Object.assign(base.style, { display: 'block' });
  for (const c of [overlay, labels])
    Object.assign(c.style, { position: 'absolute', left: '0', top: `${ReadoutLinePx}px`, pointerEvents: 'none' });
  // the readout is the part's first line, one place: its top right corner, never moving with the window or the pointer
  // (seqx-core _spec "Render", Navigators). One line; the number is the information: the system's name shrinks to an
  // ellipsis, then to nothing, the value and what follows never do; a value wider than the part runs left (flex-end)
  const nav = model.navigator;
  const side = nav?.place === 'left' || nav?.place === 'right' || nav?.axes === 'rows';
  const readout = doc.createElement('div');
  readout.className = 'seqx-navigator-readout';
  readout.setAttribute('aria-live', 'off');
  Object.assign(readout.style, {
    position: 'absolute',
    right: '0',
    top: '0',
    height: `${ReadoutLinePx}px`,
    display: 'none',
    justifyContent: 'flex-end',
    // within the part's width (a miniature's text is clipped at its start); a narrow rows overview's short form runs left
    maxWidth: side ? 'none' : '100%',
    overflow: side ? 'visible' : 'hidden',
    boxSizing: 'border-box',
    whiteSpace: 'nowrap',
    font: `var(--msui-text-xs, 11px) / ${ReadoutLinePx}px var(--msui-font-sans, sans-serif)`,
    color: 'var(--msui-text-muted, gray)',
    pointerEvents: 'none',
  });
  const readName = doc.createElement('span');
  Object.assign(readName.style, { flex: '0 1 auto', minWidth: '0', overflow: 'hidden', textOverflow: 'ellipsis' });
  const readValue = doc.createElement('span');
  Object.assign(readValue.style, { flex: 'none' });
  readout.append(readName, readValue);
  container.append(base, overlay, labels, readout);
  element.append(container);
  // a square's hint: a muted line under it (wrapped to the square's width), no event of its own; made when the view
  // resolves a square (its blocks may land after the part is mounted), said by what the square does
  const hint = doc.createElement('div');
  hint.className = 'seqx-navigator-hint';
  Object.assign(hint.style, {
    font: `var(--msui-text-xs, 11px) / ${ReadoutLinePx}px var(--msui-font-sans, sans-serif)`,
    color: 'var(--msui-text-muted, gray)',
    pointerEvents: 'none',
  });
  const syncHint = () => {
    const view = viewer.view.resolved.value;
    const nav = model.navigator;
    const square =
      !!(options.hint ?? viewer.config.navigator.squareHint) && !!nav && !!view && isSquareNavigator(nav, view);
    if (!square) return hint.remove();
    hint.textContent = nav!.switch === 'include' ? SquareHint : SquareHintBase;
    if (!hint.isConnected) element.append(hint);
  };
  const backend = new CanvasBackend(
    base as unknown as CanvasLike,
    overlay as unknown as CanvasLike,
    labels as unknown as CanvasLike,
    { fonts: new FontCache(viewer.tokens.fonts), fade: viewer.context.renderConfig.fade },
  );

  // the part's text and box written only when they change (an inline style and a text read back cost no layout): the
  // part draws on every frame of the viewer's, a hover's too, and a write there dirtied the page's layout, which the
  // canvas's next pointer move then read back (its `getBoundingClientRect`): a forced layout a move
  const text = (el: HTMLElement, t: string) => {
    if (el.textContent !== t) el.textContent = t;
  };
  const style = (el: HTMLElement, key: 'width' | 'height' | 'display' | 'cursor', v: string) => {
    if (el.style[key] !== v) el.style[key] = v;
  };
  let frame: number | undefined;
  const draw = () => {
    frame = undefined;
    const f = model.draw();
    box();
    if (!f) return;
    const dirty = backend.resize(model.width, model.height, viewer.dpr);
    backend.draw(f.lists, f.rebuilt || f.base || dirty);
    const parts = model.readoutParts.value;
    text(readName, parts?.name ? `${parts.name}: ` : '');
    text(readValue, parts ? `${parts.value}${parts.rest}` : '');
    style(readout, 'display', parts && (parts.name || parts.value) ? 'flex' : 'none');
    style(container, 'cursor', model.cursor);
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
  syncHint();
  const unsubscribe = [
    model.version.listen(schedule),
    viewer.drawn.listen(schedule),
    model.values.listen(say),
    model.focused.listen(say),
    viewer.view.resolved.listen(syncHint),
  ];

  const local = (e: MouseEvent) => {
    const r = container.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top - ReadoutLinePx] as const;
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
    ['click', (e: MouseEvent) => e.button === 0 && model.click(...local(e), e.detail, e), undefined],
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
    style(container, 'width', across() ? '' : `${model.width}px`);
    // the picture's height; the readout's line is the container's padding above it
    style(container, 'height', `${model.height || 0}px`);
    style(hint, 'width', across() ? '' : `${model.width}px`);
  };
  let observer: ResizeObserver | undefined;
  const measure = () => {
    if (!model.navigator) return;
    model.resize(
      across() ? element.clientWidth : 0,
      element.style.height ? Math.max(1, element.clientHeight - ReadoutLinePx) : 0,
    );
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
    hint.remove();
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
