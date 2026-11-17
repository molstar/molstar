/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CanvasBackend, type CanvasLike, FontCache } from '@molstar/seqx-viewer/canvas/index';
import type { SeqViewer } from './viewer.js';

// The mounted UI (_spec/README.md "Input"): the canvas stack (a base canvas and the overlay above it, which takes no
// events) in a focusable container, and the DOM events translated into the viewer's input calls. Chrome beyond the
// canvas is not built in this version: the track labels are the layout's own. Every listener, observer and node it
// adds is taken away by the function it returns.

/** The container is a block as wide as its target; the canvases are as tall as the layout's content. */
export function mountCanvas(viewer: SeqViewer, target: HTMLElement): () => void {
  const doc = target.ownerDocument;
  const container = doc.createElement('div');
  container.tabIndex = 0;
  container.setAttribute('role', 'application');
  container.setAttribute('aria-label', 'Sequence view');
  Object.assign(container.style, { position: 'relative', width: '100%', outline: 'none', touchAction: 'none' });
  const base = doc.createElement('canvas');
  const overlay = doc.createElement('canvas');
  Object.assign(base.style, { display: 'block' });
  Object.assign(overlay.style, { position: 'absolute', left: '0', top: '0', pointerEvents: 'none' });
  container.append(base, overlay);
  target.appendChild(container);

  const fonts = new FontCache(viewer.context.tokens.fonts);
  const backend = new CanvasBackend(base as unknown as CanvasLike, overlay as unknown as CanvasLike, {
    fonts,
    fade: viewer.context.renderConfig.fade,
  });
  viewer.useBackend(backend);

  const input = viewer.input;
  const local = (e: MouseEvent) => {
    const r = container.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top] as const;
  };
  const listeners: [string, (e: any) => void, AddEventListenerOptions | undefined][] = [
    [
      'pointermove',
      (e: PointerEvent) => {
        const [x, y] = local(e);
        input.move(x, y, e);
      },
      undefined,
    ],
    ['pointerleave', () => input.leave(), undefined],
    [
      'pointerdown',
      (e: PointerEvent) => {
        if (e.button !== 0) return;
        const [x, y] = local(e);
        if (input.down(x, y, e)) container.setPointerCapture?.(e.pointerId);
        container.focus({ preventScroll: true });
      },
      undefined,
    ],
    [
      'pointerup',
      (e: PointerEvent) => {
        input.up();
        if (container.hasPointerCapture?.(e.pointerId)) container.releasePointerCapture(e.pointerId);
      },
      undefined,
    ],
    ['pointercancel', () => input.up(), undefined],
    [
      'wheel',
      (e: WheelEvent) => {
        const [x, y] = local(e);
        if (input.wheel(x, y, e)) e.preventDefault();
      },
      { passive: false },
    ],
    [
      'dblclick',
      (e: MouseEvent) => {
        const [x, y] = local(e);
        input.doubleClick(x, y, e);
      },
      undefined,
    ],
    [
      'keydown',
      (e: KeyboardEvent) => {
        if (e.key === ' ') input.setSpace(true);
        else if (input.key(e)) e.preventDefault();
      },
      undefined,
    ],
    ['keyup', (e: KeyboardEvent) => e.key === ' ' && input.setSpace(false), undefined],
    ['blur', () => input.setSpace(false), undefined],
  ];
  for (const [type, fn, options] of listeners) container.addEventListener(type, fn, options);

  // the width follows the container; the height is the layout's
  let observer: ResizeObserver | undefined;
  const measure = () => viewer.resize(container.clientWidth);
  if (typeof ResizeObserver !== 'undefined') {
    observer = new ResizeObserver(measure);
    observer.observe(container);
  }
  measure();

  // a web font arrived: what was measured is measured again
  const fontFaces = doc.fonts;
  const onFonts = () => viewer.invalidateText();
  fontFaces?.addEventListener?.('loadingdone', onFonts);

  return () => {
    for (const [type, fn, options] of listeners) container.removeEventListener(type, fn, options);
    observer?.disconnect();
    fontFaces?.removeEventListener?.('loadingdone', onFonts);
    viewer.useBackend(undefined);
    container.remove();
  };
}
