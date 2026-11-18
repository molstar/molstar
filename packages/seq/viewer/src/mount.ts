/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CanvasBackend, type CanvasLike, FontCache, readCssTokens, ThemeAttribute } from '@molstar/seqx-viewer/canvas';
import type { SeqViewer } from './viewer.js';

// The mounted UI (_spec/README.md "Input"): the canvas stack (the base canvas, the overlay above it and the labels above
// that; only the container takes events) in a focusable container, and the DOM events translated into the viewer's input calls. Chrome beyond the
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
  const labels = doc.createElement('canvas');
  Object.assign(base.style, { display: 'block' });
  for (const c of [overlay, labels])
    Object.assign(c.style, { position: 'absolute', left: '0', top: '0', pointerEvents: 'none' });
  container.append(base, overlay, labels);
  target.appendChild(container);

  // the token table: the root's computed style (`theme.source: 'root'`) read before the first frame, and again when
  // `data-msui-theme` changes on the document element or the nearest ancestor that has it; else the config's
  let themeObserver: MutationObserver | undefined;
  if (viewer.config.theme.source === 'root') {
    viewer.setTokens(readCssTokens(container));
    if (typeof MutationObserver !== 'undefined') {
      themeObserver = new MutationObserver(() => {
        viewer.setTokens(readCssTokens(container));
      });
      const options = { attributes: true, attributeFilter: [ThemeAttribute] };
      const holder = target.closest(`[${ThemeAttribute}]`);
      themeObserver.observe(doc.documentElement, options);
      if (holder && holder !== doc.documentElement) themeObserver.observe(holder, options);
    }
  }

  const fonts = new FontCache(viewer.tokens.fonts);
  const backend = new CanvasBackend(
    base as unknown as CanvasLike,
    overlay as unknown as CanvasLike,
    labels as unknown as CanvasLike,
    { fonts, fade: viewer.context.renderConfig.fade },
  );
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
    // a drag ends however the pointer is lost: cancelled by the platform, or its capture taken away
    ['pointercancel', () => input.up(), undefined],
    ['lostpointercapture', () => input.up(), undefined],
    [
      'wheel',
      (e: WheelEvent) => {
        const [x, y] = local(e);
        if (input.wheel(x, y, e)) e.preventDefault();
      },
      { passive: false },
    ],
    [
      'click',
      (e: MouseEvent) => {
        if (e.button !== 0) return;
        const [x, y] = local(e);
        input.click(x, y, e);
      },
      undefined,
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
        if (e.key === ' ') {
          if (!e.repeat) input.setSpace(true);
          // Space is the viewer's while it is a modifier or a tap toggles the cursor's cell: the page does not scroll
          // with it then; else it is the page's
          if (input.spaceIsModifier || viewer.selection.cursor.value) e.preventDefault();
        } else if (input.key(e)) e.preventDefault();
      },
      undefined,
    ],
    ['keyup', (e: KeyboardEvent) => e.key === ' ' && input.spaceUp(e), undefined],
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

  // the device pixel ratio changes (a window moved to another screen, a zoom): the query for the ratio now holds is
  // listened to once, and made again for the new one
  const win = doc.defaultView;
  let query: MediaQueryList | undefined;
  const onRatio = () => {
    watchRatio();
    viewer.schedule();
  };
  const watchRatio = () => {
    query?.removeEventListener('change', onRatio);
    query = win?.matchMedia?.(`(resolution: ${win.devicePixelRatio}dppx)`);
    query?.addEventListener('change', onRatio);
  };
  watchRatio();

  // a web font arrived: what was measured is measured again
  const fontFaces = doc.fonts;
  const onFonts = () => viewer.invalidateText();
  fontFaces?.addEventListener?.('loadingdone', onFonts);

  return () => {
    for (const [type, fn, options] of listeners) container.removeEventListener(type, fn, options);
    observer?.disconnect();
    themeObserver?.disconnect();
    query?.removeEventListener('change', onRatio);
    fontFaces?.removeEventListener?.('loadingdone', onFonts);
    viewer.useBackend(undefined);
    container.remove();
  };
}
