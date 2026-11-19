/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  CanvasBackend,
  type CanvasLike,
  FontCache,
  forcedColorsActive,
  forcedColorTokens,
  readCssTokens,
  readSystemColors,
  ThemeAttribute,
} from '@molstar/seqx-viewer/canvas';
import { Announcer } from './announcer.js';
import { createDragSession } from './drag-session.js';
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
  const tip = doc.createElement('canvas');
  Object.assign(base.style, { display: 'block' });
  for (const c of [overlay, labels, tip])
    Object.assign(c.style, { position: 'absolute', left: '0', top: '0', pointerEvents: 'none' });
  // the one polite live region the mount owns: what the announcer says, visually hidden (`announcer.ts`)
  const live = doc.createElement('div');
  live.setAttribute('role', 'status');
  live.setAttribute('aria-live', 'polite');
  live.setAttribute('aria-atomic', 'true');
  Object.assign(live.style, {
    position: 'absolute',
    width: '1px',
    height: '1px',
    overflow: 'hidden',
    clip: 'rect(0 0 0 0)',
    whiteSpace: 'nowrap',
    pointerEvents: 'none',
  });
  container.append(base, overlay, labels, tip, live);
  target.appendChild(container);
  const announcer = new Announcer(viewer);
  // cleared, then set: the same words twice are heard twice
  const stopSaying = announcer.message.listen((said) => {
    live.textContent = '';
    live.textContent = said.text;
  });

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
    { fonts, fade: viewer.context.renderConfig.fade, tip: tip as unknown as CanvasLike },
  );
  viewer.useBackend(backend);

  const input = viewer.input;
  // the cursor over a window (a rail's grab zone, a ring's window): `grab`, `grabbing` while dragged
  const cursor = () => {
    container.style.cursor = input.cursor ?? '';
  };
  // a drag ends wherever the button is released: the window's own listeners while it runs (drag-session.ts)
  const drag = createDragSession(container, () => {
    input.up();
    cursor();
  });
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
        cursor();
      },
      undefined,
    ],
    ['pointerleave', () => input.leave(), undefined],
    [
      'pointerdown',
      (e: PointerEvent) => {
        // a mouse's side buttons (3: back, 4: forward) go through the bindings; taken, the page does not go back
        if (e.button === 3 || e.button === 4) {
          const [x, y] = local(e);
          if (input.sideButton(e.button === 3 ? 'back' : 'forward', x, y, e)) e.preventDefault();
          return;
        }
        if (e.button !== 0) return;
        // a release that never arrived: the earlier drag ends before this press begins
        drag.end();
        const [x, y] = local(e);
        if (input.down(x, y, e)) {
          container.setPointerCapture?.(e.pointerId);
          drag.start();
        }
        cursor();
        container.focus({ preventScroll: true });
      },
      undefined,
    ],
    [
      'pointerup',
      (e: PointerEvent) => {
        drag.end();
        input.up();
        cursor();
        if (container.hasPointerCapture?.(e.pointerId)) container.releasePointerCapture(e.pointerId);
      },
      undefined,
    ],
    // the page's own back and forward on a side button's release, once the viewer took it
    [
      'mouseup',
      (e: MouseEvent) => {
        if ((e.button === 3 || e.button === 4) && input.sideButtonTaken(e.button === 3 ? 'back' : 'forward'))
          e.preventDefault();
      },
      undefined,
    ],
    // a drag ends however the pointer is lost: cancelled by the platform, or its capture taken away
    ['pointercancel', () => (drag.end(), input.up()), undefined],
    ['lostpointercapture', () => (drag.end(), input.up()), undefined],
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

  // the platform's reduced motion: every jump and story step is a cut while it holds (`view.motion: platform`)
  const reduced = win?.matchMedia?.('(prefers-reduced-motion: reduce)');
  const onReduced = () => {
    viewer.reducedMotion = !!reduced?.matches;
  };
  onReduced();
  reduced?.addEventListener?.('change', onReduced);

  // forced colors begin or end (the reader's platform setting): the table is read again, from the system colors or not
  // whatever the theme source: the root's CSS (read again, forced inside `readCssTokens`) or the configured theme's
  // table under the system colors
  const forced = win?.matchMedia?.('(forced-colors: active)');
  const configured = viewer.tokens;
  const onForced = () =>
    viewer.setTokens(
      viewer.config.theme.source === 'root'
        ? readCssTokens(container)
        : forcedColorsActive(win)
          ? forcedColorTokens(configured, readSystemColors(container))
          : configured,
    );
  if (viewer.config.theme.source !== 'root' && forcedColorsActive(win)) onForced();
  forced?.addEventListener?.('change', onForced);

  // a web font arrived: what was measured is measured again
  const fontFaces = doc.fonts;
  const onFonts = () => viewer.invalidateText();
  fontFaces?.addEventListener?.('loadingdone', onFonts);

  return () => {
    for (const [type, fn, options] of listeners) container.removeEventListener(type, fn, options);
    drag.dispose();
    observer?.disconnect();
    themeObserver?.disconnect();
    query?.removeEventListener('change', onRatio);
    reduced?.removeEventListener?.('change', onReduced);
    forced?.removeEventListener?.('change', onForced);
    stopSaying();
    announcer.dispose();
    fontFaces?.removeEventListener?.('loadingdone', onFonts);
    viewer.useBackend(undefined);
    container.remove();
  };
}
