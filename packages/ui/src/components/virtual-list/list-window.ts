/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import {
  bind,
  type FrameState,
  flush,
  flushFrame,
  frameState,
  observeResize,
  on,
  type Restartable,
  restartable,
} from '../../dom/index.js';
import { state } from '../../state.js';
import { milliseconds } from '../../styles/transitions.js';
import { NativeMs, rowsAhead } from './window.js';

/** A row's height until the rows lay out: `--msui-row` (styles §4). */
export const DefaultRowHeight = 24;

/** What a virtual list lets its owner do (its `ref`). */
export interface VirtualListHandle {
  /** Brings the row at `index` into view (to the nearest edge) and renders it before returning. */
  scrollToIndex(index: number): void;
}

/** The rows in the DOM, [first, last). */
export interface ListRange {
  readonly first: number;
  readonly last: number;
}

/** A list's window over its rows: what both virtual lists share (the scrolling, the window, the rest). */
export interface ListWindow {
  /** The scroll position, in the frame lane (dom §3). */
  readonly top: FrameState<number>;
  readonly height: ReadonlyState<number>;
  /** A row's height: measured once rows lay out. */
  readonly rowHeight: ReadonlyState<number>;
  /** The window: the rows to have in the DOM. */
  readonly range: ReadonlyState<ListRange>;
  /** Sets the measured row height (a list measures its own rows). */
  setRowHeight(h: number): void;
  /** Attaches the scroll, wheel, touch and size listeners to the viewport (scoped, dom §6). */
  attach(viewport: HTMLElement): void;
  readonly handle: VirtualListHandle;
}

/** Before the first measure (no height yet), a screenful. */
export const viewHeight = (h: number) => h || 800;
/** After the last scroll, this long: the scrolling rests. */
const RestMs = 120;

/**
 * The window of a virtual list (react to its items, its height and its rows' height; place itself once a frame as the
 * list scrolls). The window stays while what shows (and a few rows past it) is in it, so a scroll within it costs
 * nothing; past that, it moves ahead of the scroll by a screenful and keeps `overscan` rows behind: the browser scrolls
 * before the scroll event reaches the list (a wheel or a flick moves the content on the compositor first), so what it
 * uncovers must already be there or the rows blink in a frame late. A jump (a scrollbar's drag, a far wheel step) makes
 * only what shows and a little around it, and the window grows a screenful each way once the scrolling rests (a
 * `restartable`: held while the list is paused, dom §6). `onScroll` is told at each placing scroll (a reveal ends),
 * `onRest` when it rests (the nodes into screen order).
 */
export function listWindow(options: {
  items: ReadonlyState<readonly unknown[]>;
  rowHeight?: number;
  overscan?: number;
  onScroll?(): void;
  onRest?(): void;
}): ListWindow {
  const overscan = options.overscan ?? 8;
  const top = frameState(0);
  const height = state(0);
  const rowHeight = state(options.rowHeight ?? DefaultRowHeight);
  const win = state<ListRange>({ first: 0, last: 0 });
  let viewport: HTMLElement | undefined;
  let count = 0;
  let lastTop = 0;
  // When the wheel or a touch last moved it: the compositor scrolls those before the list hears (a thumb's drag and
  // `scrollToIndex` set the position themselves, and need no margin for speed).
  let nativeAt = Number.NEGATIVE_INFINITY;
  const native = () => {
    nativeAt = performance.now();
  };
  const measure = () => {
    if (!viewport) return;
    // Both read before either is set: setting `top` makes rows, and a read after that would lay them out here, in the
    // scroll event, rather than once at the frame.
    const t = viewport.scrollTop;
    const h = viewport.clientHeight;
    top.offer(t);
    if (h !== height.value) height.set(h);
  };
  const place = (rest: boolean) => {
    const n = count;
    const rh = rowHeight.value;
    const t = top.value;
    const screenFirst = Math.min(n, Math.floor(t / rh));
    const screenLast = Math.min(n, Math.ceil((t + viewHeight(height.value)) / rh));
    const screen = screenLast - screenFirst;
    const down = t >= lastTop;
    // Its speed (px since the last frame), when the compositor moves it.
    const speed = performance.now() - nativeAt < NativeMs ? Math.abs(t - lastTop) : 0;
    lastTop = t;
    const { first: wf, last: wl } = win.value;
    // Ahead of the scroll, at least half a screen stays made, more at speed (`rowsAhead`); behind it, `overscan` rows.
    const guardAhead = rowsAhead(speed, rh, screen, Math.max(overscan, Math.ceil(screen / 2)));
    const covered =
      wf <= Math.max(0, screenFirst - (down ? overscan : guardAhead)) &&
      wl >= Math.min(n, screenLast + (down ? guardAhead : overscan)) &&
      wl <= n;
    let first: number;
    let last: number;
    if (rest) {
      first = screenFirst - screen;
      last = screenLast + screen;
    } else if (covered) return;
    else if (screenLast <= wf || screenFirst >= wl) {
      // A jump: what shows and two rows; ahead, where the scroll will be at its speed.
      const fast = rowsAhead(speed, rh, screen, 2);
      first = screenFirst - (down ? 2 : fast);
      last = screenLast + (down ? fast : 2);
    } else {
      const ahead = rowsAhead(speed, rh, screen, Math.max(overscan, screen));
      first = screenFirst - (down ? overscan : ahead);
      last = screenLast + (down ? ahead : overscan);
    }
    first = Math.max(0, first);
    last = Math.min(n, last);
    if (first !== wf || last !== wl) win.set({ first, last });
  };
  const resting: Restartable = restartable(() => {
    place(true);
    options.onRest?.();
  });
  bind(
    derive([options.items, height.asReadonly, rowHeight.asReadonly], (list, h, rh) => ({ list, h, rh })),
    ({ list }) => {
      count = list.length;
      place(true);
    },
  );
  bind(top, () => {
    options.onScroll?.();
    place(false);
    resting.start(RestMs);
  });
  return {
    top,
    height: height.asReadonly,
    rowHeight: rowHeight.asReadonly,
    range: win.asReadonly,
    setRowHeight: (h) => {
      if (h > 0 && Math.abs(h - rowHeight.value) > 0.5) rowHeight.set(h);
    },
    attach(el) {
      viewport = el;
      on(el, 'scroll', measure, { passive: true });
      on(el, 'wheel', native, { passive: true });
      on(el, 'touchmove', native, { passive: true });
      observeResize(el, measure);
      queueMicrotask(measure);
    },
    handle: {
      scrollToIndex(index) {
        const rh = rowHeight.value;
        const view = viewHeight(height.value);
        const y = Math.max(0, Math.min(index, count - 1)) * rh;
        let t = top.value;
        if (y < t) t = y;
        else if (y + rh > t + view) t = y + rh - view;
        if (t !== top.value) {
          // A jump, not a motion: no speed (the margin ahead stays what it is at rest).
          lastTop = t;
          if (viewport) viewport.scrollTop = t;
          // Now, not at the scroll event: the row is in the DOM when this returns (a tree's keyboard focuses it next).
          top.offer(t);
        }
        // And rendered: the frame's tick now (the window), then the `watch`es.
        flushFrame();
        flush();
      },
    },
  };
}

/**
 * How long a reveal's rows take to come in: the motion token on the list (`--msui-motion-base`, a registered time) and
 * a frame's margin for its end; 0 when the token is 0 (reduced motion: no animation, nothing to wait for). Where the
 * styles are not there to read (no token), the base duration.
 */
export function revealMs(el: Element | undefined): number {
  const token =
    el && typeof getComputedStyle === 'function'
      ? getComputedStyle(el).getPropertyValue('--msui-motion-base').trim()
      : '';
  const ms = token ? milliseconds(token) : 200;
  return ms > 0 ? ms + 40 : 0;
}
