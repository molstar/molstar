/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import { For } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { bind, flush, flushFrame, frameState, onCleanup } from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import { derive } from '../../reactive/state.js';
import { observeResize, on } from '../../dom/scoped.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/virtual-list.module.css.js';
import { NativeMs, rowsAhead } from './window.js';

/** What a `VirtualList` lets its owner do (its `ref`). */
export interface VirtualListHandle {
  /** Brings the row at `index` into view (to the nearest edge) and renders it before returning. */
  scrollToIndex(index: number): void;
}

/**
 * A long list that renders only what is on screen (and a screenful around it, at least `overscan` rows): rows of one height, in a
 * `ScrollArea`. A tree is a list here too: its shown nodes, flattened (`flattenTree`), so 5,000 cells cost what 40 do,
 * expanding one only recomputes the list, and scrolling swaps a few rows. Rows are keyed: one that stays in the window
 * keeps its nodes while its item is the same value (`keepRows`), so it should follow its own record (by id) for changes.
 *
 * Rows are of one height, measured: the window's height over its rows, read from a `ResizeObserver` entry (no layout
 * read). `rowHeight` holds until the first rows lay out, so rows that grow in compact mode (`--msp-touch-target`) need
 * nothing. What would make a row taller (an inline menu, a body opened under it) goes elsewhere (a popover, the
 * inspector), and a tree's extra lines (Show all, Loading…) are rows of their own.
 */
export function VirtualList<T>(props: {
  items: Observable<readonly T[]>;
  key: (item: T) => unknown;
  /** A row's height (px) until the first rows lay out; measured from then on. */
  rowHeight: number;
  overscan?: number;
  children: (item: T) => Child;
  class?: string;
  /** Reaches into its parent's side padding (`ScrollArea`'s `edge`): a list that fills a panel. */
  edge?: boolean;
  role?: string;
  'aria-label'?: string;
  /** The element that holds the rows (the one with `role`): a tree binds its keyboard there (`TreeNavigation`). */
  contentRef?: (el: HTMLElement) => void;
  ref?: (list: VirtualListHandle) => void;
}) {
  const overscan = props.overscan ?? 8;
  // The scroll position in the frame lane (rendering.md §2.3): scroll events can come several a frame; the window
  // moves once a frame, in the frame's tick (before its paint).
  const top = frameState(0);
  const height = observableState(0);
  const rowHeight = observableState(props.rowHeight);
  let viewport: HTMLElement | undefined;
  const measure = () => {
    if (!viewport) return;
    // Both read before either is set: setting `top` makes rows, and a read after that would lay them out here, in the
    // scroll event (measured: 1.7 of a jump's 4.4 ms), rather than once at the frame.
    const t = viewport.scrollTop;
    const h = viewport.clientHeight;
    top.offer(t);
    if (h !== height.value) height.set(h);
  };
  // Before the first measure (no height yet), a screenful.
  const viewHeight = (h: number) => h || 800;
  // The window: the rows in the DOM, [first, last). It stays while what shows (and a few rows past it) is in it, so a
  // scroll within it costs nothing; past that, it moves ahead of the scroll by a screenful and keeps `overscan` rows
  // behind: the browser scrolls before the scroll event reaches the list (a wheel or a flick moves the content on the
  // compositor first), so what it uncovers must already be there or the rows blink in a frame late. A jump (a
  // scrollbar's drag, a far wheel step) makes only what shows and `overscan` rows around it, about a screenful, and the
  // window grows a screenful each way once the scrolling rests: making a row costs, and a drag makes a window a frame.
  const win = observableState({ first: 0, last: 0 });
  let items: readonly unknown[] = [];
  let lastTop = 0;
  // When the wheel or a touch last moved it: the compositor scrolls those before the list hears (a thumb's drag and
  // `scrollToIndex` set the position themselves, and need no margin for speed).
  let nativeAt = Number.NEGATIVE_INFINITY;
  const native = () => {
    nativeAt = performance.now();
  };
  let resting: ReturnType<typeof setTimeout> | undefined;
  const RestMs = 120;
  const place = (rest: boolean) => {
    const n = items.length;
    const rh = rowHeight.value;
    const t = top.value;
    const screenFirst = Math.min(n, Math.floor(t / rh));
    const screenLast = Math.min(n, Math.ceil((t + viewHeight(height.value)) / rh));
    const screen = screenLast - screenFirst;
    const down = t >= lastTop;
    // Its speed (px since the last frame: the window is placed once a frame, the frame lane), when the compositor moves
    // it.
    const speed = performance.now() - nativeAt < NativeMs ? Math.abs(t - lastTop) : 0;
    lastTop = t;
    const { first: wf, last: wl } = win.value;
    // Ahead of the scroll, at least half a screen stays made (a wheel's step or a flick is about that much in one frame),
    // more at speed (`rowsAhead`, window.ts); behind it, `overscan` rows.
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
      // What shows and two rows (a drag makes a window a frame, so it makes little more than it shows); ahead, where the
      // scroll will be at its speed (a wheel's large steps).
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
  const scrolled = () => {
    place(false);
    if (resting !== undefined) clearTimeout(resting);
    resting = setTimeout(() => {
      resting = undefined;
      place(true);
    }, RestMs);
  };
  onCleanup(() => {
    if (resting !== undefined) clearTimeout(resting);
  });
  bind(combineLatest([props.items, height.asReadonly, rowHeight.asReadonly]), ([next]) => {
    items = next;
    place(true);
  });
  bind(top, scrolled);
  const range = combineLatest([props.items, win.asReadonly, rowHeight.asReadonly]).pipe(
    map(([list, w, rh]) => ({
      items: list,
      first: Math.min(w.first, list.length),
      last: Math.min(w.last, list.length),
      rh,
    })),
    distinctUntilChanged((a, b) => a.items === b.items && a.first === b.first && a.last === b.last && a.rh === b.rh),
  );
  // How many rows the window holds now (the measure's divisor).
  let count = 0;
  const shown = range.pipe(
    map((r) => {
      count = r.last - r.first;
      return r.items.slice(r.first, r.last);
    }),
  );
  props.ref?.({
    scrollToIndex(index) {
      const rh = rowHeight.value;
      const view = viewHeight(height.value);
      const y = index * rh;
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
  });
  return (
    <ScrollArea
      class={`${classes.list} ${props.class ?? ''}`}
      contentClass={classes.content}
      edge={props.edge}
      fade={false}
      viewportRef={(el) => {
        viewport = el;
        // Scoped (dom/scoped.ts): detached while hidden, gone with the list.
        on(el, 'scroll', measure, { passive: true });
        on(el, 'wheel', native, { passive: true });
        on(el, 'touchmove', native, { passive: true });
        observeResize(el, measure);
        queueMicrotask(measure);
      }}
      role={props.role}
      aria-label={props['aria-label']}
      ref={(el: HTMLElement) => props.contentRef?.(el)}
      style={{
        height: combineLatest([props.items, rowHeight.asReadonly]).pipe(map(([items, rh]) => `${items.length * rh}px`)),
        // The skeleton's row height (virtual-list.module.css): what the browser uncovers before the rows are made.
        '--msp-virtual-row': derive(rowHeight.asReadonly, (rh) => `${rh}px`),
      }}
    >
      <div
        class={classes.window}
        style={{ transform: range.pipe(map((r) => `translateY(${r.first * r.rh}px)`)) }}
        ref={(el: HTMLElement) => {
          // The rows' height: the window's over its rows, from the observer's entry (after layout, nothing read).
          observeResize(el, (entries) => {
            const h = entries[entries.length - 1].contentRect.height;
            if (count <= 0 || h <= 0) return;
            const rh = h / count;
            if (Math.abs(rh - rowHeight.value) > 0.5) rowHeight.set(rh);
          });
        }}
      >
        <For each={shown} key={props.key}>
          {(item) => props.children(item)}
        </For>
      </div>
    </ScrollArea>
  );
}

/** A shown node of a tree, flattened: its depth (0 at the roots) and its place among its siblings. */
export interface FlatNode<N> {
  readonly node: N;
  readonly depth: number;
  /** Where it is among its siblings (1-based) and how many they are: `aria-posinset`, `aria-setsize`. */
  readonly posinset: number;
  readonly setsize: number;
}

/** A tree's shown nodes as a list (depth first, children of the expanded ones), each with its depth and place. */
export function flattenTree<N>(
  roots: readonly N[],
  children: (node: N) => readonly N[] | undefined,
  expanded: (node: N) => boolean,
): FlatNode<N>[] {
  const out: FlatNode<N>[] = [];
  const walk = (nodes: readonly N[], depth: number) => {
    for (let i = 0; i < nodes.length; i++) {
      const node = nodes[i];
      out.push({ node, depth, posinset: i + 1, setsize: nodes.length });
      const kids = children(node);
      if (kids?.length && expanded(node)) walk(kids, depth + 1);
    }
  };
  walk(roots, 0);
  return out;
}

/**
 * Keeps a flattened list's rows by key: a row the same as the one before it under its key (`same`) is that object
 * again, so a keyed list (`VirtualList`, `For`) keeps its nodes and scope and only the rows that changed render anew.
 * Call it on each new list; it remembers the last.
 */
export function keepRows<R>(key: (row: R) => unknown, same: (a: R, b: R) => boolean): (rows: readonly R[]) => R[] {
  let kept = new Map<unknown, R>();
  return (rows) => {
    const next = new Map<unknown, R>();
    const out = rows.map((row) => {
      const k = key(row);
      const old = kept.get(k);
      const r = old !== undefined && same(old, row) ? old : row;
      next.set(k, r);
      return r;
    });
    kept = next;
    return out;
  };
}
