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
import { onCleanup } from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/virtual-list.module.css.js';

/** What a `VirtualList` lets its owner do (its `ref`). */
export interface VirtualListHandle {
  /** Brings the row at `index` into view (to the nearest edge) and renders it before returning. */
  scrollToIndex(index: number): void;
}

/**
 * A long list that renders only what is on screen (and `overscan` rows around it): rows of one height, in a
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
  const top = observableState(0);
  const height = observableState(0);
  const rowHeight = observableState(props.rowHeight);
  let viewport: HTMLElement | undefined;
  const measure = () => {
    if (!viewport) return;
    top.set(viewport.scrollTop);
    if (viewport.clientHeight !== height.value) height.set(viewport.clientHeight);
  };
  // Before the first measure (no height yet), a screenful.
  const viewHeight = (h: number) => h || 800;
  // The window: first and last index on screen, padded; changes only when a row enters or leaves it.
  const range = combineLatest([props.items, top.asReadonly, height.asReadonly, rowHeight.asReadonly]).pipe(
    map(([items, t, h, rh]) => {
      const first = Math.max(0, Math.floor(t / rh) - overscan);
      const shown = Math.ceil(viewHeight(h) / rh) + 2 * overscan;
      return { items, first, last: Math.min(items.length, first + shown), rh };
    }),
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
      if (t === top.value) return;
      if (viewport) viewport.scrollTop = t;
      // Now, not at the scroll event: the row is in the DOM when this returns (a tree's keyboard focuses it next).
      top.set(t);
    },
  });
  return (
    <ScrollArea
      class={`${classes.list} ${props.class ?? ''}`}
      contentClass={classes.content}
      edge={props.edge}
      viewportRef={(el) => {
        viewport = el;
        el.addEventListener('scroll', measure, { passive: true });
        const resize = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(measure) : undefined;
        resize?.observe(el);
        onCleanup(() => {
          el.removeEventListener('scroll', measure);
          resize?.disconnect();
        });
        queueMicrotask(measure);
      }}
      role={props.role}
      aria-label={props['aria-label']}
      ref={(el: HTMLElement) => props.contentRef?.(el)}
      style={{
        height: combineLatest([props.items, rowHeight.asReadonly]).pipe(map(([items, rh]) => `${items.length * rh}px`)),
      }}
    >
      <div
        class={classes.window}
        style={{ transform: range.pipe(map((r) => `translateY(${r.first * r.rh}px)`)) }}
        ref={(el: HTMLElement) => {
          // The rows' height: the window's over its rows, from the observer's entry (after layout, nothing read).
          if (typeof ResizeObserver === 'undefined') return;
          const rows = new ResizeObserver((entries) => {
            const h = entries[entries.length - 1].contentRect.height;
            if (count <= 0 || h <= 0) return;
            const rh = h / count;
            if (Math.abs(rh - rowHeight.value) > 0.5) rowHeight.set(rh);
          });
          rows.observe(el);
          onCleanup(() => rows.disconnect());
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
