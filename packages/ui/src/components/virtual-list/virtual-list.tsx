/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import { type Child, For, observeResize, onCleanup, restartable, setAttr, watch } from '../../dom/index.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/virtual-list.module.css.js';
import { listWindow, revealMs, type VirtualListHandle } from './list-window.js';
import { insertedRun } from './window.js';

export { DefaultRowHeight, type VirtualListHandle } from './list-window.js';

/**
 * A long list that renders only what is on screen (and a screenful around it, at least `overscan` rows): rows of
 * one height (a row, `--msui-row`, by default), in a `ScrollArea`. A tree is a list here too: its shown nodes,
 * flattened (`flattenTree`), so 5,000 cells cost what 40 do, expanding one only recomputes the list, and scrolling swaps
 * a few rows. Rows are keyed: one that stays in the window keeps its nodes while its item is the same value
 * (`keepRows`), so it should follow its own record (by id) for changes. The window is `listWindow`'s (list-window.ts).
 *
 * Rows are of one height, measured: the window's height over its rows, read from a `ResizeObserver` entry (no layout
 * read). `rowHeight` holds until the first rows lay out, so rows a host makes taller (a touch layout) need nothing.
 * What would make a row taller (an inline menu, a body opened under it) goes elsewhere (a popover, the inspector), and a
 * tree's extra lines (Show all, Loading…) are rows of their own.
 */
export function VirtualList<T>(props: {
  items: ReadonlyState<readonly T[]>;
  /**
   * An item's key: a row stays its item's while the key stays in the list. Not `key`: the JSX runtime keeps `key` for
   * itself (the compiler passes it apart), so a component never gets it.
   */
  itemKey: (item: T) => unknown;
  /** A row's height (px) until the first rows lay out (default 24, `--msui-row`); measured from then on. */
  rowHeight?: number;
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
  /**
   * Rows inserted in one run (a tree's node opened, its children loaded) grow in (reset.css' `data-reveal`), the rows
   * below gliding down, as a tree's rows open. Any other change, and a scroll, is instant.
   */
  reveal?: boolean;
  /**
   * Each row says its place in the whole list (`aria-posinset`, `aria-setsize` on its element), not in the window: a
   * flat list's options (a listbox). A tree's rows say their place among their siblings themselves (`flattenTree`).
   */
  positions?: boolean;
}) {
  let windowEl: HTMLElement | undefined;
  // A reveal (`reveal`): the keys of the run inserted, whose rows grow in as they are made; the rows' height is not
  // measured meanwhile (a growing row is not a row's height). It ends after its time, or at a scroll.
  const entering = new Set<unknown>();
  const grown: Element[] = [];
  let revealing = false;
  function endReveal() {
    if (!revealing) return;
    revealing = false;
    revealEnd.cancel();
    entering.clear();
    for (const el of grown) el.removeAttribute('data-reveal');
    grown.length = 0;
  }
  // Held while the list is paused (dom §6), so a hidden list ends its reveal when it shows again.
  const revealEnd = restartable(endReveal);
  const list = listWindow({
    items: props.items,
    rowHeight: props.rowHeight,
    overscan: props.overscan,
    onScroll: () => endReveal(),
  });
  const { rowHeight } = list;
  // Compared where the rows are cut (`range`, read before the rows are made), with the rows last cut.
  const keyOf = (item: T): unknown => props.itemKey(item);
  let laid: readonly T[] = [];
  const noteRun = (rows: readonly T[]) => {
    if (!props.reveal || rows === laid) return;
    const run = insertedRun(laid, rows, keyOf);
    laid = rows;
    if (!run) return;
    endReveal();
    const ms = revealMs(windowEl);
    if (!ms) return;
    for (let i = run.at; i < run.at + run.count; i++) entering.add(keyOf(rows[i]));
    revealing = true;
    revealEnd.start(ms);
  };
  onCleanup(endReveal);
  const range = derive(
    [props.items, list.range, rowHeight],
    (rows, w, rh) => {
      noteRun(rows);
      return { items: rows, first: Math.min(w.first, rows.length), last: Math.min(w.last, rows.length), rh };
    },
    (a, b) => a.items === b.items && a.first === b.first && a.last === b.last && a.rh === b.rh,
  );
  const shown = derive(range, (r) => r.items.slice(r.first, r.last));
  props.ref?.(list.handle);
  return (
    <ScrollArea
      class={`${classes.list} ${props.class ?? ''}`}
      contentClass={classes.content}
      edge={props.edge}
      fade={false}
      viewportRef={(el) => list.attach(el)}
      role={props.role}
      aria-label={props['aria-label']}
      ref={(el: HTMLElement) => props.contentRef?.(el)}
      style={{
        height: derive([props.items, rowHeight], (items, rh) => `${items.length * rh}px`),
        // The skeleton's row height (virtual-list.module.css): what the browser uncovers before the rows are made.
        '--msui-virtual-row': derive(rowHeight, (rh) => `${rh}px`),
      }}
    >
      <div
        class={classes.window}
        style={{ transform: derive(range, (r) => `translateY(${r.first * r.rh}px)`) }}
        ref={(el: HTMLElement) => {
          windowEl = el;
          // The rows' height: the window's over the rows made (not the rows planned: `For` time-slices new rows, and
          // a window with fewer than planned would read short), from the observer's entry (after layout, nothing read).
          observeResize(el, (entries) => {
            const h = entries[entries.length - 1].contentRect.height;
            const count = el.childElementCount;
            if (count <= 0 || h <= 0 || revealing) return;
            list.setRowHeight(h / count);
          });
          // Each row's place in the whole list, in window order (the rows are in it in order).
          if (props.positions)
            watch([range, shown], () => {
              const r = range.value;
              const total = String(r.items.length);
              let i = r.first;
              for (const row of el.children) {
                setAttr(row, 'aria-posinset', String(++i));
                setAttr(row, 'aria-setsize', total);
              }
            });
        }}
      >
        <For each={shown} by={props.itemKey}>
          {(item) => {
            const row = props.children(item);
            if (entering.size && row instanceof Element && entering.has(keyOf(item))) {
              row.setAttribute('data-reveal', '');
              grown.push(row);
            }
            return row;
          }}
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
