/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

/**
 * A long list that renders only what is on screen (and `overscan` rows around it): rows of one height, in a
 * `ScrollArea`. A tree is a list here too: its shown nodes, flattened (`flattenTree`), so 5,000 cells cost what 40 do,
 * expanding one only recomputes the list, and scrolling swaps a few rows. Rows are keyed: one that stays in the window
 * keeps its nodes, so it should follow its own record (by id) for changes.
 */
export function VirtualList<T>(props: {
  items: Observable<readonly T[]>;
  key: (item: T) => unknown;
  /** Every row's height (px). */
  rowHeight: number;
  overscan?: number;
  children: (item: T) => Child;
  class?: string;
  role?: string;
  'aria-label'?: string;
}) {
  const overscan = props.overscan ?? 8;
  const top = observableState(0);
  const height = observableState(0);
  let viewport: HTMLElement | undefined;
  const measure = () => {
    if (!viewport) return;
    top.set(viewport.scrollTop);
    if (viewport.clientHeight !== height.value) height.set(viewport.clientHeight);
  };
  // The window: first and last index on screen, padded; changes only when a row enters or leaves it.
  const range = combineLatest([props.items, top.asReadonly, height.asReadonly]).pipe(
    map(([items, t, h]) => {
      const first = Math.max(0, Math.floor(t / props.rowHeight) - overscan);
      // Before the first measure (no height yet), a screenful.
      const shown = Math.ceil((h || 800) / props.rowHeight) + 2 * overscan;
      return { items, first, last: Math.min(items.length, first + shown) };
    }),
    distinctUntilChanged((a, b) => a.items === b.items && a.first === b.first && a.last === b.last),
  );
  const shown = range.pipe(map((r) => r.items.slice(r.first, r.last)));
  return (
    <ScrollArea
      class={`${classes.list} ${props.class ?? ''}`}
      contentClass={classes.content}
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
      style={{ height: props.items.pipe(map((items) => `${items.length * props.rowHeight}px`)) }}
    >
      <div
        class={classes.window}
        style={{ transform: range.pipe(map((r) => `translateY(${r.first * props.rowHeight}px)`)) }}
      >
        <For each={shown} key={props.key}>
          {(item) => props.children(item)}
        </For>
      </div>
    </ScrollArea>
  );
}

/** A tree's shown nodes as a list (depth first, children of the expanded ones), each with its depth. */
export function flattenTree<N>(
  roots: readonly N[],
  children: (node: N) => readonly N[] | undefined,
  expanded: (node: N) => boolean,
): { node: N; depth: number }[] {
  const out: { node: N; depth: number }[] = [];
  const walk = (nodes: readonly N[], depth: number) => {
    for (const node of nodes) {
      out.push({ node, depth });
      const kids = children(node);
      if (kids?.length && expanded(node)) walk(kids, depth + 1);
    }
  };
  walk(roots, 0);
  return out;
}
