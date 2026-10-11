/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable, of } from 'rxjs';
import { type ObservableState, observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import type { Child, Component } from './jsx-types.js';
import { insert, insertRegion, Passthrough, removeBetween } from './render.js';
import { bind, type MaybeReactive, onCleanup, requireScope, runInScope, type Scope, toStream } from './scope.js';

function toObservable<T>(value: MaybeReactive<T>): Observable<T> {
  return toStream(value) ?? of(value as T);
}

/** A fragment holding a region that follows `source`. */
function region(source: Observable<Child>): Node {
  const fragment = document.createDocumentFragment();
  insertRegion(fragment, source, null);
  return fragment;
}

/**
 * Renders `children` with the value while it is truthy, else `fallback`. The content is rendered again (in a fresh
 * scope) when the value changes identity; for a flag, only when it flips.
 */
export function Show<T>(props: {
  when: MaybeReactive<T>;
  fallback?: () => Child;
  children: (value: NonNullable<T>) => Child;
}): Node {
  return region(
    toObservable(props.when).pipe(
      distinctUntilChanged(),
      map((v) => (v ? () => props.children(v as NonNullable<T>) : (props.fallback ?? null))),
    ),
  );
}

/** Renders a component that may change; props are passed as they are. */
export function Dynamic<P>(props: { component: MaybeReactive<Component<P> | undefined> } & P): Node {
  const { component, ...rest } = props;
  return region(
    toObservable(component).pipe(
      distinctUntilChanged(),
      map((C) => (C ? () => C(rest as P) : null)),
    ),
  );
}

interface MatchCaseValue {
  readonly [Passthrough]: 'match';
  readonly when: MaybeReactive<unknown>;
  readonly children: (value: any) => Child;
}

/** A case of a `Switch`: rendered when it is the first case whose `when` is truthy. */
export function Match<T>(props: { when: MaybeReactive<T>; children: (value: NonNullable<T>) => Child }): Child {
  return { [Passthrough]: 'match', when: props.when, children: props.children } as MatchCaseValue as unknown as Child;
}

/** Renders the first `Match` whose `when` is truthy, else `fallback`. */
export function Switch(props: { fallback?: () => Child; children: Child }): Node {
  const cases = (Array.isArray(props.children) ? props.children : [props.children]).filter(
    (c): c is MatchCaseValue => !!c && typeof c === 'object' && (c as Partial<MatchCaseValue>)[Passthrough] === 'match',
  ) as MatchCaseValue[];
  const active = combineLatest(cases.map((c) => toObservable(c.when))).pipe(
    map((values) => {
      const index = values.findIndex((v) => !!v);
      return { index, value: index < 0 ? undefined : values[index] };
    }),
    distinctUntilChanged((a, b) => a.index === b.index && Object.is(a.value, b.value)),
  );
  return region(
    active.pipe(map(({ index, value }) => (index < 0 ? (props.fallback ?? null) : () => cases[index].children(value)))),
  );
}

/**
 * Time slicing for long lists. Every `For` in a task shares one slice (`SliceMs`) and one count: once `SliceFrom` new
 * entries were made in the slice (by any list: a tree's levels are short lists, its whole is long) and the slice is
 * spent, the rest wait for the next frame (a timer when the page is hidden,
 * which gets no frames), so a 5,000-row tree never holds the main thread for long. Only new entries wait: existing ones
 * update at once, so nothing on screen flickers.
 */
const SliceMs = 12;
const SliceFrom = 64;
let sliceEnd = 0;
let sliceMade = 0;
const waiting = new Set<() => void>();
let flushing = false;

/** Counts a new entry; whether the slice still has room for it (it starts one when none is running). */
function sliceRoom(): boolean {
  if (sliceEnd === 0) {
    sliceEnd = performance.now() + SliceMs;
    sliceMade = 0;
    // The slice is this task's: it ends with it.
    queueMicrotask(() => {
      sliceEnd = 0;
    });
  }
  return ++sliceMade <= SliceFrom || performance.now() < sliceEnd;
}

function nextSlice(fn: () => void) {
  waiting.add(fn);
  if (flushing) return;
  flushing = true;
  const flush = () => {
    flushing = false;
    const fns = [...waiting];
    waiting.clear();
    for (const f of fns) f();
  };
  if (typeof requestAnimationFrame === 'function' && !document.hidden) requestAnimationFrame(flush);
  else setTimeout(flush, 0);
}

interface ForEntry<T> {
  key: unknown;
  item: T;
  scope: Scope;
  index: ObservableState<number>;
  /** The entry's first and last node (the same node for a single element). */
  first: Node;
  last: Node;
}

/**
 * Renders a list, keyed: an item keeps its nodes and scope as long as its key stays in the list (by default the item
 * itself), so moving or updating other items does not touch it. An item whose key stays but which is a different value
 * is rendered again. `children` gets the item's position as a state.
 */
export function For<T>(props: {
  each: MaybeReactive<readonly T[]>;
  key?: (item: T) => unknown;
  fallback?: () => Child;
  children: (item: T, index: ReadonlyState<number>) => Child;
}): Node {
  const owner = requireScope('For');
  const keyOf = props.key ?? ((item: T) => item);
  const start = document.createComment('for');
  const end = document.createComment('/for');
  const fragment = document.createDocumentFragment();
  fragment.append(start, end);
  let entries = new Map<unknown, ForEntry<T>>();
  let fallback: Scope | undefined;
  let warned = false;

  const create = (key: unknown, item: T, i: number): ForEntry<T> => {
    const scope = owner.child();
    const index = observableState(i);
    const content = document.createDocumentFragment();
    runInScope(scope, () => insert(content, props.children(item, index.asReadonly), null));
    let first: Node;
    let last: Node;
    if (content.childNodes.length === 1) {
      first = last = content.firstChild!;
    } else {
      // Several (or no) nodes: anchors keep the entry's place and extent.
      first = document.createComment('');
      last = document.createComment('');
      content.insertBefore(first, content.firstChild);
      content.appendChild(last);
    }
    return { key, item, scope, index, first, last };
  };

  const nodesOf = (entry: ForEntry<T>): Node[] => {
    const nodes: Node[] = [entry.first];
    let node: Node | null = entry.first;
    while (node !== entry.last && (node = node.nextSibling)) nodes.push(node);
    return nodes;
  };

  const remove = (entry: ForEntry<T>) => {
    entry.scope.dispose();
    for (const node of nodesOf(entry)) node.parentNode?.removeChild(node);
  };

  // The latest list, and whether some of it waits for the next slice (a long list renders over several frames).
  let latest: readonly T[] = [];
  let continuing = false;
  let disposed = false;
  onCleanup(() => {
    disposed = true;
  });
  const continueLater = () => {
    if (continuing) return;
    continuing = true;
    nextSlice(() => {
      continuing = false;
      if (!disposed) reconcile(latest);
    });
  };

  const reconcile = (items: readonly T[]) => {
    latest = items;
    const parent = end.parentNode!;
    const next = new Map<unknown, ForEntry<T>>();
    const order: ForEntry<T>[] = [];
    let deferred = false;
    for (let i = 0; i < items.length; i++) {
      const item = items[i];
      const key = keyOf(item);
      if (next.has(key)) {
        if (!warned) {
          warned = true;
          console.warn('For: duplicate key; later items with the same key are not rendered.', key);
        }
        continue;
      }
      let entry = entries.get(key);
      if (entry && entry.item !== item) {
        remove(entry);
        entry = undefined;
      }
      if (!entry) {
        // Many new entries: past the first ones, they wait once the slice is spent.
        if (!sliceRoom()) {
          deferred = true;
          continue;
        }
        entry = create(key, item, i);
      }
      entry.index.set(i);
      next.set(key, entry);
      order.push(entry);
    }
    for (const [key, entry] of entries) if (!next.has(key)) remove(entry);
    entries = next;

    // Place entries from the end: an entry already right before the next one stays where it is.
    let ref: Node = end;
    for (let i = order.length - 1; i >= 0; i--) {
      const entry = order[i];
      if (entry.last.nextSibling !== ref || entry.last.parentNode !== parent) {
        for (const node of nodesOf(entry)) parent.insertBefore(node, ref);
      }
      ref = entry.first;
    }

    if (order.length === 0 && props.fallback && !fallback) {
      fallback = owner.child();
      const content = document.createDocumentFragment();
      runInScope(fallback, () => insert(content, props.fallback!(), null));
      parent.insertBefore(content, end);
    } else if (order.length > 0 && fallback) {
      fallback.dispose();
      fallback = undefined;
      // The fallback sits between the start anchor and the first entry.
      removeBetween(start, order[0].first);
    }
    if (deferred) continueLater();
  };
  bind(props.each, reconcile);
  return fragment;
}
