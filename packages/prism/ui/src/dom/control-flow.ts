/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Scheduler } from '@molstar/core/task/util/scheduler';
import { combineLatest, distinctUntilChanged, map, Observable, of } from 'rxjs';
import { type ObservableState, observableState } from '../reactive/observable-state.js';
import { sameShallow } from '../reactive/same-shallow.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import type { Child, Component } from './jsx-types.js';
import { playClose, swapIn } from './motion.js';
import { insert, insertRegion, Passthrough, removeBetween } from './render.js';
import { bind, type MaybeReactive, onCleanup, requireScope, runInScope, type Scope, toStream } from './scope.js';

function toObservable<T>(value: MaybeReactive<T>): Observable<T> {
  return toStream(value) ?? of(value as T);
}

/**
 * A fragment holding a region that follows `source`. Its content is kept when the region's subscription is dropped and
 * made again (the scope was paused, `Scope.setActive`): the sources below emit the same render function while their
 * content stays, and the region does not render it twice.
 */
function region(source: Observable<Child>): Node {
  const fragment = document.createDocumentFragment();
  insertRegion(fragment, source, null, true);
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
  /**
   * Content that opens (hidden to shown after the first render; not one shown value replacing another, like a tab
   * switch) is revealed: its element gets `data-reveal`, which `reset.css` animates (it grows to its height and fades
   * in; slides in from the side with `'side'`; only fades in with `'fade'`, for what must not change size, like a band
   * drawn on a canvas). Closing plays it back (collapses, slides out or fades, then goes). None with reduced motion.
   * Overview › Foundations › Motion.
   */
  reveal?: boolean | 'side' | 'fade';
  /**
   * One shown value replaced by another (a panel by another) runs as a view transition scoped to this element, where
   * the browser has them (`swapIn`, motion.ts); opening and closing are `reveal`'s.
   */
  swap?: () => Element | null | undefined;
}): Node {
  let first = true;
  let shown = false;
  // The content's top elements, kept to play a closing back on them; and the shown value with its render function, so
  // that subscribing again (after a pause) to the same value keeps the content.
  let tops: Element[] = [];
  let shownValue: unknown;
  let shownFn: (() => Child) | undefined;
  const content = (v: NonNullable<T>, reveal: boolean) => () => {
    const child = props.children(v);
    tops = topElements(child);
    if (reveal)
      for (const el of tops) el.setAttribute('data-reveal', props.reveal === true ? '' : String(props.reveal));
    return child;
  };
  return region(
    new Observable<Child>((subscriber) => {
      let closing: Animation[] | undefined;
      // The first value of this subscription: the one a subscription made again after a pause starts with.
      let fresh = true;
      const sub = toObservable(props.when)
        .pipe(distinctUntilChanged())
        .subscribe((v) => {
          // Reopened while closing: the closing stops, the content is rendered anew.
          if (closing) for (const a of closing) a.cancel();
          closing = undefined;
          const wasShown = shown;
          const opening = !first && !wasShown && !!v;
          const closingNow = !first && wasShown && !v;
          first = false;
          shown = !!v;
          const resumed = fresh;
          fresh = false;
          if (v) {
            // Resumed with the value it was showing (a value derived anew on subscribing, like `{ section }`, is the
            // same when its fields are): the content stays.
            const same = wasShown && shownFn && (Object.is(v, shownValue) || (resumed && sameShallow(v, shownValue)));
            if (!same) {
              shownValue = v;
              shownFn = content(v as NonNullable<T>, !!props.reveal && opening);
              // A value replacing another: swapped in, if asked (the content rendered when the transition has the old
              // state, unless something newer came by then).
              if (wasShown && props.swap) {
                const fn = shownFn;
                swapIn(props.swap(), () => {
                  if (shownFn === fn) subscriber.next(fn);
                });
                return;
              }
            }
            subscriber.next(shownFn);
            return;
          }
          shownValue = shownFn = undefined;
          const kind = props.reveal === true ? 'grow' : props.reveal;
          const animations = closingNow && kind ? tops.map((el) => playClose(el, kind)) : [];
          if (animations.length === 0 || animations.some((a) => !a)) {
            subscriber.next(props.fallback ?? null);
            return;
          }
          const these = animations as Animation[];
          closing = these;
          // Still on screen for the length of the animation, but gone for the pointer, the keyboard and assistive tech.
          for (const el of tops) setInert(el, true);
          Promise.all(these.map((a) => a.finished)).then(
            () => {
              if (closing !== these || subscriber.closed) return;
              closing = undefined;
              subscriber.next(props.fallback ?? null);
            },
            () => {},
          );
        });
      return () => {
        sub.unsubscribe();
        if (closing) for (const a of closing) a.cancel();
      };
    }),
  );
}

function setInert(el: Element, inert: boolean) {
  el.toggleAttribute('inert', inert);
}

/**
 * Renders what can be hidden while it stays mounted (a panel closed behind another card, a folded section's body) so
 * that hidden UI is **inert** and **quiet**: `render` runs in a child scope that follows `active`, so while `active` is
 * false its bindings hold no subscriptions (nothing is written to the hidden DOM, nothing derived is computed) and
 * `onActive` hooks pause what is imperative; when `active` turns true again they come back, with the current values, in
 * one write each. The DOM, the scopes and with them what lives in them (a scroll position, a field's draft) are kept.
 * Its top elements get the `inert` attribute meanwhile (no pointer, no focus, no tooltips, not announced to assistive
 * tech): bind the element that hides itself (its `hidden`, a class) outside, since what `render` makes pauses.
 * A part created while `active` is false renders its static structure and fills in on the first resume.
 *
 * Content that is removed when it is hidden (`Show`) needs none of this.
 */
export function activity<T extends Child>(active: MaybeReactive<boolean>, render: () => T): T {
  const owner = requireScope('activity');
  const scope = owner.child();
  // Before rendering: a part that starts hidden binds nothing until it is shown.
  scope.follow(active);
  const result = runInScope(scope, render);
  const tops = topElements(result);
  const stream = toStream(active);
  if (!stream) {
    if (!active) for (const el of tops) setInert(el, true);
  } else if (tops.length > 0) {
    // Not `bind`: the flag must keep reaching the elements while an outer scope is paused.
    owner.own(
      stream.subscribe((a) => {
        for (const el of tops) setInert(el, !a);
      }),
    );
  }
  return result;
}

/** The elements a child renders at its top (an element, a list's elements, a fragment's children). */
function topElements(child: Child): Element[] {
  if (child instanceof Element) return [child];
  if (child instanceof DocumentFragment) return [...child.children];
  if (Array.isArray(child)) return child.flatMap((c) => topElements(c as Child));
  return [];
}

/** Renders a component that may change; props are passed as they are. */
export function Dynamic<P>(props: { component: MaybeReactive<Component<P> | undefined> } & P): Node {
  const { component, ...rest } = props;
  let shown: Component<P> | undefined;
  let render: (() => Child) | undefined;
  return region(
    toObservable(component).pipe(
      distinctUntilChanged(),
      map((C) => {
        if (!C) {
          shown = render = undefined;
          return null;
        }
        // The same component as before (a new subscription after a pause): the same function, which is not rendered again.
        if (C !== shown || !render) {
          shown = C;
          render = () => C(rest as P);
        }
        return render;
      }),
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
  let shownIndex = -1;
  let shownValue: unknown;
  let render: (() => Child) | undefined;
  return region(
    active.pipe(
      map(({ index, value }) => {
        if (index < 0) {
          render = undefined;
          shownIndex = -1;
          return props.fallback ?? null;
        }
        // The same case and value as before (a new subscription after a pause): the same function, not rendered again.
        if (index !== shownIndex || !Object.is(value, shownValue) || !render) {
          shownIndex = index;
          shownValue = value;
          render = () => cases[index].children(value);
        }
        return render;
      }),
    ),
  );
}

/**
 * Time slicing for long lists. Every `For` in a task shares one slice (`SliceMs`) and one count: once `SliceFrom` new
 * entries were made in the slice (by any list: a tree's levels are short lists, its whole is long) and the slice is
 * spent, the rest continue in a fresh task, so a 5,000-row tree never holds the main thread for long. The continuation
 * is Mol*'s `setImmediate` (a MessageChannel task: not clamped like a timer, not throttled in a hidden page), so the
 * browser takes input between any two slices and paints at its frames. 12 ms measured best (Internals › Performance,
 * 5,025 cells: 5 ms slices took 40% longer in all, 12 ms matched frame-aligned slicing). Only new entries wait: existing
 * ones update at once, so nothing on screen flickers.
 */
const SliceMs = 12;
const SliceFrom = 64;
let sliceEnd = 0;
let sliceMade = 0;
const waiting = new Set<() => void>();
let flushing = false;

let slicing: boolean | 'auto' = 'auto';

/**
 * Turns time slicing on or off (`auto`, the default: on in a browser, off in jsdom, which never renders, so there is
 * nothing to yield to and tests stay synchronous).
 */
export function timeSlicing(on: boolean | 'auto') {
  slicing = on;
}

function slicingOn(): boolean {
  if (slicing !== 'auto') return slicing;
  return typeof document !== 'undefined' && !/jsdom/i.test(globalThis.navigator?.userAgent ?? '');
}

/** Counts a new entry; whether the slice still has room for it (it starts one when none is running). */
function sliceRoom(): boolean {
  if (!slicingOn()) return true;
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
  Scheduler.setImmediate(flushWaiting);
}

function flushWaiting() {
  flushing = false;
  // A fresh task, a fresh slice (whatever ran before it; under Node the continuation is a nextTick, which runs before
  // the microtask that would end the last slice).
  sliceEnd = 0;
  const fns = [...waiting];
  waiting.clear();
  for (const f of fns) {
    // Spent: the rest wait as they are (running them would only walk their lists to defer them again).
    if (sliceSpent()) nextSlice(f);
    else f();
  }
}

/**
 * The slice, for a list that makes its entries itself (`RecycledList`'s rows): `sliceHasRoom()` counts one new entry
 * with `For`'s and tells whether the slice has room for it; `nextSlice(fn)` runs `fn` in the next one.
 */
export { sliceRoom as sliceHasRoom, nextSlice };

function sliceSpent(): boolean {
  return sliceEnd !== 0 && sliceMade > SliceFrom && performance.now() >= sliceEnd;
}

interface ForEntry<T> {
  key: unknown;
  item: T;
  scope: Scope;
  /** Only when `children` takes the index (most rows do not, and a state per row is not free). */
  index: ObservableState<number> | undefined;
  /** The entry's first and last node (the same node for a single element). */
  first: Node;
  last: Node;
  /** Its place in the list as last rendered (-1: new), for the move diff. */
  pos: number;
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

  const withIndex = props.children.length > 1;
  const create = (key: unknown, item: T, i: number): ForEntry<T> => {
    const scope = owner.child();
    const index = withIndex ? observableState(i) : undefined;
    const content = document.createDocumentFragment();
    runInScope(scope, () => insert(content, props.children(item, index?.asReadonly as ReadonlyState<number>), null));
    let first: Node;
    let last: Node;
    if (content.childNodes.length === 1 && content.firstChild!.nodeType === Node.ELEMENT_NODE) {
      first = last = content.firstChild!;
    } else {
      // Several (or no) nodes, or one that is not an element (a region's placeholder Text, which the region replaces
      // when its content arrives): anchors keep the entry's place and extent.
      first = document.createComment('');
      last = document.createComment('');
      content.insertBefore(first, content.firstChild);
      content.appendChild(last);
    }
    return { key, item, scope, index, first, last, pos: -1 };
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
      // Paused meanwhile: the rest is rendered when the scope resumes (the list subscribes again, and reconciles).
      if (!disposed && owner.active) reconcile(latest);
    });
  };

  // The entries as last placed in full: the same items in another array (a list recomputed to what it was) find them in
  // order, and nothing is done (the quick tier read 33 KB a step for 200 rows before; hot.quick.ts).
  let placed: readonly ForEntry<T>[] = [];
  const unchanged = (items: readonly T[]) => {
    const n = items.length;
    if (n === 0 || n !== placed.length) return false;
    for (let i = 0; i < n; i++) if (placed[i].item !== items[i]) return false;
    return true;
  };

  const reconcile = (items: readonly T[]) => {
    latest = items;
    if (unchanged(items)) return;
    const parent = end.parentNode!;
    // Emptied: every entry goes at once, and its nodes in one range.
    if (items.length === 0 && entries.size > 0) {
      for (const entry of entries.values()) entry.scope.dispose();
      entries = new Map();
      // One range for all their nodes (measured: 2.5× faster than emptying the parent with textContent).
      const range = document.createRange();
      range.setStartAfter(start);
      range.setEndBefore(end);
      range.deleteContents();
    }
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
      entry.index?.set(i);
      next.set(key, entry);
      order.push(entry);
    }
    for (const [key, entry] of entries) if (!next.has(key)) remove(entry);
    entries = next;

    // Place entries with as few moves as possible: those whose old places form the longest increasing run stay where
    // they are; the rest (moved or new) go before the next one, from the end (a swap is two moves, not a cascade).
    const stable = stayingPut(order);
    let ref: Node = end;
    for (let i = order.length - 1; i >= 0; i--) {
      const entry = order[i];
      if (!stable[i] || entry.last.parentNode !== parent) {
        for (const node of nodesOf(entry)) parent.insertBefore(node, ref);
      }
      ref = entry.first;
    }
    for (let i = 0; i < order.length; i++) order[i].pos = i;
    placed = deferred ? [] : order;

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

/**
 * Which entries keep their place: the longest increasing run of their old places (patience sort, O(n log n)); new
 * entries (-1) never do. Already in order (the common case): all of them, without the search.
 */
function stayingPut(order: readonly { pos: number }[]): boolean[] {
  const n = order.length;
  const stable = new Array<boolean>(n).fill(false);
  let sorted = true;
  for (let i = 0, last = -1; i < n; i++) {
    const p = order[i].pos;
    if (p < 0) continue;
    if (p < last) {
      sorted = false;
      break;
    }
    last = p;
  }
  if (sorted) {
    for (let i = 0; i < n; i++) stable[i] = order[i].pos >= 0;
    return stable;
  }
  // tails[k]: the index in `order` ending the best run of length k + 1; prev links runs back.
  const tails: number[] = [];
  const prev = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    const p = order[i].pos;
    if (p < 0) continue;
    let lo = 0;
    let hi = tails.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (order[tails[mid]].pos < p) lo = mid + 1;
      else hi = mid;
    }
    if (lo > 0) prev[i] = tails[lo - 1];
    tails[lo] = i;
  }
  for (let i = tails.length ? tails[tails.length - 1] : -1; i >= 0; i = prev[i]) stable[i] = true;
  return stable;
}
