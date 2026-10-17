/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { nextSlice, sliceHasRoom } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { insert } from '../../dom/render.js';
import {
  bind,
  flush,
  flushFrame,
  frameState,
  onCleanup,
  requireScope,
  runInScope,
  type Scope,
} from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive } from '../../reactive/state.js';
import { observeResize, on } from '../../dom/scoped.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { classes } from './_generated/virtual-list.module.css.js';
import { insertedRun, NativeMs, RevealMs, rowsAhead } from './window.js';
import type { VirtualListHandle } from './virtual-list.js';

/**
 * A recycled row of a `RecycledList` (its `slot`): its element is made once and shows one item, then another, as the
 * list scrolls (view-models.md §3, recycled rows). `show` is called with the item it shows now: another one (the slot
 * re-pointed), the same key's new value, or undefined (put away: hidden, following nothing). It writes what differs
 * into its element and follows the item's own keys (a `watch` it re-points: `WatchHandle.retarget`).
 */
export interface VirtualSlot<T> {
  readonly el: HTMLElement;
  show(item: T | undefined): void;
}

/** The rest's reorder passes, for the benches (rendering.md §1): how many, the nodes moved, the last one's ms. */
export const RecycledListStats = { settles: 0, moved: 0, lastMs: 0 };

/** A row in the window: a slot showing its item, or a keyed row rendered for it. */
interface Entry<T> {
  key: unknown;
  item: T;
  /** The positioned wrapper (`classes.row`). */
  el: HTMLElement;
  /** Its offset as last written (px). */
  y: number;
  slot?: Slot<T>;
  scope?: Scope;
}

interface Slot<T> {
  el: HTMLElement;
  view: VirtualSlot<T>;
  /** Shown (in the window) or hidden (spare). */
  shown: boolean;
}

/**
 * A long list that renders only what is on screen (and a screenful around it, at least `overscan` rows): rows of one height, in a
 * `ScrollArea`. A tree is a list here too: its shown nodes, flattened (`flattenTree`), so 5,000 cells cost what 40 do,
 * expanding one only recomputes the list, and scrolling swaps a few rows.
 *
 * Rows come in two kinds (view-models.md §3). **Recycled** (`slot`, for the items `recycle` takes, all by default): a
 * pool of slots, as many as the largest window, each with its element and scope for good; a window that moves points
 * the slots that left it at the items that came in (`VirtualSlot.show`), so scrolling creates, destroys and reorders no
 * DOM and allocates next to nothing. **Keyed** (`children`): a row rendered for its item, kept while its key stays in
 * the window and its item is the same value (`keepRows`), so it should follow its own record (by id) for changes. Each
 * row is positioned on its own (`transform`), so the DOM order is not the screen order.
 *
 * Rows are of one height, measured: a row's height, read from a `ResizeObserver` entry (no layout read). `rowHeight`
 * holds until the first rows lay out, so rows that grow in compact mode (`--msp-touch-target`) need nothing. What would
 * make a row taller (an inline menu, a body opened under it) goes elsewhere (a popover, the inspector), and a tree's
 * extra lines (Show all, Loading…) are rows of their own.
 */
export function RecycledList<T>(props: {
  items: ReadonlyState<readonly T[]>;
  /**
   * An item's key: a row stays its item's while the key stays in the window. Not `key`: the JSX runtime keeps `key`
   * for itself (the compiler passes it apart), so a component never gets it.
   */
  itemKey: (item: T) => unknown;
  /** A row's height (px) until the first rows lay out; measured from then on. */
  rowHeight: number;
  overscan?: number;
  /** A keyed row for an item (one `slot` does not take). */
  children?: (item: T) => Child;
  /** A recycled row (made in its own scope, which lives as long as the list). */
  slot?: () => VirtualSlot<T>;
  /** Which items go to slots (with `slot`; all when absent); the others are keyed rows (`children`). */
  recycle?: (item: T) => boolean;
  class?: string;
  /** Reaches into its parent's side padding (`ScrollArea`'s `edge`): a list that fills a panel. */
  edge?: boolean;
  role?: string;
  'aria-label'?: string;
  /** The element that holds the rows (the one with `role`): a tree binds its keyboard there (`TreeNavigation`). */
  contentRef?: (el: HTMLElement) => void;
  ref?: (list: VirtualListHandle) => void;
  /**
   * Rows inserted in one run (a tree's node opened, its children loaded) come in as a tree's did before it was
   * virtualized: they fade in where they belong while the rows below glide down over them. Any other change, and a
   * scroll, is instant.
   */
  reveal?: boolean;
}) {
  const owner = requireScope('RecycledList');
  const overscan = props.overscan ?? 8;
  // The scroll position in the frame lane (rendering.md §2.3): scroll events can come several a frame; the window
  // moves once a frame, in the frame's tick (before its paint).
  const top = frameState(0);
  const height = observableState(0);
  const rowHeight = observableState(props.rowHeight);
  let viewport: HTMLElement | undefined;
  const measure = () => {
    if (!viewport) return;
    // Both read before either is set: setting `top` lays out rows, and a read after would lay them out in the event.
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
  // scrollbar's drag, a far wheel step) shows only what shows and `overscan` rows around it, about a screenful, and the
  // window grows a screenful each way once the scrolling rests: a drag moves the window every frame.
  const win = observableState({ first: 0, last: 0 });
  let items: readonly T[] = [];
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
      // What shows and two rows (a drag moves the window a frame, so it shows little more than that); ahead, where the
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
    endReveal();
    place(false);
    if (resting !== undefined) clearTimeout(resting);
    resting = setTimeout(() => {
      resting = undefined;
      place(true);
      settleOrder();
    }, RestMs);
  };
  // The rows of a reveal (`reveal`): those gliding to their new place and those coming in, their classes taken off when
  // it ends (its time, or a scroll: a slot re-pointed while gliding would slide across the list).
  const revealed: HTMLElement[] = [];
  let revealEnd: ReturnType<typeof setTimeout> | undefined;
  function endReveal() {
    if (revealEnd === undefined) return;
    clearTimeout(revealEnd);
    revealEnd = undefined;
    for (const el of revealed) el.classList.remove(classes.glide, classes.enter);
    revealed.length = 0;
  }
  let animating = false;
  onCleanup(() => {
    if (resting !== undefined) clearTimeout(resting);
    endReveal();
  });
  bind(
    derive([props.items, height.asReadonly, rowHeight.asReadonly], (list, h, rh) => ({ list, h, rh })),
    ({ list: next }) => {
      items = next;
      place(true);
    },
  );
  bind(top, scrolled);

  // --- The rows ---
  let content: HTMLElement | undefined;
  let entries = new Map<unknown, Entry<T>>();
  /** Slots not in the window (hidden), and those that left it in this layout (still showing, until taken or hidden). */
  const spare: Slot<T>[] = [];
  const recycles = (item: T) => !!props.slot && (!props.recycle || props.recycle(item));
  // The rows' height: one row's (the first made, while it lives), from the observer's entry (after layout, nothing read).
  let measured: HTMLElement | undefined;
  const rows = observeResize(undefined, (list) => {
    const h = list[list.length - 1].contentRect.height;
    if (h > 0 && Math.abs(h - rowHeight.value) > 0.5) rowHeight.set(h);
  });
  const observe = (el: HTMLElement) => {
    if (measured || !rows) return;
    measured = el;
    rows.observe(el);
  };
  const wrapper = () => {
    const el = document.createElement('div');
    el.className = classes.row;
    return el;
  };
  const makeSlot = (): Slot<T> => {
    const scope = owner.child();
    const view = runInScope(scope, () => props.slot!());
    const el = wrapper();
    el.append(view.el);
    content!.append(el);
    observe(el);
    return { el, view, shown: true };
  };
  const makeKeyed = (key: unknown, item: T): Entry<T> => {
    const scope = owner.child();
    const el = wrapper();
    runInScope(scope, () => insert(el, props.children?.(item), null));
    content!.append(el);
    observe(el);
    return { key, item, el, y: Number.NaN, scope };
  };
  const dropKeyed = (entry: Entry<T>) => {
    entry.scope!.dispose();
    entry.el.remove();
    if (measured === entry.el) {
      rows?.unobserve(entry.el);
      measured = undefined;
    }
  };
  const position = (entry: Entry<T>, index: number, rh: number) => {
    const y = index * rh;
    if (entry.y === y) return;
    if (animating) {
      // Placed before (glides from there) or new to the window (comes in where it belongs).
      entry.el.classList.add(Number.isNaN(entry.y) ? classes.enter : classes.glide);
      revealed.push(entry.el);
    }
    entry.y = y;
    entry.el.style.transform = `translateY(${y}px)`;
  };
  // Adds the row for `index` (not in the window yet): a spare slot, or one made, or a keyed row; false when this task's
  // slice has no room for a new one.
  const add = (next: Map<unknown, Entry<T>>, index: number, rh: number): boolean => {
    const item = items[index];
    const key = props.itemKey(item);
    if (next.has(key)) return true;
    let entry: Entry<T>;
    if (recycles(item)) {
      let slot = spare.pop();
      if (!slot) {
        if (!sliceHasRoom()) return false;
        slot = makeSlot();
      } else if (!slot.shown) {
        slot.shown = true;
        slot.el.hidden = false;
      }
      slot.view.show(item);
      entry = { key, item, el: slot.el, y: Number.NaN, slot };
    } else {
      if (!sliceHasRoom()) return false;
      entry = makeKeyed(key, item);
    }
    next.set(key, entry);
    position(entry, index, rh);
    return true;
  };
  let continuing = false;
  let warned = false;
  const layout = (first: number, last: number, rh: number) => {
    if (!content) return;
    const next = new Map<unknown, Entry<T>>();
    // What stays in the window stays as it is (moved, or shown again when its item is a new value).
    for (let i = first; i < last; i++) {
      const item = items[i];
      const key = props.itemKey(item);
      if (next.has(key)) {
        if (!warned) {
          warned = true;
          console.warn('VirtualList: duplicate key; later items with the same key are not rendered.', key);
        }
        continue;
      }
      const entry = entries.get(key);
      if (!entry) continue;
      if (entry.item !== item) {
        if (entry.slot && recycles(item)) {
          entry.item = item;
          entry.slot.view.show(item);
        } else continue;
      }
      next.set(key, entry);
      position(entry, i, rh);
    }
    // What left: keyed rows go, slots are spare (still showing what they showed, until taken).
    const leaving = spare.length;
    for (const [key, entry] of entries) {
      if (next.get(key) === entry) continue;
      if (entry.slot) spare.push(entry.slot);
      else dropKeyed(entry);
    }
    entries = next;
    // What came in, what shows first (a slice may run out).
    const t = top.value;
    const screenFirst = Math.max(first, Math.min(last, Math.floor(t / rh)));
    const screenLast = Math.max(screenFirst, Math.min(last, Math.ceil((t + viewHeight(height.value)) / rh)));
    let deferred = false;
    for (let i = screenFirst; i < screenLast && !deferred; i++) if (!add(next, i, rh)) deferred = true;
    for (let i = first; i < screenFirst && !deferred; i++) if (!add(next, i, rh)) deferred = true;
    for (let i = screenLast; i < last && !deferred; i++) if (!add(next, i, rh)) deferred = true;
    // Slots left over are put away: hidden, following nothing.
    for (let i = spare.length - 1; i >= leaving; i--) {
      const slot = spare[i];
      if (!slot.shown) continue;
      slot.shown = false;
      slot.view.show(undefined);
      slot.el.hidden = true;
    }
    if (deferred && !continuing) {
      continuing = true;
      nextSlice(() => {
        continuing = false;
        // Paused meanwhile: the rest comes when the list resumes (its bindings replay, and lay it out again).
        if (!owner.disposed && owner.active) layout(range.first, range.last, range.rh);
      });
    }
  };
  // The window as last laid out (a continuation lays out the latest).
  let range = { first: 0, last: 0, rh: props.rowHeight };
  /** The rows last laid out (what a change of them is compared with, for a reveal). */
  let laid: readonly T[] = [];
  const ranges = derive(
    [props.items, win.asReadonly, rowHeight.asReadonly],
    (list, w, rh) => ({ list, first: Math.min(w.first, list.length), last: Math.min(w.last, list.length), rh }),
    (a, b) => a.list === b.list && a.first === b.first && a.last === b.last && a.rh === b.rh,
  );
  onCleanup(() => {
    for (const entry of entries.values()) if (!entry.slot) entry.scope!.dispose();
    entries.clear();
  });
  // At rest, the slots' nodes into screen order, in one pass (David, 2026-10-17): a screen reader browsing the list
  // reads its rows in order whenever it is still. Never while scrolling (the rest signal, as the margin's). A row
  // moved with the focus in it gets it back (moving a node drops its focus).
  const settleOrder = () => {
    if (!content) return;
    const t0 = performance.now();
    const list = [...entries.values()].sort((a, b) => a.y - b.y);
    const shown = new Set<Node>();
    for (const e of list) shown.add(e.el);
    const root = content.getRootNode() as Document | ShadowRoot;
    const focused = root.activeElement as HTMLElement | null;
    let refocus: HTMLElement | undefined;
    let ref: ChildNode | null = content.firstChild;
    let moved = 0;
    for (const e of list) {
      // Past what is not a shown row (a spare slot, put away): it stays where it is.
      while (ref && !shown.has(ref)) ref = ref.nextSibling;
      if (e.el === ref) {
        ref = ref.nextSibling;
        continue;
      }
      if (focused && !refocus && e.el.contains(focused)) refocus = focused;
      content.insertBefore(e.el, ref);
      moved++;
    }
    refocus?.focus({ preventScroll: true });
    RecycledListStats.settles++;
    RecycledListStats.moved += moved;
    RecycledListStats.lastMs = performance.now() - t0;
  };

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
      ref={(el: HTMLElement) => {
        content = el;
        props.contentRef?.(el);
        bind(ranges, (r) => {
          const run = props.reveal && r.list !== laid ? insertedRun(laid, r.list, props.itemKey) : undefined;
          laid = r.list;
          items = r.list;
          range = r;
          if (run) {
            endReveal();
            animating = true;
            revealEnd = setTimeout(endReveal, RevealMs);
          }
          layout(r.first, r.last, r.rh);
          animating = false;
        });
      }}
      style={{
        height: derive([props.items, rowHeight.asReadonly], (items, rh) => `${items.length * rh}px`),
        // The skeleton's row height (virtual-list.module.css, `.content`).
        '--msp-virtual-row': derive(rowHeight.asReadonly, (rh) => `${rh}px`),
      }}
    />
  );
}
