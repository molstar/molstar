/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { flush, onCleanup } from './scope.js';
import { on } from './scoped.js';

// The keyboard of an ARIA tree (accessibility.md, WAI-ARIA's tree pattern), over the tree's rows as a list (a
// virtualized tree has only the rows on screen in the DOM): Down and Up move, Right opens (open: enters), Left closes
// (closed: goes to the parent), Home and End jump, Enter and Space do what a click on the row does. A move scrolls the
// row into view (`scrollTo`, which renders it) and focuses it. Focus roves: the active row has `tabindex="0"`, the
// others `-1`; while the active row is not rendered (scrolled away), the tree itself is the tab stop and takes its
// keys. A row's open state is its `aria-expanded`; its open/close button is its `[data-tree-toggle]`, its click target
// `[data-tree-row]` (or `[data-tree-action]`).

const Item = '[role="treeitem"]';

export interface TreeNavigationOptions<R> {
  /** The tree's shown rows, flattened (expanded-aware), as they are now. */
  rows(): readonly R[];
  key(row: R): string;
  /** 1-based (`aria-level`). */
  level(row: R): number;
  /** Brings the row at `index` into view, rendered when it returns (`VirtualListHandle.scrollToIndex`). */
  scrollTo(index: number): void;
  /** The row Tab enters at when none was active yet (the picked one). */
  selected?(): string | undefined;
}

/**
 * A tree's keyboard and roving focus, by the rows' model rather than the DOM. The tree's element binds `tree`, and each
 * row's element `row(key)` (a treeitem, or an element holding one: an extension's row).
 */
export class TreeNavigation<R> {
  readonly #o: TreeNavigationOptions<R>;
  readonly #elements = new Map<string, HTMLElement>();
  #tree: HTMLElement | undefined;
  #active: string | undefined;
  /** The active row's index at its last move: where focus goes when its row leaves the list. */
  #index = 0;
  #parking = false;

  constructor(options: TreeNavigationOptions<R>) {
    this.#o = options;
  }

  /** The active row's key: the one with `tabindex="0"`. */
  get active(): string | undefined {
    return this.#active;
  }

  /** Ref for the tree element (`role="tree"`). */
  readonly tree = (el: HTMLElement): void => {
    this.#tree = el;
    // Scoped (rendering.md §2.4): detached while the tree's scope is paused, gone with it.
    on(el, 'keydown', this.#onKey);
    on(el, 'focusin', this.#onFocusIn);
    this.#sync();
    onCleanup(() => {
      this.#tree = undefined;
    });
  };

  /** Ref for a row's element, while it is rendered. */
  row(key: string): (el: HTMLElement) => void {
    return (el) => {
      el.dataset.treeKey = key;
      this.#elements.set(key, el);
      const item = target(el);
      if (item) item.tabIndex = key === this.#active ? 0 : -1;
      this.#sync();
      onCleanup(() => this.#leave(key, el));
    };
  }

  /**
   * For a recycled row (`RecycledList`'s `slot`): its element, which shows one row, then another. The setter moves it to
   * the row `key` (undefined: it shows none); a focus inside it goes as it goes when a row leaves the DOM (to the tree).
   */
  slot(el: HTMLElement): (key: string | undefined) => void {
    let shown: string | undefined;
    onCleanup(() => {
      if (shown !== undefined) this.#leave(shown, el);
    });
    return (key) => {
      if (key === shown) return;
      if (shown !== undefined) this.#leave(shown, el);
      shown = key;
      const item = target(el);
      if (key === undefined) {
        delete el.dataset.treeKey;
        if (item) item.tabIndex = -1;
        this.#sync();
        return;
      }
      el.dataset.treeKey = key;
      this.#elements.set(key, el);
      if (item) item.tabIndex = key === this.#active ? 0 : -1;
      this.#sync();
    };
  }

  /** Focuses the row at `index` (clamped), scrolling it into view. */
  focusIndex(index: number): void {
    const rows = this.#o.rows();
    if (!rows.length) return;
    const i = Math.max(0, Math.min(rows.length - 1, index));
    this.#index = i;
    this.#setActive(this.#o.key(rows[i]));
    this.#o.scrollTo(i);
    this.#focusActive();
  }

  #indexOf(key: string | undefined): number {
    if (key === undefined) return -1;
    const rows = this.#o.rows();
    for (let i = 0; i < rows.length; i++) if (this.#o.key(rows[i]) === key) return i;
    return -1;
  }

  #setActive(key: string | undefined): void {
    if (key === this.#active) return;
    const old = this.#active !== undefined ? this.#elements.get(this.#active) : undefined;
    const oldItem = old && target(old);
    if (oldItem) oldItem.tabIndex = -1;
    this.#active = key;
    const now = key !== undefined ? this.#elements.get(key) : undefined;
    const item = now && target(now);
    if (item) item.tabIndex = 0;
    this.#sync();
  }

  #focusActive(): void {
    const el = this.#active !== undefined ? this.#elements.get(this.#active) : undefined;
    target(el)?.focus({ preventScroll: true });
  }

  /** The tree is the tab stop while the active row is not rendered. */
  #sync(): void {
    if (!this.#tree) return;
    const rendered = this.#active !== undefined && this.#elements.has(this.#active);
    this.#tree.tabIndex = rendered ? -1 : 0;
  }

  /** A row leaves the DOM (scrolled away, closed, removed): if it had the focus, the tree keeps it. */
  #leave(key: string, el: HTMLElement): void {
    const current = this.#elements.get(key);
    if (current === el) {
      this.#elements.delete(key);
      this.#sync();
    }
    const tree = this.#tree;
    if (!tree) return;
    const root = el.getRootNode() as Document | ShadowRoot;
    const focused = root.activeElement;
    if (!focused || !el.contains(focused)) return;
    this.#parking = true;
    tree.focus({ preventScroll: true });
    this.#parking = false;
    // Rendered anew (its record replaced, the window moved back), or gone (Show all, a parent closed): the focus goes
    // back to it, or to the row now in its place, once the list has settled (the new element is in the DOM).
    queueMicrotask(() => {
      if ((tree.getRootNode() as Document | ShadowRoot).activeElement !== tree) return;
      const at = this.#indexOf(key);
      if (at >= 0) {
        target(this.#elements.get(key))?.focus({ preventScroll: true });
        return;
      }
      const rows = this.#o.rows();
      if (!rows.length) return;
      const i = Math.min(this.#index, rows.length - 1);
      const next = this.#o.key(rows[i]);
      this.#index = i;
      this.#setActive(next);
      this.#focusActive();
    });
  }

  readonly #onFocusIn = (e: FocusEvent): void => {
    const tree = e.currentTarget as HTMLElement;
    if (e.target === tree) {
      // Tab into the tree with no active row on screen: the active row, else the picked one, else the first.
      if (this.#parking || (e.relatedTarget instanceof Node && tree.contains(e.relatedTarget))) return;
      let i = this.#indexOf(this.#active);
      if (i < 0) i = this.#indexOf(this.#o.selected?.());
      this.focusIndex(Math.max(i, 0));
      return;
    }
    // Focus on a row (a click, a button in it): it is the active one.
    const owner = (e.target as HTMLElement).closest<HTMLElement>('[data-tree-key]');
    const key = owner?.dataset.treeKey;
    if (key === undefined || !tree.contains(owner)) return;
    const i = this.#indexOf(key);
    if (i >= 0) this.#index = i;
    this.#setActive(key);
  };

  readonly #onKey = (e: KeyboardEvent): void => {
    // The rows' DOM as their state is now (a row's `aria-expanded` is read below; `watch` renders in a microtask).
    flush();
    const tree = e.currentTarget as HTMLElement;
    const rows = this.#o.rows();
    if (e.target === tree) {
      // The tree itself has the focus (the active row scrolled away): the keys move from the active row.
      const at = this.#indexOf(this.#active);
      const from = at >= 0 ? at : this.#index;
      const to: Record<string, number> = {
        ArrowDown: at >= 0 ? from + 1 : from,
        ArrowUp: at >= 0 ? from - 1 : from,
        Home: 0,
        End: rows.length - 1,
      };
      if (e.key in to) {
        e.preventDefault();
        this.focusIndex(to[e.key]);
      }
      return;
    }
    const item = (e.target as HTMLElement).closest<HTMLElement>(Item);
    // Keys on a button inside a row (the eye, the menu) are the button's.
    if (!item || e.target !== item) return;
    const key = item.closest<HTMLElement>('[data-tree-key]')?.dataset.treeKey;
    const at = this.#indexOf(key);
    if (at < 0) return;
    const level = this.#o.level(rows[at]);
    const expanded = item.getAttribute('aria-expanded');
    const toggle = item.querySelector<HTMLElement>('[data-tree-toggle]');
    let handled = true;
    switch (e.key) {
      case 'ArrowDown':
        if (at + 1 < rows.length) this.focusIndex(at + 1);
        break;
      case 'ArrowUp':
        if (at > 0) this.focusIndex(at - 1);
        break;
      case 'Home':
        this.focusIndex(0);
        break;
      case 'End':
        this.focusIndex(rows.length - 1);
        break;
      case 'ArrowRight':
        if (expanded === 'false') toggle?.click();
        else if (expanded === 'true' && at + 1 < rows.length && this.#o.level(rows[at + 1]) > level)
          this.focusIndex(at + 1);
        break;
      case 'ArrowLeft':
        if (expanded === 'true') toggle?.click();
        else
          for (let i = at - 1; i >= 0; i--) {
            if (this.#o.level(rows[i]) < level) {
              this.focusIndex(i);
              break;
            }
          }
        break;
      case 'Enter':
      case ' ':
        item.querySelector<HTMLElement>('[data-tree-row], [data-tree-action]')?.click();
        break;
      default:
        handled = false;
    }
    if (handled) e.preventDefault();
  };
}

/** The focusable treeitem of a row's element: the element itself, or the one it holds. */
function target(el: HTMLElement | undefined): HTMLElement | undefined {
  if (!el) return undefined;
  return el.matches(Item) ? el : (el.querySelector<HTMLElement>(Item) ?? undefined);
}
