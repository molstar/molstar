/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { drag, on } from '../../dom/scoped.js';

/**
 * A list's rows for one drag: their midpoints read once (at the press, before any mark is written; again after a scroll
 * or a resize), so a move is a binary search over them, and the marks written only on the rows whose mark changes. Per
 * move, every row's rect read and every row's mark written was a drag across n rows n², with a layout read after
 * writes at each move (rendering.md §2.4, a move reads no layout; prism AGENTS.md, snappy).
 */
class Rows {
  readonly items: HTMLElement[];
  #mids: Float64Array | undefined;
  /** Where the row would land (an index of `items`, or its length: after the last); -1: nowhere yet. */
  target = -1;

  constructor(list: HTMLElement) {
    this.items = [...list.querySelectorAll<HTMLElement>(':scope > [data-id]')];
    this.#measure();
  }

  /** The rows moved on screen (a scroll, a resize): read again at the next move. */
  stale() {
    this.#mids = undefined;
  }

  #measure(): Float64Array {
    const { items } = this;
    const mids = new Float64Array(items.length);
    for (let i = 0; i < items.length; i++) {
      const r = items[i].getBoundingClientRect();
      mids[i] = r.top + r.height / 2;
    }
    this.#mids = mids;
    return mids;
  }

  /** The first row whose middle is below `y` (the rows go down the list); after the last: their count. */
  indexAt(y: number): number {
    const mids = this.#mids ?? this.#measure();
    let lo = 0;
    let hi = mids.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (y < mids[mid]) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  }

  /** The row the mark of `index` goes on: the row it lands before, or the last one (landing after it). */
  #marked(index: number): HTMLElement | undefined {
    const { items } = this;
    if (index < 0 || items.length === 0) return undefined;
    return index < items.length ? items[index] : items[items.length - 1];
  }

  /** Marks the row it would land before (`data-drop="before"`; the last one `"after"`); -1 clears the mark. */
  mark(index: number) {
    if (index === this.target) return;
    this.#marked(this.target)?.removeAttribute('data-drop');
    this.target = index;
    const el = this.#marked(index);
    if (el) el.setAttribute('data-drop', index < this.items.length ? 'before' : 'after');
  }

  /** Where `self` goes on a drop at `target` (an index of the list without it); undefined: it stays. */
  moveOf(self: HTMLElement | undefined): number | undefined {
    const from = self ? this.items.indexOf(self) : -1;
    const to = this.target;
    // Dropped above or below itself, it stays where it is.
    if (to < 0 || to === from || to === from + 1) return undefined;
    return to > from ? to - 1 : to;
  }
}

/**
 * Drags a row of a list by its handle (call it on the handle's pointerdown): the list's rows are its children with
 * `data-id`; the one being moved gets `data-dragging`, the row it would land before `data-drop="before"` (the last one
 * `"after"`), for the list's own styles; release moves it there. The document follows the pointer (not the handle), so
 * a release anywhere, a cancelled pointer or the window losing focus always ends the drag and clears the marks; only
 * the release moves the row.
 * Snapshots, measurement points.
 *
 * It attaches on the press: a component uses `reorderRows` (attached once, with its scope) instead.
 */
export function dragToReorder(
  e: PointerEvent,
  list: HTMLElement | undefined,
  id: string,
  moveTo: (index: number) => void,
) {
  if (!list || e.button !== 0) return;
  e.preventDefault();
  e.stopPropagation();
  const rows = new Rows(list);
  const self = rows.items.find((el) => el.dataset.id === id);
  self?.setAttribute('data-dragging', '');
  const onMove = (ev: PointerEvent) => rows.mark(rows.indexAt(ev.clientY));
  const onMoved = () => rows.stale();
  const doc = list.ownerDocument;
  const view = doc.defaultView;
  let done = false;
  const end = (released: boolean) => {
    if (done) return;
    done = true;
    doc.removeEventListener('pointermove', onMove, true);
    doc.removeEventListener('pointerup', onUp, true);
    doc.removeEventListener('pointercancel', onCancel, true);
    doc.removeEventListener('scroll', onMoved, true);
    view?.removeEventListener('blur', onCancel);
    view?.removeEventListener('resize', onMoved);
    self?.removeAttribute('data-dragging');
    const to = rows.moveOf(self);
    rows.mark(-1);
    // Only a release moves it (a cancel leaves it).
    if (released && to !== undefined) moveTo(to);
  };
  const onUp = () => end(true);
  const onCancel = () => end(false);
  doc.addEventListener('pointermove', onMove, true);
  doc.addEventListener('pointerup', onUp, true);
  doc.addEventListener('pointercancel', onCancel, true);
  // Any scroll (captured: scrolls do not bubble) or resize moves the rows on screen.
  doc.addEventListener('scroll', onMoved, { capture: true, passive: true });
  view?.addEventListener('blur', onCancel);
  view?.addEventListener('resize', onMoved, { passive: true });
}

interface Reorder {
  readonly rows: Rows;
  readonly self: HTMLElement;
}

/**
 * `dragToReorder`'s gesture attached once (prism AGENTS.md, gestures; `drag`, dom/scoped.ts), in a component: a press
 * on an element matching `handle` in one of `list`'s rows (its children with `data-id`) drags that row, with the same
 * marks; a release moves it (`moveTo` with its id and the index it goes to). Its listeners are on the list, with the
 * current scope; the pointer is captured, and every way the drag ends (a release lost, the capture lost, the scope
 * paused or ended) clears the marks; only a release moves the row.
 */
export function reorderRows(list: HTMLElement, handle: string, moveTo: (id: string, index: number) => void): void {
  let current: Reorder | undefined;
  drag<Reorder>(list, {
    start(e) {
      if (e.button !== 0) return undefined;
      const self = (e.target as Element).closest?.(handle)?.closest<HTMLElement>('[data-id]');
      if (!self || self.parentElement !== list) return undefined;
      e.preventDefault();
      e.stopPropagation();
      // The rows read before the first mark is written.
      const rows = new Rows(list);
      self.setAttribute('data-dragging', '');
      current = { rows, self };
      return current;
    },
    move(e, r) {
      r.rows.mark(r.rows.indexAt(e.clientY));
    },
    end(r, how) {
      if (current === r) current = undefined;
      r.self.removeAttribute('data-dragging');
      const to = r.rows.moveOf(r.self);
      r.rows.mark(-1);
      // Released: a cancel never moves it (rendering.md §2.4).
      if (how === 'up' && to !== undefined) moveTo(r.self.dataset.id!, to);
    },
  });
  // Attached once with the scope, as the gesture is: any scroll (captured: scrolls do not bubble) or resize moves the
  // rows on screen; the drag on reads them again.
  const moved = () => current?.rows.stale();
  const doc = list.ownerDocument;
  on(doc, 'scroll', moved, { capture: true, passive: true });
  if (doc.defaultView) on(doc.defaultView, 'resize', moved, { passive: true });
}
