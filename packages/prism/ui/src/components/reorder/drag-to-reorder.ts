/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { drag } from '../../dom/scoped.js';

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
  const items = [...list.querySelectorAll<HTMLElement>(':scope > [data-id]')];
  const self = items.find((el) => el.dataset.id === id);
  self?.setAttribute('data-dragging', '');
  let target = -1;
  const mark = (index: number) => {
    target = index;
    for (let i = 0; i < items.length; i++) {
      const el = items[i];
      if (i === index) el.setAttribute('data-drop', 'before');
      else if (index === items.length && i === items.length - 1) el.setAttribute('data-drop', 'after');
      else el.removeAttribute('data-drop');
    }
  };
  const onMove = (ev: PointerEvent) => {
    let index = items.length;
    for (let i = 0; i < items.length; i++) {
      const r = items[i].getBoundingClientRect();
      if (ev.clientY < r.top + r.height / 2) {
        index = i;
        break;
      }
    }
    mark(index);
  };
  const doc = list.ownerDocument;
  const view = doc.defaultView;
  let done = false;
  const end = (released: boolean) => {
    if (done) return;
    done = true;
    doc.removeEventListener('pointermove', onMove, true);
    doc.removeEventListener('pointerup', onUp, true);
    doc.removeEventListener('pointercancel', onCancel, true);
    view?.removeEventListener('blur', onCancel);
    self?.removeAttribute('data-dragging');
    const from = items.indexOf(self as HTMLElement);
    const to = target;
    mark(-1);
    // Only a release moves it (a cancel leaves it); dropped above or below itself, it stays where it is.
    if (released && to >= 0 && to !== from && to !== from + 1) moveTo(to > from ? to - 1 : to);
  };
  const onUp = () => end(true);
  const onCancel = () => end(false);
  doc.addEventListener('pointermove', onMove, true);
  doc.addEventListener('pointerup', onUp, true);
  doc.addEventListener('pointercancel', onCancel, true);
  view?.addEventListener('blur', onCancel);
}

interface Reorder {
  readonly items: HTMLElement[];
  readonly self: HTMLElement;
  target: number;
}

/**
 * `dragToReorder`'s gesture attached once (prism AGENTS.md, gestures; `drag`, dom/scoped.ts), in a component: a press
 * on an element matching `handle` in one of `list`'s rows (its children with `data-id`) drags that row, with the same
 * marks; a release moves it (`moveTo` with its id and the index it goes to). Its listeners are on the list, with the
 * current scope; the pointer is captured, and every way the drag ends (a release lost, the capture lost, the scope
 * paused or ended) clears the marks; only a release moves the row.
 */
export function reorderRows(list: HTMLElement, handle: string, moveTo: (id: string, index: number) => void): void {
  const mark = (r: Reorder, index: number) => {
    r.target = index;
    const { items } = r;
    for (let i = 0; i < items.length; i++) {
      const el = items[i];
      if (i === index) el.setAttribute('data-drop', 'before');
      else if (index === items.length && i === items.length - 1) el.setAttribute('data-drop', 'after');
      else el.removeAttribute('data-drop');
    }
  };
  drag<Reorder>(list, {
    start(e) {
      if (e.button !== 0) return undefined;
      const self = (e.target as Element).closest?.(handle)?.closest<HTMLElement>('[data-id]');
      if (!self || self.parentElement !== list) return undefined;
      e.preventDefault();
      e.stopPropagation();
      self.setAttribute('data-dragging', '');
      return { items: [...list.querySelectorAll<HTMLElement>(':scope > [data-id]')], self, target: -1 };
    },
    move(e, r) {
      let index = r.items.length;
      for (let i = 0; i < r.items.length; i++) {
        const b = r.items[i].getBoundingClientRect();
        if (e.clientY < b.top + b.height / 2) {
          index = i;
          break;
        }
      }
      mark(r, index);
    },
    end(r, how) {
      r.self.removeAttribute('data-dragging');
      const from = r.items.indexOf(r.self);
      const to = r.target;
      mark(r, -1);
      // Released: a cancel never moves it (rendering.md §2.4). Dropped above or below itself, it stays where it is.
      if (how === 'up' && to >= 0 && to !== from && to !== from + 1)
        moveTo(r.self.dataset.id!, to > from ? to - 1 : to);
    },
  });
}
