/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Drags a row of a list by its handle (call it on the handle's pointerdown): the list's rows are its children with
 * `data-id`; the one being moved gets `data-dragging`, the row it would land before `data-drop="before"` (the last one
 * `"after"`), for the list's own styles; release moves it there. The pointer is captured on the handle, so leaving the
 * list keeps the drag. Snapshots, measurement points.
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
  const handle = e.currentTarget as HTMLElement;
  handle.setPointerCapture?.(e.pointerId);
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
  const onUp = () => {
    handle.removeEventListener('pointermove', onMove);
    self?.removeAttribute('data-dragging');
    const from = items.indexOf(self as HTMLElement);
    const to = target;
    mark(-1);
    // Dropping above or below itself leaves it where it is.
    if (to >= 0 && to !== from && to !== from + 1) moveTo(to > from ? to - 1 : to);
  };
  handle.addEventListener('pointermove', onMove);
  handle.addEventListener('pointerup', onUp, { once: true });
  handle.addEventListener('pointercancel', onUp, { once: true });
}
