/** @jest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { dragToReorder } from '../drag-to-reorder.js';

/** A list of rows 20px tall (jsdom has no layout), each with a handle. */
function list(ids: string[]) {
  const el = document.createElement('div');
  ids.forEach((id, i) => {
    const row = document.createElement('div');
    row.dataset.id = id;
    row.getBoundingClientRect = () => ({ top: i * 20, height: 20, bottom: i * 20 + 20 }) as DOMRect;
    row.appendChild(document.createElement('span'));
    el.appendChild(row);
  });
  document.body.appendChild(el);
  return el;
}

function pointer(type: string, clientY: number) {
  const e = new MouseEvent(type, { bubbles: true, cancelable: true, clientY, button: 0 });
  Object.defineProperty(e, 'pointerId', { value: 1 });
  return e as unknown as PointerEvent;
}

describe('dragToReorder', () => {
  it('marks where the row would land and moves it there on release', () => {
    const el = list(['a', 'b', 'c']);
    const handle = el.children[0].firstElementChild as HTMLElement;
    const moves: number[] = [];
    const down = pointer('pointerdown', 5);
    Object.defineProperty(down, 'currentTarget', { value: handle });
    dragToReorder(down, el, 'a', (i) => moves.push(i));
    expect((el.children[0] as HTMLElement).dataset.dragging).toBe('');

    handle.dispatchEvent(pointer('pointermove', 55));
    expect((el.children[2] as HTMLElement).dataset.drop).toBe('after');

    handle.dispatchEvent(pointer('pointerup', 55));
    expect(moves).toEqual([2]);
    expect([...el.children].some((c) => (c as HTMLElement).dataset.drop)).toBe(false);
    expect((el.children[0] as HTMLElement).dataset.dragging).toBeUndefined();
  });

  it('dropping next to itself does not move it', () => {
    const el = list(['a', 'b', 'c']);
    const handle = el.children[1].firstElementChild as HTMLElement;
    const moves: number[] = [];
    const down = pointer('pointerdown', 25);
    Object.defineProperty(down, 'currentTarget', { value: handle });
    dragToReorder(down, el, 'b', (i) => moves.push(i));
    handle.dispatchEvent(pointer('pointermove', 29));
    handle.dispatchEvent(pointer('pointerup', 29));
    expect(moves).toEqual([]);
  });

  it('a release anywhere (or the window losing focus) ends the drag and clears the marks', () => {
    const el = list(['a', 'b', 'c']);
    const handle = el.children[0].firstElementChild as HTMLElement;
    const down = pointer('pointerdown', 5);
    Object.defineProperty(down, 'currentTarget', { value: handle });
    dragToReorder(down, el, 'a', () => {});
    document.body.dispatchEvent(pointer('pointermove', 35));
    expect((el.children[2] as HTMLElement).dataset.drop).toBe('before');
    window.dispatchEvent(new Event('blur'));
    expect([...el.children].some((c) => (c as HTMLElement).dataset.drop || (c as HTMLElement).dataset.dragging)).toBe(
      false,
    );
  });
});
