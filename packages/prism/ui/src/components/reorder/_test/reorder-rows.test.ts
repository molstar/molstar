/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { lostReleases, pointer, stubPointerCapture } from '../../../testing/gestures.js';
import { reorderRows } from '../drag-to-reorder.js';

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

describe('reorderRows', () => {
  it('a row’s drag (gestures attach once): no listener left by lost releases, no marks left, and the next drag moves it', () => {
    const box = document.createElement('div');
    document.body.append(box);
    const moves: [string, number][] = [];
    let list!: HTMLElement;
    render(() => {
      list = document.createElement('div');
      // Rows 20 px tall (jsdom lays nothing out), each with a grip.
      ['a', 'b', 'c'].forEach((id, i) => {
        const row = document.createElement('div');
        row.dataset.id = id;
        row.getBoundingClientRect = () => ({ top: i * 20, height: 20, bottom: i * 20 + 20 }) as DOMRect;
        const grip = document.createElement('span');
        grip.className = 'grip';
        row.append(grip);
        list.append(row);
      });
      reorderRows(list, '.grip', (id, index) => moves.push([id, index]));
      return list;
    }, box);
    // What the lost releases moved (none: only a release moves a row), read before the whole gesture resets it.
    let movedWhileLost = -1;
    const r = lostReleases({
      target: () => list.querySelector('[data-id=a] .grip')!,
      from: { x: 0, y: 5 },
      to: { x: 0, y: 55 },
      reset: () => {
        movedWhileLost = moves.length;
        moves.length = 0;
      },
      worked: () =>
        moves.length === 1 &&
        moves[0][0] === 'a' &&
        moves[0][1] === 2 &&
        ![...list.children].some((c) => c.hasAttribute('data-drop') || c.hasAttribute('data-dragging')),
    });
    expect(r.grown).toBe(0);
    expect(movedWhileLost).toBe(0);
    expect(r.worked).toBe(true);
  });

  it('a cancel never moves a row: only a release does', () => {
    const box = document.createElement('div');
    document.body.append(box);
    const moves: [string, number][] = [];
    let list!: HTMLElement;
    const dispose = render(() => {
      list = document.createElement('div');
      ['a', 'b', 'c'].forEach((id, i) => {
        const row = document.createElement('div');
        row.dataset.id = id;
        row.getBoundingClientRect = () => ({ top: i * 20, height: 20, bottom: i * 20 + 20 }) as DOMRect;
        const grip = document.createElement('span');
        grip.className = 'grip';
        row.append(grip);
        list.append(row);
      });
      reorderRows(list, '.grip', (id, index) => moves.push([id, index]));
      return list;
    }, box);
    const grip = () => list.querySelector('[data-id=a] .grip')!;
    for (const end of ['pointercancel', 'lostpointercapture', 'pointermove'] as const) {
      grip().dispatchEvent(pointer('pointerdown', { x: 0, y: 5 }));
      grip().dispatchEvent(pointer('pointermove', { x: 0, y: 55 }));
      grip().dispatchEvent(pointer(end, { x: 0, y: 55 }, 0));
    }
    // A press after a lost release cancels the one before.
    grip().dispatchEvent(pointer('pointerdown', { x: 0, y: 5 }));
    grip().dispatchEvent(pointer('pointermove', { x: 0, y: 55 }));
    grip().dispatchEvent(pointer('pointerdown', { x: 0, y: 5 }));
    dispose();
    expect(moves).toEqual([]);
    expect([...list.children].some((c) => c.hasAttribute('data-drop') || c.hasAttribute('data-dragging'))).toBe(false);
  });

  it('reads the rows once a drag (again after a scroll or a resize) and moves where the last move pointed', () => {
    const box = document.createElement('div');
    document.body.append(box);
    const moves: [string, number][] = [];
    let reads = 0;
    let shift = 0;
    let list!: HTMLElement;
    const dispose = render(() => {
      list = document.createElement('div');
      ['a', 'b', 'c', 'd'].forEach((id, i) => {
        const row = document.createElement('div');
        row.dataset.id = id;
        row.getBoundingClientRect = () => {
          reads++;
          return { top: i * 20 - shift, height: 20, bottom: i * 20 + 20 - shift } as DOMRect;
        };
        const grip = document.createElement('span');
        grip.className = 'grip';
        row.append(grip);
        list.append(row);
      });
      reorderRows(list, '.grip', (id, index) => moves.push([id, index]));
      return list;
    }, box);
    const grip = list.querySelector('[data-id=a] .grip')!;
    grip.dispatchEvent(pointer('pointerdown', { x: 0, y: 5 }));
    expect(reads).toBe(4);
    for (let y = 5; y < 80; y += 5) grip.dispatchEvent(pointer('pointermove', { x: 0, y }));
    expect(reads).toBe(4);
    const marks = () => [...list.children].map((c) => c.getAttribute('data-drop'));
    expect(marks()).toEqual([null, null, null, 'after']);
    // The list scrolled 40 px up: y 25 is now over row d's top half.
    shift = 40;
    list.dispatchEvent(new Event('scroll'));
    grip.dispatchEvent(pointer('pointermove', { x: 0, y: 25 }));
    expect(reads).toBe(8);
    expect(marks()).toEqual([null, null, null, 'before']);
    grip.dispatchEvent(pointer('pointerup', { x: 0, y: 25 }));
    expect(moves).toEqual([['a', 2]]);
    expect(marks()).toEqual([null, null, null, null]);
    dispose();
  });
});
