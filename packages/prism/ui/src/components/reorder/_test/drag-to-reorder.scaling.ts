// @vitest-environment jsdom

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeAll, it } from 'vitest';
import { checkLinear, scaling } from '../../../_test/scaling.js';
import { render } from '../../../dom/render.js';
import { pointer, stubPointerCapture } from '../../../testing/gestures.js';
import { dragToReorder, reorderRows } from '../drag-to-reorder.js';

// A drag over a list of n rows grows linearly (`test:scaling`, ui/src/_test/scaling.ts): the rows' places are read once
// at the press and each move finds its row by a search over them, marking only the rows whose mark changes (every
// row's rect read and every row's mark written per move: n² for a drag across the list).

beforeAll(stubPointerCapture);

/** n rows 20 px tall (jsdom lays nothing out), each with a grip. */
function rows(list: HTMLElement, n: number) {
  for (let i = 0; i < n; i++) {
    const row = document.createElement('div');
    row.dataset.id = `r${i}`;
    const rect = { top: i * 20, height: 20, bottom: i * 20 + 20 } as DOMRect;
    row.getBoundingClientRect = () => rect;
    const grip = document.createElement('span');
    grip.className = 'grip';
    row.append(grip);
    list.append(row);
  }
}

/** A move a row, top to bottom: the pointer crosses every row of the list. */
const moves = (n: number) => Array.from({ length: n }, (_, i) => pointer('pointermove', { x: 0, y: i * 20 + 15 }));

it('reorderRows: a drag across n rows', () => {
  const s = scaling(
    {
      setup: (n) => {
        const box = document.createElement('div');
        document.body.append(box);
        let list!: HTMLElement;
        const dispose = render(() => {
          list = document.createElement('div');
          rows(list, n);
          reorderRows(list, '.grip', () => {});
          return list;
        }, box);
        const grip = list.querySelector('[data-id=r0] .grip')!;
        return { box, dispose, grip, events: moves(n), down: pointer('pointerdown', { x: 0, y: 5 }) };
      },
      run: ({ grip, events, down }) => {
        grip.dispatchEvent(down);
        for (const e of events) grip.dispatchEvent(e);
        grip.dispatchEvent(pointer('pointerup', { x: 0, y: 5 }));
      },
      teardown: ({ box, dispose }) => {
        dispose();
        box.remove();
      },
    },
    1000,
  );
  checkLinear('reorderRows: a drag across n rows', s);
});

it('dragToReorder: a drag across n rows', () => {
  const s = scaling(
    {
      setup: (n) => {
        const list = document.createElement('div');
        document.body.append(list);
        rows(list, n);
        const grip = list.querySelector('[data-id=r0] .grip')!;
        return { list, grip, events: moves(n) };
      },
      run: ({ list, grip, events }) => {
        const down = pointer('pointerdown', { x: 0, y: 5 });
        Object.defineProperty(down, 'currentTarget', { value: grip });
        dragToReorder(down, list, 'r0', () => {});
        for (const e of events) grip.dispatchEvent(e);
        grip.dispatchEvent(pointer('pointerup', { x: 0, y: 5 }));
      },
      teardown: ({ list }) => list.remove(),
    },
    1000,
  );
  checkLinear('dragToReorder: a drag across n rows', s);
});
