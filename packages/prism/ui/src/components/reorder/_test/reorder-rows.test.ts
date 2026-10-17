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
import { lostReleases, stubPointerCapture } from '../../../testing/gestures.js';
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
    const r = lostReleases({
      target: () => list.querySelector('[data-id=a] .grip')!,
      from: { x: 0, y: 5 },
      to: { x: 0, y: 55 },
      reset: () => {
        moves.length = 0;
      },
      worked: () =>
        moves.length === 1 &&
        moves[0][0] === 'a' &&
        moves[0][1] === 2 &&
        ![...list.children].some((c) => c.hasAttribute('data-drop') || c.hasAttribute('data-dragging')),
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
