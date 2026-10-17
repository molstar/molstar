/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { lostReleases, stubPointerCapture } from '@molstar/prism-ui/testing/gestures';
import { afterEach, beforeAll, describe, expect, it } from 'vitest';
import type { SnapshotEntry, SnapshotsModel } from '../../../model/snapshots.js';
import { SnapshotsPanel } from '../snapshots-panel.js';

// The snapshots list: a row applies on Enter, and moves with Alt+Up and Alt+Down (the order is also the handle's drag).

beforeAll(stubPointerCapture);
afterEach(() => {
  document.body.innerHTML = '';
});

function mount() {
  const entries = observableState<readonly SnapshotEntry[]>(
    ['a', 'b', 'c'].map((id) => ({ id, name: `Snapshot ${id}` }) as SnapshotEntry),
  );
  const current = observableState<string | undefined>(undefined);
  const moves: [string, number][] = [];
  const moved: [string, number][] = [];
  const model = {
    entries: entries.asReadonly,
    current: current.asReadonly,
    saveOptions: undefined,
    entryOptions: () => undefined,
    apply: () => {},
    move: (id: string, by: number) => moves.push([id, by]),
    moveTo: (id: string, index: number) => moved.push([id, index]),
  } as unknown as SnapshotsModel;
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <SnapshotsPanel model={model} />, container);
  const rows = () => [...container.querySelectorAll<HTMLElement>('[role="list"] [role="button"][aria-keyshortcuts]')];
  const key = (el: HTMLElement, k: string, altKey: boolean) =>
    el.dispatchEvent(new KeyboardEvent('keydown', { key: k, altKey, bubbles: true, cancelable: true }));
  return { container, rows, key, moves, moved, current };
}

describe('the snapshots list', () => {
  it('shows the current snapshot in the inspector, following the current one when entries do not change', async () => {
    const t = mount();
    const inspector = () => t.container.querySelector<HTMLElement>('[aria-label="Snapshot"]')!;
    const name = () => inspector().querySelector<HTMLInputElement>('input')?.value;
    expect(name()).toBeUndefined();
    t.current.set('a');
    await Promise.resolve();
    await new Promise((r) => setTimeout(r, 0));
    expect(name()).toBe('Snapshot a');
    t.current.set('c');
    await new Promise((r) => setTimeout(r, 0));
    expect(name()).toBe('Snapshot c');
    t.current.set(undefined);
    await new Promise((r) => setTimeout(r, 0));
    expect(name()).toBeUndefined();
  });

  it('moves a row with Alt+Up and Alt+Down, and leaves the plain arrows alone', () => {
    const t = mount();
    expect(t.rows()).toHaveLength(3);
    t.key(t.rows()[1], 'ArrowUp', true);
    t.key(t.rows()[1], 'ArrowDown', true);
    t.key(t.rows()[1], 'ArrowDown', false);
    expect(t.moves).toEqual([
      ['b', -1],
      ['b', 1],
    ]);
  });

  it('drags a row by its grip: no listener left by lost releases, and the next drag moves it', () => {
    const t = mount();
    // jsdom lays nothing out: rows 20 px high, one under another.
    const items = [...t.container.querySelectorAll<HTMLElement>('[role="list"] > [data-id]')];
    items.forEach((el, i) => (el.getBoundingClientRect = () => ({ top: 20 * i, height: 20 }) as DOMRect));
    const grip = () => items[0].querySelector<HTMLElement>('[data-tip^="Drag to reorder"]')!;
    const r = lostReleases({
      target: grip,
      from: { x: 5, y: 5 },
      to: { x: 5, y: 55 },
      reset: () => (t.moved.length = 0),
      worked: () =>
        t.moved.length === 1 &&
        t.moved[0][0] === 'a' &&
        t.moved[0][1] === 2 &&
        !items.some((el) => el.hasAttribute('data-dragging') || el.hasAttribute('data-drop')),
    });
    expect(r.grown).toBe(0);
    expect(r.worked).toBe(true);
  });
});
