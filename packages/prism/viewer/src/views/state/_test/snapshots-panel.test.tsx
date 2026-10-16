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
import { afterEach, describe, expect, it } from 'vitest';
import type { SnapshotEntry, SnapshotsModel } from '../../../model/snapshots.js';
import { SnapshotsPanel } from '../snapshots-panel.js';

// The snapshots list: a row applies on Enter, and moves with Alt+Up and Alt+Down (the order is also the handle's drag).

afterEach(() => {
  document.body.innerHTML = '';
});

function mount() {
  const entries = observableState<readonly SnapshotEntry[]>(
    ['a', 'b', 'c'].map((id) => ({ id, name: `Snapshot ${id}` }) as SnapshotEntry),
  );
  const moves: [string, number][] = [];
  const model = {
    entries: entries.asReadonly,
    current: observableState<string | undefined>(undefined).asReadonly,
    saveOptions: undefined,
    entryOptions: () => undefined,
    apply: () => {},
    move: (id: string, by: number) => moves.push([id, by]),
  } as unknown as SnapshotsModel;
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <SnapshotsPanel model={model} />, container);
  const rows = () => [...container.querySelectorAll<HTMLElement>('[role="list"] [role="button"][aria-keyshortcuts]')];
  const key = (el: HTMLElement, k: string, altKey: boolean) =>
    el.dispatchEvent(new KeyboardEvent('keydown', { key: k, altKey, bubbles: true, cancelable: true }));
  return { rows, key, moves };
}

describe('the snapshots list', () => {
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
});
