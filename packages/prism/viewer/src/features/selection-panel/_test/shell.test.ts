/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import type { BuiltInShellContext } from '../../../spec.js';
import { selectionPanelShell } from '../shell.js';

// The Selection tab's shell: showing a component just made in the Structure panel waits for its row (polling).

const made: { show?: (ref: string) => void } = {};
vi.mock('../panels.js', () => ({
  selectionTab: () => ({ id: 'selection' }),
  selectionSections: (_api: unknown, _selection: unknown, _layers: unknown, show: (ref: string) => void) => {
    made.show = show;
    return [];
  },
}));
vi.mock('../../../ui/interaction.js', () => ({ selectionSummaryModel: () => ({}) }));

afterEach(() => vi.useRealTimers());

function shellWith(reveal: (ref: string) => boolean) {
  const ctx = {
    shell: { registerTabs: () => () => {}, registerSections: () => () => {}, show: vi.fn() },
    core: { sources: { interaction: { selection: {} }, layers: {} }, views: { scene: { reveal } } },
    person: {},
    own: (_kind: string, parts: unknown[]) => parts,
  } as unknown as BuiltInShellContext;
  return selectionPanelShell(ctx);
}

describe('selection panel shell', () => {
  it('reveals one component at a time, the latest, and stops with the shell', () => {
    vi.useFakeTimers();
    const reveal = vi.fn((_ref: string) => false);
    const dispose = shellWith(reveal);
    made.show!('a');
    made.show!('b');
    expect(reveal.mock.calls.map((c) => c[0])).toEqual(['a', 'b']);
    vi.advanceTimersByTime(120);
    // Only the second keeps trying: the first was replaced.
    expect(reveal.mock.calls.map((c) => c[0])).toEqual(['a', 'b', 'b', 'b']);
    dispose();
    vi.advanceTimersByTime(2000);
    expect(reveal).toHaveBeenCalledTimes(4);
  });

  it('stops when the row is there, and gives up after a second', () => {
    vi.useFakeTimers();
    let there = false;
    const found = vi.fn((_ref: string) => there);
    const dispose = shellWith(found);
    made.show!('a');
    vi.advanceTimersByTime(100);
    there = true;
    vi.advanceTimersByTime(50);
    const calls = found.mock.calls.length;
    vi.advanceTimersByTime(2000);
    expect(found).toHaveBeenCalledTimes(calls);

    there = false;
    found.mockClear();
    made.show!('b');
    vi.advanceTimersByTime(5000);
    expect(found).toHaveBeenCalledTimes(21);
    dispose();
  });
});
