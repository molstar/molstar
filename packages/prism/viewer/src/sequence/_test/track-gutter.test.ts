/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EntryKind } from '@molstar/seq-viewer';
import { describe, expect, it, vi } from 'vitest';
import { LettersRef, PlddtTrack, SecondaryTrack } from '../document.js';
import type { HostActions, StructureColoring } from '../model.js';
import { gutterRows } from '../track-gutter.js';
import { headlessViewers } from './fixtures.js';

const viewerWith = (now: () => StructureColoring | undefined, colorStructure = vi.fn()) => {
  const actions: HostActions = { colorStructure, coloring: now };
  return { viewer: headlessViewers({ actions }).viewer, colorStructure };
};

describe('the track icons', () => {
  it("every track's icon is the one paintbrush, of its first menu action, and tracks without a coloring have none", () => {
    const { viewer } = viewerWith(() => undefined);
    expect(viewer.trackColoring(PlddtTrack)).toMatchObject({ icon: 'paintbrush', coloring: 'plddt' });
    expect(viewer.trackColoring(SecondaryTrack)).toMatchObject({ icon: 'paintbrush', coloring: 'secondary' });
    expect(viewer.trackColoring(LettersRef)).toMatchObject({ icon: 'paintbrush', coloring: 'default' });
    expect(viewer.trackColoring('ruler')).toBeUndefined();
    expect(headlessViewers().viewer.trackColoring(PlddtTrack)).toBeUndefined();
  });

  it('is active where the 3D view has that coloring now (read each time), and a click applies it unless it is', () => {
    let now: StructureColoring | undefined = 'default';
    const { viewer, colorStructure } = viewerWith(() => now);
    const plddt = viewer.trackColoring(PlddtTrack);
    expect([viewer.trackColoring(LettersRef)?.active(), plddt?.active()]).toEqual([true, false]);
    plddt?.run();
    expect(colorStructure).toHaveBeenCalledExactlyOnceWith('plddt');
    now = 'plddt';
    expect([viewer.trackColoring(LettersRef)?.active(), plddt?.active()]).toEqual([false, true]);
    plddt?.run();
    expect(colorStructure).toHaveBeenCalledTimes(1);
    now = undefined;
    expect(viewer.trackColoring(SecondaryTrack)?.active()).toBe(false);
  });
});

describe('gutterRows', () => {
  it('lists the track rows that have an icon, pinned at their y and scrolled after the pinned part', () => {
    const layout = {
      entryCount: 4,
      pinnedCount: 2,
      entryKind: Uint8Array.of(EntryKind.track, EntryKind.track, EntryKind.track, EntryKind.group),
      entryTrack: Int32Array.of(0, 1, 2, -1),
      entryTop: Float64Array.of(0, 20, 5, 40),
      entryHeight: Float32Array.of(18, 22, 30, 16),
      scrolledTop: 50,
    };
    const refs = ['ruler', 'letters', 'ss'];
    const rows = gutterRows(
      layout,
      (t) => refs[t],
      (ref) => ref !== 'ruler',
    );
    expect(rows).toEqual([
      { track: 'letters', top: 20, height: 22 },
      { track: 'ss', top: 55, height: 30 },
    ]);
  });
});
