/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { flushFrame } from '@molstar/prism-ui/reactive/frame-state';
import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';

// The hover label: Mol*'s label manager reuses its labels array, so each highlight must still reach the lines.

describe('the hover lines', () => {
  it('follow every highlight, one after another', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seen: string[] = [];
    const sub = v.sources.pointer.hoverLines.subscribe((l) => seen.push(l.join(' / ')));
    for (const target of ['A:10', 'A:20', 'ligand']) {
      const r = await v.api.highlight.set({ target });
      expect(r.ok, target).toBe(true);
      // A frame between them: the lines are input, published once a frame (rendering.md §2.3).
      flushFrame();
    }
    sub.unsubscribe();
    const shown = seen.filter((l) => l.length > 0);
    expect(new Set(shown).size).toBeGreaterThanOrEqual(3);
  });

  it('ten highlights in a frame are one change of the lines (the frame lane)', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const told: string[] = [];
    const stop = v.sources.pointer.hoverLines.listen!((l) => told.push(l.join(' / ')));
    for (let i = 0; i < 10; i++) {
      const r = await v.api.highlight.set({ target: i % 2 ? 'A:10' : 'A:20' });
      expect(r.ok).toBe(true);
    }
    expect(told).toHaveLength(0);
    flushFrame();
    expect(told).toHaveLength(1);
    stop();
  });

  it('hover events come once per element under the same keys', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const structure = v.views.sequence!.source.structure.value!.structure;
    const loci = StructureElement.Loci.firstResidue(StructureElement.Loci.all(structure));
    const events: unknown[] = [];
    const sub = v.sources.pointer.hover.subscribe((e) => events.push(e));
    // The behavior replays its last value on subscribing.
    events.length = 0;
    const none = { shift: false, control: false, alt: false, meta: false };
    const hover = (modifiers: typeof none) =>
      v.plugin.behaviors.interaction.hover.next({ current: { loci }, buttons: 0, button: 0, modifiers } as never);
    hover(none);
    hover(none);
    expect(events).toHaveLength(1);
    hover({ ...none, shift: true });
    expect(events).toHaveLength(2);
    sub.unsubscribe();
  });
});
