/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';

// The hover label (G8): Mol*'s label manager reuses its labels array, so each highlight must still reach the lines.

describe('the hover lines', () => {
  it('follow every highlight, one after another', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const seen: string[] = [];
    const sub = v.sources.pointer.hoverLines.subscribe((l) => seen.push(l.join(' / ')));
    for (const target of ['A:10', 'A:20', 'ligand']) {
      const r = await v.api.highlight.set({ target });
      expect(r.ok, target).toBe(true);
    }
    sub.unsubscribe();
    const shown = seen.filter((l) => l.length > 0);
    expect(new Set(shown).size).toBeGreaterThanOrEqual(3);
  });
});
