/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure link's performance representative (seqx-viewer _spec "Structure link"): on p53 and PDB 2OCJ as the
// example records it, a hover moving along P04637 carried through `linkStructure` onto the four chains a move (EC-48),
// the side told each new image. Quick unbounded, heavy within the bound.

import { linkStructure } from '@molstar/seqx-viewer/link-structure';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, describe, expect } from 'vitest';
import { countingSide, hoverAt, open2ocj } from './link-structure-support.js';

const Protein = '/data/p53';
const Section = '/views/main/protein';

let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
});

describe('the structure link on 2OCJ', () => {
  itWithin('100 hovers on P04637 carried onto the four chains (EC-48)', async (within) => {
    viewer = await open2ocj();
    const v = viewer;
    const side = countingSide();
    const off = linkStructure(v, side, { select: 'none', focus: 'none' });
    const owner = {};
    const hit = { section: Section, system: Protein, position: hoverAt(0) };
    v.hover.set(hit, owner); // warm: the paths found once
    // measured: 100 hovers in about 2 ms (an M-series laptop); 25 ms leaves room for a parallel heavy run
    within(25, () => {
      for (let k = 1; k <= 100; k++) {
        hit.position = hoverAt(k);
        v.hover.set(hit, owner);
      }
    });
    expect(side.told).toBe(101);
    expect(side.last?.parts.length).toBe(4);
    off();
  });
});
