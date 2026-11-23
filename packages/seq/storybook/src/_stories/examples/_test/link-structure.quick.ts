/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure link on the hover's path (seqx-viewer _spec "Structure link"): on p53 and PDB 2OCJ, a hover moving
// along P04637 (a new image a move) carried onto the four chains and told to a side, against the same hover with no
// link. Not zero: the target is a new Location a move (the side keeps it). Measured 11.3 kB a move with the link (from
// 31.8 kB before the translation's cheap cuts): `toUnits` onto four units about 9.2 kB of it, the hover's Location
// about 2 kB; the bound holds it. By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { linkStructure } from '@molstar/seqx-viewer/link-structure';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, it } from 'vitest';
import { countingSide, hoverAt, open2ocj } from './link-structure-support.js';

const Owner = {};
const hit = { section: '/views/main/protein', system: '/data/p53', position: hoverAt(0) };

function mover(v: SeqViewer): (i: number) => void {
  v.hover.set(hit, Owner);
  return (i) => {
    hit.position = hoverAt(i & 7);
    v.hover.set(hit, Owner);
  };
}

describe('the structure link on the hover’s path', () => {
  it('2OCJ: a hover along P04637 carried onto the four chains, with the link and without', async ({ task }) => {
    const without = await open2ocj();
    const v = await open2ocj();
    const side = countingSide();
    const off = linkStructure(v, side, { select: 'none', focus: 'none' });
    const plainStep = mover(without);
    const step = mover(v);
    const reads = [plainStep, step, plainStep, step].map((s) => measureBytesPerStep(s).bytesPerStep);
    const plain = Math.min(reads[0], reads[2]);
    const linked = Math.min(reads[1], reads[3]);
    const { r } = quickCase(task, 'viewer hover through a structure link', step);
    expect(side.last?.parts.length).toBe(4);
    off();
    without.retire();
    v.retire();
    console.info(
      `a 2OCJ hover a move: ${linked.toFixed(0)} bytes with the link, ${plain.toFixed(0)} without; the case ${r.bytesPerStep.toFixed(0)} (${r.note})`,
    );
    // the link's own bound: the case above has no baseline file, so this holds it
    expect(linked - plain, reads.map((n) => n.toFixed(0)).join(' / ')).toBeLessThanOrEqual(12_500);
  });
});
