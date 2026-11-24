/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Axis } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { blocksRect, type RectBlock } from '../behaviors/blocks-rect.js';

// A rect dragged on a grid over blocks (_spec/README.md "Behaviors", A rect over blocks): the span between two cells in
// the blocks' display units, split at their boundaries, each cell's units on the hit strip's own block.

const axis = (system: string, n: number, dir: 1 | -1 = 1): Axis => ({
  system,
  length: n,
  ranges: [0, n],
  cyclic: false,
  dir,
  pieces: 'whole',
});

/** Blocks end to end: each its system, its length and its direction. */
function row(...bs: (readonly [string, number, (1 | -1)?])[]): RectBlock[] {
  let o = 0;
  return bs.map(([system, n, dir]) => {
    const b = { system, axis: axis(system, n, dir ?? 1), offset: o };
    o += n;
    return b;
  });
}

const text = (loc: ReturnType<typeof blocksRect>) =>
  loc!.rects.map((r) => `${r.a.system}:${[...r.a.ranges]} × ${r.b.system}:${[...r.b.ranges]}`).sort();

describe('blocksRect', () => {
  const ab = row(['a', 10], ['b', 20]);

  it('across a boundary: a rect a pair of blocks, each part its positions', () => {
    const r = blocksRect(
      ab,
      { block: 0, position: 7, yBlock: 0, yPosition: 8 },
      { block: 1, position: 2, yBlock: 1, yPosition: 4 },
    );
    expect(text(r)).toEqual(['a:7,10 × a:8,10', 'a:7,10 × b:0,5', 'b:0,3 × a:8,10', 'b:0,3 × b:0,5']);
  });

  it('a drag backwards gives the same rect', () => {
    const a = { block: 0, position: 7, yBlock: 0, yPosition: 8 };
    const b = { block: 1, position: 2, yBlock: 1, yPosition: 4 };
    expect(text(blocksRect(ab, b, a))).toEqual(text(blocksRect(ab, a, b)));
  });

  it('inside one block: one rect, the cells included', () => {
    const r = blocksRect(
      ab,
      { block: 1, position: 3, yBlock: 1, yPosition: 9 },
      { block: 1, position: 5, yBlock: 1, yPosition: 6 },
    );
    expect(text(r)).toEqual(['b:3,6 × b:6,10']);
  });

  it('copies split (one system twice): the hit strip’s block, not the first on its system', () => {
    const copies = row(['a', 10], ['b', 20], ['a', 10]);
    // from the second copy of a (units 30..39) to b's 15: units 15..33 on x, the first copy of a never touched
    const r = blocksRect(
      copies,
      { block: 2, position: 3, yBlock: 2, yPosition: 3 },
      { block: 1, position: 15, yBlock: 1, yPosition: 15 },
    );
    expect(text(r)).toEqual(['a:0,4 × a:0,4', 'a:0,4 × b:15,20', 'b:15,20 × a:0,4', 'b:15,20 × b:15,20']);
    // on the first copy the same positions are other units: a rect within it alone
    const first = blocksRect(
      copies,
      { block: 0, position: 3, yBlock: 0, yPosition: 3 },
      { block: 0, position: 5, yBlock: 0, yPosition: 5 },
    );
    expect(text(first)).toEqual(['a:3,6 × a:3,6']);
  });

  it('a reversed block (`dir: -1`): its positions from its last unit', () => {
    const rev = row(['a', 10, -1], ['b', 20]);
    // unit 0 of a is position 9: from a's 2 (unit 7) to b's 1 (unit 11)
    const r = blocksRect(
      rev,
      { block: 0, position: 2, yBlock: 1, yPosition: 0 },
      { block: 1, position: 1, yBlock: 1, yPosition: 0 },
    );
    expect(text(r)).toEqual(['a:0,3 × b:0,1', 'b:0,2 × b:0,1']);
  });

  it('a cell on no block: none', () => {
    expect(
      blocksRect(
        ab,
        { block: 5, position: 0, yBlock: 0, yPosition: 0 },
        { block: 0, position: 0, yBlock: 0, yPosition: 0 },
      ),
    ).toBeUndefined();
  });
});
