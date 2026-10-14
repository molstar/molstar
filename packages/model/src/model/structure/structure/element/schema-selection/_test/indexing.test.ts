/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { GlobalResidueIndices } from '../indexing.js';

describe('GlobalResidueIndices', () => {
  const hierarchy = (sourceIndex: number[]) =>
    ({
      residues: { _rowCount: sourceIndex.length },
      residueSourceIndex: Column.ofIntArray(sourceIndex),
    }) as any;

  it('maps residue_index to the residue when it is unique', () => {
    const map = new GlobalResidueIndices(hierarchy([10, 11, 12, 40])).residueBySourceIndex!;
    expect(map.get(12)).toBe(2);
    expect(map.get(40)).toBe(3);
    expect(map.get(13)).toBeUndefined();
  });
  it('is null when residue_index repeats (the per-chain indices are used then)', () => {
    expect(new GlobalResidueIndices(hierarchy([0, 1, 1, 2])).residueBySourceIndex).toBeNull();
  });
});
