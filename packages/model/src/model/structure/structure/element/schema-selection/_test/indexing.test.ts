/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { AtomIdIndices, GlobalResidueIndices } from '../indexing.js';
import type { NumberMap } from '@molstar/core/util/number-map';

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

describe('AtomIdIndices', () => {
  const ids = [5, 6, 5, 9, 6, 5, 100]; // (ids 5 and 6 are repeated)
  const sourceIndex = [0, 1, 2, 3, 4, 5, 6];
  const make = () => {
    const atomId = Column.ofIntArray(ids);
    const atomSourceIndex = Column.ofIntArray(sourceIndex);
    const idReads = { count: 0 };
    const counting = {
      ...atomId,
      value: (i: number) => (idReads.count++, atomId.value(i)),
    } as unknown as Column<number>;
    return { indices: new AtomIdIndices(ids.length, counting, atomSourceIndex), idReads };
  };

  it('ids which do not start near 0 (the 2nd model of a file) are in the table, not in the overflow map', () => {
    const second = Array.from({ length: 10 }, (_, i) => 11 + i); // ids 11..20 for 10 atoms
    const indices = new AtomIdIndices(
      10,
      Column.ofIntArray(second),
      Column.ofIntArray(Array.from({ length: 10 }, (_, i) => i)),
    );
    expect([11, 15, 20, 10, 21].map((id) => indices.byId.get(id))).toEqual([0, 4, 9, undefined, undefined]);
    expect((indices.byId as NumberMap<number, number>).overflowSize).toBe(0);
    // very sparse ids do not make a huge table: the table is not larger than a small multiple of the number of atoms
    const sparse = new AtomIdIndices(3, Column.ofIntArray([1, 1_000_000_000, 7]), Column.ofIntArray([0, 1, 2]));
    expect([1, 1_000_000_000, 7].map((id) => sparse.byId.get(id))).toEqual([0, 1, 2]);
  });
  it('is not built until a lookup asks for it, and is built once', () => {
    const { indices, idReads } = make();
    expect(idReads.count).toBe(0);
    indices.byId;
    const afterBuild = idReads.count;
    expect(afterBuild).toBeGreaterThanOrEqual(ids.length); // (every id is read)
    indices.bySourceIndex;
    indices.repeatedIds;
    indices.byId;
    expect(idReads.count).toBe(afterBuild);
  });
  it('the first atom of each id, all the atoms of the repeated ids, and the atom of each source index', () => {
    const { indices } = make();
    expect([5, 6, 9, 100, 7].map((id) => indices.byId.get(id))).toEqual([0, 1, 3, 6, undefined]);
    expect(indices.repeatedIds.get(5)).toEqual([0, 2, 5]);
    expect(indices.repeatedIds.get(6)).toEqual([1, 4]);
    expect(indices.repeatedIds.get(9)).toBeUndefined();
    expect([0, 3, 6, 7].map((i) => indices.bySourceIndex.get(i))).toEqual([0, 3, 6, undefined]);
  });
});
