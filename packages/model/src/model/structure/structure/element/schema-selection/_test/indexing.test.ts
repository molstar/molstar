/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { AtomIdIndices, GlobalResidueIndices } from '../indexing.js';
import { NumberMap } from '../utils.js';

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

describe('NumberMap', () => {
  it('keeps non-negative integer values for keys in and out of the limit, and has no value for the others', () => {
    const map = new NumberMap<number, number>(5);
    expect(map.get(0)).toBeUndefined();
    map.set(0, 0); // (0 is a value, not a missing one)
    map.set(3, 7);
    map.set(10, 2); // (out of the limit)
    map.set(-1, 4);
    expect([map.get(0), map.get(3), map.get(4), map.get(10), map.get(-1), map.get(11)]).toEqual([
      0,
      7,
      undefined,
      2,
      4,
      undefined,
    ]);
    map.set(3, 9);
    expect(map.get(3)).toBe(9);
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

  it('is not built until a lookup asks for it, and is built once', () => {
    const { indices, idReads } = make();
    expect(idReads.count).toBe(0);
    indices.byId;
    expect(idReads.count).toBe(ids.length);
    indices.bySourceIndex;
    indices.repeatedIds;
    expect(idReads.count).toBe(ids.length);
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
