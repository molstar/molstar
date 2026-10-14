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

describe('NumberMap with an offset and strict keys', () => {
  it('keeps the keys from the offset in the array, others in the map', () => {
    const map = new NumberMap<number, number>(5, 11); // 11..15
    for (let k = 11; k <= 15; k++) map.set(k, k - 11);
    expect([11, 13, 15].map((k) => map.get(k))).toEqual([0, 2, 4]);
    expect(map.overflowSize).toBe(0);
    map.set(10, 7);
    map.set(16, 8);
    map.set(1000, 9);
    expect([10, 16, 1000, 17].map((k) => map.get(k))).toEqual([7, 8, 9, undefined]);
    expect(map.overflowSize).toBe(3);
  });
  it('a key which is not a number has no value (an array would convert it)', () => {
    const map = new NumberMap<number, number>(5, 0);
    map.set(3, 1);
    expect(map.get('3' as any)).toBeUndefined();
    expect(map.get(3)).toBe(1);
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
