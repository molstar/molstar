/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column, Table } from '@molstar/core/data/db';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { BasicSchema } from '../schema.js';
import { sortAtomSite } from '../sort.js';

/** The array behind a column or field: its `'array'` storage. */
function arrayOf(c: { storage?(): { kind: string; array?: ArrayLike<any> } | undefined }) {
  const s = c.storage?.();
  return s?.kind === 'array' ? s.array : void 0;
}

// Model 1 (rows 0-4): entity 2's rows first and seq ids out of order, so it sorts. Model 2 (rows 5-8): already in
// order, a contiguous range. Typed and plain array columns, one column with a missing value.
const atom_site = Table.ofArrays(BasicSchema.atom_site, {
  id: Int32Array.of(1, 2, 3, 4, 5, 6, 7, 8, 9),
  pdbx_PDB_model_num: Int32Array.of(1, 1, 1, 1, 1, 2, 2, 2, 2),
  label_entity_id: ['2', '1', '1', '1', '2', '1', '1', '1', '2'],
  label_asym_id: ['B', 'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B'],
  label_seq_id: Int32Array.of(1, 3, 1, 2, 1, 1, 2, 3, 1),
  label_comp_id: ['HOH', 'VAL', 'ALA', 'GLY', 'HOH', 'ALA', 'GLY', 'VAL', 'HOH'],
  Cartn_x: Float32Array.of(0, 1, 2, 3, 4, 5, 6, 7, 8),
  B_iso_or_equiv: [10, 11, 12, 13, 14, 15, 16, 17, 18],
});

/** Every column read row by row: values and value kinds. */
function read(table: Table<BasicSchema['atom_site']>) {
  const out: Record<string, unknown> = { rowCount: table._rowCount };
  for (const k of table._columns) {
    const c = table[k] as Column<unknown>;
    out[k] = {
      isDefined: c.isDefined,
      rowCount: c.rowCount,
      values: Array.from({ length: c.rowCount }, (_, i) => c.value(i)),
      kinds: Array.from({ length: c.rowCount }, (_, i) => c.valueKind(i)),
      array: c.toArray(),
    };
  }
  return out;
}

describe('sortAtomSite', () => {
  it('sorts a model by entity, chain and seq id', async () => {
    // entities in their first rows' order, seq ids sorted within a chain
    const { atom_site: sorted, sourceIndex } = await sortAtomSite(SyncRuntimeContext, atom_site, 0, 5);
    expect(Array.from(sourceIndex.toArray())).toEqual([0, 4, 2, 3, 1]);
    expect(read(sorted)).toEqual(read(Table.view(atom_site, atom_site._schema, [0, 4, 2, 3, 1])));
  });

  it('reads a model already in order as the rows of its range, typed arrays in place', async () => {
    const { atom_site: sorted, sourceIndex } = await sortAtomSite(SyncRuntimeContext, atom_site, 5, 9);
    expect(Array.from(sourceIndex.toArray())).toEqual([5, 6, 7, 8]);
    expect(read(sorted)).toEqual(read(Table.view(atom_site, atom_site._schema, [5, 6, 7, 8])));
    expect((arrayOf(sorted.Cartn_x) as Float32Array).buffer).toBe((arrayOf(atom_site.Cartn_x) as Float32Array).buffer);
    expect(Array.isArray(arrayOf(sorted.label_comp_id))).toBe(true);
  });

  it('returns the table itself when it is one model in order', async () => {
    const one = Table.window(atom_site, atom_site._schema, 5, 9);
    const { atom_site: sorted } = await sortAtomSite(SyncRuntimeContext, one, 0, 4);
    expect(sorted).toBe(one);
  });

  it('makes a column of a sorted or windowed model on its first read', async () => {
    const reads = new Set<string>();
    const watched = Object.create(null);
    watched._rowCount = atom_site._rowCount;
    watched._columns = atom_site._columns;
    watched._schema = atom_site._schema;
    for (const k of atom_site._columns) {
      Object.defineProperty(watched, k, {
        get: () => {
          reads.add(k);
          return (atom_site as any)[k];
        },
        enumerable: true,
      });
    }
    for (const [start, end] of [
      [0, 5],
      [5, 9],
    ]) {
      const { atom_site: sorted } = await sortAtomSite(SyncRuntimeContext, watched, start, end);
      // the sort reads its keys only
      expect([...reads].sort()).toEqual(['label_asym_id', 'label_entity_id', 'label_seq_id']);
      expect(Array.from(sorted.Cartn_x.toArray())).toEqual(
        Array.from(
          Table.view(
            atom_site,
            atom_site._schema,
            Array.from(sorted.id.toArray(), (id) => id - 1),
          ).Cartn_x.toArray(),
        ),
      );
      expect(reads.has('Cartn_x')).toBe(true);
      expect(sorted.Cartn_x).toBe(sorted.Cartn_x);
      reads.clear();
    }
  });

  it('sorts by dictionary ids as by their values, repeated and empty values included', async () => {
    // entity and chain ids as dictionaries holding '' twice, as BinaryCIF writes them, and 'B' twice
    const values = ['', 'A', '', 'B', 'B'];
    const entity = Uint8Array.of(3, 1, 1, 4, 3, 0, 2, 1, 4);
    const chain = Uint8Array.of(1, 3, 4, 0, 2, 1, 3, 4, 1);
    const dictionaries = Table.ofColumns(BasicSchema.atom_site, {
      ...atom_site,
      label_entity_id: Column.ofDictionary({ indices: entity, values, schema: Column.Schema.str }),
      label_asym_id: Column.ofDictionary({ indices: chain, values, schema: Column.Schema.str }),
    });
    const arrays = Table.ofColumns(BasicSchema.atom_site, {
      ...atom_site,
      label_entity_id: Column.ofStringArray(Array.from(entity, (i) => values[i])),
      label_asym_id: Column.ofStringArray(Array.from(chain, (i) => values[i])),
    });
    for (const [start, end] of [
      [0, 5],
      [5, 9],
      [0, 9],
    ]) {
      const a = await sortAtomSite(SyncRuntimeContext, dictionaries, start, end);
      const b = await sortAtomSite(SyncRuntimeContext, arrays, start, end);
      expect(Array.from(a.sourceIndex.toArray())).toEqual(Array.from(b.sourceIndex.toArray()));
      expect(read(a.atom_site)).toEqual(read(b.atom_site));
    }
  });
});
