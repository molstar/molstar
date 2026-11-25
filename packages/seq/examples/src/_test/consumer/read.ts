/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { DataStore, Tables, Want } from '@molstar/seqx-core/data';
import type { LayerSource } from '@molstar/seqx-core/render';

// The second consumer's reading (EC-183, `_spec` Resolved view): what a layer draws is in the store, found through
// `layerSource`, the render entry's data access (a layer's source, system and shape from the model's plans: not a layout
// or a draw list), then the store's tiles.

/** Where a layer's rows are in the store: its source, the name the source knows the system by, and the shape. */
export type Access = Pick<LayerSource, 'source' | 'system' | 'shape'>;

/** The tables of the tiles covering `[start, end)` of a layer's data, finest level, once they settle; `pending` counts tiles not loaded. */
export async function readTables(
  store: DataStore,
  access: Access,
  start: number,
  end: number,
): Promise<{ tables: Tables[]; pending: number }> {
  const level = store.coverage(access.source)?.levels[0] ?? 1;
  const want: Want = { source: access.source, system: access.system, start, end, level };
  const wanter = store.wanter();
  wanter.want([want]);
  await wanter.whenSettled(2000);
  const tables: Tables[] = [];
  let pending = 0;
  for (const key of store.keysOf(want)) {
    const e = store.get(key);
    if (e?.status === 'loaded') tables.push(e.tables);
    else if (e?.status !== 'empty') pending++;
  }
  wanter.dispose();
  return { tables, pending };
}

/** A cell of a column: a number, a string, or null where missing. */
export function cell(column: Column | undefined, i: number): string | number | null {
  if (!column) return null;
  return column instanceof CategoryColumn ? column.get(i) : column[i];
}
