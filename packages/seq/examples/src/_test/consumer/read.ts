/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { DataStore, Shape, Tables, Want } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';

// The second consumer's reading (EC-183, `_spec` Resolved view): what a layer draws is in the store, found through the
// resolved view's layer and the model's access plans, with the public entries of `@molstar/seqx-core` alone. Core's
// own reader of the same is `render/layer-data.ts`; nothing of it is used here.

/** Where a layer's rows are in the store: its source, the name the source knows the system by, and the shape. */
export interface Access {
  readonly source: string;
  readonly system: string;
  readonly shape: Shape;
}

/** The store's source of a placed layer's data, from the layer's data and the model's plans (inline data: its node). */
export function accessOf(model: StaticModel, layer: ResolvedLayer): Access | undefined {
  const data = layer.data;
  if (!data || !layer.placed) return undefined;
  const plan = data.access ?? model.systems.get(data.address)?.access ?? model.annotations.get(data.address)?.access;
  if (!plan) return undefined;
  const source = plan.kind === 'inline' ? plan.source.address : plan.source;
  return { source, system: plan.system, shape: plan.shape };
}

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
