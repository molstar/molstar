/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query model's property hook over values channels (_spec/README.md "Queries", Properties; properties.md, How a
// property is read): a `values` annotation on a system read from the store's tiles, per position; pending while a tile
// is not loaded (the forms stage the system, as a numbering not arrived), undefined when the system has no such
// channel. `b-factor` and `occupancy` are the channels of those fields on the system (the mmCIF runtime writes them on
// a unit's atoms); any other name is a values annotation by its address or its id.

import type { PropertyColumn, PropertyRead } from '@molstar/seqx-core/query';
import type { AnnotationEntry, StaticModel } from './model.js';
import type { Overlay } from './overlay.js';
import type { RowsSource } from './query-model.js';

/** The values channels' fields of the built-in properties. */
const ChannelFields: Readonly<Record<string, string>> = { 'b-factor': 'b_factor', occupancy: 'occupancy' };

const idOf = (address: string) => address.slice(address.lastIndexOf('/') + 1);

/** The `values` annotations on a system, the model's then the overlay's. */
function valuesOn(model: StaticModel, overlay: Overlay, system: string): AnnotationEntry[] {
  const out: AnnotationEntry[] = [];
  for (const a of model.annotations.values()) if (a.shape === 'values' && a.target === system) out.push(a);
  for (const n of overlay.nodes) {
    const a = overlay.annotation(n.address);
    if (a?.shape === 'values' && a.target === system) out.push(a);
  }
  return out;
}

/** The channel a property name or annotation names on a system, or undefined. */
function channelOf(
  model: StaticModel,
  overlay: Overlay,
  system: string,
  name: string,
  field: string | undefined,
): AnnotationEntry | undefined {
  if (name.startsWith('/')) {
    const a = model.annotations.get(name) ?? overlay.annotation(name);
    return a?.shape === 'values' && a.target === system ? a : undefined;
  }
  const on = valuesOn(model, overlay, system);
  const channel = ChannelFields[name];
  if (channel) return on.find((a) => a.access?.kind === 'source' && a.access.source.endsWith(channel));
  return on.find((a) => idOf(a.address) === name && (field === undefined || a.access?.kind !== 'inline'));
}

/**
 * The hook over the store's tiles (`src`): the channel's tiles over the whole system, asked for; a column of its
 * values (NaN where none: the element lacks it), pending while any tile is not loaded.
 */
export function valuesProperty(
  model: StaticModel,
  overlay: Overlay,
  src: RowsSource | undefined,
): (system: string, name: string, field?: string) => PropertyRead {
  return (system, name, field) => {
    const a = channelOf(model, overlay, system, name, field);
    if (!a?.access) return undefined;
    if (!src) return 'pending';
    const source = a.access.kind === 'source' ? a.access.source : a.access.source.address;
    const n = src.data.coverage(source)?.systems.find((s) => s.name === a.access!.system)?.length;
    if (n === undefined) return 'pending';
    const want = { source, system: a.access.system, start: 0, end: n, level: 1 };
    src.want?.(want);
    const keys = src.data.keysOf(want);
    if (!keys.length && src.data.tooMany?.(want)) return 'pending';
    const values = new Float32Array(n).fill(Number.NaN);
    for (const k of keys) {
      const e = src.data.get(k);
      if (e?.status !== 'loaded') return 'pending';
      const t = e.tables.values;
      if (t?.form !== 'positions') continue;
      for (let i = 0; i < t.value.length && t.start + i < n; i++) values[t.start + i] = t.value[i];
    }
    const column: PropertyColumn = {
      length: n,
      get: (i) => (i >= 0 && i < n && !Number.isNaN(values[i]) ? values[i] : undefined),
    };
    return column;
  };
}
