/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { DataEntry, RangesList, Want } from '@molstar/seqx-core/data';
import {
  findNumbering,
  labelsNumbering,
  type Numbering,
  type NumberingOf,
  type QueryModel,
  type QuerySystem,
} from '@molstar/seqx-core/query';
import type { NumberingEntry, StaticModel } from '@molstar/seqx-core/resolve';
import type { StoreReader } from './layer-type.js';

// A source's numbering (_spec/README.md "Render", Numberings from a source): its labels are its rows' (a row a
// position: its `position` or its `start`; the label its `field`), read whole on the system at level 1, held once
// every tile of the system is in (labels need no order, so a half-read numbering would snap and print wrong). A layer
// that shows the numbering (a ruler, a text's numbers) wants its tiles; the query model reads what the store holds.

/** What a source's numbering read last: the entries it was built from (compared in place) and the numbering. */
interface Read {
  readonly entries: (DataEntry | undefined)[];
  numbering: Numbering | undefined;
}

const reads = new WeakMap<NumberingEntry, Read>();

/** The length a source numbering covers: its sequence's, else the source's own system's coverage; undefined: unknown. */
function lengthOf(model: StaticModel, entry: NumberingEntry, data: StoreReader): number | undefined {
  const own = model.systems.get(entry.system)?.length;
  if (own !== undefined) return own;
  const name = entry.access?.system;
  const c = entry.source ? data.coverage(entry.source) : undefined;
  for (const s of c?.systems ?? []) if (s.name === name) return s.length;
  return undefined;
}

/** A cell's label: a category's text, a number's spelling; null for none. */
function cellLabel(c: Column, i: number): string | null {
  if (c instanceof CategoryColumn) return c.get(i);
  const v = (c as ArrayLike<number>)[i];
  return v === undefined || Number.isNaN(v) ? null : String(v);
}

/**
 * A source's numbering from the store: undefined until every tile of its system is loaded (or when its source has no
 * system of that name). `want`, when given, asks for the tiles (a layer's `ctx.want`); without it only what is held is
 * read. Kept per entry while its tiles are the same entries.
 */
export function sourceNumbering(
  model: StaticModel,
  entry: NumberingEntry,
  data: StoreReader,
  want?: (w: Want) => void,
): Numbering | undefined {
  if (!entry.source || !entry.access) return undefined;
  const length = lengthOf(model, entry, data);
  if (length === undefined) return undefined;
  const w: Want = { source: entry.source, system: entry.access.system, start: 0, end: length, level: 1 };
  want?.(w);
  if (!data.coverage(entry.source)) return undefined;
  const keys = data.keysOf(w);
  const entries = keys.map((k) => data.get(k));
  let read = reads.get(entry);
  if (read && read.entries.length === entries.length && read.entries.every((e, i) => e === entries[i]))
    return read.numbering;
  read = { entries, numbering: undefined };
  reads.set(entry, read);
  if (!entries.every((e) => e?.status === 'loaded' || e?.status === 'empty')) return undefined;
  const labels: (string | null)[] = new Array(length).fill(null);
  const field = entry.field ?? 'label';
  for (const e of entries) {
    if (e?.status !== 'loaded') continue;
    const t = e.tables.ranges as RangesList | undefined;
    if (t?.form !== 'list') continue;
    const column = field === 'id' ? t.id : t.fields?.find((f) => f.name === field)?.column;
    if (!column) continue;
    for (let k = 0; k < t.start.length; k++) {
      const p = t.start[k];
      // a row a position; the first row of a position keeps it (a tile edge repeats none: rows are one position)
      if (p < 0 || p >= length || labels[p] !== null) continue;
      labels[p] = cellLabel(column, k);
    }
  }
  read.numbering = labelsNumbering(entry.id, { length, get: (i) => labels[i] ?? null });
  return read.numbering;
}

/**
 * A system's numbering by id as a layer shows it: the query lookup (built-ins, inline labels, a molecule's), else a
 * source's numbering of the system read from the store (its tiles wanted through `want`); undefined when none is,
 * or a source's labels are not all in yet.
 */
export function shownNumbering(
  model: StaticModel,
  system: QuerySystem,
  id: string,
  data: StoreReader,
  want?: (w: Want) => void,
): NumberingOf | undefined {
  const found = findNumbering(system, id);
  if (found) return found;
  const entry = model.numberingOf(system.address, id);
  const numbering = entry ? sourceNumbering(model, entry, data, want) : undefined;
  return numbering ? { numbering, shift: 0 } : undefined;
}

/**
 * A query model whose systems also answer their sources' numberings from what the store holds (no want): what an
 * address, a link and the copied address read once a layer has brought the labels. Everything else is `base`'s.
 */
export function withSourceNumberings(base: QueryModel, model: StaticModel, data: StoreReader): QueryModel {
  if (![...model.numberings.values()].some((n) => n.source)) return base;
  const systems = new Map<string, QuerySystem | undefined>();
  return {
    system(address) {
      if (systems.has(address)) return systems.get(address);
      const s = base.system(address);
      const out: QuerySystem | undefined = s && {
        address: s.address,
        length: s.length,
        ...(s.cyclic !== undefined ? { cyclic: s.cyclic } : {}),
        ...(s.parent ? { parent: s.parent } : {}),
        numberingIds: () => s.numberingIds(),
        numbering: (id) => {
          const own = s.numbering(id);
          if (own) return own;
          const entry = model.numberingOf(address, id);
          return entry ? sourceNumbering(model, entry, data) : undefined;
        },
      };
      systems.set(address, out);
      return out;
    },
    annotation: (address) => base.annotation(address),
    ...(base.node ? { node: (owner, key) => base.node!(owner, key) } : {}),
    ...(base.chain ? { chain: (owner, id) => base.chain!(owner, id) } : {}),
    ...(base.row ? { row: (rows, key) => base.row!(rows, key) } : {}),
  };
}
