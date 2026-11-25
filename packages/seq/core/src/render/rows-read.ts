/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column, issueOf } from '@molstar/core/repr/tree';
import type { CellsCodes, RangesList } from '@molstar/seqx-core/data';
import { RenderIssues } from './issues.js';
import { type LayerSource, readTiles } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { type HeldCells, RowMeasureFields } from './rows-model.js';
import type { RowsLay, RowsPrep } from './rows-lay.js';

// What a rows layer reads of the store (_spec/README.md "Render", The rows layer; Row columns): an alignment held
// whole (shared by the layers reading it), a lazy one's rows in view, and the records sources' fields joined by row key.

/** A column's value at `i` (`null` missing: a missing category, a NaN). */
function cellOf(c: Column, i: number): string | number | null {
  if (c instanceof CategoryColumn) return c.get(i);
  const v = (c as ArrayLike<number>)[i];
  return v === undefined || Number.isNaN(v) ? null : v;
}

/**
 * The fields of a records source whose system is the rows system, joined by row key (_spec/README.md "Render", Row
 * columns): a record's `id` names its row (a held row's label), else its position is the row. Only the fields a column
 * or the order names, and neither an own field nor a measure (those shadow it: said once, as a field two sources
 * have). Kept by the held rows, the prep and the entries' versions (compared in place: nothing made once warm),
 * `fieldsVersion` moved when they change. `joinPending` while a source that could join (its `system` the rows system,
 * a runtime bound) has not attached and the view has not settled: one that failed, or never attaches, holds nothing.
 */
export function readJoined(lay: RowsLay, prep: RowsPrep, ctx: LayoutContext, data: StoreReader) {
  const held = lay.held;
  lay.joinPending = false;
  if (!held || !prep.joinSources.length || !prep.joinFields.length) return;
  const of = lay.joinedOf;
  let same = of[0] === held && of[1] === prep;
  let at = 2;
  let nt = 0;
  for (let s = 0; s < prep.joinSources.length; s++) {
    const source = prep.joinSources[s];
    const coverage = data.coverage(source);
    if (!coverage) {
      if (!ctx.settled && prep.joinCandidates.includes(source)) lay.joinPending = true;
      continue;
    }
    let length = -1;
    for (let k = 0; k < coverage.systems.length; k++)
      if (coverage.systems[k].name === prep.rowsSystem) length = coverage.systems[k].length ?? held.rows;
    if (length < 0) continue;
    JoinRegion[1] = length;
    const read = readTiles(ctx, data, joinSource(lay, source, prep.rowsSystem, length), JoinRegion, 1);
    if (read.pending.length) lay.joinPending = true;
    for (let k = 0; k < read.entries.length; k++) {
      const e = read.entries[k];
      const t = e.tables.ranges as RangesList | undefined;
      if (t?.form !== 'list') continue;
      if (of[at] !== source || of[at + 1] !== t || of[at + 2] !== e.version) same = false;
      of[at] = source;
      of[at + 1] = t;
      of[at + 2] = e.version;
      at += 3;
      nt++;
    }
  }
  if (same && of.length === at) return;
  of[0] = held;
  of[1] = prep;
  of.length = at;
  lay.fieldsVersion++;
  lay.joinIssue = undefined;
  if (!nt) {
    lay.joined = undefined;
    return;
  }
  const joined = new Map<string, (string | number | null)[]>();
  const from = new Map<string, string>();
  const said = new Set<string>();
  const byLabel = new Map<string, number>();
  held.labels.forEach((l, i) => {
    if (!byLabel.has(l)) byLabel.set(l, i);
  });
  for (let q = 2; q < at; q += 3) {
    const source = of[q] as string;
    const t = of[q + 1] as RangesList;
    for (const f of t.fields ?? []) {
      if (!prep.joinFields.includes(f.name)) continue;
      if (held.fields?.has(f.name) || RowMeasureFields.includes(f.name)) {
        said.add(`${f.name} of ${source} is the alignment’s own ${held.fields?.has(f.name) ? 'field' : 'measure'}`);
        continue;
      }
      const first = from.get(f.name);
      if (first !== undefined && first !== source) {
        said.add(`${f.name} is a field of ${first} and of ${source}: ${first}’s read`);
        continue;
      }
      from.set(f.name, source);
      let a = joined.get(f.name);
      if (!a) joined.set(f.name, (a = new Array(held.rows).fill(null)));
      for (let k = 0; k < t.start.length; k++) {
        const id = t.id?.get(k);
        const row = id != null ? byLabel.get(id) : t.start[k];
        if (row === undefined || row < 0 || row >= held.rows) continue;
        a[row] = cellOf(f.column, k);
      }
    }
  }
  lay.joined = joined;
  if (said.size)
    lay.joinIssue = issueOf(
      RenderIssues,
      'seq.render.rows.field',
      `A joined field is not read: ${[...said].join('; ')}`,
      {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.rows.field'].hint,
      },
    );
}

/** The region a join reads: the whole rows system (its end set before each read: a frame makes no array). */
const JoinRegion = [0, 0];

/** The layer source a joined source is read through, kept by the layer (one a source and length). */
function joinSource(lay: RowsLay, source: string, system: string, length: number): LayerSource {
  let s = lay.joinSrc.get(source);
  if (!s || s.length !== length) {
    s = { source, system, dataSystem: system, shape: 'ranges', length, ratio: 1 };
    lay.joinSrc.set(source, s);
  }
  return s;
}

/** Joined held cells by store and source, shared by the layers reading one alignment (the rows, its logo). */
const Joined = new WeakMap<StoreReader, Map<string, { key: string; held: HeldCells }>>();

/**
 * An alignment held whole (_spec/README.md "Render", The rows layer): every tile of its rows wanted at level 1, joined
 * once they are all in, kept by the entries' versions and shared by the layers reading it on one store; undefined while
 * some are pending or failed.
 */
export function readHeld(ctx: LayoutContext, data: StoreReader, src: LayerSource, rows: number, columns: number) {
  const read = readTiles(ctx, data, src, [0, rows], 1);
  if (read.pending.length || read.failed.length) return undefined;
  let key = `${rows}|${columns}`;
  for (const e of read.entries) key += `,${e.version}`;
  let bySource = Joined.get(data);
  if (!bySource) Joined.set(data, (bySource = new Map()));
  const kept = bySource.get(src.source);
  if (kept?.key === key) return kept.held;
  const C = columns;
  const codes = new Uint8Array(rows * C);
  const labels: string[] = new Array(rows).fill('');
  const first = new Int32Array(rows).fill(1);
  let fields: Map<string, (string | number | null)[]> | undefined;
  for (const e of read.entries) {
    const t = e.tables.cells;
    if (!t || t.form !== 'codes' || t.columns !== C) continue;
    const n = Math.floor(t.codes.length / Math.max(C, 1));
    codes.set(t.codes.subarray(0, Math.min(n, rows - t.start) * C), t.start * C);
    for (let i = 0; i < n && t.start + i < rows; i++) {
      labels[t.start + i] = t.label?.get(i) ?? `row ${t.start + i + 1}`;
      if (t.first) first[t.start + i] = t.first[i];
    }
    // the rows' own fields (Row columns), a value a row
    for (const f of t.fields ?? []) {
      fields ??= new Map();
      let a = fields.get(f.name);
      if (!a) fields.set(f.name, (a = new Array(rows).fill(null)));
      for (let i = 0; i < n && t.start + i < rows; i++) a[t.start + i] = cellOf(f.column, i);
    }
  }
  const held: HeldCells = {
    rows,
    columns: C,
    codes,
    labels,
    first,
    version: (kept?.held.version ?? 0) + 1,
    ...(fields ? { fields } : {}),
  };
  bySource.set(src.source, { key, held });
  return held;
}

/** The rows system's source of a layer (its data's columns path dropped) and the alignment's size, when attached. */
export function rowsSource(src: LayerSource, data: StoreReader) {
  const coverage = data.coverage(src.source);
  const rows = coverage?.systems.find((s) => s.name === src.system)?.length;
  const columns = coverage?.systems.find((s) => s.name !== src.system)?.length;
  const { path: _p, back: _b, ...rest } = src;
  return { src: { ...rest, ...(rows !== undefined ? { length: rows } : {}) } as LayerSource, rows, columns };
}

/** A lazy alignment: the tiles of the rows in view at level 1 (their labels and cells), the rest pending. */
export function readLazy(lay: RowsLay, ctx: LayoutContext, data: StoreReader, src: LayerSource) {
  const { first, last } = lay;
  if (last <= first) return;
  // bounded by the screen: a band of more rows than device pixels reads none; the band a neutral fill (a source's
  // coarse row levels, later)
  if ((last - first) * lay.rowScale * ctx.dpr < (last - first) * 0.99) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.rows.reduced', 'The rows are below a pixel and the source serves no row bins', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.rows.reduced'].hint,
      }),
    );
    lay.neutral = true;
    return;
  }
  const read = readTiles(ctx, data, src, [first, last], 1);
  for (const e of read.entries) {
    const t = e.tables.cells as CellsCodes | undefined;
    if (!t || t.form !== 'codes') continue;
    const C = t.columns;
    const rows = Math.floor(t.codes.length / Math.max(C, 1));
    const a = Math.max(first, t.start);
    const b = Math.min(last, t.start + rows);
    for (let r = a; r < b; r++) {
      lay.lazyLabels.set(r, t.label?.get(r - t.start) ?? `row ${r + 1}`);
      lay.lazyCodes.set(r, t.codes.subarray((r - t.start) * C, (r - t.start + 1) * C));
      if (t.first) lay.lazyFirst.set(r, t.first[r - t.start]);
    }
  }
  for (let i = 0; i + 1 < read.pending.length; i += 2)
    for (let r = read.pending[i]; r < read.pending[i + 1]; r++) {
      lay.pending.push(lay.rowY(r), lay.rowScale);
    }
}
