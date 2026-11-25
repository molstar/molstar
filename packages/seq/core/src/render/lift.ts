/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover (_spec/README.md "Render", A liftover; EC-212): data over one assembly drawn on a sequence of another
// through a map read by region (`map.source`). The map's runs are read on the section's sequence by the end there
// (`RegionRequest.side`: `b`, the map's `to`), the section's window carried back along them onto each far sequence,
// the data read there, and each row carried forward along the runs: a part that maps nowhere is not drawn (EC-21), a
// part with several images is drawn at each (EC-20), its thick part carried with it and its strand turned by an
// inverted run. Hover and selection cross the map both ways, and between images of one item through its own system.
// Loaded with the genome chunk: a liftover is over assemblies.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { CategoryColumn, type Column, categoryColumn, issueOf } from '@molstar/core/repr/tree';
import type { DataEntry, PairsList, RangesList, Want } from '@molstar/seqx-core/data';
import {
  type AccessPlan,
  type AnnotationEntry,
  type Run,
  RunsLinear,
  runTargets,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import { regionMaps } from '@molstar/seqx-core/view/lift';
import { spanText } from './accessibility.js';
import { compileRowFilter } from './filter.js';
import { RenderIssues } from './issues.js';
import {
  dataRegions,
  type LayerSource,
  type LiftSource,
  mergeNear,
  pairsLevels,
  placeholderSpans,
  readTiles,
} from './layer-data.js';
import { collectItems, type DensityTile, type ItemRows, type ItemsBase, newItemRows } from './items.js';
import type { LayoutContext, StoreReader } from './layer-type.js';

export { regionMaps };

type Loaded = Extract<DataEntry, { status: 'loaded' }>;
type Fields = RangesList['fields'];

/** The items level of a read of `span` positions: the finest whose tiles a want's bound holds. */
function itemsLevel(data: StoreReader, source: string, span: number): number {
  return pairsLevels(data.coverage(source)?.levels ?? [1], 1, span, data).items;
}

/** A pairs list's rows as runs from their far end to their near end on `system` into `out`; the unequal rows' count. */
function runsOf(t: PairsList, system: string, onB: boolean, seen: Set<number>, out: Map<string, Run[]>): number {
  let unequal = 0;
  for (let i = 0; i < t.aStart.length; i++) {
    const row = t.row?.[i];
    if (row !== undefined) {
      if (seen.has(row)) continue;
      seen.add(row);
    }
    const aName = t.aSystem?.get(i) ?? system;
    const bName = t.bSystem?.get(i) ?? system;
    if ((onB ? bName : aName) !== system) continue;
    const length = t.aEnd[i] - t.aStart[i];
    if (length <= 0 || t.bEnd[i] - t.bStart[i] !== length) {
      unequal++;
      continue;
    }
    const far = onB ? aName : bName;
    const farStart = onB ? t.aStart[i] : t.bStart[i];
    const nearStart = onB ? t.bStart[i] : t.aStart[i];
    const dir: 1 | -1 = (t.aStrand?.[i] || 1) * (t.bStrand?.[i] || 1) < 0 ? -1 : 1;
    const run: Run = { from: farStart, to: dir === 1 ? nearStart : nearStart + length - 1, length, dir };
    let list = out.get(far ?? '');
    if (!list) out.set(far ?? '', (list = []));
    list.push(run);
  }
  return unequal;
}

/** Whether two reads are the same entries (by identity, in order): what was built from them holds. */
function sameEntries(a: readonly DataEntry[] | undefined, b: readonly DataEntry[]): boolean {
  if (!a || a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** A table's field by name. */
function fieldOf(fields: Fields, name: string): Column | undefined {
  if (!fields) return undefined;
  for (let k = 0; k < fields.length; k++) if (fields[k].name === name) return fields[k].column;
  return undefined;
}

/** A column's number at `i` (a strand's `+`/`-` as 1/−1), NaN when none. */
function numberAt(c: Column | undefined, i: number): number {
  if (!c) return Number.NaN;
  if (c instanceof CategoryColumn) {
    const v = c.get(i);
    return v === '+' ? 1 : v === '-' ? -1 : Number.NaN;
  }
  const v = (c as ArrayLike<number>)[i];
  return typeof v === 'number' ? v : Number.NaN;
}

/** The coordinate fields a liftover carries through the runs (the rest ride along as they are). */
const ThickStart = 'thick_start';
const ThickEnd = 'thick_end';
const Strand = 'strand';

/**
 * What a lifted layer keeps between layouts (_spec/README.md "Render", A liftover: kept): the runs' tiles and the runs
 * built from them, each far system's tiles, and the rows carried from them, so a layout whose reads hold rebuilds
 * nothing (a pan within the tiles copies the rows).
 */
interface LiftState {
  runEntries: readonly DataEntry[] | undefined;
  byFar: Map<string, RunsLinear>;
  unequal: number;
  farEntries: Map<string, readonly DataEntry[]>;
  built: boolean;
  /** The rows carried, on the near system: copied into the lay's each layout. */
  readonly rows: ItemRows;
  readonly scratch: ItemRows;
  readonly buf: number[];
  readonly back: number[];
  readonly occ: Map<string | number, number>;
  /** Each carried row's own table and row, its thick part and strand as carried. */
  readonly srcFields: Fields[];
  readonly srcRow: number[];
  readonly thickS: number[];
  readonly thickE: number[];
  readonly strand: number[];
}

const States = new WeakMap<ItemsBase, LiftState>();
function stateOf(lay: ItemsBase): LiftState {
  let s = States.get(lay);
  if (!s)
    States.set(
      lay,
      (s = {
        runEntries: undefined,
        byFar: new Map(),
        unequal: 0,
        farEntries: new Map(),
        built: false,
        rows: newItemRows(),
        scratch: newItemRows(),
        buf: [],
        back: [],
        occ: new Map(),
        srcFields: [],
        srcRow: [],
        thickS: [],
        thickE: [],
        strand: [],
      }),
    );
  return s;
}

/**
 * A lifted layer's items into `lay.rows`, on the section's system (`readItems`' for `src.lift`): each far row carried
 * along the runs, a row an image (EC-20; a part with none is not drawn, EC-21), the placeholders where the runs or the
 * far data are not read. No density tables: the far data is read as items.
 */
export function readLifted(
  lay: ItemsBase,
  ctx: LayoutContext,
  data: StoreReader,
  src: LayerSource,
  fields: readonly string[],
): DensityTile[] {
  const lift = src.lift!;
  const axis = ctx.axis;
  const st = stateOf(lay);
  const regions = dataRegions(src, axis, true);
  let span = 0;
  for (let i = 0; i + 1 < regions.length; i += 2) span += regions[i + 1] - regions[i];
  const near = lift.near;
  const read = readTiles(ctx, data, near, regions, itemsLevel(data, near.source, span), { items: true });
  placeholderSpans(src, axis, read.pending, lay.pending);
  placeholderSpans(src, axis, read.failed, lay.failed);
  if (!sameEntries(st.runEntries, read.entries)) {
    const byName = new Map<string, Run[]>();
    const seen = new Set<number>();
    st.unequal = 0;
    for (const e of read.entries as Loaded[]) {
      const t = e.tables.pairs;
      if (t?.form === 'list') st.unequal += runsOf(t, near.system, lift.side === 'b', seen, byName);
    }
    st.byFar.clear();
    for (const [name, runs] of byName) st.byFar.set(name, new RunsLinear(runs));
    st.runEntries = read.entries;
    st.built = false;
  }
  if (st.unequal)
    ctx.report(
      issueOf(
        RenderIssues,
        'seq.render.lift.run',
        `${st.unequal} rows of the map have ends of unequal length: no runs`,
        {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.lift.run'].hint,
        },
      ),
    );
  const { buf, back } = st;
  for (const [name, linear] of st.byFar) {
    const far = lift.far.get(name);
    if (!far) continue;
    // the near window carried back onto the far system
    back.length = 0;
    const inv = linear.inverse();
    for (let i = 0; i + 1 < regions.length; i += 2) inv.targets(regions[i], regions[i + 1], back);
    const farRegions = mergeNear(HalfOpenRanges.normalize(back));
    let farSpan = 0;
    for (let i = 0; i + 1 < farRegions.length; i += 2) farSpan += farRegions[i + 1] - farRegions[i];
    const farRead = readTiles(ctx, data, far, farRegions, itemsLevel(data, far.source, farSpan), { items: true });
    for (let k = 0; k < 2; k++) {
      const list = k === 0 ? farRead.pending : farRead.failed;
      if (!list.length) continue;
      buf.length = 0;
      for (let i = 0; i + 1 < list.length; i += 2) linear.targets(list[i], list[i + 1], buf);
      placeholderSpans(src, axis, buf, k === 0 ? lay.pending : lay.failed);
    }
    if (!sameEntries(st.farEntries.get(name), farRead.entries)) {
      st.farEntries.set(name, farRead.entries);
      st.built = false;
    }
  }
  if (!st.built) build(st, lift, fields);
  copyRows(st.rows, lay.rows);
  return [];
}

/** The rows carried from the reads held (`st.byFar`, `st.farEntries`), their fields carried too. */
function build(st: LiftState, lift: LiftSource, fields: readonly string[]) {
  const rows = st.rows;
  rows.n = 0;
  rows.filtered = 0;
  rows.clipped = 0;
  rows.leftOut = 0;
  rows.clamped = 0;
  st.occ.clear();
  for (const [name, linear] of st.byFar) {
    const far = lift.far.get(name);
    const entries = st.farEntries.get(name);
    if (!far || !entries) continue;
    const wrap = far.length !== undefined && far.wraps ? { length: far.length, cyclic: false } : undefined;
    collectItems(st.scratch, entries as Loaded[], fields, far.filter, far.base, wrap);
    carry(st, linear);
    rows.filtered += st.scratch.filtered;
    rows.clipped += st.scratch.clipped;
    rows.leftOut += st.scratch.leftOut;
  }
  carriedFields(st);
  st.built = true;
}

/**
 * Each far row's images along `linear` appended to the carried rows (an image a row: its runs' targets joined where
 * they meet in one direction), its thick part carried and clipped to the image (thin where none of it lands there) and
 * its strand turned by an inverted run.
 */
function carry(st: LiftState, linear: RunsLinear) {
  const from = st.scratch;
  const out = st.rows;
  const pieces = st.buf;
  for (let j = 0; j < from.n; j++) {
    const runs = linear.forward(from.starts[j], from.ends[j]);
    if (!runs.length) continue;
    // `[start, end, dir]` a piece, joined where two runs' targets meet in one direction
    pieces.length = 0;
    for (const r of runs) {
      const [t0, t1] = runTargets(r);
      pieces.push(t0, t1, r.dir);
    }
    sortTriples(pieces);
    const f = from.fields[j];
    const row = from.row[j];
    const ts = numberAt(fieldOf(f, ThickStart), row);
    const te = numberAt(fieldOf(f, ThickEnd), row);
    const fieldStrand = numberAt(fieldOf(f, Strand), row);
    const strand = Number.isNaN(fieldStrand) ? from.strand[j] : fieldStrand;
    // the thick part's images, when it has one
    const thick = !Number.isNaN(ts) && !Number.isNaN(te) && te > ts ? linear.forward(ts, te) : undefined;
    const key = from.keys[j];
    for (let k = 0; k < pieces.length; ) {
      const s = pieces[k];
      let e = pieces[k + 1];
      const dir = pieces[k + 2];
      for (k += 3; k < pieces.length && pieces[k] <= e && pieces[k + 2] === dir; k += 3) e = Math.max(e, pieces[k + 1]);
      const n = out.n++;
      const o = st.occ.get(key) ?? 0;
      st.occ.set(key, o + 1);
      out.keys[n] = key;
      out.starts[n] = s;
      out.ends[n] = e;
      out.values[n] = from.values[j];
      out.labels[n] = from.labels[j];
      out.occ[n] = o;
      out.strand[n] = strand * dir;
      out.itemStarts[n] = -1;
      out.itemEnds[n] = -1;
      st.srcFields[n] = f;
      st.srcRow[n] = row;
      st.strand[n] = strand * dir;
      if (Number.isNaN(ts) || Number.isNaN(te)) st.thickS[n] = st.thickE[n] = Number.NaN;
      else {
        // the thick images within this piece, their hull; none: the piece is thin (an empty thick part)
        let lo = Number.POSITIVE_INFINITY;
        let hi = Number.NEGATIVE_INFINITY;
        for (const r of thick ?? []) {
          const [t0, t1] = runTargets(r);
          const a = Math.max(t0, s);
          const b = Math.min(t1, e);
          if (b > a) {
            lo = Math.min(lo, a);
            hi = Math.max(hi, b);
          }
        }
        st.thickS[n] = hi > lo ? lo : s;
        st.thickE[n] = hi > lo ? hi : s;
      }
    }
  }
}

/** `[start, end, dir, …]` sorted by start (an insertion sort: a row's runs are a few). */
function sortTriples(a: number[]) {
  for (let i = 3; i < a.length; i += 3) {
    const s = a[i];
    const e = a[i + 1];
    const d = a[i + 2];
    let j = i - 3;
    for (; j >= 0 && a[j] > s; j -= 3) {
      a[j + 3] = a[j];
      a[j + 4] = a[j + 1];
      a[j + 5] = a[j + 2];
    }
    a[j + 3] = s;
    a[j + 4] = e;
    a[j + 5] = d;
  }
}

/**
 * The carried rows' fields: one table over them (a row its own), each far table's columns gathered by row, the thick
 * part and the strand as carried (what the gene line and the themes read), so no coordinate of the far system remains.
 */
function carriedFields(st: LiftState) {
  const rows = st.rows;
  const n = rows.n;
  const names = new Map<string, boolean>();
  const tables = new Set<Fields>();
  for (let j = 0; j < n; j++) tables.add(st.srcFields[j]);
  for (const f of tables)
    for (const c of f ?? []) names.set(c.name, (names.get(c.name) ?? false) || c.column instanceof CategoryColumn);
  const out: { name: string; column: Column }[] = [];
  for (const [name, category] of names) {
    if (name === ThickStart || name === ThickEnd || name === Strand) continue;
    if (category) {
      const values: (string | null)[] = [];
      for (let j = 0; j < n; j++) {
        const c = fieldOf(st.srcFields[j], name);
        const v =
          c instanceof CategoryColumn ? c.get(st.srcRow[j]) : c ? String((c as ArrayLike<number>)[st.srcRow[j]]) : null;
        values.push(v ?? null);
      }
      out.push({ name, column: categoryColumn(values) });
    } else {
      const values = new Float64Array(n);
      for (let j = 0; j < n; j++) values[j] = numberAt(fieldOf(st.srcFields[j], name), st.srcRow[j]);
      out.push({ name, column: values });
    }
  }
  if (names.has(ThickStart) || names.has(ThickEnd)) {
    out.push({ name: ThickStart, column: Float64Array.from(st.thickS.slice(0, n)) });
    out.push({ name: ThickEnd, column: Float64Array.from(st.thickE.slice(0, n)) });
  }
  out.push({ name: Strand, column: Int8Array.from(st.strand.slice(0, n), (v) => (v < 0 ? -1 : v > 0 ? 1 : 0)) });
  const fields: Fields = out;
  for (let j = 0; j < n; j++) {
    rows.fields[j] = fields;
    rows.row[j] = j;
  }
}

/** `from`'s rows copied into `to` (its arrays kept, written over). */
function copyRows(from: ItemRows, to: ItemRows) {
  const n = from.n;
  to.n = n;
  to.filtered = from.filtered;
  to.clipped = from.clipped;
  to.leftOut = from.leftOut;
  to.clamped = from.clamped;
  for (let j = 0; j < n; j++) {
    to.keys[j] = from.keys[j];
    to.starts[j] = from.starts[j];
    to.ends[j] = from.ends[j];
    to.values[j] = from.values[j];
    to.labels[j] = from.labels[j];
    to.occ[j] = from.occ[j];
    to.fields[j] = from.fields[j];
    to.strand[j] = from.strand[j];
    to.row[j] = from.row[j];
    to.itemStarts[j] = from.itemStarts[j];
    to.itemEnds[j] = from.itemEnds[j];
  }
  to.link.length = n;
  to.link.fill(-1);
}

// ---- hover and selection through the map ----

/** What a projection through a region map reads of the store: the tiles held now, never a want. */
export interface RegionReader {
  coverage(source: string): { readonly levels?: readonly number[] } | undefined;
  keysOf(want: Want): string[];
  get(key: string): DataEntry | undefined;
}

/** The side a sequence is on in a region map's annotation: `a` its `from` assembly's, `b` its `to`'s. */
function sideOf(model: StaticModel, m: AnnotationEntry, system: string): 'a' | 'b' | undefined {
  if (model.assemblies.get(m.target)?.systems.includes(system)) return 'a';
  if (m.ySystem !== undefined && model.assemblies.get(m.ySystem)?.systems.includes(system)) return 'b';
  return undefined;
}

/** The region map joining `system` to the other assembly, and the side `system` is on. */
function mapOf(model: StaticModel, system: string): { m: AnnotationEntry; side: 'a' | 'b' } | undefined {
  for (const address of regionMaps(model.spec)) {
    const m = model.annotations.get(address);
    const side = m?.bySystem?.has(system) ? sideOf(model, m, system) : undefined;
    if (m && side) return { m, side };
  }
  return undefined;
}

/** A map's systems on one side by the names its source lists them, kept per annotation. */
const Named = new WeakMap<AnnotationEntry, Map<string, string>[]>();
function systemNamed(model: StaticModel, m: AnnotationEntry, side: 'a' | 'b', name: string): string | undefined {
  let by = Named.get(m);
  if (!by) {
    by = [new Map(), new Map()];
    for (const [system, plan] of m.bySystem ?? []) {
      const s = sideOf(model, m, system);
      if (s) by[s === 'a' ? 0 : 1].set(plan.system, system);
    }
    Named.set(m, by);
  }
  return by[side === 'a' ? 0 : 1].get(name);
}

/**
 * The runs from `system` (on `side`) to each system of the other side with an end over `ranges`, from the items tiles
 * the store holds at any level (read by `system`'s end), by the far system: each a linear map from `system` to it.
 */
function runsFrom(
  model: StaticModel,
  data: RegionReader,
  m: AnnotationEntry,
  side: 'a' | 'b',
  system: string,
  ranges: ArrayLike<number>,
): Map<string, RunsLinear> {
  const out = new Map<string, RunsLinear>();
  const near = m.bySystem?.get(system);
  if (!near || near.kind !== 'source') return out;
  const byName = new Map<string, Run[]>();
  const seen = new Set<number>();
  for (const level of data.coverage(near.source)?.levels ?? [1])
    for (let i = 0; i + 1 < ranges.length; i += 2) {
      const want: Want = {
        source: near.source,
        system: near.system,
        start: ranges[i],
        end: ranges[i + 1],
        level,
        items: true,
        ...(near.side ? { side: near.side } : {}),
      };
      for (const key of data.keysOf(want)) {
        const e = data.get(key);
        const t = e?.status === 'loaded' ? e.tables.pairs : undefined;
        if (t?.form === 'list') runsOf(t, near.system, side === 'b', seen, byName);
      }
    }
  const other = side === 'a' ? 'b' : 'a';
  for (const [name, runs] of byName) {
    const far = systemNamed(model, m, other, name);
    // the runs go from the far end to the near: read back, from `system` to the far one
    if (far !== undefined) out.set(far, new RunsLinear(runs).inverse() as RunsLinear);
  }
  return out;
}

/**
 * `ranges` on `from` carried onto `to` through a map read by region (`RegionProjector.project`): across the map when
 * they are on its two sides, through the other side and back when on one (the other images of one item, EC-20), every
 * image appended to `out`; false when no such map joins them.
 */
export function imagesThrough(
  model: StaticModel,
  data: RegionReader,
  from: string,
  ranges: ArrayLike<number>,
  to: string,
  out: number[],
): boolean {
  const a = mapOf(model, from);
  if (!a || !a.m.bySystem?.has(to)) return false;
  const toSide = sideOf(model, a.m, to);
  if (!toSide) return false;
  const runs = runsFrom(model, data, a.m, a.side, from, ranges);
  if (toSide !== a.side) {
    const l = runs.get(to);
    if (l) for (let i = 0; i + 1 < ranges.length; i += 2) l.targets(ranges[i], ranges[i + 1], out);
    return true;
  }
  // one side: through each system of the other the ranges reach, and back
  const mid: number[] = [];
  for (const [system, l] of runs) {
    mid.length = 0;
    for (let i = 0; i + 1 < ranges.length; i += 2) l.targets(ranges[i], ranges[i + 1], mid);
    const back = runsFrom(model, data, a.m, toSide === 'a' ? 'b' : 'a', system, mid).get(to);
    if (back) for (let i = 0; i + 1 < mid.length; i += 2) back.targets(mid[i], mid[i + 1], out);
  }
  return true;
}

/** The most places a card names: past it, the count alone. */
const MaxNamed = 4;

/**
 * `[start, end)` on `system` as places on the map's `to` side (`RegionProjector.places`): from the `from` side across
 * the map, from the `to` side through the `from` side and back (an image's item's places, itself among them). A place
 * is a colinear chain of the item's runs (EC-20): its pieces in the source's order, on one sequence, in one direction,
 * each after the last on both; a run over a source position already placed (a duplication), another sequence, a turn
 * or a step back starts another. So an insertion on the new side splits no place. Named when at least two (one is the
 * projected hover itself), else undefined.
 */
export function liftPlaces(
  model: StaticModel,
  data: RegionReader,
  view: ResolvedView,
  system: string,
  start: number,
  end: number,
): string | undefined {
  const a = mapOf(model, system);
  if (!a) return undefined;
  // the item's extent on the `from` side, by system
  const sources = new Map<string, number[]>();
  if (a.side === 'a') sources.set(system, [start, end]);
  else
    for (const [s, l] of runsFrom(model, data, a.m, 'b', system, [start, end])) {
      const r: number[] = [];
      l.targets(start, end, r);
      sources.set(s, Array.from(HalfOpenRanges.normalize(r)));
    }
  const shown = new Set<string>();
  for (const s of view.sections) for (const b of s.blocks ?? []) shown.add(b.axis.system);
  const named: string[] = [];
  let n = 0;
  for (const [source, ranges] of sources)
    for (const [to, l] of runsFrom(model, data, a.m, 'a', source, ranges)) {
      if (!shown.has(to)) continue;
      const name = model.systems.get(to)?.name ?? to;
      for (let i = 0; i + 1 < ranges.length; i += 2) {
        const runs = l.forward(ranges[i], ranges[i + 1]);
        let pEnd = Number.NEGATIVE_INFINITY;
        let pDir = 0;
        let lo = 0;
        let hi = 0;
        let open = false;
        const close = () => {
          if (!open) return;
          n++;
          if (named.length < MaxNamed) named.push(`${name} ${spanText(lo, hi)}`);
        };
        for (const r of runs) {
          const [t0, t1] = runTargets(r);
          const on = open && r.dir === pDir && r.from >= pEnd && (r.dir === 1 ? t0 >= hi : t1 <= lo);
          if (!on) {
            close();
            open = true;
            lo = t0;
            hi = t1;
          } else {
            lo = Math.min(lo, t0);
            hi = Math.max(hi, t1);
          }
          pDir = r.dir;
          pEnd = r.from + r.length;
        }
        close();
      }
    }
  if (n < 2) return undefined;
  return `${n} places: ${named.join(', ')}${n > named.length ? ', …' : ''}`;
}

/** A projection's or a card's places kept for a version of the store (`RegionProjector`: a warm hover reads them). */
interface Kept {
  version: number;
  model: unknown;
  view: unknown;
  from: string;
  to: string;
  start: number;
  end: number;
  ok: boolean;
  readonly out: number[];
  places: string | undefined;
}

/** The kept answers of one projector: the last few hovers' (a ring). */
export interface LiftCache {
  readonly kept: Kept[];
  next: number;
}

const CacheSize = 32;

export function newLiftCache(): LiftCache {
  return { kept: [], next: 0 };
}

function slot(cache: LiftCache): Kept {
  let k = cache.kept[cache.next];
  if (!k)
    cache.kept[cache.next] = k = {
      version: -1,
      model: undefined,
      view: undefined,
      from: '',
      to: '',
      start: 0,
      end: 0,
      ok: false,
      out: [],
      places: undefined,
    };
  cache.next = (cache.next + 1) % CacheSize;
  return k;
}

/** `imagesThrough`, a single range's answer kept for the store's version (0 B a warm hover). */
export function projectLifted(
  cache: LiftCache,
  version: number,
  model: StaticModel,
  data: RegionReader,
  from: string,
  ranges: ArrayLike<number>,
  to: string,
  out: number[],
): boolean {
  if (ranges.length !== 2) return imagesThrough(model, data, from, ranges, to, out);
  const s = ranges[0];
  const e = ranges[1];
  for (const k of cache.kept)
    if (k.version === version && k.model === model && k.view === undefined && k.from === from && k.to === to)
      if (k.start === s && k.end === e) {
        for (let i = 0; i < k.out.length; i++) out.push(k.out[i]);
        return k.ok;
      }
  const k = slot(cache);
  k.out.length = 0;
  k.ok = imagesThrough(model, data, from, ranges, to, k.out);
  k.version = version;
  k.model = model;
  k.view = undefined;
  k.from = from;
  k.to = to;
  k.start = s;
  k.end = e;
  for (let i = 0; i < k.out.length; i++) out.push(k.out[i]);
  return k.ok;
}

/** `liftPlaces`, kept for the store's version (0 B a warm hover's card). */
export function placesLifted(
  cache: LiftCache,
  version: number,
  model: StaticModel,
  data: RegionReader,
  view: ResolvedView,
  system: string,
  start: number,
  end: number,
): string | undefined {
  for (const k of cache.kept)
    if (k.version === version && k.model === model && k.view === view && k.from === system)
      if (k.start === start && k.end === end) return k.places;
  const k = slot(cache);
  k.places = liftPlaces(model, data, view, system, start, end);
  k.version = version;
  k.model = model;
  k.view = view;
  k.from = system;
  k.to = '';
  k.start = start;
  k.end = end;
  k.out.length = 0;
  return k.places;
}

// ---- a lifted layer's source ----

/** The annotation a layer's data reads: its address's, else the one a query starts at. */
function annotationOf(model: StaticModel, address: string) {
  const a = model.annotations.get(address);
  if (a) return a;
  for (const [at, entry] of model.annotations) if (address.startsWith(`${at}/`)) return entry;
  return undefined;
}

/** A plan's store source and system (an inline node's own address). */
const planOf = (p: AccessPlan) =>
  p.kind === 'inline'
    ? { source: p.source.address, system: p.system, shape: p.shape }
    : { source: p.source, system: p.system, shape: p.shape, ...(p.side ? { side: p.side } : {}) };

/**
 * A lifted layer's source (_spec/README.md "Render", A liftover): its rows on the section's system, read as the map's
 * runs there and the data along them on each far system the map's source names and the data's annotation reads.
 */
export function liftSource(model: StaticModel, data: NonNullable<ResolvedLayer['data']>): LayerSource | undefined {
  const lift = data.lift!;
  const m = model.annotations.get(lift.map);
  const a = annotationOf(model, data.address);
  if (!m?.bySystem || !a?.bySystem) return undefined;
  const filter = data.filter ? compileRowFilter(data.filter) : undefined;
  const base = data.base ? compileRowFilter(data.base) : undefined;
  const length = model.systems.get(data.system)?.length;
  const len = (n: number | undefined) => (n !== undefined ? { length: n } : {});
  const far = new Map<string, LayerSource>();
  for (const [system, plan] of a.bySystem) {
    const named = m.bySystem.get(system);
    if (!named) continue;
    const n = model.systems.get(system)?.length;
    far.set(named.system, {
      ...planOf(plan),
      dataSystem: system,
      ...len(n),
      ...(n !== undefined && plan.kind !== 'inline' ? { wraps: true } : {}),
      ratio: 1,
      ...(filter ? { filter } : {}),
      ...(base ? { base } : {}),
    });
  }
  const near: LayerSource = { ...planOf(lift.near), dataSystem: data.system, ...len(length), ratio: 1 };
  return {
    source: near.source,
    system: near.system,
    dataSystem: data.system,
    shape: 'ranges',
    ...len(length),
    ratio: 1,
    ...(filter ? { filter } : {}),
    ...(base ? { base } : {}),
    ...((data.filter && !filter) || (data.base && !base) || data.unread ? { unreadFilter: true } : {}),
    lift: { map: lift.map, side: lift.side, near, far },
  };
}
