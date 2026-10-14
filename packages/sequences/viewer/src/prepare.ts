/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the kernel renderers prepare (seqspec:viewer#extension-points; seqspec:sources#levels-and-the-zoom-gate):
// headless and pure, typed arrays over each
// block's coordinate system, segment tables in its indices and per-position arrays by `i − base`
// (seqspec:model#pieces-and-layout-maps), never shown positions or pixels, so a filter view, a zoom or a piece never
// re-prepares. Each placement's kept rows go through its path's map (`Coords`, runs in to order: a range takes their
// hull, or through a codon step from the protein each run as a fragment); through a codon step a per-position mark
// covers its residue's larger part. Prepare-time code allocates; the frame and the hit run only `featureAt`, `ofPartAt`
// and `positionKey`.

import { endsOf, letterMapTable, reach, splitRows } from "@molstar/seq-core/format";
import {
  type Alphabet,
  type Annotation,
  Annotations,
  type ComposedPath,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  type Numbering,
  Paths,
  type SeqSet,
  Systems,
} from "@molstar/seq-core/model";
import { intervalTree } from "./intervals.js";
import { predecessor, prefixReach, type ValuePyramid, valuePyramid } from "./pyramids.js";
import type { LayerInput, LayerRenderer, PrepareContext } from "./registry.js";
import { prepareMap, prepareMatrix, prepareRows, RowsHeader, rowsStream } from "./rows.js";
import { laneLineHeights } from "./runs.js";
import type { StreamLayer } from "./stream-contract.js";
import { lettersStream, trackStream } from "./stream-sources.js";
import type { FieldTable } from "./themes.js";

/** One placement's kept rows on its block, with its path. */
interface Placed {
  readonly block: number;
  readonly ann: Annotation;
  readonly slot: number;
  readonly rows: readonly number[];
  readonly path: ComposedPath | null;
  /** On a rows track, the row's lane offset (its entry draws it); -1: every row (column data). */
  readonly on: number;
}
/** Each placement's kept rows; on a rows track's lanes only where `lanes` (features: drawn on their row's entry). */
function placed(input: LayerInput, lanes = false): Placed[] {
  const { layer, set, section } = input;
  const out: Placed[] = [];
  const first = section.tracks[layer.track]?.laneStart ?? 0;
  for (const pl of layer.placements) {
    const ann = set.annotations[pl.annotation];
    const end = ann ? endsOf(ann).find((e) => e.slot === pl.slot) : undefined;
    if (!ann || !end || (pl.lane >= 0 && !lanes)) continue;
    const kept = layer.kept[layer.data.indexOf(ann.id)];
    const rows = kept ? end.rows.filter((r) => kept[r]) : end.rows;
    const path = pl.path >= 0 ? layer.paths[pl.path] : null;
    out.push({ block: pl.block, ann, slot: pl.slot, rows, path, on: pl.lane >= 0 ? pl.lane - first : -1 });
  }
  return out;
}

/**
 * The larger contiguous part [out[0], out[1]) of the cells a source position lands on: its own, its mapped one, or
 * its codon's three (a codon split by an intron: the larger part, the first on a tie). False: none.
 */
export function partOf(path: ComposedPath | null, i: number, out: Int32Array): boolean {
  out[0] = out[1] = -1;
  if (!path) {
    out[0] = i;
    out[1] = i + 1;
    return i >= 0;
  }
  if (!path.heldSpans) {
    const j = Paths.at(path, i);
    out[0] = j;
    out[1] = j + 1;
    return j >= 0;
  }
  let best = 0;
  let run = 0;
  let lo = 0;
  let hi = 0;
  let prev = -2;
  for (let k = 0; k < 3; k++) {
    const j = Paths.spanAt(path, i, k);
    if (j < 0) {
      prev = -2;
      continue;
    }
    if (run > 0 && Math.abs(j - prev) === 1) {
      run++;
      lo = Math.min(lo, j);
      hi = Math.max(hi, j + 1);
    } else {
      run = 1;
      lo = j;
      hi = j + 1;
    }
    prev = j;
    if (run > best) {
      best = run;
      out[0] = lo;
      out[1] = hi;
    }
  }
  return best > 0;
}

/** A ligand in the `position` form is one position, its domain's start: whatever lands on its atoms lands there. */
function ligandPart(input: LayerInput, b: number, part: Int32Array) {
  const block = input.section.blocks[b];
  if (!block.ligand || part[0] < 0) return;
  part[0] = block.base;
  part[1] = block.base + 1;
}

/** Every cell a source position lands on (a codon's three, split or not); prepare-time. */
function eachCell(path: ComposedPath | null, i: number, visit: (j: number) => void) {
  if (!path) visit(i);
  else
    for (let k = 0; k < path.span; k++) {
      const j = Paths.spanAt(path, i, k);
      if (j >= 0) visit(j);
    }
}

const runs: number[] = [];
const pushRun = (a: number, b: number, out: number[]) => {
  out.push(a, b);
};
/**
 * For [s, t) between elements the path pairs (one before s, one at t or after; its runs are sorted): the first paired
 * element at t or after, whose boundary is the insert's; -1: not between them.
 */
function pairedAfter(path: ComposedPath, s: number, t: number): number {
  const r = path.runs;
  if (r.length === 0 || r[0] >= s) return -1;
  for (let k = 0; k < r.length; k += 2) if (r[k] >= t) return r[k];
  return -1;
}
/** Source indices [s, t) (s === t: a boundary) through a path's map, as block index runs in `runs`, in to order. */
function runsOf(path: ComposedPath, s: number, t: number) {
  runs.length = 0;
  if (Coords.sameSystem(path.map.from, path.from)) Coords.forEachTo(path.map, s, t, pushRun, runs);
  else Coords.forEachFrom(path.map, s, s === t ? s + 1 : t, pushRun, runs);
}

/**
 * Features over block coordinate system indices, sorted by block, lane and start; `seg` indexes each (block, lane).
 * `LaneSegments` (runs.ts) is its lanes' part, `SplitFeatures` (layout.ts) its split rows'.
 */
export interface Features {
  readonly count: number;
  readonly block: Int32Array;
  /** [start, stop) on the block; start === stop: a boundary. */
  readonly start: Int32Array;
  readonly stop: Int32Array;
  /** The theme's location: the annotation and its row. */
  readonly annotation: Int32Array;
  readonly row: Int32Array;
  /** A type's own number (secondary's class, a variant's place in its stack, its `alt` letter). */
  readonly kind: Int32Array;
  readonly code: Int32Array;
  /** An item's fragments share its lane, so they may overlap there (an exon and its CDS: `featureAt` walks back). */
  readonly lane: Int32Array;
  readonly lanes: number;
  /** On a rows track: the row's lane offset whose entry draws it; -1: every row. */
  readonly on: Int32Array;
  /** Some feature is on a row (`on` ≥ 0): no density over the union of rows. */
  readonly hasRows: boolean;
  /** (block × lanes + lane) → its first feature; one more at the end. */
  readonly seg: Int32Array;
  /** Per feature, the largest stop in its segment up to it (`prefixReach`): a window is two searches (`segWindow`). */
  readonly reach: Int32Array;
  /**
   * The fragment of its item a hairline joins it from: the one reaching furthest before it, when every earlier
   * fragment of the item ends before its start (one hairline across a gap); -1: none.
   */
  readonly joinFrom: Int32Array;
  readonly label: readonly (string | null)[];
  /** A split layer's rows (seqspec:format#layer-types), null: none (one shape either way). */
  readonly split: SplitRows | null;
  /** Its streamed items' description (ST3a: the stream manager loads them by region); none: not streamed. */
  readonly stream?: StreamLayer | null;
  /** A streamed region's table: its items' fields, for field themes (`Theme.raw`); none: rows in the set. */
  readonly fields?: StreamFields | null;
}
/** A streamed region's items' fields (ST3b): its patch's columns, and its items' values of the last field theme drawn. */
export interface StreamFields {
  readonly table: FieldTable;
  theme: unknown;
  keys: Int32Array | null;
  /** The last glyph field theme drawn, and the items' glyphs by it (ST4b). */
  glyphTheme: unknown;
  glyphs: Int32Array | null;
}
/** A split layer's values in row order, each row's own features, and per data entry each row's value (-1: none). */
export interface SplitRows {
  readonly labels: readonly string[];
  readonly rows: readonly Features[];
  readonly rowOf: readonly Int32Array[];
  /** Per row, the unit every one of its data rows sits on; -1: none or several. */
  readonly unit: Int32Array;
}
/** Per split row, the one unit its data rows target (`SplitRows.unit`); prepare-time. */
function rowUnits(anns: readonly Annotation[], rowOf: readonly Int32Array[], n: number): Int32Array {
  const unit = new Int32Array(n).fill(-2);
  for (let d = 0; d < anns.length; d++)
    for (let r = 0; r < rowOf[d].length; r++) {
      const [v, u] = [rowOf[d][r], anns[d].slots.unit[anns[d].slot[r]]];
      if (v >= 0) unit[v] = unit[v] === -2 || unit[v] === u ? u : -1;
    }
  return unit.map((u) => Math.max(u, -1));
}
/**
 * The features whole (the compact form) and, with a `split` param, each value's prepared by `prep` alone, through a
 * filter of its rows; not on a rows track's lanes. Prepare-time.
 */
const withSplit =
  (prep: (input: LayerInput) => Features) =>
  (input: LayerInput, rc?: PrepareContext): Features => {
    const { layer, set } = input;
    const sp = input.params.split as { by: string; order?: string[] | null; empty?: "show" | "hide" } | null;
    const anns = layer.data.map((id) => set.annotations[id] as Annotation);
    const on =
      sp && !layer.placements.some((p) => p.lane >= 0) ? splitRows(set, anns, layer.kept, sp, rc?.registry.svs) : null;
    const rows = on?.values.map((_, r) => {
      const kept = on.rowOf.map((of) => Uint8Array.from(of, (x) => +(x === r)));
      return prep({ ...input, layer: { ...layer, kept } });
    });
    const unit = on ? rowUnits(anns, on.rowOf, on.values.length) : null;
    return Object.assign(prep(input), {
      split: on && rows && unit ? { labels: on.values, rows, rowOf: on.rowOf, unit } : null,
    });
  };
interface F {
  b: number;
  s: number;
  e: number;
  a: number;
  r: number;
  k: number;
  c: number;
  l: string | null;
  lane: number;
  on: number;
  /** An item across a cyclic origin (`arcs`): its arc's start and its tail's stop, [ws, L) then [0, wt); -1: none. */
  ws?: number;
  wt?: number;
}

/**
 * Ranges as runs (`runs`: a range across an exon junction in two fragments), or every position a row covers as one
 * (its residue's part through a codon step). `kind` and `code` per row; lanes packed per block when `pack`.
 */
export function features(
  input: LayerInput,
  o: {
    runs: boolean;
    pack: boolean;
    kind?: (a: Annotation, r: number) => number;
    code?: (a: Annotation, r: number) => number;
  },
): Features {
  const list: F[] = [];
  const part = new Int32Array(2);
  for (const p of placed(input, true)) {
    const a = p.ann;
    if (a.geometry !== "ranges") continue;
    for (const r of p.rows) {
      const lab = a.rowLabel && a.rowLabel[r] >= 0 ? input.set.strings.get(a.rowLabel[r]) : null;
      const base = {
        b: p.block,
        a: a.id,
        r,
        k: o.kind?.(a, r) ?? 0,
        c: o.code?.(a, r) ?? 0,
        l: lab,
        lane: 0,
        on: p.on,
      };
      const boundary = a.start[r] === a.stop[r];
      const push = (s: number, e: number) => {
        part[0] = s;
        part[1] = e;
        ligandPart(input, p.block, part);
        list.push({ ...base, s: part[0], e: boundary ? part[0] : part[1] });
      };
      // A boundary stays one (an insertion's mark between positions), whatever the type's mode.
      const whole = o.runs || boundary;
      const before = list.length;
      if (whole && !p.path) push(a.start[r], Math.max(a.stop[r], a.start[r]));
      else if (whole && p.path) {
        runsOf(p.path, a.start[r], a.stop[r]);
        if (p.path.heldSpans) for (let k = 0; k < runs.length; k += 2) push(runs[k], runs[k + 1]);
        else if (runs.length > 0) {
          let lo = runs[0];
          let hi = runs[1];
          for (let k = 2; k < runs.length; k += 2) [lo, hi] = [Math.min(lo, runs[k]), Math.max(hi, runs[k + 1])];
          push(lo, hi);
        }
      } else
        for (let i = a.start[r]; i < Math.max(a.stop[r], a.start[r] + 1); i++)
          if (partOf(p.path, i, part)) {
            ligandPart(input, p.block, part);
            list.push({ ...base, s: part[0], e: part[1] });
          }
      // An insert: a range wholly inside what the path leaves unmapped between paired elements (a design's grafted
      // loop on the wild type) lands on the boundary after it (seqspec:format#placement); outside the pairs, nothing.
      const next = list.length === before && !boundary && p.path ? pairedAfter(p.path, a.start[r], a.stop[r]) : -1;
      if (next >= 0 && p.path) {
        runsOf(p.path, next, next);
        if (runs.length > 0) list.push({ ...base, s: runs[0], e: runs[0] });
      }
    }
  }
  // Greedy lanes per block (and row), by item over its whole extent (a transcript's exons on one lane): an item takes
  // the first lane free at its start (boundaries one position wide); its fragments take its lane.
  const items = new Map<string, F>();
  const itemOf = (f: F) => `${f.b}:${f.on}:${f.a}:${input.set.annotations[f.a]?.item[f.r] ?? f.r}`;
  for (const f of list) {
    const x = items.get(itemOf(f));
    if (!x) items.set(itemOf(f), { ...f, ws: -1, wt: -1 });
    else [x.s, x.e] = [Math.min(x.s, f.s), Math.max(x.e, f.e)];
  }
  arcs(input, list, items, itemOf);
  // An item across a cyclic origin packs by its arc: first, a lane each, its tail [0, wt) taken and the lane free only
  // up to ws; the others then as on a line.
  const wraps = (x: F) => ((x.ws as number) >= 0 ? 0 : 1);
  const spans = [...items.values()].sort(
    (x, y) => x.b - y.b || x.on - y.on || wraps(x) - wraps(y) || x.s - y.s || y.e - x.e,
  );
  let lanes = 1;
  const ends: number[] = [];
  const limit: number[] = [];
  for (let k = 0; o.pack && k < spans.length; k++) {
    const f = spans[k];
    if (k === 0 || spans[k - 1].b !== f.b || spans[k - 1].on !== f.on) ends.length = limit.length = 0;
    const wrap = (f.ws as number) >= 0;
    const [s, e] = wrap ? [0, f.wt as number] : [f.s, Math.max(f.e, f.s + 1)];
    let l = 0;
    while (l < ends.length && (wrap || ends[l] > s || e > limit[l])) l++;
    ends[l] = e;
    limit[l] = wrap ? (f.ws as number) : Number.POSITIVE_INFINITY;
    f.lane = l;
    lanes = Math.max(lanes, l + 1);
  }
  for (const f of list) f.lane = (items.get(itemOf(f)) as F).lane;
  list.sort((x, y) => x.b - y.b || x.lane - y.lane || x.s - y.s);
  const n = list.length;
  const I = () => new Int32Array(n);
  const out = {
    ...{ count: n, block: I(), start: I(), stop: I(), annotation: I(), row: I(), kind: I(), code: I(), lane: I() },
    on: Int32Array.from(list, (f) => f.on),
    lanes,
    hasRows: list.some((f) => f.on >= 0),
    seg: new Int32Array(input.section.blocks.length * lanes + 1),
    reach: new Int32Array(n),
    joinFrom: new Int32Array(n).fill(-1),
    label: list.map((f) => f.l),
    split: null,
  };
  list.forEach((f, k) => {
    out.block[k] = f.b;
    out.start[k] = f.s;
    out.stop[k] = f.e;
    out.annotation[k] = f.a;
    out.row[k] = f.r;
    out.kind[k] = f.k;
    out.code[k] = f.c;
    out.lane[k] = f.lane;
  });
  for (let g = 0, k = 0; g < out.seg.length; g++) {
    while (k < n && out.block[k] * lanes + out.lane[k] < g) k++;
    out.seg[g] = k;
  }
  for (let g = 0; g + 1 < out.seg.length; g++) prefixReach(out.start, out.stop, out.seg[g], out.seg[g + 1], out.reach);
  // Per item (block, row, annotation, item) its fragment reaching furthest so far: a gap after it is one hairline.
  // An item across a cyclic origin is one arc: the gap outside it, [wt, ws), is never joined.
  const furthest = new Map<string, number>();
  for (let k = 0; o.pack && k < n; k++) {
    const key = itemOf(list[k]);
    const j = furthest.get(key);
    const x = items.get(key) as F;
    const outside = j !== undefined && out.stop[j] === x.wt && out.start[k] === x.ws;
    if (j !== undefined && out.stop[j] < out.start[k] && !outside) out.joinFrom[k] = j;
    if (j === undefined || out.stop[k] >= out.stop[j]) furthest.set(key, k);
  }
  return out;
}

/**
 * Items across the origin of a whole cyclic unit (seqspec:model#cyclic-systems: their fragments reach both ends) get
 * their arc
 * (`ws`, `wt`): the circle less the widest gap between their fragments (`Coords.cyclicHull`). Prepare-time.
 */
function arcs(input: LayerInput, list: readonly F[], items: Map<string, F>, itemOf: (f: F) => string) {
  const runs = new Map<string, number[]>();
  for (const f of list) {
    const x = items.get(itemOf(f)) as F;
    const blk = input.section.blocks[x.b];
    const mol = blk.alignment < 0 ? input.set.molecules[input.set.units.molecule[blk.unit]] : null;
    const whole = mol?.cyclic && blk.base === 0 && blk.length === mol.length;
    if (whole && x.s === 0 && x.e === mol.length)
      runs.get(itemOf(f))?.push(f.s, f.e) ?? runs.set(itemOf(f), [f.s, f.e]);
  }
  const out = new Int32Array(4);
  for (const [key, r] of runs) {
    const x = items.get(key) as F;
    const order = Array.from({ length: r.length / 2 }, (_, k) => k).sort((p, q) => r[2 * p] - r[2 * q]);
    const sorted = Int32Array.from(order.flatMap((k) => [r[2 * k], r[2 * k + 1]]));
    const L = input.section.blocks[x.b].length;
    if (Coords.cyclicHull(sorted, order.length, L, out) < 2) continue;
    x.ws = out[0];
    x.wt = out[3];
  }
}

/**
 * The feature of a block's lane covering index i, the one drawn on top (the last in its segment; boundaries: at it),
 * walking back while the lane's `reach` passes i, so an earlier, longer fragment under a later one is found;
 * -1: none. Logarithmic plus the overlap at i, allocation-free.
 */
export function featureAt(f: Features, block: number, lane: number, i: number): number {
  if (lane < 0 || lane >= f.lanes) return -1;
  const g = block * f.lanes + lane;
  const lo = f.seg[g];
  for (let k = predecessor(f.start, lo, f.seg[g + 1], i); k >= lo && f.reach[k] > i; k--)
    if (f.stop[k] > i || (f.stop[k] === f.start[k] && f.start[k] === i)) return k;
  return -1;
}

/**
 * Per block, by `i − base`: a value (NaN: none) and the annotation row that set it; the last covering row wins. A
 * block with one dense placement over its own positions shares its annotation's array (a `subarray`, no copy);
 * `annotation[b]` is null where one annotation covers the block (`only[b]`), `row[b]` null where no row does (-1) and
 * for `values` (nothing reads it), `value[b]` null where the layer's data has no values (a `color` layer over ranges:
 * its rows say what it covers), so nothing per position is held that the data does not have.
 */
export interface PerPosition {
  readonly value: (Float64Array | null)[];
  readonly annotation: (Int32Array | null)[];
  readonly row: (Int32Array | null)[];
  /** Per block, its one annotation where `annotation[b]` is null; -1: none. */
  readonly only: Int32Array;
  /** The domain drawn: the layer's (its `domain`, else over its data), from zero when the data are not negative. */
  readonly lo: number;
  readonly hi: number;
  /** Per block its values' pyramid: `values` min, max, sum and count; `color` sums and counts. */
  readonly pyramid: readonly (ValuePyramid | null)[];
  /** `values`: its streamed data's description (the stream manager loads it); null: none streamed. */
  readonly stream?: StreamLayer | null;
}
/** `values`: per position, with each block's pyramid; streamed blocks by region (ST1: `stream`). */
const perPositionValues = (input: LayerInput): PerPosition => ({
  ...perPosition(input, false, true),
  stream: trackStream(input),
});
export function perPosition(input: LayerInput, rows = true, minmax = false): PerPosition {
  const blocks = input.section.blocks;
  const value: (Float64Array | null)[] = blocks.map(() => null);
  const annotation: (Int32Array | null)[] = blocks.map(() => null);
  const row: (Int32Array | null)[] = blocks.map(() => null);
  const list = placed(input);
  const valued = list.some((p) => p.ann.sparse || p.ann.dense?.[p.slot]);
  const only = new Int32Array(blocks.length).fill(-1);
  const count = new Int32Array(blocks.length);
  for (const p of list) {
    only[p.block] = only[p.block] === -1 || only[p.block] === p.ann.id ? p.ann.id : -2;
    count[p.block]++;
  }
  for (const p of list) {
    const b = p.block;
    const { length: n, base } = blocks[b];
    const a = p.ann;
    const dense = a.dense?.[p.slot];
    if (dense && !p.path && count[b] === 1 && !blocks[b].ligand && dense.length >= base + n) {
      value[b] = dense.subarray(base, base + n);
      continue;
    }
    if (valued && !value[b]) value[b] = new Float64Array(n).fill(Number.NaN);
    if (only[b] === -2 && !annotation[b]) annotation[b] = new Int32Array(n).fill(-1);
    const lig = blocks[b].ligand === 1;
    const set = (i: number, v: number, r: number) => {
      const j = lig && i >= 0 ? 0 : i - base;
      if (i < 0 || j < 0 || j >= n) return;
      const vb = value[b];
      if (vb) vb[j] = v;
      const at = annotation[b];
      if (at) at[j] = a.id;
      if (rows && r >= 0 && !row[b]) row[b] = new Int32Array(n).fill(-1);
      const rb = row[b];
      if (rb) rb[j] = r;
    };
    if (dense) for (let i = 0; i < dense.length; i++) eachCell(p.path, i, (j) => set(j, dense[i], -1));
    else
      for (const r of p.rows) {
        const v = a.sparse ? a.sparse[r] : Number.NaN;
        if (a.geometry === "values")
          for (let i = a.start[r]; i < Math.max(a.stop[r], a.start[r] + 1); i++)
            eachCell(p.path, i, (j) => set(j, v, r));
        else if (!p.path) for (let j = a.start[r]; j < Math.max(a.stop[r], a.start[r] + 1); j++) set(j, v, r);
        else {
          runsOf(p.path, a.start[r], a.stop[r]);
          for (let k = 0; k < runs.length; k += 2) for (let j = runs[k]; j < runs[k + 1]; j++) set(j, v, r);
        }
      }
  }
  const l = input.layer;
  const given = l.domainGiven;
  const lo = Number.isFinite(l.domainMin) ? (given ? l.domainMin : Math.min(0, l.domainMin)) : 0;
  const hi = Number.isFinite(l.domainMax) ? l.domainMax : 1;
  const pyramid = value.map((v) => (v ? valuePyramid(v, minmax) : null));
  return { value, annotation, row, only, lo, hi: hi > lo ? hi : lo + 1, pyramid };
}

/**
 * A cell's key for its class (reductions.ts `cellClasses`), at `j = i − base`: where the theme colors per row or
 * coarser (`byRow`), the row's annotation and row, so a row's run takes one theme call; -1 where no row covers it, or
 * the theme varies finer (a key never stands for a value: a dense layer is colored per position); 0: nothing.
 */
export function positionKey(s: PerPosition, b: number, j: number, byRow: boolean): number {
  const rows = s.row[b];
  const values = s.value[b];
  const r = rows ? rows[j] : -1;
  if (r < 0 && (!values || Number.isNaN(values[j]))) return 0;
  if (!byRow || r < 0) return -1;
  const ann = s.annotation[b];
  return ((ann ? ann[j] : s.only[b]) + 1) * 16777216 + r + 2;
}
/**
 * A `mean` theme over a row field: block b's values per position (`i − base`), each covering row's through
 * `read(annotation, row)`, NaN where none; the caller builds its sum-and-count pyramid once per theme.
 */
export function rowValues(s: PerPosition, b: number, read: (annotation: number, row: number) => number): Float64Array {
  const rows = s.row[b];
  const ann = s.annotation[b];
  const out = new Float64Array(rows ? rows.length : 0).fill(Number.NaN);
  for (let j = 0; rows && j < rows.length; j++) if (rows[j] >= 0) out[j] = read(ann ? ann[j] : s.only[b], rows[j]);
  return out;
}

export const LetterForm = { residues: 0, ligand: 1, atoms: 2 } as const;
/**
 * Other polymers' residues on the blocks they reach (`of`: one, a list, or every codon mapping's), as a segment table
 * of parts: each residue's larger part (`main`, where its letter is drawn) and the cells of a codon an intron splits
 * off, sorted by block and start (stable: placement order where polymers overlap), with `seg` per block and `reach`
 * (`segWindow`). Arrays per part only: nothing the length of a block.
 */
export interface OfLetters {
  readonly count: number;
  readonly block: Int32Array;
  readonly start: Int32Array;
  readonly stop: Int32Array;
  readonly code: Uint8Array;
  /** The residue's index in its own sequence, and its polymer's unit (a sequence: its first; -1: none). */
  readonly source: Int32Array;
  readonly unit: Int32Array;
  readonly main: Uint8Array;
  /** Per block its first part; one more at the end. */
  readonly seg: Int32Array;
  readonly reach: Int32Array;
  /** The last lookup's block and first part reaching past it: this layer's own. */
  readonly cursor: Int32Array;
}
/**
 * The first part of block b over index i (the one drawn), -1: none. From the layer's cursor, galloping forward
 * (amortized constant for increasing i, logarithmic for a jump), else a binary search.
 */
export function ofPartAt(o: OfLetters, b: number, i: number): number {
  const lo = o.seg[b];
  const hi = o.seg[b + 1];
  const c = o.cursor;
  let k = c[1];
  if (c[0] !== b || k < lo || k > hi || (k > lo && o.reach[k - 1] > i)) k = predecessor(o.reach, lo, hi, i) + 1;
  else if (k < hi && o.reach[k] <= i) {
    let step = 1;
    while (k + step < hi && o.reach[k + step] <= i) {
      k += step;
      step *= 2;
    }
    k = predecessor(o.reach, k, Math.min(hi, k + step), i) + 1;
  }
  c[0] = b;
  c[1] = k;
  for (let j = k; j < hi && o.start[j] <= i; j++) if (o.stop[j] > i) return j;
  return -1;
}
export interface LettersState {
  /** Per block, by `i − base`: ASCII residue codes, an alignment row's (gaps included), or atomic numbers. */
  readonly codes: readonly (Uint8Array | null)[];
  /** `LetterForm` per block. */
  readonly form: Int8Array;
  /** A ligand position's label (its chemical component). */
  readonly names: readonly (string | null)[];
  readonly unit: Int32Array;
  readonly alignment: Int32Array;
  /** `row`: an alignment row's letters; -1: the block's. */
  readonly row: number;
  readonly of: OfLetters | null;
  /** `letter_map`: the code shown per stored code (128 entries); null: as stored. Never for atoms. */
  readonly map: Uint8Array | null;
  /** Per block, its letters' alphabet (the letters LOD's widest letter); null: atoms or none. */
  readonly alphabet: readonly (Alphabet | null)[];
  /** Blocks over a lazy molecule (null codes): their letters by region (ST3a); null: none streamed. */
  readonly stream: StreamLayer | null;
}

function prepareLetters(input: LayerInput): LettersState {
  const { set, section, params } = input;
  const n = section.blocks.length;
  const form = new Int8Array(n);
  const names: (string | null)[] = [];
  const unit = new Int32Array(n);
  const alignment = new Int32Array(n);
  const row = typeof params.row === "number" ? params.row : -1;
  const alphabet: (Alphabet | null)[] = [];
  const codes = section.blocks.map((b, i) => {
    unit[i] = b.unit;
    alignment[i] = b.alignment;
    names.push(null);
    const aln = set.alignments[b.alignment];
    alphabet.push(
      aln ? aln.alphabet : b.unit >= 0 ? (set.molecules[set.units.molecule[b.unit]]?.alphabet ?? null) : null,
    );
    if (aln) {
      const s = aln.lazy ? -1 : aln.rows.slotOf((row >= 0 ? row : Math.max(0, aln.referenceRow)) as never);
      return s < 0 ? null : aln.rows.codes.subarray(s * aln.columnCount, (s + 1) * aln.columnCount);
    }
    const mol = b.unit >= 0 ? set.molecules[set.units.molecule[b.unit]] : null;
    if (mol?.element && b.ligand) {
      form[i] = LetterForm.ligand;
      names[i] = mol.chemCompId ? set.strings.get(mol.chemCompId[0]) : "";
    } else if (mol?.element) form[i] = LetterForm.atoms;
    if (mol?.lazy) return null;
    const all = mol ? (mol.element ?? mol.code) : null;
    return all && !b.ligand ? all.subarray(b.base, b.base + b.length) : all;
  });
  const of = params.of != null ? ofLetters(input) : null;
  const map = letterMapTable(params.letter_map as never);
  return { codes, form, names, unit, alignment, row, of, map, alphabet, stream: of ? null : lettersStream(input) };
}

function ofLetters(input: LayerInput): OfLetters {
  const { set, section, layer } = input;
  const K = CoordinateSystemKind;
  const list: { b: number; s: number; e: number; c: number; i: number; u: number; m: number }[] = [];
  const part = new Int32Array(2);
  const first = new Map<number, number>();
  for (const pl of layer.placements) {
    const path = pl.path >= 0 ? layer.paths[pl.path] : null;
    const from = path?.from ?? { kind: K.unit, id: section.blocks[pl.block].unit, row: -1 };
    const mol = set.molecules[from.kind === K.unit ? set.units.molecule[from.id] : from.id];
    let unit = from.kind === K.unit ? from.id : (first.get(from.id) ?? -1);
    for (let u = 0; from.kind === K.molecule && unit < 0 && u < set.units.count; u++)
      if (set.units.alive[u] && set.units.molecule[u] === from.id) unit = u;
    if (from.kind === K.molecule) first.set(from.id, unit);
    for (let i = 0; mol && i < mol.length; i++) {
      if (!partOf(path, i, part)) continue;
      const r = { b: pl.block, c: mol.code[i], i, u: unit };
      list.push({ ...r, s: part[0], e: part[1], m: 1 });
      // A codon an intron splits: its other cells, each a part of its own.
      eachCell(path, i, (j) => {
        if (j < part[0] || j >= part[1]) list.push({ ...r, s: j, e: j + 1, m: 0 });
      });
    }
  }
  list.sort((x, y) => x.b - y.b || x.s - y.s);
  const n = list.length;
  const seg = new Int32Array(section.blocks.length + 1);
  for (const x of list) seg[x.b + 1]++;
  for (let b = 0; b < section.blocks.length; b++) seg[b + 1] += seg[b];
  const o = {
    count: n,
    block: Int32Array.from(list, (x) => x.b),
    start: Int32Array.from(list, (x) => x.s),
    stop: Int32Array.from(list, (x) => x.e),
    code: Uint8Array.from(list, (x) => x.c),
    source: Int32Array.from(list, (x) => x.i),
    unit: Int32Array.from(list, (x) => x.u),
    main: Uint8Array.from(list, (x) => x.m),
    seg,
    reach: new Int32Array(n),
    cursor: Int32Array.of(-1, 0),
  };
  for (let b = 0; b < section.blocks.length; b++) prefixReach(o.start, o.stop, seg[b], seg[b + 1], o.reach);
  return o;
}

/**
 * Where a ruler's labels come from per block, in order, the first that labels a position winning: with `of`, each
 * placement's coordinate system through its path (a list `of`: prefixed with its label); else the block's own system in
 * its `numbering` (a list: the first that labels it). Null: the section's numbering, read live by the draw, as is a
 * named numbering the block does not have.
 */
export interface RulerSource {
  /** From the source's coordinate system to the block's; null: the block's own. */
  readonly path: ComposedPath | null;
  /** The first that labels a position wins; none: ordinals. */
  readonly numberings: readonly Numbering[];
  readonly prefix: string;
}
export interface RulerState {
  readonly sources: readonly (readonly RulerSource[] | null)[];
  readonly names: readonly string[];
}
/** A coordinate system's numberings by name, in the order asked; none named: its default (none: ordinals). */
const numberingsOf = (set: SeqSet, s: CoordinateSystem, named: readonly string[]): Numbering[] =>
  Systems.numberingsOf(set, s.kind, s.id, s.row, named);
function prepareRuler(input: LayerInput): RulerState {
  const { set, section, params, layer } = input;
  const own = Array.isArray(params.numbering) ? params.numbering : params.numbering != null ? [params.numbering] : [];
  const named = own.filter((n): n is string => typeof n === "string");
  const list = Array.isArray(params.of);
  const sources = section.blocks.map((b, i): RulerSource[] | null => {
    if (params.of == null) {
      const mine =
        b.unit >= 0 ? numberingsOf(set, { kind: CoordinateSystemKind.unit, id: b.unit, row: -1 }, named) : [];
      return named.length > 0 && mine.length > 0 ? [{ path: null, numberings: mine, prefix: "" }] : null;
    }
    const out: RulerSource[] = [];
    for (const pl of layer.placements) {
      if (pl.block !== i) continue;
      const path = pl.path >= 0 ? layer.paths[pl.path] : null;
      const from = path ? path.from : { kind: CoordinateSystemKind.unit, id: b.unit, row: -1 };
      const label = from.kind === CoordinateSystemKind.unit ? set.units.label[from.id] : set.molecules[from.id]?.label;
      const prefix = list && label !== undefined && label >= 0 ? set.strings.get(label) : "";
      out.push({ path, numberings: numberingsOf(set, from, named), prefix });
    }
    return out;
  });
  const names = section.blocks.map((b) => (b.label >= 0 ? set.strings.get(b.label) : ""));
  return { sources, names };
}

/** A position's key on the axis: block-major, so keys order as the axis does (side by side). Exact in a double. */
export const keyOf = (block: number, index: number) => block * 4294967296 + index;

/**
 * Pairs with both ends placed, as parts on their blocks, and those with one end in a block the section's `positions`
 * filter left out (`leftOut`): that end is block -1, the other's glyph drawn and its arc a tick toward the partner
 * (`pushArc`), as for an end filtered out within its block. An interval index over their keys (arcs' hits).
 */
export interface PairsState {
  readonly set: LayerInput["set"];
  readonly count: number;
  /** Per pair, per end (a, b): block, part start, stop; -1, -1, -1 for an end in a block left out. */
  readonly block: Int32Array;
  readonly start: Int32Array;
  readonly stop: Int32Array;
  readonly annotation: Int32Array;
  readonly row: Int32Array;
  /** Pairs by their first key; per pair its last key; the interval tree over the last keys in that order. */
  readonly order: Int32Array;
  readonly lo: Float64Array;
  readonly hi: Float64Array;
  readonly tree: Float64Array;
  /**
   * The placed ends by key: their pair × 2 + side; the longest end part (an end starting that far before a run may
   * meet it).
   */
  readonly endKey: Float64Array;
  readonly endOf: Int32Array;
  readonly maxEnd: number;
  /**
   * The kinds of contact at a cell (several ends on one block, start and stop): a pair's kind is its `category` value
   * interned per state (0 when there is none). `endAt` is an end's place in `endOf` (-1: not placed), `endCell` the
   * cell of each place; a cell's places are `cellLo[c]` to `cellLo[c + 1]`, its distinct kinds `cellKinds[c]` of them
   * from `cellKindAt[c]` in `kindList` (ascending), each with the place of its first end in `kindRep`.
   */
  readonly kind: Int32Array;
  readonly endAt: Int32Array;
  readonly endCell: Int32Array;
  readonly cellLo: Int32Array;
  readonly cellKinds: Int32Array;
  readonly cellKindAt: Int32Array;
  readonly kindList: Int32Array;
  readonly kindRep: Int32Array;
  /**
   * Per pair its place among the overlay's last arcs, from 1 (0: not drawn there; the last on top): with
   * `arcs: "hover"` only they hit.
   */
  readonly drawn: Int32Array;
  /** The geometry in the last view (Arcs.frame): base, direction, max height, widest, and its key. */
  readonly frame: Float64Array;
  /** Two ends' x (Arcs.ends), then their virtual x, the lower first (Arcs.virtual). */
  readonly xs: Float64Array;
  /** The axis the frame was computed on. */
  readonly at: { axis: unknown };
}
const endPart = new Int32Array(2);
/** An end's parts on each placement where its position lands, as (block, start, stop) triples; prepare-time. */
function parts(input: LayerInput, at: readonly Placed[] | undefined, i: number): number[] {
  const out: number[] = [];
  for (const P of at ?? []) {
    if (!partOf(P.path, i, endPart)) continue;
    ligandPart(input, P.block, endPart);
    out.push(P.block, endPart[0], endPart[1]);
  }
  return out;
}
let reached = false;
const onReach = () => {
  reached = true;
};
/** Whether end `slot` of annotation `id` reaches a block the section's `positions` filter left out; prepare-time. */
function inLeftOut(input: LayerInput, id: number, slot: number): boolean {
  const a = input.set.annotations[id];
  const end = a ? endsOf(a).find((e) => e.slot === slot) : undefined;
  reached = false;
  if (end) for (const b of input.section.leftOut) if (!reached) reach(input.set, end.system, b, onReach);
  return reached;
}
/**
 * The pair kinds and the cells they share (see `PairsState`): ends are in key order, so the ends of one block, start
 * and stop are a run. Prepare-time.
 */
function pairCells(
  set: LayerInput["set"],
  rows: readonly number[][],
  s: { block: Int32Array; start: Int32Array; stop: Int32Array },
  ends: Int32Array,
) {
  const names = new Map<string, number>();
  const kind = new Int32Array(rows.length);
  rows.forEach((x, k) => {
    const a = set.annotations[x[6]];
    const v = a ? Annotations.rowValue(a, set, x[7], "category") : undefined;
    const name = v == null ? "" : String(v);
    if (!names.has(name)) names.set(name, names.size);
    kind[k] = names.get(name) as number;
  });
  const endAt = new Int32Array(2 * rows.length).fill(-1);
  ends.forEach((e, q) => {
    endAt[e] = q;
  });
  const endCell = new Int32Array(ends.length);
  const cellLo: number[] = [];
  const cellKinds: number[] = [];
  const cellKindAt: number[] = [];
  const kindList: number[] = [];
  const kindRep: number[] = [];
  const same = (a: number, b: number) =>
    s.block[a] === s.block[b] && s.start[a] === s.start[b] && s.stop[a] === s.stop[b];
  for (let q = 0; q < ends.length; ) {
    let r = q + 1;
    while (r < ends.length && same(ends[r], ends[q])) r++;
    const first = new Map<number, number>();
    for (let i = q; i < r; i++) {
      endCell[i] = cellLo.length;
      const kd = kind[ends[i] >> 1];
      if (!first.has(kd)) first.set(kd, i);
    }
    cellLo.push(q);
    cellKindAt.push(kindList.length);
    cellKinds.push(first.size);
    for (const kd of [...first.keys()].sort((a, b) => a - b)) {
      kindList.push(kd);
      kindRep.push(first.get(kd) as number);
    }
    q = r;
  }
  cellLo.push(ends.length);
  return {
    kind,
    endAt,
    endCell,
    cellLo: Int32Array.from(cellLo),
    cellKinds: Int32Array.from(cellKinds),
    cellKindAt: Int32Array.from(cellKindAt),
    kindList: Int32Array.from(kindList),
    kindRep: Int32Array.from(kindRep),
  };
}
function preparePairs(input: LayerInput): PairsState {
  const byEnd = new Map<string, Placed[]>();
  for (const p of placed(input)) {
    const k = `${p.ann.id}:${p.slot}`;
    byEnd.set(k, [...(byEnd.get(k) ?? []), p]);
  }
  const rows: number[][] = [];
  // Per end (annotation:slot), whether it reaches a block the section's filter left out.
  const outs = new Map<string, boolean>();
  const leftOut = (id: number, slot: number) => {
    if (input.section.leftOut.length === 0) return false;
    const k = `${id}:${slot}`;
    if (!outs.has(k)) outs.set(k, inLeftOut(input, id, slot));
    return outs.get(k) as boolean;
  };
  for (const id of input.layer.data) {
    const a = input.set.annotations[id];
    if (!a?.pairs) continue;
    const bSlot = new Int32Array(a.count).fill(-1);
    for (const e of endsOf(a)) if (e.side === 1) for (const r of e.rows) bSlot[r] = e.slot;
    const kept = input.layer.kept[input.layer.data.indexOf(id)];
    for (let r = 0; r < a.count; r++) {
      if (kept && !kept[r]) continue;
      const As = parts(input, byEnd.get(`${id}:${a.slot[r]}`), a.start[r]);
      const Bs = parts(input, byEnd.get(`${id}:${bSlot[r]}`), a.pairs.bStart[r]);
      // An end in a block the filter left out keeps its partner's glyph (a pocket's filter leaves its ligand out).
      if (Bs.length === 0 && leftOut(id, bSlot[r]))
        for (let i = 0; i < As.length; i += 3) rows.push([As[i], As[i + 1], As[i + 2], -1, -1, -1, id, r]);
      else if (As.length === 0 && leftOut(id, a.slot[r]))
        for (let j = 0; j < Bs.length; j += 3) rows.push([-1, -1, -1, Bs[j], Bs[j + 1], Bs[j + 2], id, r]);
      for (let i = 0; i < As.length; i += 3)
        for (let j = 0; j < Bs.length; j += 3)
          rows.push([As[i], As[i + 1], As[i + 2], Bs[j], Bs[j + 1], Bs[j + 2], id, r]);
    }
  }
  const n = rows.length;
  const I = (m: number) => new Int32Array(m);
  const s = { block: I(2 * n), start: I(2 * n), stop: I(2 * n), annotation: I(n), row: I(n) };
  const lo = new Float64Array(n);
  const hi = new Float64Array(n);
  rows.forEach((x, k) => {
    for (let e = 0; e < 2; e++) {
      s.block[2 * k + e] = x[3 * e];
      s.start[2 * k + e] = x[3 * e + 1];
      s.stop[2 * k + e] = x[3 * e + 2];
    }
    s.annotation[k] = x[6];
    s.row[k] = x[7];
    // A half-placed pair's interval is its placed end.
    const ka = x[0] < 0 ? keyOf(x[3], x[4]) : keyOf(x[0], x[1]);
    const kb = x[3] < 0 ? ka : keyOf(x[3], x[4]);
    lo[k] = Math.min(ka, kb);
    hi[k] = Math.max(ka, kb);
  });
  const order = Int32Array.from(lo.keys()).sort((x, y) => lo[x] - lo[y]);
  const ends = Int32Array.from({ length: 2 * n }, (_, e) => e)
    .filter((e) => s.block[e] >= 0)
    .sort((x, y) => keyOf(s.block[x], s.start[x]) - keyOf(s.block[y], s.start[y]));
  return {
    set: input.set,
    count: n,
    ...s,
    ...pairCells(input.set, rows, s, ends),
    order,
    lo: Float64Array.from(order, (k) => lo[k]),
    hi,
    tree: intervalTree(Float64Array.from(order, (k) => hi[k])),
    endKey: Float64Array.from(ends, (e) => keyOf(s.block[e], s.start[e])),
    endOf: ends,
    maxEnd: ends.reduce((m, e) => Math.max(m, s.stop[e] - s.start[e]), 1),
    drawn: new Int32Array(n),
    frame: new Float64Array(12).fill(Number.NaN),
    xs: new Float64Array(4),
    at: { axis: null },
  };
}

/**
 * A row's `alt` as its first character's code (a variant's substitution); a boundary's (an insertion's inserted
 * letters) as minus its length; 0: none.
 */
function firstCode(a: Annotation, set: LayerInput["set"], r: number): number {
  const id = a.alt ? a.alt[r] : -1;
  const v = id >= 0 ? set.strings.get(id) : Annotations.rowValue(a, set, r, "alt");
  if (typeof v !== "string" || v.length === 0) return 0;
  return a.start[r] === a.stop[r] && /^[A-Za-z]+$/.test(v) ? -v.length : v.charCodeAt(0);
}

const SECONDARY: Readonly<Record<string, number>> = {
  helix: 0,
  H: 0,
  G: 0,
  I: 0,
  strand: 1,
  E: 1,
  B: 1,
  turn: 2,
  T: 2,
  S: 2,
};
/** Secondary structure classes: helix, strand, turn; anything else is drawn as a range. */
export const SecondaryKind = { helix: 0, strand: 1, turn: 2, other: 3 } as const;

// `display: "union"`: every item in one lane, drawn over each other (their union, as coverage).
const ranges = (input: LayerInput) => features(input, { runs: true, pack: input.params.display !== "union" });
const secondary = (input: LayerInput) =>
  features(input, {
    runs: true,
    pack: false,
    kind: (a, r) =>
      SECONDARY[a.category && a.category[r] >= 0 ? input.set.strings.get(a.category[r]) : ""] ?? SecondaryKind.other,
  });
const positions = (input: LayerInput) => features(input, { runs: false, pack: false });
/** Variants whole (a deletion over several positions is one), each at its place in its position's stack. */
const variants = (input: LayerInput) => {
  const f = features(input, { runs: true, pack: false, code: (a, r) => firstCode(a, input.set, r) });
  // Boundaries (insertions, between positions) stack apart from the position's own variants.
  let at = 0;
  let marks = f.count > 0 && f.stop[0] <= f.start[0] ? 1 : 0;
  for (let k = 1; k < f.count; k++) {
    if (f.block[k] !== f.block[at] || f.start[k] !== f.start[at] || f.lane[k] !== f.lane[at]) {
      at = k;
      marks = f.stop[k] <= f.start[k] ? 1 : 0;
      continue;
    }
    const mark = f.stop[k] <= f.start[k];
    f.kind[k] = mark ? marks++ : k - at - marks;
  }
  return f;
};
/** Emphasis: its positions, and the block's letters to draw them with (through its track's `letter_map`, as resolved). */
export interface EmphasisState {
  readonly f: Features;
  readonly letters: LettersState;
}
const emphasis = (input: LayerInput): EmphasisState => ({
  f: positions(input),
  letters: prepareLetters({ ...input, params: { letter_map: input.params.letter_map } }),
});

/** A feature layer's streamed items (ST3a): their description beside the features in memory. */
const withStream =
  (prep: (input: LayerInput, rc?: PrepareContext) => Features) =>
  (input: LayerInput, rc?: PrepareContext): Features => {
    const stream = trackStream(input);
    return stream ? { ...prep(input, rc), stream } : prep(input, rc);
  };
const featureStream = (f: Features) => f.stream ?? null;

/** A kernel type's prepare and what goes with its state: whether it reads content, wrapped bands, the tree header. */
export type PreparedPart = Pick<LayerRenderer, "content" | "lineHeights" | "header" | "stream"> &
  Required<Pick<LayerRenderer, "prepare">>;
const lanes = { content: false, lineHeights: laneLineHeights };
/**
 * Every kernel type's prepare (LY3's entries spread them): `ranges`, `points` and `secondary` banded in wrapped lines
 * by their lanes, `rows` with the tree's header and its stream; only `letters`, `emphasis` and `rows` read residues
 * (seqspec:viewer#updates-and-status).
 */
export const PreparedParts: Readonly<Record<string, PreparedPart>> = {
  letters: { prepare: prepareLetters, stream: (s: LettersState) => s.stream },
  ranges: { prepare: withStream(withSplit(ranges)), stream: featureStream, ...lanes },
  points: { prepare: withStream(withSplit(positions)), stream: featureStream, ...lanes },
  emphasis: { prepare: emphasis },
  variants: { prepare: withStream(withSplit(variants)), content: false, stream: featureStream },
  secondary: { prepare: secondary, ...lanes },
  values: { prepare: perPositionValues, content: false, stream: (s: PerPosition) => s.stream ?? null },
  color: { prepare: (input) => perPosition(input), content: false },
  ruler: { prepare: prepareRuler, content: false },
  pairs: { prepare: preparePairs, content: false },
  rows: { prepare: prepareRows, header: RowsHeader, stream: rowsStream },
  logo: { prepare: (input) => prepareMatrix(input, true), content: false },
  matrix: { prepare: (input) => prepareMatrix(input), content: false },
  "seqv:map": { prepare: prepareMap, content: false },
};
