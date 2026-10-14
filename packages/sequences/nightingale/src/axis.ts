/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Nightingale's axis (seqspec:libraries#nightingale): one per `nightingale-manager`, 1-based, one length. Side by side,
// a section's blocks laid end to end on one manager; stacked, a manager per block. Only shown positions are on it. Also
// the edges `mount` uses: block-system runs onto the axis (`AxisSink`, one converter for `materialize`'s placements and
// for a host's highlight through the shared projection), an axis position back to a block and index, an event's extent
// as runs of block indices, and a host's loci as each manager's highlight.

import {
  Arrangement,
  type Closure,
  type Projection,
  projectionOf,
  type ResolvedView,
  type RunSink,
  sectionHighlights,
} from '@molstar/seq-core/format';
import {
  type CoordinateMap,
  Coords,
  Extents,
  type Loci,
  type MapHit,
  OrderedSet,
  type SeqSet,
} from '@molstar/seq-core/model';

/** One block on a manager's axis; plain data (it travels in the element description). */
export interface AxisBlock {
  /** The view's block. */
  readonly block: number;
  /** Its first axis index (0-based). */
  readonly offset: number;
  /** Positions it puts on the axis. */
  readonly count: number;
  /** Its first index (`view.blocks.first`: a piece's start; 0 for columns). */
  readonly first: number;
  /** A ligand in the position form: every index is its one position. */
  readonly ligand: boolean;
  /** The shown subset (`positions`), within the piece; null: all of the piece. */
  readonly shown: OrderedSet | null;
  /** A reversed piece (`dir` -1): `shown`'s ranks run from its last. */
  readonly reversed: boolean;
  /** The block's layout map where indices are not `first + k` (a reversed piece, one across a cyclic origin); else null. */
  readonly layout: CoordinateMap | null;
  readonly unit: number;
  readonly members: readonly number[] | null;
  readonly alignment: number;
}
export interface AxisMeta {
  readonly section: number;
  readonly scrollGroup: string | null;
  readonly length: number;
  readonly blocks: readonly AxisBlock[];
  /** The section's `highlights` on this axis (seqspec:format#sections): kept lit under every host highlight. */
  readonly highlights?: AxisRuns;
}

/** The managers of one section: one side by side, one per block stacked (each the longest block's length). */
export function axesOf(view: ResolvedView, set: SeqSet, section: number): AxisMeta[] {
  const S = view.sections;
  const B = view.blocks;
  const scrollGroup = S.scrollGroup[section] >= 0 ? set.strings.get(S.scrollGroup[section]) : null;
  const blocks: AxisBlock[] = [];
  let offset = 0;
  let longest = 1;
  for (let b = S.blockStart[section]; b < S.blockStart[section] + S.blockCount[section]; b++) {
    const shown = B.positions[b];
    const ligand = B.ligand[b] === 1;
    const layout = B.dir[b] !== 1 || B.layout[b].count !== 1 ? B.layout[b] : null;
    let count = ligand ? 1 : shown ? OrderedSet.size(shown) : 0;
    for (let s = 0; !ligand && !shown && s < B.layout[b].count; s++) count += B.layout[b].length[s];
    const members = B.members[b] ? Array.from(B.members[b] ?? []) : null;
    const [unit, alignment, first] = [B.unit[b], B.alignment[b], B.first[b]];
    blocks.push({
      block: b,
      offset,
      count,
      first,
      ligand,
      shown: ligand ? null : shown,
      reversed: B.dir[b] < 0,
      layout,
      unit,
      members,
      alignment,
    });
    offset += count;
    if (count > longest) longest = count;
  }
  if (S.arrangement[section] !== Arrangement.stacked) return [{ section, scrollGroup, length: offset, blocks }];
  return blocks.map((x) => ({ section, scrollGroup, length: longest, blocks: [{ ...x, offset: 0 }] }));
}

/** Runs of axis positions, half-open and 0-based, reused: `start`/`stop` double when full and are kept. */
export interface AxisRuns {
  count: number;
  start: Int32Array;
  stop: Int32Array;
  /** False once a run was added below the last one: `finish` sorts. */
  sorted: boolean;
}
const SHIFT = 2 ** 26;
/** Sort keys (`start * 2^26 + length`), reused by every `finish` of unsorted runs. */
let keys = new Float64Array(64);
/** A copy at twice the length (a reused buffer that is full). */
function doubled(a: Int32Array): Int32Array {
  const b = new Int32Array(2 * a.length);
  b.set(a);
  return b;
}

function createRuns(capacity = 8): AxisRuns {
  return { count: 0, start: new Int32Array(capacity), stop: new Int32Array(capacity), sorted: true };
}
/** Appends [s, t); a run that touches or overlaps the last (in order) joins it. */
function addRun(r: AxisRuns, s: number, t: number): void {
  if (s >= t) return;
  const n = r.count;
  if (n > 0) {
    if (s < r.start[n - 1]) r.sorted = false;
    else if (s <= r.stop[n - 1]) {
      if (t > r.stop[n - 1]) r.stop[n - 1] = t;
      return;
    }
  }
  if (n === r.start.length) {
    r.start = doubled(r.start);
    r.stop = doubled(r.stop);
  }
  r.start[n] = s;
  r.stop[n] = t;
  r.count = n + 1;
}
/** Heapsort of the first `n` keys, in place (the typed array's own sort would need a view, or sort them all). */
function heapsort(a: Float64Array, n: number): void {
  for (let i = (n >> 1) - 1; i >= 0; i--) sift(a, i, n);
  for (let e = n - 1; e > 0; e--) {
    const top = a[0];
    a[0] = a[e];
    a[e] = top;
    sift(a, 0, e);
  }
}
function sift(a: Float64Array, from: number, n: number): void {
  const x = a[from];
  let i = from;
  for (let c = 2 * i + 1; c < n; c = 2 * i + 1) {
    if (c + 1 < n && a[c + 1] > a[c]) c++;
    if (a[c] <= x) break;
    a[i] = a[c];
    i = c;
  }
  a[i] = x;
}
/** Sorts unsorted runs in place (packed keys, a heapsort), then joins runs that touch. Returns `r`. */
function finishRuns(r: AxisRuns): AxisRuns {
  const n = r.count;
  if (!r.sorted) {
    if (keys.length < n) keys = new Float64Array(2 * n);
    for (let i = 0; i < n; i++) keys[i] = r.start[i] * SHIFT + (r.stop[i] - r.start[i]);
    heapsort(keys, n);
    for (let i = 0; i < n; i++) {
      const s = Math.floor(keys[i] / SHIFT);
      r.start[i] = s;
      r.stop[i] = s + (keys[i] - s * SHIFT);
    }
    r.sorted = true;
  }
  let w = 0;
  for (let i = 1; i < n; i++) {
    if (r.start[i] <= r.stop[w]) {
      if (r.stop[i] > r.stop[w]) r.stop[w] = r.stop[i];
    } else {
      w++;
      r.start[w] = r.start[i];
      r.stop[w] = r.stop[i];
    }
  }
  r.count = n > 0 ? w + 1 : 0;
  return r;
}
function sameRuns(a: AxisRuns, b: AxisRuns): boolean {
  if (a.count !== b.count) return false;
  for (let i = 0; i < a.count; i++) if (a.start[i] !== b.start[i] || a.stop[i] !== b.stop[i]) return false;
  return true;
}
function copyRuns(to: AxisRuns, from: AxisRuns): void {
  if (to.start.length < from.count) {
    to.start = new Int32Array(from.start.length);
    to.stop = new Int32Array(from.start.length);
  }
  for (let i = 0; i < from.count; i++) {
    to.start[i] = from.start[i];
    to.stop[i] = from.stop[i];
  }
  to.count = from.count;
  to.sorted = true;
}
/** Nightingale's ranges (`a:b,…`, 1-based, inclusive); the string is its boundary, built once per changed manager. */
function runsText(r: AxisRuns): string {
  let out = '';
  for (let i = 0; i < r.count; i++) out += `${i > 0 ? ',' : ''}${r.start[i] + 1}:${r.stop[i]}`;
  return out;
}
/** Parses Nightingale's `a:b,…` (1-based, inclusive) into `out` with a char-code scan; ranges without a colon are skipped. */
function parseRuns(text: string, out: AxisRuns): AxisRuns {
  let a = 0;
  let b = 0;
  let side = 0;
  let n = 0;
  for (let i = 0; i <= text.length; i++) {
    const c = i < text.length ? text.charCodeAt(i) : 44;
    if (c >= 48 && c <= 57) {
      if (side === 0) a = a * 10 + c - 48;
      else b = b * 10 + c - 48;
      n++;
    } else if (c === 58) side = 1;
    else if (c === 44) {
      if (side === 1 && n > 0 && a >= 1 && b >= a) addRun(out, a - 1, b);
      a = b = side = n = 0;
    }
  }
  return finishRuns(out);
}
export const AxisRuns = {
  create: createRuns,
  add: addRun,
  finish: finishRuns,
  equal: sameRuns,
  copy: copyRuns,
  text: runsText,
  parse: parseRuns,
};

const hit: MapHit = { element: -1, phase: 0 };
function toAxisRun(a: number, z: number, sink: AxisSink): void {
  addRun(sink.runs, sink.block.offset + a, sink.block.offset + z);
}
/**
 * Block-system runs onto a block's axis, into `runs` (the converter of `materialize`'s placements and of highlights):
 * a ligand block one position; a block shown whole `first + k` (clipped); a `positions` subset by its ranks (cut-out
 * stretches abut, so a run stays one run); a reversed piece or one across a cyclic origin through its layout map.
 */
export class AxisSink implements RunSink {
  block!: AxisBlock;
  runs: AxisRuns;
  constructor(runs: AxisRuns = createRuns()) {
    this.runs = runs;
  }
  add(s: number, t: number): void {
    const b = this.block;
    if (b.ligand) addRun(this.runs, b.offset, b.offset + 1);
    else if (b.shown) {
      const k0 = OrderedSet.lowerBound(b.shown, s);
      const k1 = OrderedSet.lowerBound(b.shown, t);
      if (k0 < k1 && b.reversed) addRun(this.runs, b.offset + b.count - k1, b.offset + b.count - k0);
      else if (k0 < k1) addRun(this.runs, b.offset + k0, b.offset + k1);
    } else if (b.layout) Coords.forEachTo(b.layout, s, t, toAxisRun, this);
    else {
      const lo = s > b.first ? s : b.first;
      const hi = t < b.first + b.count ? t : b.first + b.count;
      if (lo < hi) addRun(this.runs, b.offset + lo - b.first, b.offset + hi - b.first);
    }
  }
}

let found = -1;
function pick(h: MapHit): void {
  if (found < 0) found = h.element;
}
/** The first `from` element whose image under `map` is `x`; -1: none. Allocates nothing (NG-20). */
export function firstFrom(map: CoordinateMap, x: number): number {
  found = -1;
  Coords.forEachFromAt(map, x, hit, pick, 0);
  return found;
}
/** The block index of a block's k-th position on the axis; -1 when the layout holds none. Allocates nothing. */
export function indexAt(b: AxisBlock, k: number): number {
  if (b.ligand) return b.first;
  if (b.shown) return OrderedSet.getAt(b.shown, b.reversed ? b.count - 1 - k : k);
  if (!b.layout) return b.first + k;
  return firstFrom(b.layout, k);
}

/** Where an axis position lies, written by `locate` (reused by the caller). */
export interface Located {
  block: number;
  index: number;
}
/** The block and index at an axis index (0-based) by binary search over the blocks' offsets; false: none. */
export function locate(meta: AxisMeta, at: number, out: Located): boolean {
  const bs = meta.blocks;
  let lo = 0;
  let hi = bs.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (bs[mid].offset <= at) lo = mid;
    else hi = mid - 1;
  }
  const b = bs[lo];
  if (!b || at < b.offset || at >= b.offset + b.count) return false;
  out.block = b.block;
  out.index = indexAt(b, at - b.offset);
  return out.index >= 0;
}

/** A run of one block's indices, [start, stop). */
export interface BlockRun {
  readonly block: number;
  readonly start: number;
  readonly stop: number;
}
/** Block runs, reused (an event's extent made from it once per change). */
export interface BlockRuns {
  count: number;
  block: Int32Array;
  start: Int32Array;
  stop: Int32Array;
}
function addBlockRun(r: BlockRuns, block: number, start: number, stop: number): void {
  const n = r.count;
  if (n > 0 && r.block[n - 1] === block && r.stop[n - 1] >= start && start >= r.start[n - 1]) {
    if (stop > r.stop[n - 1]) r.stop[n - 1] = stop;
    return;
  }
  // A reversed piece's runs descend: one that ends where the last starts joins it.
  if (n > 0 && r.block[n - 1] === block && stop === r.start[n - 1]) {
    r.start[n - 1] = start;
    return;
  }
  if (n === r.block.length) {
    r.block = doubled(r.block);
    r.start = doubled(r.start);
    r.stop = doubled(r.stop);
  }
  r.block[n] = block;
  r.start[n] = start;
  r.stop[n] = stop;
  r.count = n + 1;
}
export const BlockRuns = {
  create: (capacity = 8): BlockRuns => ({
    count: 0,
    block: new Int32Array(capacity),
    start: new Int32Array(capacity),
    stop: new Int32Array(capacity),
  }),
  add: addBlockRun,
  /** The runs as event data (made once per change). */
  toArray: (r: BlockRuns): BlockRun[] =>
    Array.from({ length: r.count }, (_, i) => ({ block: r.block[i], start: r.start[i], stop: r.stop[i] })),
};
const ctx = { runs: null as BlockRuns | null, block: 0 };
function fromLayout(s: number, t: number, c: typeof ctx): void {
  if (c.runs) addBlockRun(c.runs, c.block, s, t);
}

/** Axis runs as runs of block indices into `out`, in order; a `positions` subset splits a run where its indices jump. */
export function extentOf(meta: AxisMeta, axis: AxisRuns, out: BlockRuns): BlockRuns {
  out.count = 0;
  for (let i = 0; i < axis.count; i++) {
    for (let j = 0; j < meta.blocks.length; j++) {
      const x = meta.blocks[j];
      const k0 = Math.max(axis.start[i], x.offset) - x.offset;
      const k1 = Math.min(axis.stop[i], x.offset + x.count) - x.offset;
      if (k0 >= k1) continue;
      if (x.ligand) addBlockRun(out, x.block, x.first, x.first + 1);
      else if (x.shown) for (let k = k0; k < k1; k++) addBlockRun(out, x.block, indexAt(x, k), indexAt(x, k) + 1);
      else if (x.layout) {
        ctx.runs = out;
        ctx.block = x.block;
        Coords.forEachFrom(x.layout, k0, k1, fromLayout, ctx);
      } else addBlockRun(out, x.block, x.first + k0, x.first + k1);
    }
  }
  return out;
}

/** Runs of one block, as an ordered set: sorted, joined, one pass (an Interval when it is one run). */
function setOf(runs: readonly BlockRun[]): OrderedSet {
  const sorted = runs.length > 1 ? runs.slice().sort((a, b) => a.start - b.start) : runs;
  const [starts, stops]: number[][] = [[], []];
  for (const r of sorted) {
    const n = stops.length - 1;
    if (n >= 0 && r.start <= stops[n]) stops[n] = Math.max(stops[n], r.stop);
    else {
      starts.push(r.start);
      stops.push(r.stop);
    }
  }
  if (starts.length === 1) return OrderedSet.ofRange(starts[0], stops[0]);
  const out: number[] = [];
  for (let i = 0; i < starts.length; i++) for (let x = starts[i]; x < stops[i]; x++) out.push(x);
  return OrderedSet.ofSortedArray(Int32Array.from(out));
}
/**
 * An extent as loci, so every manager lights its counterpart (`highlightOf`): residues of a block's unit (each
 * member's, merged copies), else the columns of the first block's alignment. Null when empty.
 */
export function extentLoci(view: ResolvedView, set: SeqSet, extent: readonly BlockRun[]): Loci | null {
  if (extent.length === 0) return null;
  const B = view.blocks;
  const alignment = B.alignment[extent[0].block];
  if (alignment >= 0) {
    const columns = setOf(extent.filter((r) => B.alignment[r.block] === alignment)) as never;
    return { kind: 'column-loci', set, alignment: alignment as never, columns };
  }
  const units = new Map<number, BlockRun[]>();
  for (const r of extent) {
    if (B.alignment[r.block] >= 0) continue;
    for (const u of B.members[r.block] ?? [B.unit[r.block]]) (units.get(u) ?? units.set(u, []).get(u))?.push(r);
  }
  const elements = Array.from(units, ([unit, runs]) => ({ unit: unit as never, indices: setOf(runs) as never }));
  return { kind: 'seq-loci', set, elements };
}

/**
 * A manager's highlight runs for the entries of `ext`: each of its blocks' projection (the shared one: the same routes,
 * closure and coordinate-system rules as the viewer), converted by `sink` into `out`, sorted and joined.
 */
export function highlightRuns(p: Projection, ext: Extents, meta: AxisMeta, sink: AxisSink, out: AxisRuns): AxisRuns {
  out.count = 0;
  out.sorted = true;
  sink.runs = out;
  for (let i = 0; i < meta.blocks.length; i++) {
    sink.block = meta.blocks[i];
    p.project(ext, meta.blocks[i].block, sink);
  }
  return finishRuns(out);
}
/** The section's `highlights` on a manager's axis: every entry's items, whatever its style. */
export function sectionRuns(view: ResolvedView, set: SeqSet, meta: AxisMeta): AxisRuns {
  const ext = Extents.create(16);
  for (const h of sectionHighlights(view, set, meta.section)) if (h.loci) Extents.addLoci(ext, h.loci);
  return highlightRuns(projectionOf(view, set), ext, meta, new AxisSink(), createRuns());
}
/** `kept` (a section's `highlights`) added to a manager's runs `out`, sorted and joined: allocation-free once grown. */
export function keepLit(out: AxisRuns, kept: AxisRuns | undefined): AxisRuns {
  if (!kept || kept.count === 0) return out;
  for (let k = 0; k < kept.count; k++) addRun(out, kept.start[k], kept.stop[k]);
  return finishRuns(out);
}
/** `ext` made from a host's loci (cleared first), closed when asked; false: none, or loci of another set. */
export function prepareLoci(p: Projection, ext: Extents, loci: Loci | null, closure: Closure): boolean {
  ext.count = 0;
  if (!loci || loci.set !== p.set) return false;
  Extents.addLoci(ext, loci);
  if (closure === 'containing') p.close(ext);
  return ext.count > 0;
}
