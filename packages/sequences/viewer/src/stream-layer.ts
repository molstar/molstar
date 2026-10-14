/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// One layer's held regions and every per-frame read of them (seqspec:sources#the-viewers-side): `Region`, `Layer` (the
// `StreamView` drawing reads), and the 1-D routes through a block's path. Allocation-free on its read paths.

import { type StreamedData, TILE, type Tile, Tiles } from "@molstar/seq-core/format";
import { type ComposedPath, Paths } from "@molstar/seq-core/model";
import type { Hint } from "./prefetch.ts";
import { LETTERS_PER_PX, RegionState, type StreamItems, type StreamLayer, type StreamView } from "./stream-contract.ts";
import * as Levels from "./stream-levels.ts";

// Read once into module constants: a read through an import is a property read of the module per use, which
// the test runner's module wrapper turns into an allocation on per-frame paths.
const { D1, L, LABELS, LV, servedLevels } = Levels;
const LETTERS_GATE = LETTERS_PER_PX;
const PENDING = RegionState.pending;
const LOADED = RegionState.loaded;
const EMPTY = RegionState.empty;
const FAILED = RegionState.failed;

/** A 1-D region of data (level -1): this many positions (items: `ITEMS_SPAN`); of a level, `TILE` bins. */
export const DATA_SPAN = TILE * 16;
const ITEMS_SPAN = TILE * 256;
/**
 * The zoom gate's finer data regions (ST4b): the level slot of an items layer's regions of `DATA_SPAN` positions, asked
 * inside a data region its provider answered with bins, where a device px holds at most `FINE_PER_PX` positions.
 */
export const FINE = 30;
/** Finer regions where a device px holds at most this many positions: about a device px's 64, 16 such regions at 1,024 px. */
export const FINE_PER_PX = 64;
/** An items layer's regions in view at most, a coarser level taken past it (ST4b). */
export const MAX_REGIONS = 16;
/** What `items` gives for a region not loaded (its state in `out[4]`). */
const NO_ITEMS: StreamItems = {
  features: null,
  first: new Int32Array(0),
  bins: null,
  width: 1,
  table: null,
  item: null,
  ends: null,
  hint: 0,
};
/** An item's lane, its number (the hover's `part`) and how many held regions hold it (ST3b: forgotten at 0). */
interface ItemLane {
  lane: number;
  readonly id: number;
  held: number;
}
/** Past this, a region's age reads as this (ms; a Smi). */
const AGE_MAX = 1 << 30;

export class Region {
  readonly owner: Layer;
  /** ly·16 + lx, or `LABELS` (a chunk of TILE rows at `ty`). */
  readonly level: number;
  readonly ty: number;
  readonly tx: number;
  readonly index: number;
  state: number = PENDING;
  tile: Tile | null = null;
  /** 1-D values: per bin `min, max, mean` (stride 3) or per position (stride 1), NaN none. */
  vals: Float64Array | null = null;
  /** 1-D letters: the data's upper-case codes from the region's start. */
  codes: Uint8Array | null = null;
  /** 1-D items: their table or bins; and the block indices its change redraws (-1: its extent). */
  items: StreamItems | null = null;
  /** The same items folded by the layer's `collapse` field (built beside `items`); null: none. */
  folded: StreamItems | null = null;
  lo = -1;
  hi = -1;
  stride = 1;
  bytes = 0;
  /** The pass that last needed it: the LRU's stamp; the current pass's pins it. */
  need = 0;
  /** 0 coarse, 1 the view's level, 2 the vertical overview, 3 labels, 4 the margin. */
  group = 0;
  /** The pass that last needed it for the view (a group below the margin's): `settled` counts it. */
  viewed = 0;
  flight: AbortController | null = null;
  /** Failures in this outage, and when it may be asked again (`NotServed`: never). */
  tries = 0;
  retryAt = 0;
  /** In a run of failures (its layer's outage). */
  failing = false;
  /** In its layer's `all`. */
  at = -1;
  /** The dynamic nodes its patch carried (held while it is). */
  carried: string[] | null = null;
  /** The frame time it was first drawn at after it landed (-1: not yet; the overlay's fade-in); items: their keys. */
  landed = -1;
  /** An items data region its provider answered with this level's bins (the zoom gate); -1: items. */
  gate = -1;
  /** A gated data region clicked: its part in view is asked in finer regions whatever the zoom, each forced past the
   * gate (a "Load anyway"; ST4b). */
  force = false;
  keys: string[] | null = null;
  foldedKeys: string[] | null = null;
  constructor(owner: Layer, level: number, ty: number, tx: number, index: number) {
    this.owner = owner;
    this.level = level;
    this.ty = ty;
    this.tx = tx;
    this.index = index;
  }
}

/** Whether block name `name` is `target`, or begins with it and a "/" (`chr17/unit` is `chr17`'s). No allocation. */
export function names(name: string, target: string): boolean {
  return name.startsWith(target) && (name.length === target.length || name.charCodeAt(target.length) === 47);
}

/** A tile's index in its level; -1 outside the source's columns. */
export function indexOf(columns: number, lx: number, ty: number, tx: number): number {
  const across = Tiles.count(columns, lx);
  return ty < 0 || tx < 0 || tx >= across ? -1 : (ty * across + tx) | 0;
}

/**
 * Region r's age into `out[at]` at the frame time `out[0]` (ms; negative: no time, nothing stamped), stamped at the first
 * frame that asks after it landed (so regions landing together fade together); none: `AGE_MAX`. Through the array: a
 * double returned is boxed.
 */
function aged(r: Region | undefined, out: Float64Array, at: number) {
  const held = r !== undefined && (r.state === LOADED || r.state === EMPTY);
  if (held && r.landed < 0 && out[0] >= 0) r.landed = out[0];
  out[at] = held && r.landed >= 0 ? out[0] - r.landed : AGE_MAX;
}

/** A range walk's runs as their hull in `h` (`Paths`' visitor: no closure a call). */
function hullOf(s: number, t: number, h: Int32Array) {
  if (s < h[0]) h[0] = s;
  if (t > h[1]) h[1] = t;
}
/** A path's runs a pass asks for, appended (grown at edge rate). */
export const routes = { runs: new Int32Array(64), n: 0 };
export function routeOf(s: number, t: number, r: typeof routes) {
  if (2 * r.n + 2 > r.runs.length) {
    const g = new Int32Array(2 * r.runs.length);
    g.set(r.runs);
    r.runs = g;
  }
  r.runs[2 * r.n] = s;
  r.runs[2 * r.n + 1] = t;
  r.n++;
}
/** `routes`' runs sorted by start (an insertion sort: a few runs, no allocation). */
export function sortRuns(r: typeof routes) {
  const x = r.runs;
  for (let i = 1; i < r.n; i++) {
    const a = x[2 * i];
    const e = x[2 * i + 1];
    let j = i - 1;
    for (; j >= 0 && x[2 * j] > a; j--) {
      x[2 * j + 2] = x[2 * j];
      x[2 * j + 3] = x[2 * j + 1];
    }
    x[2 * j + 2] = a;
    x[2 * j + 3] = e;
  }
}
/** The distinct regions of `span` positions sorted runs [0, n) of `runs` touch. */
export function regionsOver(runs: Int32Array, n: number, span: number): number {
  let count = 0;
  let last = -1;
  for (let q = 0; q < n; q++) {
    const j0 = Math.max(last + 1, Math.floor(runs[2 * q] / span));
    const j1 = Math.floor((runs[2 * q + 1] - 1) / span);
    if (j1 >= j0) count += j1 - j0 + 1;
    if (j1 > last) last = j1;
  }
  return count;
}
function valuesRun(s: number, t: number, c: Layer) {
  c.valuesRun(s, t);
}
const NO_VALUES = new Float64Array(4);
export const NO_HINTS: readonly Hint[] = [];
/** The letter at offset o of a letters region held; 0: none. */
export function letterIn(r: Region | undefined, o: number): number {
  const codes = r !== undefined && r.state === LOADED ? r.codes : null;
  return codes !== null && o < codes.length ? codes[o] : 0;
}

/** One layer's regions, by level and index (small-integer keys), and in a list for sweeps without an iterator. */
/** Changes a layer logs for `changed`. */
const LOG = 16;
export class Layer implements StreamView {
  readonly ref: string;
  source: StreamLayer;
  version = 0;
  /** Per change, by its version mod `LOG`: the version and the region's columns [c0, c1) (a labels chunk: -1, -1). */
  private readonly log = new Int32Array(3 * LOG);
  budget = 0;
  /** The layer's index in the current view; -1: paused. */
  layer = -1;
  bytes = 0;
  /** Its regions failing (an outage while above 0; `NotServed` ones aside). */
  failing = 0;
  closed = false;
  /**
   * A layer of no view (ST7's codon label): a lazy molecule's letters, asked a few positions at a time
   * (`lettersAt`), the regions wanted (system, index) until held.
   */
  aside = false;
  wanted = new Int32Array(8);
  wantedN = 0;
  /** The local tile in progress, and its next sample row. */
  job: Region | null = null;
  row = 0;
  readonly levels: (Map<number, Region> | undefined)[] = [];
  readonly all: Region[] = [];
  /** 1-D: per block the view's level (-1 data; -2 none), and the manifest's widths; 2-D: grid level → listed level. */
  chosen = new Int16Array(0);
  /**
   * 1-D, where a miniature draws the layer: per block the map's level (-2: none) and the level a map px gives before
   * `MAX_REGIONS` made it coarser (a read at least that coarse reads the map's).
   */
  mapChosen = new Int16Array(0);
  mapBase = new Int16Array(0);
  widths = new Int32Array(0);
  grid = new Int16Array(L * L);
  /** 2-D: per grid level the level served in its place (`servedLevels`); null: a local stream's, every level. */
  served: Int16Array | null = null;
  /** 2-D: the view's fine level and its first lane tile, for `state`. */
  readonly fine = new Int32Array(3);
  /** A letters stream: its data only, within `gate` positions a device px; an items stream; a data region's span. */
  seq = false;
  isItems = false;
  span = DATA_SPAN;
  gate = LETTERS_GATE;
  /** Items: the lanes used, and each item's lane and number by its key, while a held region holds it. */
  lanes = 1;
  readonly laneOf = new Map<string, ItemLane>();
  /** The folded items' (`collapse`): their lanes, and each group's lane by its key. */
  foldedLanes = 1;
  readonly foldedLaneOf = new Map<string, ItemLane>();
  /** The last `fold`: the items given are the folded ones. */
  folding = false;
  /**
   * How its section's view last moved, from every change of it (not only the passes: a pan before the first pass has an
   * earlier view too): the scroll and cell it had, the direction it last panned (-1, 0, 1) and zoomed (-1: in).
   */
  readonly moves = new Float64Array(4);
  /** The layer's track, and whether the person expanded it (set when the view is). */
  track = "";
  isExpanded: (() => boolean) | null = null;
  ids = 0;
  /** The items cursor's walk: per region its level key, index, start and stop, and how many (grown only). */
  walks = new Int32Array(64);
  walked = 0;
  /** Per block, its system's positions a block index covers (through a path: 3 a residue; else 1). */
  ratio = new Float64Array(0);
  /** `back`'s and `toBlock`'s hull. */
  readonly hull = new Int32Array(2);
  /** The walk's last region and finer region listed (a path's runs may share one: listed once). */
  private lastJ = -1;
  private lastQ = -1;
  /** `region` through a path: the regions its runs reach at the view's level, listed once per read (k 0). */
  private rlist = new Int32Array(16);
  private rcount = 0;
  /** `values` through a path: per px, then the minimum, maximum, sum and count over its runs; its system and `out`. */
  private readonly acc = new Float64Array(5);
  private sys1 = -1;
  private out1: Float64Array = NO_VALUES;
  private held1 = false;
  constructor(ref: string, source: StreamLayer) {
    this.ref = ref;
    this.source = source;
    this.listed();
    this.routed();
  }
  /** The source's paths read (edge-rate): each block's ratio. */
  routed() {
    const p = this.source.paths;
    const n = this.source.systems?.length ?? 0;
    this.ratio = Float64Array.from({ length: n }, (_, b) => (p && b < p.length && p[b] ? Paths.ratio(p[b]) : 1));
  }
  /** Block b's path from its system (ST7); null: shown directly. */
  pathOf(b: number): ComposedPath | null {
    const p = this.source.paths;
    return p != null && b < p.length ? p[b] : null;
  }
  /** Block b's indices [i0, i1) as the hull of its system's positions that reach them, into `hull`; false: none. */
  back(b: number, i0: number, i1: number): boolean {
    const p = this.pathOf(b);
    const h = this.hull;
    if (p === null) {
      h[0] = i0;
      h[1] = i1;
      return i1 > i0;
    }
    h[0] = 0x7fffffff;
    h[1] = -1;
    Paths.forEachFrom(p, i0, i1, hullOf, h);
    return h[1] > h[0];
  }
  toBlock(b: number, s0: number, s1: number, out: Int32Array): boolean {
    const p = this.pathOf(b);
    if (p === null) {
      out[0] = s0;
      out[1] = s1;
      return s1 > s0;
    }
    out[0] = 0x7fffffff;
    out[1] = -1;
    Paths.forEachTo(p, s0, s1, hullOf, out);
    return out[1] > out[0];
  }
  /** The manifest's levels read once per source: 1-D widths; an alignment's [2^ly, 2^lx] by grid level (-1: not served). */
  listed() {
    const lv = this.source.remote?.manifest.levels ?? [];
    this.seq = this.source.remote?.source.data_kind === "letters";
    this.isItems = this.source.remote?.source.data_kind === "items";
    this.span = this.isItems ? ITEMS_SPAN : DATA_SPAN;
    // Positions a bin, integers: the reads' index arithmetic stays in int32 (ST3b: a double phi there was boxed).
    this.widths = Int32Array.from(lv, (w) => (typeof w === "number" ? w : 0));
    this.grid.fill(-1);
    lv.forEach((w, i) => {
      if (typeof w === "number") return;
      const [ly, lx] = [Math.log2(w[0]), Math.log2(w[1])];
      if (Number.isInteger(ly) && Number.isInteger(lx) && ly < L && lx < L) this.grid[ly * L + lx] = i;
    });
    this.served = null;
    if (this.source.remote === null || this.source.columns <= 0) return;
    this.served = new Int16Array(L * L);
    servedLevels(this.grid, this.served);
  }
  /** Positions a region of level `lvl` covers (-1: the data; `FINE`: the zoom gate's finer data). */
  spanOf(lvl: number): number {
    return lvl === FINE ? DATA_SPAN : lvl < 0 ? this.span : TILE * this.widths[lvl];
  }
  find(level: number, index: number): Region | undefined {
    const m = this.levels[level];
    return m === undefined ? undefined : m.get(index);
  }
  tile(ly: number, lx: number, ty: number, tx: number): Tile | null {
    const r = this.find(ly * L + lx, indexOf(this.source.columns, lx, ty, tx));
    return r !== undefined && r.state === LOADED ? r.tile : null;
  }
  /** Outside the source: empty. */
  state(ly: number, lx: number, ty: number, tx: number): number {
    const index = indexOf(this.source.columns, lx, ty, tx);
    const r = this.find(ly * L + lx, index);
    return r !== undefined ? r.state : index < 0 || ty * TILE * 2 ** ly >= this.source.rows ? EMPTY : PENDING;
  }
  label(row: number): string | null {
    return this.source.remote?.labels.get(row) ?? null;
  }
  values(b: number, i0: number, i1: number, out: Float64Array): boolean {
    const sys = this.source.systems;
    const s = sys !== null && b < sys.length ? sys[b] : -1;
    if (s < 0) return false;
    const p = this.pathOf(b);
    if (p === null) return this.valuesOf(s, i0, i1, out);
    // Through a path: each run of its system the indices reach (a codon split by an intron: two), merged.
    const a = this.acc;
    a[0] = out[0] * this.ratio[b];
    a[1] = a[2] = Number.NaN;
    a[3] = a[4] = 0;
    this.sys1 = s;
    this.out1 = out;
    this.held1 = false;
    Paths.forEachFrom(p, i0, i1, valuesRun, this);
    out[0] = a[1];
    out[1] = a[2];
    out[2] = a[3];
    out[3] = a[4];
    this.out1 = NO_VALUES;
    return this.held1;
  }
  /** One run of `values` through a path, merged into `acc`. */
  valuesRun(s0: number, s1: number) {
    const o = this.out1;
    const a = this.acc;
    o[0] = a[0];
    if (!this.valuesOf(this.sys1, s0, s1, o)) return;
    this.held1 = true;
    if (o[3] <= 0) return;
    if (!(o[0] >= a[1])) a[1] = o[0];
    if (!(o[1] <= a[2])) a[2] = o[1];
    a[3] += o[2];
    a[4] += o[3];
  }
  /** `values` over system s's positions [i0, i1). */
  private valuesOf(s: number, i0: number, i1: number, out: Float64Array): boolean {
    const W = this.widths;
    let lvl = -1;
    for (let k = 0; k < W.length; k++) if (W[k] <= out[0]) lvl = k;
    for (let l = lvl; l < W.length; l++) {
      // In int32 throughout (indices are non-negative): a double here was boxed, about 31 bytes a frame at 100 kb with
      // an items layer beside it (ST3b; ST3a-2 measured 23).
      const w = l < 0 ? 1 : W[l];
      let j0 = ((i0 + w - 1) / w) | 0;
      let j1 = ((i1 + w - 1) / w) | 0;
      if (j1 <= j0) {
        j0 = (i0 / w) | 0;
        j1 = j0 + 1;
      }
      if (this.bins(s, l, j0, j1, out)) return true;
    }
    return false;
  }
  letters(b: number, i0: number, i1: number, out: Uint8Array): number {
    const sys = this.source.systems;
    const s = sys !== null && b < sys.length ? sys[b] : -1;
    const m = s < 0 ? undefined : this.levels[D1 + s * LV];
    const n = Math.min(i1 - i0, out.length);
    let held = 0;
    for (let i = 0; i < n; ) {
      const ri = Math.floor((i0 + i) / this.span);
      const e = Math.min(n, (ri + 1) * this.span - i0);
      const r = m === undefined ? undefined : m.get(ri);
      const codes = r !== undefined && r.state === LOADED ? r.codes : null;
      // `o`: the index in `out` of the region's first code.
      const o = ri * this.span - i0;
      const to = codes === null ? i : Math.max(i, Math.min(e, o + codes.length));
      for (let j = i; j < to; j++) out[j] = (codes as Uint8Array)[j - o];
      out.fill(0, to, e);
      held += to - i;
      i = e;
    }
    return held;
  }
  items(b: number, i0: number, i1: number, k: number, out: Int32Array): StreamItems | null {
    const sys = this.source.systems;
    const s = sys !== null && b < sys.length && this.isItems ? sys[b] : -1;
    if (s < 0) return null;
    // The regions over [i0, i1) in order, listed once a walk (k 0): a gated data region's finer ones in its place.
    // Through a path: over the hull of the positions the indices reach, at its system's positions a device px.
    // Through a path: only the regions its runs reach, at its system's positions a device px (-1: a data read).
    const p = this.pathOf(b);
    if (k === 0) this.walk(b, s, i0, i1, p === null || out[0] < 0 ? out[0] : (out[0] * this.ratio[b]) | 0);
    if (k >= this.walked) return null;
    const w = this.walks;
    const r = this.find(w[4 * k], w[4 * k + 1]);
    out[1] = w[4 * k + 2];
    out[2] = w[4 * k + 3];
    // Below the previous region's stop its items are an earlier region's (a region between two listed ones, never
    // listed through a path, draws nothing: its items are the next one's).
    out[3] = k === 0 ? -1 : w[4 * k - 1];
    if (out.length > 6) out[6] = out[1];
    if (p !== null && this.toBlock(b, out[1], out[2], this.hull)) {
      out[1] = this.hull[0];
      out[2] = this.hull[1];
    } else if (p !== null) out[1] = out[2] = 0;
    out[4] = r === undefined ? PENDING : r.state;
    out[5] = this.folding ? this.foldedLanes : this.lanes;
    if (r === undefined || r.state !== LOADED) return NO_ITEMS;
    return this.folding && r.folded !== null ? r.folded : (r.items ?? NO_ITEMS);
  }
  fold(fit: number) {
    this.folding = this.source.collapse != null && fit < this.lanes && !(this.isExpanded?.() ?? false);
  }
  /**
   * The regions over block b's [i0, i1) at `per` positions a device px (-1: a data read) into `walks` (level key,
   * index, start, stop): its system's positions there, or through its path the runs that reach them (ST7).
   */
  private walk(b: number, s: number, i0: number, i1: number, per: number) {
    const W = this.widths;
    let lvl = -1;
    for (let q = 0; q < W.length; q++) if (W[q] <= per) lvl = q;
    // The level the pass chose for the block where it chose one at this zoom (`MAX_REGIONS`); a data read asks for
    // the data.
    const chosen = b < this.chosen.length ? this.chosen[b] : -2;
    if (per >= 0 && chosen >= lvl) lvl = chosen;
    // A miniature's read (as coarse as the map's px, or coarser): the level the pass asked for the map.
    const mapped = b < this.mapChosen.length ? this.mapChosen[b] : -2;
    if (per >= 0 && mapped > lvl && lvl >= this.mapBase[b]) lvl = mapped;
    this.walked = 0;
    this.lastJ = this.lastQ = -1;
    const p = this.pathOf(b);
    if (p === null) {
      this.walkOver(s, lvl, i0, i1, per);
      return;
    }
    routes.n = 0;
    Paths.forEachFrom(p, i0, i1, routeOf, routes);
    sortRuns(routes);
    for (let q = 0; q < routes.n; q++) this.walkOver(s, lvl, routes.runs[2 * q], routes.runs[2 * q + 1], per);
  }
  /** `walk`'s regions over system s's [i0, i1) at level `lvl`, each listed once a walk. */
  private walkOver(s: number, lvl: number, i0: number, i1: number, per: number) {
    const span = this.spanOf(lvl);
    const ext = (this.source.remote as StreamedData).manifest.systems[s].extent;
    const end = Math.min(i1, ext[1]);
    for (let j = Math.floor(Math.max(i0, ext[0]) / span) | 0; j * span < end; j++) {
      const r = this.find(D1 + s * LV + lvl + 1, j);
      const fine = lvl < 0 && r !== undefined && r.state === LOADED && r.gate >= 0 && (per <= FINE_PER_PX || r.force);
      if (!fine) {
        if (j > this.lastJ)
          this.walkTo(D1 + s * LV + lvl + 1, j, Math.max(ext[0], j * span), Math.min(ext[1], (j + 1) * span));
        if (j > this.lastJ) this.lastJ = j;
        continue;
      }
      if (j > this.lastJ) this.lastJ = j;
      const a = Math.max(i0, ext[0], j * span);
      const e = Math.min(end, (j + 1) * span);
      for (let q = Math.floor(a / DATA_SPAN) | 0; q * DATA_SPAN < e; q++) {
        if (q <= this.lastQ) continue;
        this.lastQ = q;
        this.walkTo(D1 + s * LV + FINE + 1, q, Math.max(ext[0], q * DATA_SPAN), Math.min(ext[1], (q + 1) * DATA_SPAN));
      }
    }
  }
  private walkTo(key: number, index: number, start: number, stop: number) {
    if (4 * this.walked === this.walks.length) {
      const w = new Int32Array(2 * this.walks.length);
      w.set(this.walks);
      this.walks = w;
    }
    const o = 4 * this.walked++;
    this.walks[o] = key;
    this.walks[o + 1] = index;
    this.walks[o + 2] = start;
    this.walks[o + 3] = stop;
  }
  stateAt(b: number, index: number): number {
    // Through a path: the first position of its system that reaches the index; none: nothing is there (empty).
    if (this.pathOf(b) !== null && !this.back(b, index, index + 1)) return EMPTY;
    const i = this.pathOf(b) !== null ? this.hull[0] : index;
    const sys = this.source.systems;
    const s = sys !== null && b < sys.length ? sys[b] : -1;
    const l = b < this.chosen.length ? this.chosen[b] : -2;
    const ext = s < 0 ? null : (this.source.remote as StreamedData).manifest.systems[s].extent;
    if (ext === null || i < ext[0] || i >= ext[1] || l === -3) return EMPTY;
    // Out of view: no level is needed there, nothing known (pending, never empty).
    if (l < -1) return PENDING;
    const r = this.find(D1 + s * LV + l + 1, Math.floor(i / this.spanOf(l)));
    return r === undefined ? PENDING : r.state;
  }
  region(b: number, i0: number, i1: number, k: number, out: Float64Array): boolean {
    const p = this.pathOf(b);
    if (p === null) return this.regionOf(b, i0, i1, k, out);
    // Through a path: the regions its runs reach, listed once a read (k 0); their extents back on the block.
    const sys = this.source.systems;
    const s = sys !== null && b < sys.length ? sys[b] : -1;
    const l = b < this.chosen.length ? this.chosen[b] : -2;
    if (s < 0 || l < -1) return false;
    if (k === 0) this.listRegions(p, s, l, i0, i1);
    if (k >= this.rcount) return false;
    const j = this.rlist[k];
    this.regionAt(s, l, j, out);
    const ok = this.toBlock(b, out[1], out[2], this.hull);
    out[1] = ok ? this.hull[0] : 0;
    out[2] = ok ? this.hull[1] : 0;
    return true;
  }
  /** The regions of system s at level l the path's runs over [i0, i1) reach, in order, into `rlist`. */
  private listRegions(p: ComposedPath, s: number, l: number, i0: number, i1: number) {
    routes.n = 0;
    Paths.forEachFrom(p, i0, i1, routeOf, routes);
    sortRuns(routes);
    const span = this.spanOf(l);
    const ext = (this.source.remote as StreamedData).manifest.systems[s].extent;
    let last = -1;
    this.rcount = 0;
    for (let q = 0; q < routes.n; q++) {
      const a = Math.max(ext[0], routes.runs[2 * q]);
      const e = Math.min(ext[1], routes.runs[2 * q + 1]);
      for (let j = Math.max(last + 1, Math.floor(a / span) | 0); j * span < e; j++) {
        if (this.rcount === this.rlist.length) {
          const g = new Int32Array(2 * this.rlist.length);
          g.set(this.rlist);
          this.rlist = g;
        }
        this.rlist[this.rcount++] = j;
        last = j;
      }
    }
  }
  private regionOf(b: number, i0: number, i1: number, k: number, out: Float64Array): boolean {
    const sys = this.source.systems;
    const s = sys !== null && b < sys.length ? sys[b] : -1;
    const l = b < this.chosen.length ? this.chosen[b] : -2;
    if (s < 0 || l < -1) return false;
    const span = this.spanOf(l);
    const ext = (this.source.remote as StreamedData).manifest.systems[s].extent;
    const j = (Math.floor(Math.max(i0, ext[0]) / span) + k) | 0;
    if (j * span >= Math.min(i1, ext[1])) return false;
    this.regionAt(s, l, j, out);
    return true;
  }
  /** Region j of system s at level l into `out[1..5]`, as `region` gives it. */
  private regionAt(s: number, l: number, j: number, out: Float64Array) {
    const span = this.spanOf(l);
    const ext = (this.source.remote as StreamedData).manifest.systems[s].extent;
    out[1] = Math.max(ext[0], j * span);
    out[2] = Math.min(ext[1], (j + 1) * span);
    const r = this.find(D1 + s * LV + l + 1, j);
    out[3] = r === undefined ? PENDING : r.state;
    aged(r, out, 4);
    // Values draw the next coarser level loaded in a pending or failed region's place (`values`): no skeleton, no
    // fade after, and a failure there no band (`out[5]` 1: something is drawn in its place; ST4b).
    out[5] = 0;
    for (let c = l + 1; !this.seq && !this.isItems && c < this.widths.length && out[5] === 0; c++) {
      if (out[3] === LOADED) break;
      const q = this.find(D1 + s * LV + c + 1, Math.floor(out[1] / (TILE * this.widths[c])) | 0);
      if (q === undefined || q.state !== LOADED) continue;
      if (out[3] === PENDING) out[3] = LOADED;
      out[5] = 1;
      aged(q, out, 4);
    }
  }
  tileAt(ly: number, lx: number, ty: number, tx: number, out: Float64Array): number {
    const s = this.state(ly, lx, ty, tx);
    aged(this.find(ly * L + lx, indexOf(this.source.columns, lx, ty, tx)), out, 1);
    return s;
  }
  /** Bins [j0, j1) of system s at level l into `out` (`values`); false when a region they lie in is not loaded. */
  private bins(s: number, l: number, j0: number, j1: number, out: Float64Array): boolean {
    const m = this.levels[D1 + s * LV + l + 1];
    if (m === undefined) return false;
    const n = l < 0 ? DATA_SPAN : TILE;
    out[0] = out[1] = Number.NaN;
    out[2] = out[3] = 0;
    // The region's values and stride, read once a region (ST3b: read per bin, beside an items layer's regions, the
    // loads went polymorphic and boxed about 30 bytes a frame at 100 kb).
    let vals: Float64Array | null = null;
    let stride = 1;
    let at = -1;
    for (let j = j0; j < j1; j++) {
      const ri = Math.floor(j / n) | 0;
      if (ri !== at) {
        at = ri;
        const r = m.get(ri);
        if (r === undefined || r.state === PENDING || r.state === FAILED) return false;
        vals = r.vals;
        stride = r.stride;
      }
      const o = ((j - ri * n) * stride) | 0;
      if (vals === null || o >= vals.length) continue;
      const mean = vals[stride === 3 ? o + 2 : o];
      if (Number.isNaN(mean)) continue;
      const lo = vals[o];
      const hi = vals[stride === 3 ? o + 1 : o];
      if (!(lo >= out[0])) out[0] = lo;
      if (!(hi <= out[1])) out[1] = hi;
      out[2] += mean;
      out[3]++;
    }
    return true;
  }
  /** Region r arrived, failed or left: the version bumped and its columns logged (labels and 1-D: unknown, all). */
  bump(r: Region) {
    const v = ++this.version;
    const o = 3 * (v % LOG);
    this.log[o] = v;
    if (r.level >= D1) {
      // A 1-D region: its extent, or what its items reach (an item drawn past it), in its system's positions; logged as
      // the hull of the block indices it reaches on the blocks that show the system (through a path: `toBlock`).
      const span = this.spanOf(((r.level - D1) % LV) - 1);
      const a = r.lo >= 0 ? r.lo : r.tx * span;
      const e = r.lo >= 0 ? r.hi : (r.tx + 1) * span;
      const s = ((r.level - D1) / LV) | 0;
      const sys = this.source.systems;
      let lo = 0x7fffffff;
      let hi = -1;
      for (let b = 0; sys !== null && b < sys.length; b++) {
        if (sys[b] !== s || !this.toBlock(b, a, e, this.hull)) continue;
        lo = Math.min(lo, this.hull[0]);
        hi = Math.max(hi, this.hull[1]);
      }
      this.log[o + 1] = hi < lo ? 0 : lo;
      this.log[o + 2] = hi < lo ? 0 : hi;
      return;
    }
    const labels = r.level >= LABELS;
    this.log[o + 1] = labels ? -1 : (r.tx * TILE) << (r.level & 15);
    this.log[o + 2] = labels ? -1 : ((r.tx + 1) * TILE) << (r.level & 15);
  }
  changed(since: number, out: Int32Array): number {
    if (since < 0 || this.version - since > LOG || 2 * (this.version - since) > out.length) return -1;
    for (let v = since + 1; v <= this.version; v++) {
      const o = 3 * (v % LOG);
      if (this.log[o] !== v || this.log[o + 1] < 0) return -1;
      out[2 * (v - since - 1)] = this.log[o + 1];
      out[2 * (v - since) - 1] = this.log[o + 2];
    }
    return this.version - since;
  }
  add(r: Region) {
    let m = this.levels[r.level];
    if (m === undefined) this.levels[r.level] = m = new Map();
    m.set(r.index, r);
    r.at = this.all.length;
    this.all.push(r);
  }
  remove(r: Region) {
    this.levels[r.level]?.delete(r.index);
    const last = this.all.pop() as Region;
    if (last !== r) {
      this.all[r.at] = last;
      last.at = r.at;
    }
    r.at = -1;
    this.bytes -= r.bytes;
    if (this.job === r) this.job = null;
  }
}
