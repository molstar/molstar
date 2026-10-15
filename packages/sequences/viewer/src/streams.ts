/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LY4: `tileLevels` and `reduce` ported (`reduce` in one pass, by sample rows), the cache and the manager rewritten.
// ST1: the port (`fetch(request, signal)` and data patches) in place of the tile client, 1-D regions per coordinate system
// and level, the margin, `state` and `settled`, eviction in one sorted sweep.
// The stream manager's region cache (seqspec:sources#the-viewers-side; seqspec:sources#levels-and-the-zoom-gate): what
// drawing needs of a layer whose data is held by region, loaded by view, a data layer under the document (patches never
// reach it). An alignment's region is one tile of a level: [ly, lx] (2^ly rows × 2^lx columns a sample), tile (ty, tx),
// in the source's coordinates (rows in lane order, columns the alignment's), in `Tiles`' addressing; a 1-D layer's
// (letters, items, values) is `TILE` bins of a level of its manifest, or `DATA_SPAN` positions of its data, in one
// coordinate system of the manifest, the one a block shows. Its layer's description (`renderer.stream(s)`) names a
// streamed source (its provider, the port) or an alignment in memory, reduced in place in the same wire format, a
// sample row at a time under the frame's budget, one tile in progress per layer. Per pass (one a frame, from viewport
// and layout versions, never from draws, so it runs headlessly): the regions in view per block-system run of every
// block, coarse before fine from the center out, then a margin of `Stream.Margin` view widths (1-D layers), each asked
// once, what an older view still lacks cancelled, at most `Stream.MaxInFlight` in flight. Past `Stream.CacheBytes` the
// least recently needed regions are evicted, the pass's pinned, down to 7/8 of it in one sorted sweep. A failure
// retries at 1, 2, 4 s, then cools down
// 30 s; each outage is reported once; `StreamNotServed` is never retried. Reads (`StreamView`) and a pass over a cached
// view allocate nothing. ST1b-a: a patch's dynamic nodes are held beside its region (`DynamicNodes`, carried by it,
// counted in the cache) and reach the set through `#assemble.stream`, one graph update a frame. ST3a: a letters
// stream's data held as codes per region, asked only where a device px holds at most `LETTERS_PER_PX` positions, and
// read by the letters draw (`StreamView.letters`). ST3a-2: an items stream's regions (65,536 positions of data) held as
// feature tables built at arrival, lanes the layer's by item, or a level's `count` bins (`StreamView.items`); a 1-D
// region's change logs its columns; `stateAt` for the draws' placeholders; the letters gate `Stream.LettersPerPx`.
// ST3b: per region the frame it was first drawn in (`region`, `tileAt`: the overlay's skeletons, error bands and
// fade-in), the patch's table and an item number per feature (hover, labels, field themes), lanes held by the regions
// that hold their items (forgotten with the last), a miniature's overview tiles without an overview column. ST4b: the
// zoom gate: an items data region answered with a level's bins (`gate`) is drawn as density, and where a device px
// holds at most `FINE_PER_PX` positions its part in view is asked again in regions of `DATA_SPAN` (`FINE`). ST7: a 1-D
// layer placed on a block through a path (a protein's block over a genome's data, by its codon map) asks the regions its
// block's runs reach back through the path (its exons, never its introns), and its reads take block indices back
// through the path and give extents in block indices; its items are placed on the block at arrival.

import {
  type DataPatch,
  type StreamedData,
  type StreamNode,
  StreamNotServed,
  StreamProp,
  type StreamRequest,
  type StreamTable,
  streamRequest,
  TILE,
  type Tile,
  Tiles,
} from '@molstar/seq-core/format';
import { Paths } from '@molstar/seq-core/model';
import type { AxisView, ShownIndex } from './axis.js';
import { type ConfigManager, SeqViewerConfig } from './config.js';
import { DynamicNodes, heldNodeRefs } from './dynamic.js';
import type { ViewModel } from './graph.js';
import { itemSystemNodes } from './item-systems.js';
import { type Layouts, type RowModel, type SectionLayout, SectionLayouts } from './layout.js';
import type { Hint } from './prefetch.js';
import type { SeqViewerRegistry } from './registry.js';
import type { Scheduler } from './scheduler.js';
import { live, state, type Unsubscribable } from './state.js';
import {
  type LocalCells,
  RegionState,
  type StreamManager,
  type StreamStats,
  type StreamView,
} from './stream-contract.js';
import type { Layer, Region } from './stream-layer.js';
import * as Store from './stream-layer.js';
import * as Levels from './stream-levels.js';
import { held1, tileOf } from './stream-patches.js';
import type { Viewports } from './viewport.js';

export {
  LETTERS_PER_PX,
  type LocalCells,
  RegionState,
  type StreamItems,
  type StreamLayer,
  type StreamManager,
  type StreamStats,
  type StreamView,
} from './stream-contract.js';
export { blockLevels, listedLevels, MAP_WIDTH, overviewLevels, reduce, tileLevels } from './stream-levels.js';
export { lettersStream, streamedBlock, streamStats, trackStream } from './stream-sources.js';

// Read once into module constants: a read through an import is a property read of the module per use, which
// the test runner's module wrapper turns into an allocation on per-frame paths.
const {
  DATA_SPAN,
  FINE,
  FINE_PER_PX,
  indexOf,
  letterIn,
  MAX_REGIONS,
  NO_HINTS,
  names,
  regionsOver,
  routeOf,
  routes,
  sortRuns,
} = Store;
const { blockLevels, D1, L, LABELS, LV, overviewLevels, reduce } = Levels;
const PENDING = RegionState.pending;
const LOADED = RegionState.loaded;
const EMPTY = RegionState.empty;
const FAILED = RegionState.failed;

const PRIORITY = ['view', 'view', 'view', 'view', 'margin', 'prefetch'] as const;
/** The groups from this one on are never waited for (`settled`): the margin's, and the prefetch's. */
const MARGIN = 4;
const PREFETCH = 5;
/** A prefetch's lead past the margin, in view widths. */
const LEAD = 1;
/** A frame's local work: at most this many cells or ms, at least one sample row. */
const LOCAL_CELLS = 1 << 21;
const LOCAL_MS = 6;
const RETRIES = 3;
const COOLDOWN_MS = 30_000;

/** A display row's lane offset, or the next row's in `dir` past wedges. */
function laneAt(r: RowModel, d: number, dir: number): number {
  for (let k = d; k >= 0 && k < r.count; k += dir) {
    const lane = SectionLayouts.laneOf(r, k);
    if (lane >= 0) return lane - r.laneStart;
  }
  return dir > 0 ? r.laneCount : -1;
}

/**
 * A rows layer's lanes in view, the manager's and the draws' alike: into `out` [k0, k1) lane offsets, display rows
 * [d0, d1], the surface y of display row 0 and a display row's height (CSS px). False: none.
 */
export function lanesInView(l: SectionLayout, v: AxisView, out: Float64Array): boolean {
  const r = l.rows;
  if (r.layer < 0 || r.count === 0 || r.height <= 0) return false;
  const top = l.scrolledTop + l.entryTop[r.entry] - v.scrollY;
  const d0 = Math.max(0, Math.floor((l.scrolledTop - top) / r.height));
  const d1 = Math.min(r.count - 1, Math.floor((l.height - top) / r.height));
  if (d1 < d0) return false;
  out[0] = laneAt(r, d0, 1);
  out[1] = laneAt(r, d1, -1) + 1;
  out[2] = d0;
  out[3] = d1;
  out[4] = top;
  out[5] = r.height;
  return out[1] > out[0];
}

export interface StreamHost {
  readonly registry: SeqViewerRegistry;
  readonly config: ConfigManager;
  readonly scheduler: Scheduler;
  readonly layout: Layouts;
  readonly viewports: Viewports;
  active(): boolean;
  /** The shown document still holds this layer (of another view: paused, not closed). */
  has(ref: string): boolean;
  /** The person expanded this track's rows (`UI.expanded`): its folded items show every member. */
  expanded?(track: string): boolean;
  report(code: string, message: string): void;
  /** The dynamic nodes changed (at most once a frame): `#assemble.stream` is to run. */
  nodesChanged?(): void;
  /** Adds the keys the selection holds: a node it names stays. */
  held?(out: Set<string>): void;
}

/** A queue entry: the layer's index in `order`, level, ty, tx, group. */
const Q = 5;
export class Streams implements StreamManager {
  readonly version = state(0);
  readonly stats = live<StreamStats>({
    requests: 0,
    hits: 0,
    inFlight: 0,
    cancelled: 0,
    coarse: 0,
    fine: 0,
    bytes: 0,
    errors: 0,
    computed: 0,
  });
  private readonly host: StreamHost;
  private readonly clients = new Map<string, Layer>();
  /** Every open layer, and the current view's in its order. */
  private readonly open: Layer[] = [];
  /** `lettersAt`'s layers, by stream, and in a list (a pass reads it without an iterator). */
  private readonly aside = new Map<StreamedData, Layer>();
  private readonly asides: Layer[] = [];
  /** Item map faults reported, by source and code (`itemSystems`). */
  private readonly reported = new Set<string>();
  private readonly order: Layer[] = [];
  private readonly flying: Region[] = [];
  private model: ViewModel | null = null;
  /** The snapshots' focus regions to prefetch (`hint`). */
  private hints: readonly Hint[] = [];
  private subs: Unsubscribable[] = [];
  private frame = -1;
  private wake = -1;
  private wakeDue = Number.POSITIVE_INFINITY;
  private clock = 0;
  private now = 0;
  private bytes = 0;
  private disposed = false;
  private queue = new Int32Array(64 * Q);
  private prio = new Float64Array(64);
  private queued = 0;
  /** This frame's local work: its start and the cells read. */
  private t0 = 0;
  private cells = 0;
  private readonly lv = new Int32Array(4);
  /** Each pass's time by its clock mod 64, for the grace on cancellation (`graced`). */
  private readonly passAt = new Float64Array(64);
  private grace = 0;
  private readonly lanes = new Float64Array(6);
  private readonly center = new Float64Array(2);
  private runs = new Int32Array(64);
  /** `routesOf`'s runs: `runs`, or the module's `routes`. */
  private routed = this.runs;
  /** This pass's regions in view not yet loaded or failed; `settled()`'s waiters; eviction's sweep. */
  private missing = 0;
  private stale = true;
  private waiters: (() => void)[] = [];
  private readonly sweep: Region[] = [];
  /** The dynamic nodes held, and the frame that tells the graph they changed. */
  readonly nodes: DynamicNodes;
  private nodeFrame = -1;
  private readonly onNodes = () => {
    this.nodeFrame = -1;
    this.host.nodesChanged?.();
  };

  constructor(host: StreamHost) {
    this.host = host;
    this.nodes = new DynamicNodes(
      () => {
        if (this.nodeFrame < 0 && !this.disposed) this.nodeFrame = host.scheduler.requestFrame(this.onNodes);
      },
      (code, message) => host.report(code, message),
      heldNodeRefs(host.registry.svs),
    );
  }

  get(layer: string): StreamView | null {
    return this.clients.get(layer) ?? null;
  }

  state(layer: string, block: number, index: number, lane = -1): number {
    const c = this.clients.get(layer);
    if (c === undefined) return PENDING;
    const sys = c.source.systems;
    if (sys === null) {
      // 2-D: the tile row of `lane` (the layer's own, from 0) at the view's fine level, else the view's first.
      const ty = lane >= 0 ? Math.floor(lane / (TILE * 2 ** c.fine[0])) : c.fine[2];
      return c.state(c.fine[0], c.fine[1], ty, Math.floor(index / (TILE * 2 ** c.fine[1])));
    }
    // Letters zoomed out past the gate: nothing is drawn there, nothing needed (empty).
    return c.stateAt(block, index);
  }

  settled(): Promise<void> {
    if (this.calm()) return Promise.resolve();
    return new Promise((ok) => {
      this.waiters.push(ok);
      this.onChange();
    });
  }
  /** Nothing to wait for: no pass can run, or the last pass found every region in view loaded or failed. */
  private calm(): boolean {
    return this.disposed || !this.host.active() || this.order.length === 0 || (!this.stale && this.missing === 0);
  }
  private release() {
    // Checked first: `splice` makes an array even when there is nothing to take (a pass's 32 bytes).
    if (this.waiters.length > 0) for (const ok of this.waiters.splice(0)) ok();
  }

  hasFailed(layer: string): boolean {
    return (this.clients.get(layer)?.failing ?? 0) > 0;
  }

  force(layer: string, block: number, index: number): boolean {
    const c = this.clients.get(layer);
    const sys = c?.source.systems;
    if (!c?.isItems || !sys || block >= sys.length) return false;
    const s = sys[block];
    // Through a path: the first position of the system that reaches the index.
    const at = c.back(block, index, index + 1) ? c.hull[0] : -1;
    const r = s >= 0 && at >= 0 ? c.find(D1 + s * LV, Math.floor(at / c.span)) : undefined;
    if (r === undefined || r.gate < 0 || r.force) return false;
    // Its part in view asked in finer regions at the next pass, forced past the gate (not the whole 65 kb region).
    r.force = true;
    (c as Layer).bump(r);
    this.onChange();
    return true;
  }
  /** `setActive(true)`: a pass in the next frame, since the passes while inactive asked for nothing. */
  resume() {
    this.onChange();
  }

  retry(layer: string): number {
    const c = this.clients.get(layer);
    let n = 0;
    // A label's letters (`lettersAt`) too: they have no track of their own.
    for (const x of [c, ...this.asides])
      for (let k = 0; x && k < x.all.length; k++) {
        const r = x.all[k];
        if (r.state !== FAILED || r.retryAt === Number.POSITIVE_INFINITY) continue;
        r.tries = 0;
        r.retryAt = 0;
        n++;
      }
    if (n > 0) this.onChange();
    return n;
  }

  /** The graph's new view: layers for its streamed layers; others pause (still in the document) or close. Edge-rate. */
  update(model: ViewModel | null) {
    this.model = model;
    for (const s of this.subs) s.unsubscribe();
    this.subs = [];
    this.order.length = 0;
    for (const c of this.open) c.layer = -1;
    const rv = model?.resolved;
    for (let k = 0; model && rv && k < rv.layers.count; k++) {
      const type = rv.layers.type[k];
      const r = type >= 0 ? this.host.registry.get('renderers', model.set.strings.get(type)) : undefined;
      const src = r?.stream?.(model.prepared[k]?.state) ?? null;
      if (!src) continue;
      const ref = model.set.strings.get(rv.layers.ref[k]);
      let c = this.clients.get(ref);
      if (c && c.source.key !== src.key) this.close(c);
      c = this.clients.get(ref);
      if (!c) {
        c = new Store.Layer(ref, src);
        this.clients.set(ref, c);
        this.open.push(c);
      }
      const listed = c.source.remote?.manifest !== src.remote?.manifest;
      c.source = src;
      if (listed) c.listed();
      c.routed();
      c.layer = k;
      c.track = model.set.strings.get(rv.tracks.ref[rv.layers.track[k]]);
      c.isExpanded ??= this.expandedOf(c);
      this.order.push(c);
      const sec = this.host.layout.sectionRef(rv.tracks.section[rv.layers.track[k]]);
      const vp = sec ? this.host.viewports.get(sec) : null;
      const l = sec ? this.host.layout.get(sec) : null;
      if (vp) {
        // Its motion from here on, the view it has now the one to compare with (the first change is a pan too).
        this.moved(c, vp);
        this.subs.push(
          vp.live.version.subscribe(() => {
            this.moved(c, vp);
            this.onChange();
          }),
        );
      }
      if (l) this.subs.push(l.version.subscribe(this.onChange));
    }
    // A label's letters layers of streams no lazy molecule of the set reads any longer (another document) close; then
    // the layers no view shows. Each closes once, outside the walk over `open` it shrinks.
    const live = new Set<unknown>();
    for (const m of model?.set.molecules ?? []) if (m?.lazy) live.add(m.custom.props.get(StreamProp));
    for (const c of [...this.asides]) if (!live.has(c.source.remote)) this.close(c);
    for (let j = this.open.length - 1; j >= 0; j--) {
      const c = this.open[j];
      if (c.layer < 0 && !c.aside && !this.host.has(c.ref)) this.close(c);
    }
    this.rehold();
    this.onChange();
  }

  /** Whether the person expanded layer c's track (made once a layer). */
  private expandedOf(c: Layer): () => boolean {
    return () => this.host.expanded?.(c.track) ?? false;
  }

  lettersAt(data: StreamedData, s: number, at: Int32Array, n: number, out: Uint8Array): number {
    // This read's positions replace what an earlier one wanted (a hover moved on: its regions are cancelled).
    let c = this.aside.get(data);
    if (c !== undefined) c.wantedN = 0;
    let held = 0;
    let lost = 0;
    for (let k = 0; k < n; k++) {
      const x = at[k];
      out[k] = 0;
      if (x < 0) continue;
      // A view's letters layer of the same stream first: its regions are not asked twice.
      out[k] = this.viewLetter(data, s, x);
      if (out[k] > 0) {
        held++;
        continue;
      }
      if (c === undefined) {
        const systems = Int32Array.of(s);
        c = new Store.Layer(`${data.source.ref}#letters`, {
          key: '',
          rows: 0,
          columns: 0,
          remote: data,
          local: null,
          systems,
        });
        c.aside = true;
        this.aside.set(data, c);
        this.asides.push(c);
        this.open.push(c);
      }
      const j = Math.floor(x / c.span) | 0;
      const r = c.find(D1 + s * LV, j);
      const v = letterIn(r, x - j * c.span);
      if (v > 0) {
        out[k] = v;
        held++;
      } else if (r === undefined || r.state === PENDING) this.wantOf(c, s, j);
      else if (r.state === FAILED && r.retryAt !== Number.POSITIVE_INFINITY) {
        // Failed: asked again after its backoff (or a Retry); said meanwhile.
        this.wantOf(c, s, j);
        if (lost === 0) lost = -1;
      } else lost = -2;
    }
    return lost < 0 ? lost : held;
  }
  /** Position x of system s from a view's letters layer of stream `data`; 0: not held there. */
  private viewLetter(data: StreamedData, s: number, x: number): number {
    for (let q = 0; q < this.order.length; q++) {
      const c = this.order[q];
      if (!c.seq || c.source.remote !== data) continue;
      const j = Math.floor(x / c.span) | 0;
      const v = letterIn(c.find(D1 + s * LV, j), x - j * c.span);
      if (v > 0) return v;
    }
    return 0;
  }
  /** Region j of system s wanted by `lettersAt` (once). */
  private wantOf(c: Layer, s: number, j: number) {
    for (let q = 0; q < c.wantedN; q++) if (c.wanted[2 * q] === s && c.wanted[2 * q + 1] === j) return;
    if (2 * c.wantedN + 2 > c.wanted.length) {
      const g = new Int32Array(2 * c.wanted.length);
      g.set(c.wanted);
      c.wanted = g;
    }
    c.wanted[2 * c.wantedN++] = s;
    c.wanted[2 * c.wantedN - 1] = j;
    this.onChange();
  }
  /** An aside layer's wanted regions asked (this pass's order), those held dropped from the list. */
  private wantLetters(c: Layer) {
    const i = this.order.length;
    this.order.push(c);
    let kept = 0;
    for (let q = 0; q < c.wantedN; q++) {
      const s = c.wanted[2 * q];
      const j = c.wanted[2 * q + 1];
      const r = c.find(D1 + s * LV, j);
      if (r !== undefined && (r.state === LOADED || r.state === EMPTY)) continue;
      // The view's group (`settled` waits: a label's two regions at most), after every region in view: as far as 10^4
      // regions from the center.
      this.center[0] = j - 1e4;
      this.want(i, c, D1 + s * LV, 0, j, j, 1);
      c.wanted[2 * kept] = s;
      c.wanted[2 * kept++ + 1] = j;
    }
    c.wantedN = kept;
  }

  /** The holders changed (the view, the selection): nodes nobody holds any longer are evicted. Edge-rate. */
  rehold() {
    if (this.nodes.entries.size === 0) return;
    const keys = shownKeys(this.model);
    this.host.held?.(keys);
    this.nodes.show(keys);
  }

  private readonly onChange = () => {
    this.stale = true;
    if (this.frame < 0 && !this.disposed && this.order.length > 0)
      this.frame = this.host.scheduler.requestFrame(this.onFrame);
  };
  private readonly onFrame = () => {
    this.frame = -1;
    this.pass();
  };
  private readonly onWake = () => {
    this.wake = -1;
    this.wakeDue = Number.POSITIVE_INFINITY;
    this.onChange();
  };

  /** One pass now (a frame runs it; tests call it): what the views need, older requests cancelled, requests started. */
  pass() {
    if (this.disposed || !this.host.active() || !this.model) {
      this.release();
      return;
    }
    this.queued = 0;
    this.missing = 0;
    this.clock++;
    this.now = this.host.scheduler.now();
    this.passAt[this.clock & 63] = this.now;
    this.grace = this.host.config.get(SeqViewerConfig.Stream.CancelGraceMs);
    for (let i = 0; i < this.order.length; i++)
      if (this.order[i].source.systems !== null) this.need1(i);
      else this.need(i);
    // The letters a label asked for (`lettersAt`), the view's priority, in the order for this pass only.
    const n = this.order.length;
    for (let q = 0; q < this.asides.length; q++) if (this.asides[q].wantedN > 0) this.wantLetters(this.asides[q]);
    this.cancel();
    this.dispatch();
    this.order.length = n;
    this.stale = this.frame >= 0;
    if (this.missing === 0) this.release();
  }

  /**
   * 1-D layer i's regions: per block in view, its system's coarsest level where it is not the view's, the view's level
   * (the coarsest listed whose bin is at most the positions a device px; finer: the data), then the margin.
   */
  private need1(i: number) {
    const c = this.order[i];
    const rv = (this.model as ViewModel).resolved;
    const ref = this.host.layout.sectionRef(rv.tracks.section[rv.layers.track[c.layer]]);
    const v = ref ? this.host.layout.axis(ref) : null;
    const ax = v?.axis;
    if (!v || !ax) return;
    const sys = c.source.systems as Int32Array;
    const W = c.widths;
    const top = W.length - 1;
    const margin = this.host.config.get(SeqViewerConfig.Stream.Margin);
    c.gate = this.host.config.get(SeqViewerConfig.Stream.LettersPerPx);
    const ahead = this.host.config.get(SeqViewerConfig.Stream.Prefetch);
    if (c.chosen.length !== sys.length) c.chosen = new Int16Array(sys.length);
    for (let b = 0; b < sys.length && b < ax.count; b++) {
      const s = sys[b];
      const lo = v.visLo[b];
      const hi = v.visHi[b];
      c.chosen[b] = -2;
      if (s < 0 || hi <= lo) continue;
      // Its system's positions a device px (through a path, a residue's codon is 3).
      const per = c.ratio[b] / (v.cell * v.span[b] * v.dpr);
      let lvl = -1;
      for (let k = 0; k < W.length && !c.seq; k++) if (W[k] <= per) lvl = k;
      // Items: a coarser level where the view would take more than `MAX_REGIONS` of this one's regions (a file's
      // levels can lie far apart: a BigBed's counted 256, then its zoom records; ST4b), never their items.
      // Through a path the regions its runs in view reach are counted, the data's too (a gene's exons over 2 Mb).
      const routed = c.pathOf(b) !== null;
      if (c.isItems && !routed) while (lvl >= 0 && lvl < top && (hi - lo) / (TILE * W[lvl]) > MAX_REGIONS) lvl++;
      else if (c.isItems) while (lvl < top && this.regionsIn(c, v, b, lo, hi, c.spanOf(lvl)) > MAX_REGIONS) lvl++;
      // Letters: the data only, and only where a device px holds at most `LETTERS_PER_PX` positions.
      c.chosen[b] = c.seq && per > c.gate ? -3 : lvl;
      if (c.chosen[b] === -3) continue;
      // The center in the system's positions: through a path the center residue's first, else the middle of what the
      // view reaches; nothing reached, nothing to ask.
      const at = v.indexOf(b, ((lo + hi) >> 1) as ShownIndex);
      let mid = at;
      if (routed && c.back(b, at, at + 1)) mid = c.hull[0];
      else if (routed && this.routesOf(c, v, b, lo, hi) > 0) mid = this.routed[0];
      else if (routed) continue;
      // Values: the coarsest first, drawn while finer loads; items show placeholders instead.
      if (lvl !== top && top >= 0 && !c.seq && !c.isItems) this.runs1(i, v, b, s, top, lo, hi, mid, 0);
      this.runs1(i, v, b, s, lvl, lo, hi, mid, 1);
      // The zoom gate: inside a data region answered with bins, its part in view again in finer regions.
      if (c.isItems && lvl < 0) this.fine1(i, v, b, s, lo, hi, mid, per <= FINE_PER_PX);
      const w = Math.round(margin * (hi - lo));
      if (w <= 0) continue;
      if (lo > 0) this.runs1(i, v, b, s, lvl, Math.max(0, lo - w), lo, mid, 4);
      if (hi < ax.shown[b]) this.runs1(i, v, b, s, lvl, hi, Math.min(ax.shown[b], hi + w), mid, 4);
      if (ahead) this.ahead1(i, v, b, s, lvl, lo, hi, w, mid);
    }
    this.map1(i, v, this.mapped(c, true));
  }

  /**
   * A miniature draws 1-D layer i (`mapped`): every block's whole shown range at the level whose bin is about a map
   * device px, at most `MAX_REGIONS` regions a block (else coarser, up to the coarsest), in the overview's group (after
   * the view's, before the margin). Else none: the map's levels go.
   */
  private map1(i: number, v: AxisView, on: boolean) {
    const c = this.order[i];
    const ax = v.axis as NonNullable<AxisView['axis']>;
    const sys = c.source.systems as Int32Array;
    if (c.mapChosen.length !== sys.length) {
      c.mapChosen = new Int16Array(sys.length);
      c.mapBase = new Int16Array(sys.length);
    }
    c.mapChosen.fill(-2);
    c.mapBase.fill(-2);
    if (!on || c.seq || v.navCell <= 0) return;
    const W = c.widths;
    const top = W.length - 1;
    for (let b = 0; b < sys.length && b < ax.count; b++) {
      const s = sys[b];
      const shown = ax.shown[b];
      if (s < 0 || shown <= 0) continue;
      const per = c.ratio[b] / (v.navCell * v.span[b] * v.dpr);
      let lvl = -1;
      for (let k = 0; k < W.length; k++) if (W[k] <= per) lvl = k;
      c.mapBase[b] = lvl;
      const routed = c.pathOf(b) !== null;
      if (!routed) while (lvl < top && shown / c.spanOf(lvl) > MAX_REGIONS) lvl++;
      else while (lvl < top && this.regionsIn(c, v, b, 0, shown, c.spanOf(lvl)) > MAX_REGIONS) lvl++;
      c.mapChosen[b] = lvl;
      const at = v.indexOf(b, (shown >> 1) as ShownIndex);
      let mid = at;
      if (routed && c.back(b, at, at + 1)) mid = c.hull[0];
      else if (routed && this.routesOf(c, v, b, 0, shown) > 0) mid = this.routed[0];
      else if (routed) continue;
      this.runs1(i, v, b, s, lvl, 0, shown, mid, 2);
    }
  }

  /**
   * Prefetch (ST6), after the view's and the margin's: past the margin, a view width (`LEAD`) on the side the view last
   * panned toward; a zoom's next level (finer where it zoomed in, coarser out) over the view, where that is at most
   * `MAX_REGIONS` of its regions; and the regions the snapshots next to the one shown focus on (`hint`), at the level a
   * view of them would need. The direction is kept until the view moves the other way, so the lead is not dropped (and
   * aborted) when the view stops.
   */
  private ahead1(
    i: number,
    v: AxisView,
    b: number,
    s: number,
    lvl: number,
    lo: number,
    hi: number,
    w: number,
    mid: number,
  ) {
    const c = this.order[i];
    const ax = v.axis as NonNullable<AxisView['axis']>;
    const dir = c.moves[2];
    const zoom = c.moves[3];
    const shown = ax.shown[b];
    if (dir > 0 && hi + w < shown)
      this.runs1(i, v, b, s, lvl, hi + w, Math.min(shown, hi + w + LEAD * (hi - lo)), mid, PREFETCH);
    else if (dir < 0 && lo - w > 0)
      this.runs1(i, v, b, s, lvl, Math.max(0, lo - w - LEAD * (hi - lo)), lo - w, mid, PREFETCH);
    // The next level a zoom would need: only where the layer has levels to move between.
    const next = zoom < 0 ? lvl - 1 : zoom > 0 ? lvl + 1 : lvl;
    const fits =
      next !== lvl &&
      next >= -1 &&
      next < c.widths.length &&
      !c.seq &&
      (c.pathOf(b) === null
        ? (hi - lo) / c.spanOf(next) <= MAX_REGIONS
        : this.regionsIn(c, v, b, lo, hi, c.spanOf(next)) <= MAX_REGIONS);
    if (fits)
      this.runs1(
        i,
        v,
        b,
        s,
        next,
        Math.max(0, lo - (zoom > 0 ? w : 0)),
        Math.min(shown, hi + (zoom > 0 ? w : 0)),
        mid,
        PREFETCH,
      );
    // The snapshots' focus regions, in the data's positions, where this block shows the hint's target (directly).
    const hints = c.pathOf(b) === null ? this.hints : NO_HINTS;
    const name = ax.names[b];
    for (let k = 0; k < hints.length; k++) {
      const h = hints[k];
      if (!names(name, h.target)) continue;
      // The level a view of the hint's range across this view's width needs (device px: this view's positions over those a
      // device px holds).
      const px = (hi - lo) * v.cell * v.span[b] * v.dpr;
      let at = -1;
      for (let q = 0; q < c.widths.length && !c.seq; q++) if (c.widths[q] <= (h.stop - h.start) / px) at = q;
      if (c.seq && (h.stop - h.start) / px > c.gate) continue;
      if ((h.stop - h.start) / c.spanOf(at) > MAX_REGIONS) continue;
      const ext = (c.source.remote as StreamedData).manifest.systems[s].extent;
      this.span1(i, c, s, at, Math.max(ext[0], h.start), Math.min(ext[1], h.stop), (h.start + h.stop) / 2, PREFETCH);
    }
  }

  /** A change of layer c's section's view: the direction it panned or zoomed, kept until it moves the other way. */
  private moved(
    c: Layer,
    vp: { readonly live: { readonly current: { readonly scrollX: number; readonly cell: number } } },
  ) {
    const m = c.moves;
    const { scrollX, cell } = vp.live.current;
    // A cell that changed is a zoom (it moves the scroll too); else a scroll that changed is a pan.
    if (m[1] > 0 && cell !== m[1]) {
      m[2] = 0;
      m[3] = cell > m[1] ? -1 : 1;
    } else if (m[1] > 0 && scrollX !== m[0]) {
      m[2] = scrollX > m[0] ? 1 : -1;
      m[3] = 0;
    }
    m[0] = scrollX;
    m[1] = cell;
  }

  /** The prefetch hints (`snapshotHints`): replaced as a whole when the document or the snapshot changes. */
  hint(hints: readonly Hint[]) {
    this.hints = hints;
    this.onChange();
  }
  /** System s's regions at level `lvl` over block b's shown [lo, hi), by distance from index `cx`. */
  private runs1(
    i: number,
    v: AxisView,
    b: number,
    s: number,
    lvl: number,
    lo: number,
    hi: number,
    cx: number,
    g: number,
  ) {
    const c = this.order[i];
    const ext = (c.source.remote as StreamedData).manifest.systems[s].extent;
    const n = this.routesOf(c, v, b, lo, hi);
    const runs = this.routed;
    for (let j = 0; j < n; j++)
      this.span1(i, c, s, lvl, Math.max(ext[0], runs[2 * j]), Math.min(ext[1], runs[2 * j + 1]), cx, g);
  }
  /** System s's regions at level `lvl` over the positions [c0, c1), ordered by distance from position `cx`. */
  private span1(i: number, c: Layer, s: number, lvl: number, c0: number, c1: number, cx: number, g: number) {
    const span = c.spanOf(lvl);
    const level = D1 + s * LV + lvl + 1;
    this.center[0] = cx / span - 0.5;
    for (let tx = Math.floor(c0 / span); tx <= Math.floor((c1 - 1) / span) && c1 > c0; tx++)
      this.want(i, c, level, 0, tx, tx, g);
  }

  /** System s's finer data regions over block b's shown [lo, hi) inside its gated data regions (the zoom gate). */
  private fine1(i: number, v: AxisView, b: number, s: number, lo: number, hi: number, cx: number, near: boolean) {
    const c = this.order[i];
    const ext = (c.source.remote as StreamedData).manifest.systems[s].extent;
    const level = D1 + s * LV + FINE + 1;
    this.center[0] = cx / DATA_SPAN - 0.5;
    const n = this.routesOf(c, v, b, lo, hi);
    const runs = this.routed;
    for (let j = 0; j < n; j++) {
      const c0 = Math.max(ext[0], runs[2 * j]);
      const c1 = Math.min(ext[1], runs[2 * j + 1]);
      for (let tx = Math.floor(c0 / c.span); tx * c.span < c1; tx++) {
        const r = c.find(D1 + s * LV, tx);
        if (r === undefined || r.state !== LOADED || r.gate < 0 || !(near || r.force)) continue;
        const e = Math.min(c1, (tx + 1) * c.span);
        for (let q = Math.floor(Math.max(c0, tx * c.span) / DATA_SPAN); q * DATA_SPAN < e; q++)
          this.want(i, c, level, 0, q, q, 1);
      }
    }
  }

  /** Layer i's regions for its section's view, every block's runs, and the vertical overview; its labels in view. */
  private need(i: number) {
    const c = this.order[i];
    const rv = (this.model as ViewModel).resolved;
    const ref = this.host.layout.sectionRef(rv.tracks.section[rv.layers.track[c.layer]]);
    const live = ref ? this.host.layout.get(ref) : null;
    const v = ref ? this.host.layout.axis(ref) : null;
    const ax = v?.axis;
    if (!live || !v || !ax) return;
    const l = live.current;
    const r = l.rows;
    if (r.layer !== c.layer || r.count === 0 || r.height <= 0) return;
    const src = c.source;
    const lv = this.lv;
    const g = this.lanes;
    const budget = this.host.config.get(SeqViewerConfig.Stream.SampleBudget) ?? 0;
    const margin = this.host.config.get(SeqViewerConfig.Stream.Margin);
    c.budget = budget;
    if (lanesInView(l, v, g)) {
      for (let b = 0; b < ax.count; b++) {
        const lo = v.visLo[b];
        const hi = v.visHi[b];
        if (hi <= lo) continue;
        blockLevels(v, b, budget, g, lv, c.served);
        if (src.remote === null && lv[0] === 0 && lv[1] === 0) continue;
        const n = this.runsOf(v, b, lo, hi);
        const mid = v.indexOf(b, ((lo + hi) >> 1) as ShownIndex);
        for (let j = 0; j < n; j++) {
          const c0 = this.runs[2 * j];
          const c1 = this.runs[2 * j + 1];
          if (lv[2] !== lv[0] || lv[3] !== lv[1]) this.over(i, lv[2], lv[3], c0, c1, mid, 0);
          this.over(i, lv[0], lv[1], c0, c1, mid, 1);
          c.fine[0] = lv[0];
          c.fine[1] = lv[1];
          c.fine[2] = Math.floor(g[0] / (TILE * 2 ** lv[0]));
        }
        if (src.remote !== null && margin > 0) this.ring(i, v, b, lv, lo, hi, mid);
        if (src.remote !== null && this.host.config.get(SeqViewerConfig.Stream.Prefetch))
          this.lead(i, v, b, lv, lo, hi, mid);
      }
      if (src.remote !== null)
        for (let ch = Math.floor(g[0] / TILE); ch * TILE < g[1]; ch++) this.want(i, c, LABELS, ch, 0, ch, 3);
    }
    if ((l.vnavWidth > 0 || this.mapped(c)) && src.columns > 0) {
      overviewLevels(v, r.laneCount, src.columns, budget, lv, c.served);
      g[0] = 0;
      g[1] = r.laneCount;
      if (src.remote !== null || lv[2] > 0 || lv[3] > 0)
        this.over(i, lv[2], lv[3], 0, src.columns, src.columns >> 1, 2);
    }
  }

  /**
   * A miniature names layer c's rows track (the drawn navigator's `tracks`, a `seqv:map` layer's) where the section has
   * no overview column: its overview's tiles are needed all the same (ST1's routed item). `all` (a 1-D layer): a
   * miniature naming no tracks draws it too, as every layer whose renderer is not `miniature: false` (the map's rule,
   * draw/miniature.ts). Per pass, allocation-free.
   */
  private mapped(c: Layer, all = false): boolean {
    const m = this.model as ViewModel;
    const rv = m.resolved;
    const t = rv.layers.track[c.layer];
    if (all && rv.tracks.hidden[t]) return false;
    const ref = m.set.strings.get(rv.tracks.ref[t]);
    const r = all ? this.host.registry.get('renderers', m.prepared[c.layer]?.renderer ?? '') : undefined;
    const any = all && r?.miniature !== false;
    const p = this.host.layout.navigatorOf(rv.tracks.section[t])?.parts;
    const maps = p && !p.onRuler && (p.axis.includes('miniature') || p.axis.includes('overview'));
    if (maps && (p.tracks?.length ? p.tracks.includes(ref) : any)) return true;
    for (let k = 0; k < rv.layers.count; k++) {
      const s = m.prepared[k]?.state as { tracks?: readonly string[] } | undefined;
      const same = rv.tracks.section[rv.layers.track[k]] === rv.tracks.section[t];
      if (same && m.prepared[k].renderer === 'seqv:map' && (s?.tracks?.length ? s.tracks.includes(ref) : any))
        return true;
    }
    return false;
  }

  /**
   * Block b's margin (`Stream.Margin` view widths and heights each side, the view's tiles included: asked twice is
   * asked once): the tiles of the view's level around it, for the next screen of a pan. A source's only: tiles of an
   * alignment in memory are computed when a view needs them.
   */
  private ring(i: number, v: AxisView, b: number, lv: Int32Array, lo: number, hi: number, mid: number) {
    const g = this.lanes;
    const margin = this.host.config.get(SeqViewerConfig.Stream.Margin);
    const k0 = g[0];
    const k1 = g[1];
    const rows = this.order[i].source.rows;
    const dk = Math.round(margin * (k1 - k0));
    const dc = Math.round(margin * (hi - lo));
    g[0] = Math.max(0, k0 - dk);
    g[1] = Math.min(rows, k1 + dk);
    const n = this.runsOf(v, b, Math.max(0, lo - dc), Math.min(v.axis?.shown[b] ?? hi, hi + dc));
    for (let j = 0; j < n; j++) this.over(i, lv[0], lv[1], this.runs[2 * j], this.runs[2 * j + 1], mid, 4);
    g[0] = k0;
    g[1] = k1;
  }

  /**
   * Block b's prefetch (ST6): past the ring, a view width of columns (`LEAD`) toward the side the view last panned to,
   * the view's tiles of the view's level in the view's lanes (a pan's next screen; the direction kept as the 1-D's).
   */
  private lead(i: number, v: AxisView, b: number, lv: Int32Array, lo: number, hi: number, mid: number) {
    const c = this.order[i];
    const dir = c.moves[2];
    if (dir === 0) return;
    const shown = v.axis?.shown[b] ?? hi;
    const dc = Math.round(this.host.config.get(SeqViewerConfig.Stream.Margin) * (hi - lo));
    const a = dir > 0 ? hi + dc : Math.max(0, lo - dc - LEAD * (hi - lo));
    const e = dir > 0 ? Math.min(shown, hi + dc + LEAD * (hi - lo)) : lo - dc;
    if (e <= a) return;
    const n = this.runsOf(v, b, a, e);
    for (let j = 0; j < n; j++) this.over(i, lv[0], lv[1], this.runs[2 * j], this.runs[2 * j + 1], mid, PREFETCH);
  }

  /**
   * Block b's shown [lo, hi) as runs of layer c's system in `routed`, their count: the block's runs, or through its path
   * the runs of its system that reach them (ST7: a protein's residues in view, their exons). Edge-rate growth.
   */
  private routesOf(c: Layer, v: AxisView, b: number, lo: number, hi: number): number {
    const n = this.runsOf(v, b, lo, hi);
    const p = c.pathOf(b);
    this.routed = this.runs;
    if (p === null) return n;
    routes.n = 0;
    for (let j = 0; j < n; j++) Paths.forEachFrom(p, this.runs[2 * j], this.runs[2 * j + 1], routeOf, routes);
    sortRuns(routes);
    this.routed = routes.runs;
    return routes.n;
  }
  /** The regions of `span` positions block b's shown [lo, hi) reaches through its path (`routesOf`'s runs). */
  private regionsIn(c: Layer, v: AxisView, b: number, lo: number, hi: number, span: number): number {
    const n = this.routesOf(c, v, b, lo, hi);
    return regionsOver(this.routed, n, span);
  }
  /** Block b's shown [lo, hi) as block-system runs into `runs` (the layout map's); their count. Edge-rate growth. */
  private runsOf(v: AxisView, b: number, lo: number, hi: number): number {
    let n = v.shownToRuns(b, lo, hi, this.runs);
    if (2 * n > this.runs.length) {
      this.runs = new Int32Array(4 * n);
      n = v.shownToRuns(b, lo, hi, this.runs);
    }
    return n;
  }

  /** The tiles over columns [c0, c1) × the lanes in `lanes` at a level, by distance from column `cx` and the middle. */
  private over(i: number, ly: number, lx: number, c0: number, c1: number, cx: number, group: number) {
    const c = this.order[i];
    const sx = TILE * 2 ** lx;
    const sy = TILE * 2 ** ly;
    const k0 = this.lanes[0];
    const k1 = this.lanes[1];
    // The center in tiles, for the queue's order (a double passed to a call V8 does not inline is boxed).
    this.center[0] = cx / sx - 0.5;
    this.center[1] = (k0 + k1) / 2 / sy - 0.5;
    const across = Tiles.count(c.source.columns, lx);
    const level = ly * L + lx;
    for (let ty = Math.floor(k0 / sy); ty <= Math.floor((k1 - 1) / sy); ty++)
      for (let tx = Math.floor(c0 / sx); tx <= Math.floor((c1 - 1) / sx); tx++)
        this.want(i, c, level, ty, tx, (ty * across + tx) | 0, group);
  }

  /** A region the pass needs: pinned, counted a hit if held, queued unless held, in flight, computing or waiting. */
  private want(i: number, c: Layer, level: number, ty: number, tx: number, index: number, group: number) {
    const r = c.find(level, index);
    if (r !== undefined) {
      // Once a pass; but one the margin asked first and a view asks after (another block's, the overview) is the
      // view's too: `settled` waits for it.
      const again = r.need === this.clock;
      if (again && (group >= MARGIN || r.viewed === this.clock)) return;
      r.need = this.clock;
      if (group < MARGIN) r.viewed = this.clock;
      if (r.state === LOADED || r.state === EMPTY) {
        if (!again) this.stats.current.hits++;
        return;
      }
      if (r.state === FAILED && this.now < r.retryAt && r.flight === null) {
        if (r.tries > 0) this.wakeAt(r.retryAt);
        return;
      }
      if (group < MARGIN) this.missing++;
      if (r.flight !== null || r === c.job) return;
    } else if (group < MARGIN) this.missing++;
    // An alignment level its manifest does not list is not asked (the grid's levels).
    if (level < LABELS && c.source.remote !== null && c.grid[level] < 0) {
      if (group < MARGIN) this.missing--;
      return;
    }
    this.enqueue(i, level, ty, tx, group);
  }

  /** By group, then distance from the center (labels: by chunk). */
  private enqueue(i: number, level: number, ty: number, tx: number, group: number) {
    const p = group * 1e9 + (level === LABELS ? ty : (tx - this.center[0]) ** 2 + (ty - this.center[1]) ** 2);
    if (this.queued === this.prio.length) {
      const q = new Int32Array(this.queue.length * 2);
      q.set(this.queue);
      this.queue = q;
      const pr = new Float64Array(this.prio.length * 2);
      pr.set(this.prio);
      this.prio = pr;
    }
    // In priority order, by a binary search and one shift (a queue of tens of regions).
    let lo = 0;
    let hi = this.queued;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.prio[mid] <= p) lo = mid + 1;
      else hi = mid;
    }
    this.prio.copyWithin(lo + 1, lo, this.queued);
    this.queue.copyWithin((lo + 1) * Q, lo * Q, this.queued * Q);
    this.queued++;
    this.prio[lo] = p;
    const q = this.queue;
    q[lo * Q] = i;
    q[lo * Q + 1] = level;
    q[lo * Q + 2] = ty;
    q[lo * Q + 3] = tx;
    q[lo * Q + 4] = group;
  }

  /**
   * What an older view asked for and still lacks, once no pass has needed it for `Stream.CancelGraceMs` (`graced`):
   * requests aborted, local tiles in progress dropped. A pan jittering across a region's edge keeps it in flight.
   */
  private cancel() {
    const st = this.stats.current;
    for (let j = this.flying.length - 1; j >= 0; j--) {
      const r = this.flying[j];
      if (r.need === this.clock || this.graced(r)) continue;
      this.discard(r);
      st.cancelled++;
    }
    for (let j = 0; j < this.open.length; j++) {
      const job = this.open[j].job;
      if (job === null || job.need === this.clock || this.graced(job)) continue;
      this.discard(job);
      st.cancelled++;
    }
  }
  /** Needed by a pass less than `Stream.CancelGraceMs` ago (the last 64 passes'). */
  private graced(r: Region): boolean {
    return this.clock - r.need < 64 && this.now - this.passAt[r.need & 63] < this.grace;
  }
  /** The in-flight region this pass does not need that a pass needed longest ago aborted, its slot free; false: none. */
  private preemptStale(): boolean {
    let oldest: Region | null = null;
    for (let j = 0; j < this.flying.length; j++) {
      const r = this.flying[j];
      if (r.need !== this.clock && (oldest === null || r.need < oldest.need)) oldest = r;
    }
    if (oldest === null) return false;
    this.discard(oldest);
    this.stats.current.cancelled++;
    return true;
  }

  /** Remote requests up to the in-flight bound; local tiles a sample row at a time within the frame's budget. */
  private dispatch() {
    const max = this.host.config.get(SeqViewerConfig.Stream.MaxInFlight);
    const st = this.stats.current;
    const q = this.queue;
    this.t0 = this.now;
    this.cells = 0;
    // The tiles in progress this pass needs first: each was the top of its layer's queue when it started. A graced one
    // no pass now needs pauses, its rows kept, until a pass needs it again or its layer's next tile replaces it.
    let spent = false;
    for (let j = 0; j < this.order.length && !spent; j++) {
      const job = this.order[j].job;
      if (job !== null && job.need === this.clock) spent = !this.work(this.order[j]);
    }
    for (let k = 0; k < this.queued; k++) {
      const c = this.order[q[k * Q]];
      const level = q[k * Q + 1];
      const ty = q[k * Q + 2];
      const tx = q[k * Q + 3];
      const index = level === LABELS ? ty : level >= D1 ? tx : indexOf(c.source.columns, level % L, ty, tx);
      let r = c.find(level, index);
      // Asked twice in this pass: the first started it.
      if (r !== undefined && (r.flight !== null || r.state === LOADED || r === c.job)) continue;
      if (c.source.remote !== null) {
        // Full: a graced region no pass now needs gives its slot first, the oldest first; then a view's region takes
        // the slot of a prefetch that is not in the view (aborted: a hint only).
        while (st.inFlight >= max && this.preemptStale());
        if (st.inFlight >= max && (q[k * Q + 4] >= MARGIN || !this.preempt())) continue;
        r = r ?? this.add(c, level, ty, tx, index);
        r.group = q[k * Q + 4];
        this.request(c, r);
      } else if (c.source.local !== null) {
        if (spent) continue;
        r = r ?? this.add(c, level, ty, tx, index);
        r.group = q[k * Q + 4];
        r.tile = { cls: new Uint8Array(TILE * TILE), share: new Uint8Array(TILE * TILE) };
        // One job a layer: a graced one the pass does not need gives way to this one.
        if (c.job !== null && c.job !== r) {
          this.discard(c.job);
          st.cancelled++;
        }
        c.job = r;
        c.row = 0;
        if (!this.work(c)) spent = true;
      }
    }
    this.queued = 0;
    if (spent) this.onChange();
  }

  /** The newest prefetch in flight that the view does not need aborted, its slot free; false: none. */
  private preempt(): boolean {
    for (let j = this.flying.length - 1; j >= 0; j--) {
      const r = this.flying[j];
      if (r.group !== PREFETCH || r.viewed === this.clock) continue;
      this.discard(r);
      this.stats.current.cancelled++;
      return true;
    }
    return false;
  }

  private add(c: Layer, level: number, ty: number, tx: number, index: number): Region {
    const r = new Store.Region(c, level, ty, tx, index);
    r.need = this.clock;
    c.add(r);
    return r;
  }

  /** The layer's local tile in progress, a sample row at a time; true when it finished. */
  private work(c: Layer): boolean {
    const r = c.job as Region;
    const t = r.tile as Tile;
    const s = c.source;
    const cells = s.local as LocalCells;
    const ly = r.level >> 4;
    const lx = r.level & (L - 1);
    while (c.row < TILE) {
      if (this.cells >= LOCAL_CELLS || this.host.scheduler.now() - this.t0 > LOCAL_MS) return false;
      this.cells += reduce(cells.codes, cells.slot, s.columns, s.rows, ly, lx, r.ty, r.tx, c.row, c.row + 1, t);
      c.row++;
    }
    c.job = null;
    this.stats.current.computed++;
    this.store(r, t);
    this.bump();
    return true;
  }

  /** Region r asked of its layer's provider: a label chunk, an alignment tile, or a 1-D region. Request-rate. */
  private request(c: Layer, r: Region) {
    const st = this.stats.current;
    const ac = new AbortController();
    r.flight = ac;
    r.state = PENDING;
    this.flying.push(r);
    st.requests++;
    st.inFlight++;
    const rem = c.source.remote as StreamedData;
    // A provider that throws instead of rejecting fails the region as a rejection does, and one that answers with no
    // promise (a patch, or nothing) is taken as resolved with it (edge rate: one try a request).
    let answer: Promise<DataPatch>;
    try {
      answer = Promise.resolve(rem.provider.fetch(this.requestOf(c, r), ac.signal));
    } catch (e) {
      answer = Promise.reject(e);
    }
    answer.then(
      (p) => this.arrived(r, ac, p),
      (e) => this.failed(r, ac, e),
    );
  }
  private requestOf(c: Layer, r: Region): StreamRequest {
    const rem = c.source.remote as StreamedData;
    const priority = PRIORITY[r.group];
    if (r.level >= D1) {
      const s = ((r.level - D1) / LV) | 0;
      const lvl = ((r.level - D1) % LV) - 1;
      const span = c.spanOf(lvl);
      const sys = rem.manifest.systems[s];
      const start = r.tx * span;
      const data = lvl < 0 || lvl === FINE;
      // A finer region inside a clicked gated region is asked past the gate.
      const forced = r.force || (lvl === FINE && c.find(D1 + s * LV, Math.floor(start / c.span))?.force === true);
      return streamRequest(rem.source.ref, sys.key, {
        start,
        stop: Math.min(start + span, sys.extent[1]),
        level: data ? -1 : lvl,
        priority,
        // The zoom gate by size, and a region asked whole (ST4b).
        ...(data && c.isItems ? { max_bytes: this.host.config.get(SeqViewerConfig.Stream.RegionBytes) } : {}),
        ...(forced ? { force: true } : {}),
      });
    }
    const system = rem.manifest.systems[0]?.key ?? rem.source.ref;
    const { rows, columns } = c.source;
    const [ly, lx] = r.level === LABELS ? [0, 0] : [r.level >> 4, r.level & (L - 1)];
    const r0 = r.ty * TILE * 2 ** ly;
    const c0 = r.level === LABELS ? 0 : r.tx * TILE * 2 ** lx;
    const span = {
      start: c0,
      stop: Math.min(columns, c0 + TILE * 2 ** lx),
      row_start: r0,
      row_stop: Math.min(rows, r0 + TILE * 2 ** ly),
    };
    if (r.level === LABELS)
      return streamRequest(rem.source.ref, system, { ...span, stop: columns, labels: true, priority });
    return streamRequest(rem.source.ref, system, { ...span, level: c.grid[r.level], priority });
  }

  /** A request answered: false when it no longer counts (cancelled, evicted, closed, disposed). */
  private landed(r: Region, ac: AbortController): boolean {
    if (r.flight !== ac) return false;
    r.flight = null;
    this.unfly(r);
    return !this.disposed && !r.owner.closed;
  }
  private unfly(r: Region) {
    const j = this.flying.indexOf(r);
    if (j >= 0) {
      this.flying[j] = this.flying[this.flying.length - 1];
      this.flying.pop();
    }
    this.stats.current.inFlight--;
  }

  /** A patch for region r: held once (a patch for a region held changes nothing), never in the document. */
  private arrived(r: Region, ac: AbortController, p: DataPatch) {
    if (!this.landed(r, ac)) return;
    // A patch the region cannot take (a provider's answer of another shape) is its failure, retried and reported as
    // one, never an unhandled rejection that leaves it pending (GI-13).
    try {
      this.take(r, p);
    } catch (e) {
      // A failed region holds nothing: the nodes its patch carried go (a retry adds them again).
      if (r.carried !== null) {
        this.nodes.release(r.carried);
        r.carried = null;
      }
      this.fail(r, e);
    }
  }
  /**
   * Region r's patch held. The steps that throw on a patch of another shape run first: only a patch taken ends the
   * region's failures and its layer's outage, so a malformed answer backs off as any failure does.
   */
  private take(r: Region, p: DataPatch) {
    const rem = r.owner.source.remote as StreamedData;
    // An answer of another shape than asked is a failure (seqspec:sources#the-viewers-side): row labels are `rows`,
    // a tile `bins`, a 1-D region its data's (`letters`; `bins` or `items` for values and items).
    const kind = rem.source.data_kind;
    const asked =
      r.level === LABELS
        ? p.kind === 'rows'
        : r.level < D1
          ? p.kind === 'bins'
          : kind === 'values' || r.owner.isItems
            ? p.kind === 'bins' || p.kind === 'items'
            : p.kind === kind;
    if (!asked)
      throw new Error(
        `a "${p.kind}" answer where "${r.level === LABELS ? 'rows' : r.level < D1 ? 'bins' : kind}" was asked`,
      );
    const rows = r.level === LABELS && p.kind === 'rows' ? p : null;
    const n = rows === null ? 0 : rows.labels.length;
    if (p.nodes !== undefined && p.nodes.length > 0) r.carried = this.nodes.add(rem.source.ref, p.nodes);
    // Items with a `map`: their coordinate systems, carried by the region as its patch's own nodes are (ST7).
    if (p.kind === 'items' && rem.source.data_kind === 'items' && r.level >= D1 && p.items.map !== undefined)
      this.itemSystems(r, rem, p.items);
    const bytes1 = r.level !== LABELS && r.level >= D1 ? held1(r, p, rem.source.data_kind === 'values') : 0;
    const tile = r.level < D1 && r.level !== LABELS ? tileOf(p, r.level & (L - 1)) : null;
    r.tries = 0;
    this.recovered(r);
    this.outageOver(r.owner);
    if (r.level === LABELS) {
      // Held as a region is: its strings' bytes counted, and evicted with it (`forgetLabels`).
      let bytes = 0;
      for (let k = 0; rows !== null && k < n; k++) {
        rem.labels.set(rows.row_start + k, rows.labels[k]);
        bytes += 2 * (rows.labels[k]?.length ?? 0);
      }
      this.hold(r, bytes);
    } else if (tile !== null) this.store(r, tile);
    else this.hold(r, bytes1);
    this.bump();
    this.onChange();
  }

  private failed(r: Region, ac: AbortController, e: unknown) {
    if (this.landed(r, ac)) this.fail(r, e);
  }
  private fail(r: Region, e: unknown) {
    const c = r.owner;
    this.stats.current.errors++;
    r.state = FAILED;
    c.bump(r);
    // An outage starts when a region starts failing and none of its layer's is (NotServed: each region, once).
    const notServed = e instanceof StreamNotServed;
    if (notServed || (!r.failing && c.failing++ === 0)) {
      const q = r.level >= D1 ? this.requestOf(c, r) : null;
      const what =
        r.level === LABELS
          ? `labels of rows from ${r.ty * TILE}`
          : q
            ? `${q.system} [${q.start}, ${q.stop}) at level ${q.level}`
            : `tile ${r.level >> 4}/${r.level & 15}/${r.ty}/${r.tx}`;
      const why = e instanceof Error ? e.message : String(e);
      const src = c.source.remote?.source;
      this.host.report('svs.stream.load', `${src?.uri ?? src?.stream ?? c.ref}: ${what} failed (${why})`);
    }
    r.failing ||= !notServed;
    const now = this.host.scheduler.now();
    if (notServed) r.retryAt = Number.POSITIVE_INFINITY;
    else if (++r.tries <= RETRIES) {
      r.retryAt = now + 1000 * 2 ** (r.tries - 1);
      this.wakeAt(r.retryAt);
    } else {
      r.tries = 0;
      r.retryAt = now + COOLDOWN_MS;
    }
    this.bump();
    // A failure settles its region: a pass decides whether the view has (`settled`).
    if (this.waiters.length > 0) this.onChange();
  }

  /**
   * A region of layer c loaded: its source is back, so its outage ends with the failed regions no view needs now (they
   * were never asked again); those in view stay, failing until they load. A later failure is a new outage, reported.
   * Arrival rate.
   */
  private outageOver(c: Layer) {
    for (let k = c.all.length - 1; c.failing > 0 && k >= 0; k--) {
      const x = c.all[k];
      if (x.failing && x.state === FAILED && x.need !== this.clock) this.discard(x);
    }
  }
  /** A region's failures ended (it loaded or left): the layer's outage ends with its last. */
  private recovered(r: Region) {
    if (!r.failing) return;
    r.failing = false;
    r.owner.failing--;
  }

  /** A pass at `t` (one timer for every backoff, the earliest). */
  private wakeAt(t: number) {
    if (this.wake >= 0 && this.wakeDue <= t) return;
    if (this.wake >= 0) this.host.scheduler.clearTimer(this.wake);
    this.wakeDue = t;
    this.wake = this.host.scheduler.setTimer(this.onWake, Math.max(0, t - this.host.scheduler.now()));
  }

  /** A region's data, held once (a replacement's bytes replace its own); then eviction past the bound. */
  private store(r: Region, t: Tile) {
    r.tile = t;
    this.hold(r, t.cls.byteLength + t.share.byteLength);
  }
  /** Region r loaded with `bytes` (-1: loaded empty), counted and the cache bounded. */
  private hold(r: Region, bytes: number) {
    const c = r.owner;
    const empty = bytes < 0;
    if (empty) bytes = 0;
    c.bytes += bytes - r.bytes;
    this.bytes += bytes - r.bytes;
    r.bytes = bytes;
    r.state = empty ? EMPTY : LOADED;
    r.landed = -1;
    c.bump(r);
    const st = this.stats.current;
    // A labels chunk is neither level's.
    if (r.level !== LABELS) {
      if (r.group === 1) st.fine++;
      else st.coarse++;
    }
    this.evict();
    st.bytes = this.bytes + this.nodes.bytes;
  }

  /**
   * Over `Stream.CacheBytes`: the least recently needed regions dropped, the current pass's pinned, down to 7/8 of the
   * bound, so a sweep (one scan and a sort) is paid per eighth of the cache filled, not per region. Arrival-rate.
   */
  private evict() {
    const cap = this.host.config.get(SeqViewerConfig.Stream.CacheBytes);
    // Dynamic nodes count, and go with the regions that carried them unless something else holds them.
    if (this.bytes + this.nodes.bytes <= cap) return;
    const xs = this.sweep;
    for (let j = 0; j < this.open.length; j++) {
      const all = this.open[j].all;
      for (let k = 0; k < all.length; k++) if (all[k].bytes > 0 && all[k].need < this.clock) xs.push(all[k]);
    }
    xs.sort(byNeed);
    const low = cap - cap / 8;
    for (let k = 0; k < xs.length && this.bytes + this.nodes.bytes > low; k++) this.discard(xs[k]);
    xs.length = 0;
  }

  /**
   * A labels chunk dropped: its rows' labels go with it, unless another layer over the same source holds that chunk (a
   * label read off the view reads none until a view asks the chunk again). Eviction-rate.
   */
  private forgetLabels(r: Region) {
    const rem = r.owner.source.remote;
    if (rem === null) return;
    for (const c of this.open)
      if (c !== r.owner && c.source.remote === rem && c.find(LABELS, r.index)?.state === LOADED) return;
    for (let row = r.ty * TILE, end = row + TILE; row < end; row++) rem.labels.delete(row);
  }

  /** Region r's items' coordinate systems held (`itemSystemNodes`), carried by it. Arrival-rate. */
  private itemSystems(r: Region, rem: StreamedData, items: StreamTable) {
    const s = ((r.level - D1) / LV) | 0;
    let to = '';
    for (const [key, i] of rem.systems) if (i === s) to = key;
    const nodes: StreamNode[] = [];
    // Each fault once per source and code (a transcripts file repeats it per transcript and region).
    const once = (code: string, message: string) => {
      const k = `${rem.source.ref}|${code}`;
      if (this.reported.has(k)) return;
      this.reported.add(k);
      this.host.report(code, message);
    };
    if (to !== '') itemSystemNodes(rem.source.ref, to, items, nodes, once);
    const more = nodes.length > 0 ? this.nodes.add(rem.source.ref, nodes, true) : null;
    if (more !== null) r.carried = r.carried === null ? more : [...r.carried, ...more];
  }

  /** A region dropped: its request aborted, its bytes released. */
  private discard(r: Region) {
    const c = r.owner;
    const ac = r.flight;
    if (ac !== null) {
      r.flight = null;
      this.unfly(r);
      ac.abort();
    }
    this.bytes -= r.bytes;
    if (r.carried !== null) {
      this.nodes.release(r.carried);
      r.carried = null;
    }
    if (r.level === LABELS && r.state === LOADED) this.forgetLabels(r);
    this.stats.current.bytes = this.bytes + this.nodes.bytes;
    this.recovered(r);
    // Items' lanes are held by the regions holding them: forgotten with the last (ST3b).
    for (let k = 0; r.keys !== null && k < r.keys.length; k++) {
      const x = c.laneOf.get(r.keys[k]);
      if (x !== undefined && --x.held <= 0) c.laneOf.delete(r.keys[k]);
    }
    r.keys = null;
    for (let k = 0; r.foldedKeys !== null && k < r.foldedKeys.length; k++) {
      const x = c.foldedLaneOf.get(r.foldedKeys[k]);
      if (x !== undefined && --x.held <= 0) c.foldedLaneOf.delete(r.foldedKeys[k]);
    }
    r.foldedKeys = null;
    c.remove(r);
    c.bump(r);
  }

  private bump() {
    this.stats.bump();
    this.version.set(this.version.value + 1);
  }

  private close(c: Layer) {
    c.closed = true;
    if (c.aside) this.aside.delete(c.source.remote as StreamedData);
    if (c.aside) this.asides.splice(this.asides.indexOf(c), 1);
    for (let k = c.all.length - 1; k >= 0; k--) this.discard(c.all[k]);
    this.clients.delete(c.ref);
    const j = this.open.indexOf(c);
    if (j >= 0) this.open.splice(j, 1);
  }

  dispose() {
    if (this.disposed) return;
    this.disposed = true;
    this.release();
    if (this.frame >= 0) this.host.scheduler.cancelFrame(this.frame);
    if (this.nodeFrame >= 0) this.host.scheduler.cancelFrame(this.nodeFrame);
    if (this.wake >= 0) this.host.scheduler.clearTimer(this.wake);
    for (const s of this.subs) s.unsubscribe();
    for (let j = this.open.length - 1; j >= 0; j--) this.close(this.open[j]);
    this.order.length = 0;
    this.version.complete();
    this.stats.complete();
  }
}

const byNeed = (a: Region, b: Region) => a.need - b.need;

/**
 * What a shown view places (its shown tracks' layers' data, its blocks' units and molecules): what holds a node. A
 * template's instances hold nothing (ST1b-b): they show what regions and manifests carry, so a template over a
 * streamed source stays under `Stream.CacheBytes`, an instance leaving with its region.
 */
function shownKeys(model: ViewModel | null): Set<string> {
  const keys = new Set<string>();
  const rv = model?.resolved;
  if (!model || !rv) return keys;
  const { set } = model;
  for (let k = 0; k < rv.layers.count; k++) {
    if (rv.tracks.hidden[rv.layers.track[k]] || rv.layers.node[k] !== null) continue;
    for (let d = rv.layers.dataStart[k], e = d + rv.layers.dataCount[k]; d < e; d++) {
      const a = set.annotations[rv.data[d]];
      if (a) keys.add(set.strings.get(a.key));
    }
  }
  for (let b = 0; b < rv.blocks.count; b++) {
    const u = rv.blocks.unit[b];
    if (u < 0) continue;
    keys.add(set.strings.get(set.units.key[u]));
    const m = set.molecules[set.units.molecule[u]];
    if (m) keys.add(set.strings.get(m.key));
  }
  return keys;
}
