/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type DataEntry, fitLevel, type Shape, type Want } from '@molstar/seqx-core/data';
import { inversePath, type Path, projectRanges, type StaticModel, withGenome } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { compileRowFilter, type RowFilter } from './filter.js';
import type { LayoutContext, StoreReader } from './layer-type.js';

// How a layer reads its data (_spec/README.md "Render", Data): the access plan of its data's node names the store's
// source and system; the region in view is carried back along the layer's path onto the data's system; the level fits
// the positions a device pixel covers there; the tiles' entries are read, and those not loaded are placeholders.

/**
 * A liftover's module (`lift.ts`: a lifted layer's source and items, hover and selection through its map), loaded with
 * the genome chunk: a liftover is over assemblies, so a model holding one waits for both. Why lazy: a reader that never
 * opens a genome pays nothing for it.
 */
export const Lift: { module?: typeof import('./lift.js') } = {};
withGenome(() =>
  import('./lift.js').then((m) => {
    Lift.module = m;
  }),
);

/** Where a layer's data comes from in the store. */
export interface LayerSource {
  /** The store's source (an inline node's own address, or a source node's). */
  readonly source: string;
  /** The name the source knows the system by. */
  readonly system: string;
  /** The data's system in the spec. */
  readonly dataSystem: string;
  readonly shape: Shape;
  /** A sided pairs source's side the system is read on (its wants carry it). */
  readonly side?: 'a' | 'b';
  /** The data's system's length, when known. */
  readonly length?: number;
  /** Whether the data's system wraps, and whether a source's rows may run past its length (then wrapped or clipped). */
  readonly cyclic?: boolean;
  readonly wraps?: boolean;
  /** The path from the data's system to the section's (none: the same) and back. */
  readonly path?: Path;
  readonly back?: Path;
  /** Data positions per section position along the path: 1/3 from residues onto nucleotides, 3 back, else 1. */
  readonly ratio: number;
  /** The filter the layer's query puts on its rows (an item layer's), compiled once; `unread`: a form it cannot read. */
  readonly filter?: RowFilter;
  /** A leading item form's predicate: its rejected rows are left out and not counted as filtered. */
  readonly base?: RowFilter;
  readonly unreadFilter?: boolean;
  /** Data over another assembly carried onto the section's sequence through a map read by region (`lift.ts`). */
  readonly lift?: LiftSource;
}

/**
 * A lifted layer's reads (_spec/README.md "Render", A liftover): the map's runs on the section's system, read by the
 * end on it (`side`), and the data on each far system, by the name the map's source lists that system by.
 */
export interface LiftSource {
  readonly map: string;
  readonly side: 'a' | 'b';
  /** The map's runs: a pairs source read on the section's system by `side`. */
  readonly near: LayerSource;
  /** The data on each far system (no path: its rows are carried along the runs), by the far name the runs give. */
  readonly far: ReadonlyMap<string, LayerSource>;
}

/** The store's source of a placed layer's data; undefined when its node has no access plan. */
export function layerSource(model: StaticModel, layer: ResolvedLayer): LayerSource | undefined {
  const data = layer.data;
  if (!data || !layer.placed) return undefined;
  // carried through a map read by region (a liftover): its source, loaded with the genome chunk (`lift.ts`)
  if (data.lift) return Lift.module?.liftSource(model, data);
  // a template instance's data is the overlay's: its plan and length travel with it
  let access = data.access ?? model.systems.get(data.address)?.access;
  if (!access) {
    let a = model.annotations.get(data.address);
    // a query reads the annotation it starts at; its form after it filters the rows (`filter`)
    if (!a) for (const [address, entry] of model.annotations) if (data.address.startsWith(`${address}/`)) a = entry;
    // a ranges layer over links reads their ends on its system (_spec/README.md "Render", Layer types: ranges)
    access = a?.shape === 'pairs' && layer.type !== 'pairs' ? a.ends?.plan(data.system) : a?.access;
  }
  if (!access) return undefined;
  const plan =
    access.kind === 'inline'
      ? { source: access.source.address, system: access.system, shape: access.shape }
      : access.side
        ? { source: access.source, system: access.system, shape: access.shape, side: access.side }
        : { source: access.source, system: access.system, shape: access.shape };
  let path = layer.path;
  let dataSystem = data.system;
  // a unit's letters from a source are its sequence's, read on the sequence's positions through the unit's step
  // (_spec/README.md "Resolve", Molecules and units)
  const own = model.systems.get(data.system);
  if (own?.kind === 'unit' && own.parent && access.kind === 'source' && !data.access) {
    const through = model.path(own.parent.system, path?.to ?? data.system).path;
    if (through) {
      path = through;
      dataSystem = own.parent.system;
    }
  }
  const filter = data.filter ? compileRowFilter(data.filter) : undefined;
  const base = data.base ? compileRowFilter(data.base) : undefined;
  const length = dataSystem === data.system ? (data.length ?? own?.length) : model.systems.get(dataSystem)?.length;
  return {
    ...plan,
    dataSystem,
    ...(length !== undefined ? { length } : {}),
    ...(model.systems.get(dataSystem)?.cyclic ? { cyclic: true } : {}),
    // a source's rows past the system's length are the reader's to wrap or clip; inline rows were wrapped by the model
    ...(length !== undefined && access.kind !== 'inline' ? { wraps: true } : {}),
    ...(path ? { path, back: inversePath(path) } : {}),
    ratio: path?.step === 1 ? 1 / 3 : path?.step === -1 ? 3 : 1,
    ...(filter ? { filter } : {}),
    ...(base ? { base } : {}),
    ...((data.filter && !filter) || (data.base && !base) || data.unread ? { unreadFilter: true } : {}),
  };
}

/**
 * The data's positions in view (`[start, end, …]`, sorted, merged): the section's carried back along the path; with
 * `hidden`, the stretches the breaks in view hide too (what a layer draws as their ticks). Near ones merged (`near`,
 * `mergeNear`: what a layer reads); without, exactly those in view (what a layer walks position by position).
 */
export function dataRegions(src: LayerSource, axis: AxisTransform, hidden = false, near = true): number[] {
  const vis = axis.visible([]);
  if (hidden) {
    const b = axis.breaks([]);
    for (let i = 0; i + 2 < b.length; i += 3)
      if (b[i + 2] >= b[i + 1]) vis.push(b[i + 1], b[i + 2]);
      else vis.push(b[i + 1], axis.length ?? b[i + 1], 0, b[i + 2]);
  }
  const exact = HalfOpenRanges.normalize(src.back ? projectRanges(src.back, vis) : vis);
  return near ? mergeNear(exact) : Array.from(exact);
}

/**
 * Sorted, disjoint `[start, end, …]` widened outward to multiples of `step` and merged (`inward`: shrunk to them, the
 * empty dropped), pushed onto `out` (a bridge's `R`: a window on the far level's grid, _spec/README.md "Render",
 * Bridges).
 */
export function widenedToGrid(ranges: readonly number[], step: number, out: number[], inward = false): number[] {
  for (let j = 0; j + 1 < ranges.length; j += 2) {
    const s = (inward ? Math.ceil(ranges[j] / step) : Math.floor(ranges[j] / step)) * step;
    const e = (inward ? Math.floor(ranges[j + 1] / step) : Math.ceil(ranges[j + 1] / step)) * step;
    if (e <= s) continue;
    if (out.length && s <= out[out.length - 1]) out[out.length - 1] = Math.max(e, out[out.length - 1]);
    else out.push(s, e);
  }
  return out;
}

/** Gaps of at most this many positions between regions are read with them (small: an intron of an exon pair stays apart). */
export const ReadGap = 16;

/**
 * Sorted, merged regions with the gaps of at most `ReadGap` positions between them closed: a collapsed axis's pieces
 * in view (150,000 single positions, EC-33) are read as a few regions, not a want each.
 */
export function mergeNear(ranges: ArrayLike<number>, gap = ReadGap): number[] {
  const out: number[] = [];
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const n = out.length;
    if (n && ranges[i] - out[n - 1] <= gap) out[n - 1] = Math.max(out[n - 1], ranges[i + 1]);
    else out.push(ranges[i], ranges[i + 1]);
  }
  return out;
}

/** A region on the data's system carried onto the section's as pieces (`[start, end, …]`), into `out`. */
export function toSection(src: LayerSource, s: number, e: number, out: number[]): number[] {
  out.length = 0;
  if (!src.path) {
    out.push(s, e);
    return out;
  }
  return projectRanges(src.path, [s, e], out);
}

/**
 * The level that fits the view on the data's system: the coarsest bin no wider than `coarser` device pixels' positions
 * (a values layer's `values.coarser`, else 1), and coarser still while the view's span at it would want more tiles
 * than the store's bound (`maxTilesPerWant`, where the reader says it: a whole chromosome at a high device pixel ratio
 * reads its coarse level, never a refused one); with `finest`, the finest level.
 */
export function layerLevel(
  src: LayerSource,
  data: StoreReader,
  axis: AxisTransform,
  finest = false,
  coarser = 1,
): number {
  const levels = data.coverage(src.source)?.levels ?? [1];
  if (finest) return levels[0];
  return boundedLevel(levels, (coarser * src.ratio) / axis.devicePx, (axis.u1 - axis.u0) * src.ratio, data);
}

/**
 * The coarsest of `levels` no wider than `positionsPerPixel`, and coarser still while `span` positions at it would
 * want more tiles than the store's bound (`maxTilesPerWant`, where the reader says it): what a layer reads and what a
 * template wants of a nested source, so neither asks for a long system a tile at level 1.
 */
export function boundedLevel(
  levels: readonly number[],
  positionsPerPixel: number,
  span: number,
  data: Pick<StoreReader, 'settings'>,
): number {
  let level = fitLevel(levels, positionsPerPixel);
  const s = data.settings;
  if (s)
    for (let i = levels.indexOf(level); i + 1 < levels.length && span / (s.tileBins * level) + 1 > s.maxTilesPerWant; )
      level = levels[++i];
  return level;
}

/**
 * A grid's level (_spec/README.md "Render", Data: the level rules by shape): the finest of `levels` whose bin is at least
 * `binPx` device pixels on both axes (the larger of the two axes' positions a pixel), then coarser while the view's x
 * tiles × y tiles pass the store's bound. The other way from values: a cell under a pixel is reduced again for nothing.
 * Why the coarser axis's level, not the finer's: a 2D bin is square, so a level fine enough for the finer axis puts
 * many bins in a device pixel on the coarser one (a block pair of a zoomed chromosome against one at the fit), reading
 * tiles of the coarser axis's whole window to draw them a pixel apart; a bin at least `binPx` on both axes reads no bin
 * the plot cannot show. The finer axis then shows bins wider than its pixels.
 */
export function gridLevel(
  levels: readonly number[],
  positionsPerPixel: { readonly x: number; readonly y: number },
  span: { readonly x: number; readonly y: number },
  data: Pick<StoreReader, 'settings'>,
  binPx = 1,
): number {
  const need = binPx * Math.max(positionsPerPixel.x, positionsPerPixel.y);
  let i = levels.findIndex((l) => l >= need);
  if (i < 0) i = levels.length - 1;
  const s = data.settings;
  const bins = s?.gridTileBins ?? s?.tileBins ?? 256;
  const tiles = (l: number) => (Math.ceil(span.x / (bins * l)) + 1) * (Math.ceil(span.y / (bins * l)) + 1);
  if (s) while (i + 1 < levels.length && tiles(levels[i]) > s.maxTilesPerWant) i++;
  return levels[i];
}

/** A pairs layer's two reads: bins at a pixel's level, items at a tile size (`RegionRequest.items`). */
export interface PairsLevels {
  readonly bins: number;
  readonly items: number;
}

/**
 * A pairs layer's levels: **bins** at the coarsest level no wider than `binPx` device pixels' positions (bounded as
 * `boundedLevel`); **items** at the finest level whose tiles over the view fit the store's bound (no pixel floor: a
 * 10 Mb view reads its items in a few large tiles, never thousands of small ones).
 */
export function pairsLevels(
  levels: readonly number[],
  positionsPerPixel: number,
  span: number,
  data: Pick<StoreReader, 'settings'>,
  binPx = 2,
): PairsLevels {
  const bins = boundedLevel(levels, binPx * positionsPerPixel, span, data);
  const s = data.settings;
  let items = levels[levels.length - 1];
  if (!s) items = levels[0];
  else
    for (const l of levels)
      if (span / (s.tileBins * l) + 1 <= s.maxTilesPerWant) {
        items = l;
        break;
      }
  return { bins, items };
}

/**
 * Whether a pairs layer reads its items (the gate, with hysteresis): asked when the bin pairs with an end in view sum
 * to at most `gate`; once shown, kept while the items in view (when loaded, else the bins' sum) stay at most
 * `hysteresis` × `gate`. Why: a pan at the threshold would otherwise flip bins and arcs every frame.
 */
export function itemsGate(
  showing: boolean,
  binPairs: number,
  items: number | undefined,
  gate: number,
  hysteresis: number,
): boolean {
  if (!showing) return binPairs <= gate;
  return (items ?? binPairs) <= gate * Math.max(1, hysteresis);
}

/** A 2D tile read: its entry and whether it is the asked tile's mirror (a symmetric grid's, drawn transposed). */
export interface GridTileRead {
  readonly entry: Extract<DataEntry, { status: 'loaded' }>;
  readonly transposed: boolean;
}

/** The 2D tiles' state over a rect read: loaded entries, and the rects (`x0, x1, y0, y1` positions) still pending or failed. */
export interface GridRead {
  readonly entries: GridTileRead[];
  readonly pending: number[];
  readonly failed: number[];
  tooMany?: boolean;
}

/**
 * Wants a 2D rect (the view at `view` priority, a ring of `marginTiles` tiles round it at `margin`) and reads its tiles:
 * a tile not loaded is a pending rect (EC-169: never drawn or hit as empty), a failed one a failed rect.
 */
export function readGridTiles(
  ctx: Pick<LayoutContext, 'want'> & { readonly key?: string },
  data: StoreReader,
  src: LayerSource,
  x: readonly [number, number],
  y: { readonly system: string; readonly start: number; readonly end: number },
  level: number,
  marginTiles = 1,
): GridRead {
  const out: GridRead = { entries: [], pending: [], failed: [] };
  // the layer's key on its wants: a layer's wants over every block pair are bounded together (`maxTilesPerLayer`)
  const want: Want = {
    source: src.source,
    system: src.system,
    start: x[0],
    end: x[1],
    level,
    y,
    priority: 'view',
    layer: ctx.key,
  };
  ctx.want(want);
  if (marginTiles > 0) {
    const m = marginTiles * level * (data.settings?.gridTileBins ?? 256);
    const yr = { ...y, start: Math.max(0, y.start - m), end: y.end + m };
    ctx.want({ ...want, start: Math.max(0, x[0] - m), end: x[1] + m, y: yr, priority: 'margin' });
  }
  if (data.coverage(src.source) === undefined || !data.tilesOf) {
    out.pending.push(x[0], x[1], y.start, y.end);
    return out;
  }
  // the layer's wants together past the store's bounds: none asked, its band "zoom in" (drawn by the track)
  if (ctx.key !== undefined && data.layerRefused?.(ctx.key)) {
    out.pending.push(x[0], x[1], y.start, y.end);
    out.tooMany = true;
    return out;
  }
  const tiles = data.tilesOf(want);
  if (!tiles.length && data.tooMany?.(want)) {
    out.pending.push(x[0], x[1], y.start, y.end);
    out.tooMany = true;
    return out;
  }
  for (const t of tiles) {
    const e = data.get(t.key);
    const r = t.region;
    // the asked tile's rect: the served one's, swapped for a mirror
    const rect = t.transposed ? [r.y!.start, r.y!.end, r.start, r.end] : [r.start, r.end, r.y!.start, r.y!.end];
    if (e?.status === 'loaded') out.entries.push({ entry: e, transposed: t.transposed });
    else if (e?.status === 'empty') continue;
    else if (e?.status === 'failed') out.failed.push(...rect);
    else out.pending.push(...rect);
  }
  return out;
}

/**
 * A layer's level as last fitted to the store's bounds, by the tiles it spans (`fittedLevel`), by reader, layer, and
 * the block (its data system, a grid's y block's too): a layer over blocks is fitted a block at a time, each kept.
 */
type FitEntry = { readonly sig: Float64Array; readonly level: number };
const Fitted = new WeakMap<object, Map<string, Map<string, Map<string, FitEntry>>>>();
/** The signature being made (its first `SigN` slots): grown, never shrunk (`length = 0` would free its store). */
let Sig = new Float64Array(64);
let SigN = 0;
function sigPut(v: number) {
  if (SigN === Sig.length) {
    const g = new Float64Array(2 * Sig.length);
    g.set(Sig);
    Sig = g;
  }
  Sig[SigN++] = v;
}

/**
 * The layer's level fitted to the store's bounds on a layer's wants together (`maxTilesPerLayer`, `maxBytesPerLayer`;
 * the design's "past it coarser, then zoom in"): `level`, else the first coarser listed level whose view wants fit
 * (`StoreReader.layerLoad`), else `level` (the store refuses it: the band's "zoom in"). Asked again only when the view's
 * tiles change (`sig` writes them), so a pan within its tiles asks nothing. 1D: `regions` on the data's system; 2D: `y`.
 */
export function fittedLevel(
  ctx: { readonly key?: string },
  data: StoreReader,
  src: LayerSource,
  regions: readonly number[],
  level: number,
  y?: { readonly system: string; readonly start: number; readonly end: number },
): number {
  const key = ctx.key;
  if (!data.layerLoad || key === undefined) return level;
  const levels = data.coverage(src.source)?.levels ?? [level];
  const bins = (y ? data.settings?.gridTileBins : data.settings?.tileBins) ?? 1024;
  SigN = 0;
  sigPut(level);
  for (let i = 0; i + 1 < regions.length; i++) sigPut(Math.floor(regions[i] / (bins * level)));
  if (y) {
    sigPut(Math.floor(y.start / (bins * level)));
    sigPut(Math.floor(y.end / (bins * level)));
  }
  let byKey = Fitted.get(data);
  if (!byKey) Fitted.set(data, (byKey = new Map()));
  let byX = byKey.get(key);
  if (!byX) byKey.set(key, (byX = new Map()));
  let byY = byX.get(src.dataSystem);
  if (!byY) byX.set(src.dataSystem, (byY = new Map()));
  const yKey = y?.system ?? '';
  const at = byY.get(yKey);
  if (at && sameSig(at.sig)) return at.level;
  let fitted = level;
  for (let i = Math.max(0, levels.indexOf(level)); i < levels.length; i++) {
    const l = levels[i];
    if (l < level) continue;
    const wants: Want[] = [];
    for (let k = 0; k + 1 < regions.length; k += 2)
      wants.push({
        source: src.source,
        system: src.system,
        start: regions[k],
        end: regions[k + 1],
        level: l,
        layer: key,
        ...(y ? { y } : {}),
      });
    if (data.layerLoad(wants).fits) {
      fitted = l;
      break;
    }
  }
  byY.set(yKey, { sig: Sig.slice(0, SigN), level: fitted });
  return fitted;
}

/** Whether a kept fit's signature is `Sig` (a loop: no closure made a call). */
function sameSig(sig: Float64Array): boolean {
  if (sig.length !== SigN) return false;
  for (let i = 0; i < sig.length; i++) if (sig[i] !== Sig[i]) return false;
  return true;
}

/** The tiles' state over the regions read: loaded entries, and the data positions still pending or failed. */
export interface TileRead {
  readonly entries: Extract<DataEntry, { status: 'loaded' }>[];
  /** `[start, end, …]` on the data's system. */
  readonly pending: number[];
  readonly failed: number[];
  /** Whether a region was past the store's bound on a want's tiles (its positions pending, none listed). */
  tooMany?: boolean;
}

/** What reading tiles reads of a layout: its wants and the margin it asks beside them. */
export type TileContext = Pick<LayoutContext, 'want'> & {
  readonly config: Pick<LayoutContext['config'], 'margin'>;
  /** The layer's key: its wants carry it, so the store bounds them together (`maxTilesPerLayer`, `maxBytesPerLayer`). */
  readonly key?: string;
};

/**
 * Wants the regions (the view at `view` priority, a margin each side) and reads the tiles of the view's. A tile not
 * asked for, pending or evicted is pending (EC-169: never drawn as empty); a failed one is failed; an unknown source
 * (not attached yet) is pending over the regions. `extra` marks the wants (items tiles, the layer they are for).
 */
export function readTiles(
  ctx: TileContext,
  data: StoreReader,
  src: LayerSource,
  regions: readonly number[],
  level: number,
  extra?: { readonly items?: boolean; readonly layer?: string; readonly priority?: 'view' | 'margin' },
): TileRead {
  const out: TileRead = { entries: [], pending: [], failed: [] };
  const known = data.coverage(src.source) !== undefined;
  // two regions can share a tile (a cyclic view wrapping a short system): each entry is read once
  const seen = new Set<string>();
  for (let i = 0; i + 1 < regions.length; i += 2) {
    const start = regions[i];
    const end = regions[i + 1];
    // a layout makes these every frame: the view's and its margins' each of one small shape (no spread, no absent
    // member: an items want is the rare other shape)
    const layer = extra?.layer ?? ctx.key;
    const items = extra?.items === true;
    const priority = extra?.priority ?? 'view';
    const want = wantOf(src, start, end, level, priority, layer, items);
    ctx.want(want);
    // a margin read has none of its own
    const m = priority === 'view' ? Math.ceil((end - start) * ctx.config.margin) : 0;
    if (m > 0) {
      const max = src.length ?? Number.MAX_SAFE_INTEGER;
      if (start > 0) ctx.want(wantOf(src, Math.max(0, start - m), start, level, 'margin', layer, items));
      if (end < max) ctx.want(wantOf(src, end, Math.min(max, end + m), level, 'margin', layer, items));
    }
    if (!known) {
      out.pending.push(start, end);
      continue;
    }
    // the layer's wants together past the store's bounds: none asked, its band "zoom in" (drawn by the track)
    if (ctx.key !== undefined && data.layerRefused?.(ctx.key)) {
      out.pending.push(start, end);
      out.tooMany = true;
      continue;
    }
    const keys = data.keysOf(want);
    // past the store's bound on a want's tiles none is listed or asked (`seq.store.too_many_tiles`): pending, never
    // drawn as empty, at the cost of no key a tile (past a source's end none is listed either: nothing there)
    if (!keys.length && data.tooMany?.(want)) {
      out.pending.push(start, end);
      out.tooMany = true;
      continue;
    }
    for (const key of keys) {
      const e = data.get(key);
      if (e?.status === 'loaded') {
        if (!seen.has(key)) out.entries.push(e);
        seen.add(key);
      } else if (e?.status === 'empty') continue;
      else if (e?.status === 'failed') out.failed.push(Math.max(start, e.region.start), Math.min(end, e.region.end));
      else if (e) out.pending.push(Math.max(start, e.region.start), Math.min(end, e.region.end));
      else out.pending.push(start, end);
    }
  }
  return out;
}

/** A 1D want of a layer: one shape a kind (an items want's the other). */
function wantOf(
  src: LayerSource,
  start: number,
  end: number,
  level: number,
  priority: 'view' | 'margin',
  layer: string | undefined,
  items: boolean,
): Want {
  // a sided source's want names its side (the rare third shape)
  if (src.side)
    return { source: src.source, system: src.system, start, end, level, priority, layer, items, side: src.side };
  return items
    ? { source: src.source, system: src.system, start, end, level, priority, layer, items: true }
    : { source: src.source, system: src.system, start, end, level, priority, layer };
}

/** Placeholders' x extents: data positions carried onto the section and to x (`x0, x1, cut` triples), into `out`. */
export function placeholderSpans(src: LayerSource, axis: AxisTransform, ranges: readonly number[], out: number[]) {
  const pieces: number[] = [];
  const merged = HalfOpenRanges.normalize(ranges);
  for (let i = 0; i + 1 < merged.length; i += 2) {
    toSection(src, merged[i], merged[i + 1], pieces);
    for (let k = 0; k + 1 < pieces.length; k += 2) axis.spans(pieces[k], pieces[k + 1], out);
  }
  return out;
}
