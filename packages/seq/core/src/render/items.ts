/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column, issueOf } from '@molstar/core/repr/tree';
import { RenderIssues } from './issues.js';
import type { DataEntry, RangesDensity, RangesList } from '@molstar/seqx-core/data';
import { type QueryItem, rangesLocation } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { axisBreaks, breakHolding, type ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type DensityStep, densityStep, type RenderConfig } from './config.js';
import { NoColor, type PassWriter } from './draw-lists.js';
import {
  dataRegions,
  type LayerSource,
  layerLevel,
  Lift,
  placeholderSpans,
  readTiles,
  toSection,
} from './layer-data.js';
import type { FilterRow, RowFilter } from './filter.js';
import type { Box, Hit, ItemDescription, LayoutContext, StoreReader } from './layer-type.js';
import { breakHit } from './marks.js';
import { addToStrip, paintStrip, resetStrip, type StripBuffers } from './ranges-density.js';
import { sectionId } from './layer-type.js';
import type { SeqTokens } from './tokens.js';

// What the layers of a `ranges` annotation's items share (`ranges`, `secondary`, `variants`; _spec/README.md "Render",
// Layer types, Data): the items of the tiles read (a key and start seen in an earlier tile is that range again, over
// the tile's edge; a key repeated with another start is another piece of the item), the layout's prologue, the
// placeholders, the ticks of the items a break hides, the count strip, an item's hit and where an item's ends are.

export interface ItemRows {
  /** The rows in view (`n` of them; the arrays are kept longer and written over frame to frame). */
  n: number;
  keys: (string | number)[];
  starts: number[];
  ends: number[];
  /** The value of the first of the fields asked for the table has, as text; the label (the id when none). */
  values: (string | undefined)[];
  labels: (string | undefined)[];
  /** Which occurrence of its key (a piece of the item) each row is in the order read. */
  occ: number[];
  /** The rows the layer's filter left out of the tiles read (the track's "n filtered"). */
  filtered: number;
  /** Rows past the end of their system (said once, `seq.render.item.range`): clipped, left out, or longer than a lap. */
  clipped: number;
  leftOut: number;
  clamped: number;
  /** Each row's table fields and its row in that table: what a theme's channels read (a field, a measure). */
  fields: (RangesList['fields'] | undefined)[];
  /** Each row's strand as its table gives it (−1, 0, 1; 0 when none): what a liftover carries through a run. */
  strand: number[];
  row: number[];
  /** Each row's item's whole extent when its table says it (`RangesList.itemStart`, `itemEnd`), else −1. */
  itemStarts: number[];
  itemEnds: number[];
  /**
   * Over a cyclic axis, a row of an item across the origin (one ending at the system's length) and the row of the same
   * key that starts at 0 (its other piece) hold each other's index; else −1. Made once as the rows are read: a hit
   * takes the item's own pieces from it, whatever is in view.
   */
  link: number[];
}

export function newItemRows(): ItemRows {
  return {
    n: 0,
    filtered: 0,
    clipped: 0,
    leftOut: 0,
    clamped: 0,
    keys: [],
    starts: [],
    ends: [],
    values: [],
    labels: [],
    occ: [],
    fields: [],
    strand: [],
    row: [],
    itemStarts: [],
    itemEnds: [],
    link: [],
  };
}

/** A table's field by name. */
function fieldOf(fields: RangesList['fields'], name: string): Column | undefined {
  if (!fields) return undefined;
  for (let k = 0; k < fields.length; k++) if (fields[k].name === name) return fields[k].column;
  return undefined;
}

function textAt(c: Column | undefined, i: number): string | undefined {
  if (!c) return undefined;
  const v = c instanceof CategoryColumn ? c.get(i) : c[i];
  return v === null || v === undefined ? undefined : String(v);
}

/** A row as the filter reads it: its key and its table's fields, a column looked up once a table. */
class TableRow implements FilterRow {
  key: string | number = '';
  i = 0;
  private fields: RangesList['fields'];
  private readonly cols = new Map<string, Column | undefined>();
  bind(fields: RangesList['fields']) {
    if (fields === this.fields) return;
    this.fields = fields;
    this.cols.clear();
  }
  field(name: string): string | number | undefined {
    let c = this.cols.get(name);
    if (c === undefined && !this.cols.has(name)) this.cols.set(name, (c = fieldOf(this.fields, name)));
    if (!c) return undefined;
    const v = c instanceof CategoryColumn ? c.get(this.i) : (c as ArrayLike<unknown>)[this.i];
    return typeof v === 'number' || typeof v === 'string' ? v : v === null || v === undefined ? undefined : String(v);
  }
}
const filterRow = new TableRow();

export interface DensityTile {
  readonly table: RangesDensity;
  readonly level: number;
}

/**
 * Reads the list tables of `entries` into `rows` (cleared first), each row's value of the first of `fields` the table
 * has and its label; the density tables are returned.
 */
export function collectItems(
  rows: ItemRows,
  entries: readonly Extract<DataEntry, { status: 'loaded' }>[],
  fields: readonly string[],
  filter?: RowFilter,
  base?: RowFilter,
  wrap?: { readonly length: number; readonly cyclic: boolean },
): DensityTile[] {
  rows.n = 0;
  rows.filtered = 0;
  rows.clipped = 0;
  rows.leftOut = 0;
  rows.clamped = 0;
  const densities: DensityTile[] = [];
  const earlier = new Map<string | number, number[]>();
  const occurrences = new Map<string | number, number>();
  for (const e of entries) {
    const t = e.tables.ranges;
    if (!t) continue;
    if (t.form === 'density') {
      densities.push({ table: t, level: e.region.level });
      continue;
    }
    const from = rows.n;
    let col: Column | undefined;
    for (let f = 0; f < fields.length && !col; f++) col = fieldOf(t.fields, fields[f]);
    const labels = fieldOf(t.fields, 'label');
    if (filter || base) filterRow.bind(t.fields);
    for (let i = 0; i < t.start.length; i++) {
      const id = t.id?.get(i);
      const key: string | number =
        id !== null && id !== undefined ? id : t.row ? t.row[i] : `${e.region.level}/${e.region.tile}:${i}`;
      if (earlier.get(key)?.includes(t.start[i])) continue;
      if (filter || base) {
        filterRow.key = key;
        filterRow.i = i;
        // the track's own rows: left out and not counted, as the user's filter's are
        if (base && !base(filterRow)) {
          const starts = earlier.get(key);
          if (starts) starts.push(t.start[i]);
          else earlier.set(key, [t.start[i]]);
          continue;
        }
        if (filter && !filter(filterRow)) {
          // left out, and its repeat in a later tile with it
          rows.filtered++;
          const starts = earlier.get(key);
          if (starts) starts.push(t.start[i]);
          else earlier.set(key, [t.start[i]]);
          continue;
        }
      }
      // a row past the system's length (a feature across a ring's origin, as GFF3 writes it): on a cyclic system a row
      // starting past it is moved back by laps, one longer than a lap is one lap, and one ending past it is two rows of
      // one item; on a linear system it is clipped, or left out when it starts past the end (EC-9)
      let start = t.start[i];
      let end = t.end[i];
      let tail = 0;
      if (wrap !== undefined && end > wrap.length) {
        if (!wrap.cyclic) {
          if (start >= wrap.length) {
            rows.leftOut++;
            continue;
          }
          rows.clipped++;
          end = wrap.length;
        } else {
          if (end - start > wrap.length) {
            rows.clamped++;
            end = start + wrap.length;
          }
          if (start >= wrap.length) {
            const back = Math.floor(start / wrap.length) * wrap.length;
            // its repeat in another tile, by the start as read
            const starts = earlier.get(key);
            if (starts) starts.push(start);
            else earlier.set(key, [start]);
            start -= back;
            end -= back;
          }
          if (end > wrap.length) tail = end - wrap.length;
          end = Math.min(end, wrap.length);
        }
      }
      const k = occurrences.get(key) ?? 0;
      occurrences.set(key, k + (tail > 0 ? 2 : 1));
      for (let piece = 0; piece <= (tail > 0 ? 1 : 0); piece++) {
        const j = rows.n++;
        rows.keys[j] = key;
        rows.starts[j] = piece ? 0 : start;
        rows.ends[j] = piece ? tail : end;
        rows.values[j] = textAt(col, i);
        rows.labels[j] = textAt(labels, i) ?? id ?? undefined;
        rows.occ[j] = k + piece;
        rows.fields[j] = t.fields;
        rows.strand[j] = t.strand ? t.strand[i] : 0;
        rows.row[j] = i;
        rows.itemStarts[j] = t.itemStart ? t.itemStart[i] : -1;
        rows.itemEnds[j] = t.itemEnd ? t.itemEnd[i] : -1;
      }
    }
    for (let j = from; j < rows.n; j++) {
      const starts = earlier.get(rows.keys[j]);
      if (starts) starts.push(rows.starts[j]);
      else earlier.set(rows.keys[j], [rows.starts[j]]);
    }
  }
  return densities;
}

/** Links the rows of an item across a cyclic origin: a row ending at the length and a row of its key starting at 0. */
export function linkPieces(rows: ItemRows, wrap: { readonly length: number; readonly cyclic: boolean } | undefined) {
  rows.link.length = rows.n;
  rows.link.fill(-1);
  if (!wrap?.cyclic) return;
  const heads = new Map<string | number, number>();
  for (let j = 0; j < rows.n; j++) if (rows.starts[j] === 0 && rows.ends[j] < wrap.length) heads.set(rows.keys[j], j);
  if (!heads.size) return;
  for (let j = 0; j < rows.n; j++) {
    if (rows.ends[j] !== wrap.length || rows.starts[j] === 0) continue;
    const h = heads.get(rows.keys[j]);
    if (h === undefined) continue;
    rows.link[j] = h;
    rows.link[h] = j;
  }
}

/** What every items layer's layout holds. */
/**
 * What an items layer used on its strip (_spec/README.md "Render", Layer types: `used`): its band's height when it drew
 * a mark, a hidden items' tick, a pending or failed span, a strip (the coarse steps) or a count strip with a column
 * counted (items past the lanes: overlapping ranges in a thin band), else 0 (text collapses it).
 */
export function itemsUsed(lay: ItemsBase & { readonly m: number }): number {
  const any =
    lay.m > 0 ||
    lay.tn > 0 ||
    lay.pending.length > 0 ||
    lay.failed.length > 0 ||
    lay.step === 'bars' ||
    lay.step === 'density';
  if (any) return lay.box.height;
  if (lay.density) {
    const counts = lay.strip.counts;
    for (let c = 0; c < lay.density.width; c++) if (counts[c] > 0) return lay.box.height;
  }
  return 0;
}

export interface ItemsBase {
  key: string;
  system: string;
  section: string;
  /** The section's address: a tick's hit names it as a break's does. */
  sectionAddress: string;
  owner: string;
  /** The data's system: an element's extent is on it. */
  dataSystem: string;
  /** The layer's source: a hit carries its element onto the section along its path. */
  src?: LayerSource;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  tokens: SeqTokens;
  rows: ItemRows;
  /** The data not here yet (pending) or never (failed), as x spans. */
  pending: number[];
  failed: number[];
  /** The hidden items' ticks: one a break and lane whatever their number, `x, lane, count, hiddenStart, hiddenEnd`. */
  tn: number;
  ticks: number[];
  /**
   * A lane's height and the gap within it (a layer of one lane: the box's, none), and where its first lane starts (the
   * band's top, or below its padding and label lane, its lanes centred: _spec/README.md "Render", Track heights).
   */
  laneHeight: number;
  gap: number;
  laneY: number;
  hoverWidth: number;
  hoverOffset: number;
  /** What a hover's squircle reads (its rule: `squircleRatio`, `squircleMax`, `squircleMin`). */
  config: RenderConfig;
  /** The count strip (the whole box, or a row of it) and its buffers, kept while the box's columns do not change. */
  density?: StripBuffers['image'];
  strip: StripBuffers;
}

export interface ItemsPrep {
  readonly src?: LayerSource;
  /** The annotation's address: the owner of its items' element keys. */
  readonly owner: string;
}

/** The fields of a new layout. */
export function newItemsBase(ctx: LayoutContext, prep: ItemsPrep): ItemsBase {
  return {
    key: ctx.key,
    system: ctx.axis.system,
    section: sectionId(ctx.section),
    sectionAddress: ctx.section.address,
    owner: prep.owner,
    dataSystem: prep.src?.dataSystem ?? '',
    box: ctx.box,
    axis: ctx.axis,
    step: 'cells',
    tokens: ctx.tokens,
    rows: newItemRows(),
    pending: [],
    failed: [],
    tn: 0,
    ticks: [],
    laneHeight: ctx.box.height,
    gap: 0,
    laneY: ctx.box.y,
    hoverWidth: ctx.config.outline.hover,
    hoverOffset: ctx.config.emphasis.hoverOffset,
    config: ctx.config,
    strip: { diff: new Float64Array(0), counts: new Float64Array(0), v: 1 },
  };
}

/** A layout's start: what the frame gives, and what the last one left emptied. */
export function beginItems(lay: ItemsBase, ctx: LayoutContext, prep: ItemsPrep) {
  lay.key = ctx.key;
  lay.system = ctx.axis.system;
  lay.section = sectionId(ctx.section);
  lay.sectionAddress = ctx.section.address;
  lay.owner = prep.owner;
  lay.dataSystem = prep.src?.dataSystem ?? '';
  lay.src = prep.src;
  lay.box = ctx.box;
  lay.axis = ctx.axis;
  lay.tokens = ctx.tokens;
  lay.hoverWidth = ctx.config.outline.hover;
  lay.hoverOffset = ctx.config.emphasis.hoverOffset;
  lay.config = ctx.config;
  lay.laneHeight = ctx.box.height;
  lay.gap = 0;
  lay.laneY = ctx.box.y;
  lay.rows.n = 0;
  lay.pending.length = 0;
  lay.failed.length = 0;
  lay.tn = 0;
  lay.ticks.length = 0;
  lay.density = undefined;
}

/**
 * Reads the items in view and the stretches its breaks hide (an item there is drawn as its break's tick) into `lay.rows`,
 * the placeholders into `pending` and `failed`; the density tables it brought.
 */
export function readItems(
  lay: ItemsBase,
  ctx: LayoutContext,
  data: StoreReader,
  src: LayerSource,
  fields: readonly string[],
): DensityTile[] {
  const axis = ctx.axis;
  lay.step = densityStep(axis.scale, axis.dpr, ctx.config);
  if (src.unreadFilter)
    ctx.report(
      issueOf(RenderIssues, 'seq.render.filter', 'The layer’s data query is not read row by row: no row is left out', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.filter'].hint,
      }),
    );
  // carried through a map read by region (a liftover): its reader, loaded with the genome chunk (`lift.ts`)
  if (src.lift) return Lift.module?.readLifted(lay, ctx, data, src, fields) ?? [];
  const regions = dataRegions(src, axis, true);
  // a row of a cyclic system that runs past its end also overlaps the last position: its tile gives the origin's piece
  const wrap = src.length !== undefined && src.wraps ? { length: src.length, cyclic: src.cyclic === true } : undefined;
  if (wrap?.cyclic) regions.push(wrap.length - 1, wrap.length);
  const read = readTiles(ctx, data, src, regions, layerLevel(src, data, axis));
  placeholderSpans(src, axis, read.pending, lay.pending);
  placeholderSpans(src, axis, read.failed, lay.failed);
  const densities = collectItems(lay.rows, read.entries, fields, src.filter, src.base, wrap);
  linkPieces(lay.rows, src.length !== undefined ? { length: src.length, cyclic: src.cyclic === true } : undefined);
  const { clipped, leftOut, clamped } = lay.rows;
  if (clipped || leftOut || clamped) {
    const said = [
      clipped && `${clipped} clipped to its end`,
      leftOut && `${leftOut} left out (starting past its end)`,
      clamped && `${clamped} longer than one lap, clamped to one`,
    ].filter(Boolean);
    ctx.report(
      issueOf(RenderIssues, 'seq.render.item.range', `Rows past the end of ${src.dataSystem}: ${said.join(', ')}`, {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.item.range'].hint,
      }),
    );
  }
  return densities;
}

/** The tick of a hidden item at section position `p` in `lane`: one a break and lane, its count kept. */
export function tickOn(lay: ItemsBase, axis: AxisTransform, p: number, lane: number) {
  const b = axisBreaks(axis);
  const i = b.length ? breakHolding(axis, p) : -1;
  if (i < 0 || b[i] < axis.u0 || b[i] > axis.u1) return;
  const x = axis.xOfUnit(b[i]);
  for (let t = 0; t < lay.tn; t++)
    if (lay.ticks[5 * t] === x && lay.ticks[5 * t + 1] === lane) {
      lay.ticks[5 * t + 2]++;
      return;
    }
  lay.ticks.push(x, lane, 1, b[i + 1], b[i + 2]);
  lay.tn++;
}

/** The hidden items' ticks: one a break and lane, in the break ink. */
export function drawTicks(lay: ItemsBase, out: PassWriter) {
  const px = 1 / lay.axis.dpr;
  const g = out.args;
  for (let t = 0; t < lay.tn; t++) {
    g[0] = Math.round(lay.ticks[5 * t] * lay.axis.dpr) / lay.axis.dpr - 1.5 * px;
    g[1] = lay.laneY + lay.ticks[5 * t + 1] * lay.laneHeight + lay.gap / 2;
    g[2] = 3 * px;
    g[3] = lay.laneHeight - lay.gap;
    g[4] = lay.tokens.colors['seq-break'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

/** A tick in `lane` within the break's slop of `x`, as the break it stands on (its count the items it hides). */
export function tickHit(lay: ItemsBase, x: number, lane: number): Hit | undefined {
  for (let t = 0; t < lay.tn; t++) {
    const at = 5 * t;
    if (lay.ticks[at + 1] !== lane || Math.abs(x - lay.ticks[at]) > 3) continue;
    return breakHit(
      { section: { address: lay.sectionAddress }, id: lay.section },
      lay.axis,
      lay.ticks[at],
      lay.ticks[at + 3],
      lay.ticks[at + 4],
      lay.ticks[at + 2],
    );
  }
  return undefined;
}

/**
 * The count strip of the items and the density bins of a read, one a device column, in the density tokens: its image
 * (the layers' strip at the `bars` and `density` steps and from a source's bins).
 */
export function countStrip(
  strip: StripBuffers,
  src: LayerSource,
  axis: AxisTransform,
  rows: ItemRows,
  densities: readonly DensityTile[],
  tokens: SeqTokens,
): StripBuffers['image'] {
  resetStrip(strip, axis.columns);
  strip.v = 1;
  for (let i = 0; i < rows.n; i++)
    addToStrip(strip, src, axis, rows.starts[i], Math.max(rows.ends[i], rows.starts[i] + 1));
  for (const { table, level } of densities)
    for (let i = 0; i < table.bin.length; i++) {
      strip.v = table.count[i];
      addToStrip(strip, src, axis, table.bin[i] * level, (table.bin[i] + 1) * level);
    }
  return paintStrip(strip, axis.columns, tokens, 1, axis);
}

/**
 * An item's hit: its element with its extent on the data's system and its hull on the section's (every piece carried
 * along the layer's path), the position under the pointer.
 */
export function itemHit(
  base: Pick<Hit, 'layer' | 'section' | 'system' | 'position'>,
  lay: ItemsBase,
  key: string | number,
  s: number,
  e: number,
  row?: number,
): Hit {
  let lo = s;
  let hi = e;
  if (lay.src?.path) {
    const pieces = toSection(lay.src, s, e, []);
    lo = Number.POSITIVE_INFINITY;
    hi = Number.NEGATIVE_INFINITY;
    for (let k = 0; k + 1 < pieces.length; k += 2) {
      lo = Math.min(lo, pieces[k]);
      hi = Math.max(hi, pieces[k + 1]);
    }
  }
  return {
    ...base,
    element: { owner: lay.owner, key },
    extent: rangesLocation(lay.dataSystem, piecesOf(lay.rows, row, s, e)),
    ...(lo < hi ? { hull: [lo, hi] as const } : {}),
  };
}

/**
 * An item's ranges for a hit on `row`: its own, and over a cyclic axis its other piece across the origin (`link`: the
 * click selects the item's extent, not the piece under the pointer). No scan, nothing for an item of one piece.
 */
function piecesOf(rows: ItemRows, row: number | undefined, s: number, e: number): number[] {
  const p = row === undefined ? -1 : (rows.link[row] ?? -1);
  return p < 0 ? [s, e] : [s, e, rows.starts[p], rows.ends[p]];
}

/** The rows of every item of an inline annotation, by key: what an end is decided from (undefined: a source's). */
export type ItemExtents = ReadonlyMap<string | number, readonly QueryItem[]>;

const extents = new WeakMap<readonly QueryItem[], ItemExtents>();

/** The extents of the layer's data when it is inline (every row is known), else undefined. */
export function itemExtents(model: StaticModel, layer: ResolvedLayer): ItemExtents | undefined {
  const address = layer.data?.address;
  if (address === undefined) return undefined;
  let items = model.annotations.get(address)?.items;
  if (!items) for (const [a, e] of model.annotations) if (address.startsWith(`${a}/`)) items = e.items;
  if (!items) return undefined;
  let m = extents.get(items);
  if (!m) {
    const by = new Map<string | number, QueryItem[]>();
    for (const it of items) {
      const g = by.get(it.key);
      if (g) g.push(it);
      else by.set(it.key, [it]);
    }
    extents.set(items, (m = by));
  }
  return m;
}

/**
 * Which rows in view hold the ends of their item, in reading order (rows sharing a key are the pieces of one item):
 * `last` the row holding its end, `first` the one holding its start, found among all the item's rows (`known`). Across
 * the origin of a cyclic system the row that starts at it comes last (the longest, if several) and the one that ends
 * there first; otherwise the row ending last and the row starting first. An item whose rows are not known (a tiled
 * source's: the rest may be outside what was read) has its ends from its rows' item extent when its table gives one
 * (the row ending at the item's end, the one starting at its start), else none. `first` and `last` are written for rows
 * `0..n`.
 */
export function markEnds(
  rows: ItemRows,
  known: ItemExtents | undefined,
  cyclic: boolean,
  length: number,
  first: Uint8Array,
  last: Uint8Array,
) {
  const { keys, starts, ends, n, itemStarts, itemEnds } = rows;
  const decided = new Map<string | number, [number, number, number, number]>();
  for (let i = 0; i < n; i++) {
    first[i] = 0;
    last[i] = 0;
    const g = known?.get(keys[i]);
    if (!g) {
      if (itemEnds[i] >= 0) {
        if (ends[i] === itemEnds[i]) last[i] = 1;
        if (starts[i] === itemStarts[i]) first[i] = 1;
      }
      continue;
    }
    let d = decided.get(keys[i]);
    if (!d) {
      const across = cyclic && g.length > 1 && g.some((r) => r.end === length) && g.some((r) => r.start === 0);
      let c = across ? -1 : 0;
      let f = across ? -1 : 0;
      for (let k = 0; k < g.length; k++) {
        if (across) {
          if (g[k].start === 0 && (c < 0 || g[k].end > g[c].end)) c = k;
          if (g[k].end === length && (f < 0 || g[k].start < g[f].start)) f = k;
        } else {
          if (g[k].end > g[c].end) c = k;
          if (g[k].start < g[f].start) f = k;
        }
      }
      d = [g[c].start, g[c].end, g[f].start, g[f].end];
      decided.set(keys[i], d);
    }
    if (starts[i] === d[0] && ends[i] === d[1]) last[i] = 1;
    if (starts[i] === d[2] && ends[i] === d[3]) first[i] = 1;
  }
}

const imageEnd: number[] = [];
const imageStart: number[] = [];

/**
 * Where the first and last positions of a row `[s, e)` land on the section, and whether the row runs against it (its
 * start image above its end image): what an arrowhead at the row's end or start needs. On the section's own system
 * the positions themselves. No number made per row: read the fields, never a call with them.
 */
export class RowImages {
  /** The last position's image `[ea, eb)` and the first's `[sa, sb)`; NaN when a position maps nowhere. */
  ea = Number.NaN;
  eb = Number.NaN;
  sa = Number.NaN;
  sb = Number.NaN;
  flip = false;

  set(src: LayerSource | undefined, s: number, e: number) {
    this.flip = false;
    if (!src?.path) {
      this.ea = e - 1;
      this.eb = e;
      this.sa = s;
      this.sb = s + 1;
      return;
    }
    toSection(src, e - 1, e, imageEnd);
    toSection(src, s, s + 1, imageStart);
    if (imageEnd.length < 2 || imageStart.length < 2) {
      this.ea = this.eb = this.sa = this.sb = Number.NaN;
      return;
    }
    this.ea = imageEnd[0];
    this.eb = imageEnd[imageEnd.length - 1];
    this.sa = imageStart[0];
    this.sb = imageStart[imageStart.length - 1];
    this.flip = e - s > 1 && this.sa > this.ea;
  }

  /**
   * The side an arrowhead stands on at the end (`atEnd`) or the start of the row, for the piece `[pa, pb)` of the
   * section it is drawn on, as it is drawn: 0 none (the end is not in this piece, or the view clips it), 1 the tip on
   * the right, 2 on the left. The end's side is the screen's right where the residues run to the right.
   */
  tip(axis: AxisTransform, box: Box, pa: number, pb: number, atEnd: boolean): number {
    const lo = atEnd ? this.ea : this.sa;
    const hi = atEnd ? this.eb : this.sb;
    if (!(lo >= pa && hi <= pb)) return 0;
    const at = atEnd === this.flip ? lo : hi;
    const x = axis.x(at);
    if (!(x >= box.x - 0.5 && x <= box.x + box.width + 0.5)) return 0;
    return ((axis.dir === 1) !== this.flip) === atEnd ? 1 : 2;
  }
}

/**
 * An item's extent from its pieces: their hull, or on a cyclic system of `length` the circle less its widest gap
 * (an item over the origin, two rows meeting there, is one arc: its start past its end).
 */
export function itemArc(pieces: number[], length: number | undefined): [number, number] {
  const n = pieces.length >> 1;
  const order = Array.from({ length: n }, (_, k) => k).sort((a, b) => pieces[2 * a] - pieces[2 * b]);
  const first = pieces[2 * order[0]];
  let last = 0;
  for (let k = 0; k < n; k++) last = Math.max(last, pieces[2 * k + 1]);
  if (length === undefined) return [first, last];
  let best = length - last + first;
  let arc: [number, number] = [first, last];
  let reach = pieces[2 * order[0] + 1];
  for (let k = 1; k < n; k++) {
    const s = pieces[2 * order[k]];
    if (s - reach > best) {
      best = s - reach;
      arc = [s, reach];
    }
    reach = Math.max(reach, pieces[2 * order[k] + 1]);
  }
  return arc;
}

/** Whether an item's pieces meet at a cyclic system's origin: one ends at its length, another starts at 0. */
function crossesOrigin(pieces: readonly number[], length: number): boolean {
  let end = false;
  let start = false;
  for (let k = 0; k + 1 < pieces.length; k += 2) {
    if (pieces[k + 1] === length) end = true;
    if (pieces[k] === 0) start = true;
  }
  return end && start;
}

/**
 * The items in view, each once: a row per piece makes one item (its extent the hull of the pieces on the data's
 * system, its place the hull on the section clamped to the box). Allocates; for an export's accessible layer and the
 * viewer's announcements, never a frame.
 */
export function describeItems(lay: ItemsBase): ItemDescription[] {
  const rows = lay.rows;
  const by = new Map<
    string | number,
    { label: string; s: number; e: number; lo: number; hi: number; pieces: number[] }
  >();
  for (let i = 0; i < rows.n; i++) {
    const key = rows.keys[i];
    let lo = rows.starts[i];
    let hi = rows.ends[i];
    if (lay.src?.path) {
      const pieces = toSection(lay.src, rows.starts[i], rows.ends[i], []);
      lo = Number.POSITIVE_INFINITY;
      hi = Number.NEGATIVE_INFINITY;
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        lo = Math.min(lo, pieces[k]);
        hi = Math.max(hi, pieces[k + 1]);
      }
    }
    const g = by.get(key);
    if (g) {
      g.s = Math.min(g.s, rows.starts[i]);
      g.e = Math.max(g.e, rows.ends[i]);
      g.lo = Math.min(g.lo, lo);
      g.hi = Math.max(g.hi, hi);
      g.pieces.push(rows.starts[i], rows.ends[i]);
    } else
      by.set(key, {
        label: rows.labels[i] || String(key),
        s: rows.starts[i],
        e: rows.ends[i],
        lo,
        hi,
        pieces: [rows.starts[i], rows.ends[i]],
      });
  }
  // an item across a cyclic origin (its rows meeting there): its two runs from its start (EC-9)
  const length = lay.src?.cyclic ? lay.src.length : undefined;
  const x0 = lay.box.x;
  const x1 = lay.box.x + lay.box.width;
  const out: ItemDescription[] = [];
  for (const [key, g] of by) {
    const a = g.lo < g.hi ? lay.axis.x(g.lo) : Number.NaN;
    const b = g.lo < g.hi ? lay.axis.x(g.hi) : Number.NaN;
    const l = Number.isFinite(a) && Number.isFinite(b) ? Math.min(a, b) : x0;
    const r = Number.isFinite(a) && Number.isFinite(b) ? Math.max(a, b) : x1;
    // across the origin only when a piece ends at the length and another starts at 0 (a spliced item's intron wider than
    // its outer gap is no crossing)
    const arc = length !== undefined && crossesOrigin(g.pieces, length) ? itemArc(g.pieces, length) : undefined;
    out.push({
      key,
      label: g.label,
      start: g.s,
      end: g.e,
      ...(arc && arc[0] > arc[1] ? { runs: [arc[0], length!, 0, arc[1]] } : {}),
      x0: Math.max(x0, Math.min(x1, l)),
      x1: Math.max(x0, Math.min(x1, r)),
    });
  }
  return out;
}
