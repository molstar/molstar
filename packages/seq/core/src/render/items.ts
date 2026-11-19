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
import { type DensityStep, densityStep } from './config.js';
import { NoColor, type PassWriter } from './draw-lists.js';
import { dataRegions, type LayerSource, layerLevel, placeholderSpans, readTiles, toSection } from './layer-data.js';
import type { FilterRow, RowFilter } from './filter.js';
import type { Box, Hit, LayoutContext, StoreReader } from './layer-type.js';
import { breakHit } from './marks.js';
import { addToStrip, paintStrip, resetStrip, type StripBuffers } from './ranges.js';
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
  /** Each row's table fields and its row in that table: what a theme's channels read (a field, a measure). */
  fields: (RangesList['fields'] | undefined)[];
  row: number[];
}

export function newItemRows(): ItemRows {
  return { n: 0, filtered: 0, keys: [], starts: [], ends: [], values: [], labels: [], occ: [], fields: [], row: [] };
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
): DensityTile[] {
  rows.n = 0;
  rows.filtered = 0;
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
    if (filter) filterRow.bind(t.fields);
    for (let i = 0; i < t.start.length; i++) {
      const id = t.id?.get(i);
      const key: string | number =
        id !== null && id !== undefined ? id : t.row ? t.row[i] : `${e.region.level}/${e.region.tile}:${i}`;
      if (earlier.get(key)?.includes(t.start[i])) continue;
      if (filter) {
        filterRow.key = key;
        filterRow.i = i;
        if (!filter(filterRow)) {
          // left out, and its repeat in a later tile with it
          rows.filtered++;
          const starts = earlier.get(key);
          if (starts) starts.push(t.start[i]);
          else earlier.set(key, [t.start[i]]);
          continue;
        }
      }
      const k = occurrences.get(key) ?? 0;
      occurrences.set(key, k + 1);
      const j = rows.n++;
      rows.keys[j] = key;
      rows.starts[j] = t.start[i];
      rows.ends[j] = t.end[i];
      rows.values[j] = textAt(col, i);
      rows.labels[j] = textAt(labels, i) ?? id ?? undefined;
      rows.occ[j] = k;
      rows.fields[j] = t.fields;
      rows.row[j] = i;
    }
    for (let j = from; j < rows.n; j++) {
      const starts = earlier.get(rows.keys[j]);
      if (starts) starts.push(rows.starts[j]);
      else earlier.set(rows.keys[j], [rows.starts[j]]);
    }
  }
  return densities;
}

/** What every items layer's layout holds. */
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
  /** A lane's height and the gap within it (a layer of one lane: the box's, none). */
  laneHeight: number;
  gap: number;
  hoverWidth: number;
  hoverOffset: number;
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
    hoverWidth: ctx.config.outline.hover,
    hoverOffset: ctx.config.emphasis.hoverOffset,
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
  lay.laneHeight = ctx.box.height;
  lay.gap = 0;
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
  const read = readTiles(ctx, data, src, dataRegions(src, axis, true), layerLevel(src, data, axis));
  placeholderSpans(src, axis, read.pending, lay.pending);
  placeholderSpans(src, axis, read.failed, lay.failed);
  return collectItems(lay.rows, read.entries, fields, src.filter);
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
    g[1] = lay.box.y + lay.ticks[5 * t + 1] * lay.laneHeight + lay.gap / 2;
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
  return paintStrip(strip, axis.columns, tokens);
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
    extent: rangesLocation(lay.dataSystem, [s, e]),
    ...(lo < hi ? { hull: [lo, hi] as const } : {}),
  };
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
 * source's: the rest may be outside what was read) has no end marked. `first` and `last` are written for rows `0..n`.
 */
export function markEnds(
  rows: ItemRows,
  known: ItemExtents | undefined,
  cyclic: boolean,
  length: number,
  first: Uint8Array,
  last: Uint8Array,
) {
  const { keys, starts, ends, n } = rows;
  const decided = new Map<string | number, [number, number, number, number]>();
  for (let i = 0; i < n; i++) {
    first[i] = 0;
    last[i] = 0;
    const g = known?.get(keys[i]);
    if (!g) continue;
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
