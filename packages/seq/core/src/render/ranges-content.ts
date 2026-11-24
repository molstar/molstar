/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { QueryItem } from '@molstar/seqx-core/query';
import { toSection } from './layer-data.js';
import type { ContentContext } from './layer-type.js';
import type { LayerContent } from './sizing.js';
import { LabelInset, LabelMinPx, minGlyphPx } from './range-labels.js';
import { knownEmpty, knownLanes } from './known.js';
import { mappedPieces, type RangesPrep } from './ranges-lay.js';

// A ranges layer's content (_spec/README.md "Render", Track heights): the lanes its whole inline data packs into and
// whether a label goes beside at the fit, kept by its prep, read before the layout sizes its track.

/** A filter's row over an inline item: its key and its fields (`base`, the track's own rows). */
class ItemRow {
  key: string | number = '';
  item: QueryItem | undefined;
  field(name: string): string | number | undefined {
    const v = this.item?.fields?.[name];
    return typeof v === 'string' || typeof v === 'number' ? v : undefined;
  }
}
const itemRow = new ItemRow();

/** A layer's lanes over its whole data, kept by its prep (its data does not change under it). */
const Lanes = new WeakMap<RangesPrep, number>();
/** Whether a label goes beside at a fit, kept by its prep with what it was found for. */
const Beside = new WeakMap<
  RangesPrep,
  {
    fit: number;
    font: string;
    text: unknown;
    version: number;
    beside: boolean;
    yes: LayerContent;
    no: LayerContent;
  }
>();

/** A layer's content with no label lane, kept by its prep (a layout makes no object a frame). */
const LanesOnly = new WeakMap<RangesPrep, LayerContent>();
function lanesOnly(prep: RangesPrep, lanes: number, partial = false): LayerContent {
  let c = LanesOnly.get(prep);
  if (!c || c.lanes !== lanes || (c.partial ?? false) !== partial)
    LanesOnly.set(prep, (c = partial ? { lanes, labels: false, partial } : { lanes, labels: false }));
  return c;
}

/** The inline items a layer draws: every row of its annotation, less those its track's own filter leaves out. */
function contentRows(prep: RangesPrep): QueryItem[] | undefined {
  const ext = prep.extents;
  const src = prep.src;
  if (!ext || !src) return undefined;
  const out: QueryItem[] = [];
  for (const [key, rows] of ext)
    for (const r of rows) {
      if (src.base) {
        itemRow.key = key;
        itemRow.item = r;
        if (!src.base(itemRow)) continue;
      }
      out.push(r);
    }
  itemRow.item = undefined;
  return out;
}

/** An item a lane (the gene line's): each key's rows as one row over their hull. */
export function itemHulls(rows: readonly QueryItem[]): QueryItem[] {
  const hulls = new Map<string | number, QueryItem>();
  for (const r of rows) {
    const h = hulls.get(r.key);
    hulls.set(r.key, h ? { ...h, start: Math.min(h.start, r.start), end: Math.max(h.end, r.end) } : r);
  }
  return [...hulls.values()];
}

/** The lanes first-fit by start needs (as few as the most rows over one position: intervals pack greedily). */
export function packedLanes(rows: readonly QueryItem[]): number {
  const order = rows.map((_, k) => k).sort((a, b) => rows[a].start - rows[b].start);
  const ends: number[] = [];
  for (const k of order) {
    const s = rows[k].start;
    const e = Math.max(rows[k].end, s + 1);
    let lane = 0;
    while (lane < ends.length && ends[lane] > s) lane++;
    ends[lane] = e;
  }
  return ends.length;
}

/** An inline item's label, as the layout reads it: its `label` field, else its id. */
function itemLabel(r: QueryItem): string | undefined {
  const l = r.fields?.label;
  if (typeof l === 'string' && l) return l;
  return typeof r.key === 'string' ? r.key : undefined;
}

/** Whether some item's label fits none of its marks with the whole axis in view (`fit` px a position). */
function labelsBeside(prep: RangesPrep, rows: readonly QueryItem[], ctx: ContentContext): boolean {
  const src = prep.src!;
  const font = ctx.config.itemFont;
  const glyph = minGlyphPx(font);
  const heads = (((prep.arrow & 1) + ((prep.arrow >> 1) & 1)) * (ctx.config.laneHeight - ctx.config.laneGap) * 9) / 10;
  for (const r of rows) {
    const label = itemLabel(r);
    if (!label) continue;
    let widest = 0;
    if (src.path) {
      const pieces = toSection(src, r.start, r.end, mappedPieces);
      for (let k = 0; k + 1 < pieces.length; k += 2) widest = Math.max(widest, pieces[k + 1] - pieces[k]);
    } else widest = r.end - r.start;
    // its arrowheads take no text (as the layout's `tipRoom`)
    const w = widest * ctx.fit - heads;
    if (w < LabelMinPx || label.length * glyph + LabelInset > w) return true;
    if (ctx.text.measure(font, label) + LabelInset > w) return true;
  }
  return false;
}

/**
 * A layer's content (`LayerType.content`): its lanes over its whole inline data, and with a label lane at a fit whether
 * a label goes beside; undefined for data not known whole (a source's).
 */
export function rangesContent(prep: RangesPrep, ctx: ContentContext): LayerContent | undefined {
  // `lanes: one`: one lane, no label lane, whatever the data
  if (prep.one) return lanesOnly(prep, 1);
  let lanes = Lanes.get(prep);
  let rows: QueryItem[] | undefined;
  if (lanes === undefined) {
    rows = contentRows(prep);
    // a source's data is not known whole: the lanes it has packed, one while it is pending (the loading lane)
    if (!rows)
      return prep.src
        ? knownEmpty(prep.known ?? prep)
          ? lanesOnly(prep, 0)
          : lanesOnly(prep, knownLanes(prep.known ?? prep), true)
        : undefined;
    // the gene line lays an item's pieces in one lane, by its extent
    lanes = packedLanes(prep.line === 'gene' ? itemHulls(rows) : rows);
    Lanes.set(prep, lanes);
  }
  // a gene model's label is beside it in its own lane (`ranges-gene.ts`): no label lane
  if (!ctx.labelLane || !(ctx.fit > 0) || prep.line === 'gene') return lanesOnly(prep, lanes);
  const font = ctx.config.itemFont;
  const kept = Beside.get(prep);
  // kept while the fit, the font and the metrics are (a web font arriving measures again)
  if (kept && kept.fit === ctx.fit && kept.font === font && kept.text === ctx.text && kept.version === ctx.text.version)
    return kept.beside ? kept.yes : kept.no;
  rows ??= contentRows(prep);
  const beside = rows ? labelsBeside(prep, rows, ctx) : false;
  const answer = {
    fit: ctx.fit,
    font,
    text: ctx.text,
    version: ctx.text.version,
    beside,
    yes: { lanes, labels: true },
    no: { lanes, labels: false },
  };
  Beside.set(prep, answer);
  return beside ? answer.yes : answer.no;
}
