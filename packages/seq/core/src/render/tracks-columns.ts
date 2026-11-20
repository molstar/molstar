/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ResolvedLayer, ResolvedTrack, ResolvedView } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';
import { elide, type TextMetrics } from './text.js';
import { SlimRailPx } from './rail.js';
import { columnsRoom, fitRowColumns, type RowColumn, rowColumnsOf } from './rows-model.js';
import type { LayoutInput, LayoutState } from './tracks-types.js';
import { sectionLayoutName } from './tracks-axes.js';

// The view's label column (_spec/README.md "Render", Configuration: `labelColumn`; The rows layer; Row columns): its
// measured width, what a rows layer's tree, numbers and row columns add to it and to the right gutter, and the track
// labels cut to it.

/**
 * The rows' header in the view's one label column (_spec/README.md "Render", The rows layer): a rows layer's tree in
 * tree order and its first residue numbers widen the label column by their width (the labels keep theirs), its last
 * numbers take a right gutter. The widest of a view's tracks sections, so every section's data starts at one x. Its
 * row columns (Row columns, in a text figure's too) widen the label column by the left ones and the right gutter by
 * the right ones, by layer.
 */
export function rowsHeader(
  view: ResolvedView,
  rails: LayoutInput['rails'],
): {
  readonly tree: number;
  readonly numbers: boolean;
  readonly rail: boolean;
  /** The view has a tracks section (a label column a narrow frame shrinks). */
  readonly tracks: boolean;
} {
  let tree = 0;
  let numbers = false;
  let rail = false;
  let tracks = false;
  for (const s of view.sections) {
    if (sectionLayoutName(s, view) !== 'tracks') continue;
    tracks = true;
    for (const t of s.tracks)
      for (const l of t.layers) {
        if (l.type !== 'rows' || !l.placed || l.params.hidden === true) continue;
        if (l.params.numbers !== 'none') numbers = true;
        if (rails?.[t.address] !== undefined) rail = true;
        const w = l.params.tree_width;
        if (l.params.order === 'tree' && l.params.tree !== undefined)
          tree = Math.max(tree, typeof w === 'number' ? w : 120);
      }
  }
  return { tree, numbers, rail, tracks };
}

const NoColumns: ReadonlyMap<string, readonly RowColumn[]> = new Map();
const ViewColumns = new WeakMap<ResolvedView, ReadonlyMap<string, readonly RowColumn[]>>();

/** The row columns of a view's rows layers (tracks and text sections), by layer: read once a view. */
export function viewColumns(view: ResolvedView): ReadonlyMap<string, readonly RowColumn[]> {
  let out = ViewColumns.get(view);
  if (out) return out;
  let columns: Map<string, readonly RowColumn[]> | undefined;
  for (const s of view.sections) {
    const name = sectionLayoutName(s, view);
    if (name !== 'tracks' && name !== 'text') continue;
    for (const t of s.tracks)
      for (const l of t.layers) {
        if (l.type !== 'rows' || !l.placed || l.params.hidden === true) continue;
        const cols = l.params.columns === undefined ? undefined : rowColumnsOf(l.params.columns);
        if (cols?.length) (columns ??= new Map()).set(layerKey(t, l), cols);
      }
  }
  out = columns ?? NoColumns;
  ViewColumns.set(view, out);
  return out;
}

/** The title row a rows track's layer with columns takes above its rows (_spec/README.md "Render", Row columns). */
export function titleRowOf(track: ResolvedTrack, config: RenderConfig): number {
  for (const l of track.layers)
    if (l.type === 'rows' && Array.isArray(l.params.columns) && l.params.columns.length)
      return config.rows.columns.titlePx;
  return 0;
}

/**
 * The row columns' widths a view grants by layer (_spec/README.md "Render", Row columns: EC-107): what the view lacks
 * for the data's least width at every column's own width is taken from the columns, right ones first, each down to a
 * floor and then hidden; and the room the widest layer's left and right columns take.
 */
export function columnsFit(
  state: LayoutState,
  columns: ReadonlyMap<string, readonly RowColumn[]>,
  width: number,
  fixed: number,
  config: RenderConfig,
  rail: boolean,
): ColumnsFit | undefined {
  if (!columns.size) return undefined;
  const kept = state.columnsFit as ColumnsFit | undefined;
  if (kept && kept.columns === columns && kept.width === width && kept.fixed === fixed && kept.config === config)
    return kept;
  const { floorPx, gapPx, dataMinPx } = config.rows.columns;
  const widths = new Map<string, number[]>();
  let left = 0;
  let right = 0;
  for (const [key, cols] of columns) {
    const own = cols.map((c) => c.width);
    const need = fixed + dataMinPx + columnsRoom(cols, own, 'left', gapPx) + columnsRoom(cols, own, 'right', gapPx);
    const w = fitRowColumns(cols, need - width, floorPx, gapPx, []);
    widths.set(key, w);
    left = Math.max(left, columnsRoom(cols, w, 'left', gapPx));
    right = Math.max(right, columnsRoom(cols, w, 'right', gapPx));
  }
  // a rows rail keeps the gutter's gap from the last right column
  if (rail && right) right += gapPx;
  const rightX = width - config.inset - (rail ? SlimRailPx : 0) - right;
  const grants = new Map<string, { widths: readonly number[]; rightX: number }>();
  for (const [key, w] of widths) grants.set(key, { widths: w, rightX });
  const fit: ColumnsFit = { columns, width, fixed, config, left, right, grants };
  state.columnsFit = fit;
  return fit;
}

/** The row columns' fit a view's layout keeps (what it was made for, the rooms, each layer's grant). */
interface ColumnsFit {
  readonly columns: ReadonlyMap<string, readonly RowColumn[]>;
  readonly width: number;
  readonly fixed: number;
  readonly config: RenderConfig;
  readonly left: number;
  readonly right: number;
  readonly grants: ReadonlyMap<string, { readonly widths: readonly number[]; readonly rightX: number }>;
}

/** A config with its label column widened by `extra` px, one a base and width. */
const Widened = new WeakMap<RenderConfig, Map<number, RenderConfig>>();
export function widened(base: RenderConfig, extra: number): RenderConfig {
  if (!extra) return base;
  let m = Widened.get(base);
  if (!m) Widened.set(base, (m = new Map()));
  let c = m.get(extra);
  if (!c) m.set(extra, (c = { ...base, labelWidth: base.labelWidth + extra }));
  return c;
}

/**
 * The label column's width a view lays out with (_spec/README.md "Render", Configuration: `labelColumn`): the widest
 * track label, between the column's least and most; `labelWidth` when the
 * column is not measured. What the viewer's data box starts after, before any layout (its subpath, off the barrel).
 */
export function labelColumnWidth(view: ResolvedView, config: RenderConfig, text: TextMetrics): number {
  const c = config.labelColumn;
  if (!c.measure) return config.labelWidth;
  const kept = Columns.get(view);
  if (kept && kept.config === config && kept.text === text && kept.version === text.version) return kept.width;
  const width = measureColumn(view, config, text);
  Columns.set(view, { config, text, version: text.version, width });
  return width;
}

/** The label column's width by view, with the config and metrics it was measured for. */
const Columns = new WeakMap<
  ResolvedView,
  { readonly config: RenderConfig; readonly text: TextMetrics; readonly version: number; readonly width: number }
>();

function measureColumn(view: ResolvedView, config: RenderConfig, text: TextMetrics): number {
  const c = config.labelColumn;
  let w = 0;
  // the labels alone: a direction mark (a section reversed) fades at the column's end, so a reverse never moves the
  // data box (its windows keep their px a position)
  for (const s of view.sections)
    for (const t of s.tracks) {
      // an alignment's rows label themselves in the column: their names, not the tracks', are what it holds (the
      // column keeps its fixed width, its row labels their own ellipsis)
      for (const l of t.layers) if (l.type === 'rows') return config.labelWidth;
      const m = text.measure(config.labelFont, t.label);
      if (m > w) w = m;
    }
  return Math.min(c.max, Math.max(c.min, Math.ceil(w)));
}

/** A config with its label column `width` px wide, one a base and width. */
const Measured = new WeakMap<RenderConfig, Map<number, RenderConfig>>();
export function measured(base: RenderConfig, width: number): RenderConfig {
  if (width === base.labelWidth) return base;
  let m = Measured.get(base);
  if (!m) Measured.set(base, (m = new Map()));
  let c = m.get(width);
  if (!c) m.set(width, (c = { ...base, labelWidth: width }));
  return c;
}

/** A track label cut to `max` px with an ellipsis, kept by label and width (a layout makes no string a track). */
export function elidedLabel(state: LayoutState, text: TextMetrics, font: string, label: string, max: number): string {
  const kept = state.elided;
  if (state.elidedVersion !== text.version) {
    kept.clear();
    state.elidedVersion = text.version;
  }
  let m = kept.get(label);
  if (!m) {
    if (kept.size > 4096) kept.clear();
    kept.set(label, (m = new Map()));
  }
  let t = m.get(max);
  if (t === undefined) m.set(max, (t = elide(text, font, label, max)));
  return t;
}

/** The right gutter's width when a rows layer draws residue numbers: five digits and a gap. */
export function numbersGutter(text: TextMetrics, config: RenderConfig): number {
  return text.measure(config.itemFont, '00000') + config.labelGap;
}

/** A layer's key: its address, else its track's and its index. */
export function layerKey(track: ResolvedTrack, layer: ResolvedLayer): string {
  return layer.address ?? `${track.address}:${layer.index}`;
}
