/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue, Value } from '@molstar/core/repr/tree';
import { matchLeaves, readNewick } from '@molstar/seqx-core/data/alignment-readers';
import type { Key, Location } from '@molstar/seqx-core/query';
import type { Channel, MeasureEntry } from '@molstar/seqx-core/theme/index';
import type { Axis, ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, type LayerChannels, PlainDatum } from './channels.js';
import { type RenderConfig } from './config.js';
import { type ImageData, NoColor } from './draw-lists.js';
import { ImageCache } from './images.js';
import { type Alphabet, LettersColorFallback } from './letters.js';
import type { LayerSource } from './layer-data.js';
import type { Box, Hit, LayoutContext } from './layer-type.js';
import { Pyramid } from './pyramid.js';
import { ColumnsLay } from './rows-columns.js';
import {
  buildClades,
  type Clades,
  type ColumnMeasures,
  columnMeasures,
  DisplayRows,
  type HeldCells,
  isGap,
  measureEntries,
  numericValues,
  type RowColumn,
  RowMeasureFields,
  upper,
} from './rows-model.js';
import type { SeqTokens } from './tokens.js';

// What a rows layer keeps frame to frame (_spec/README.md "Render", The rows layer): its prepared form, its held
// cells, order, folds and display rows, the rows' transform (a row's y, the row under a y, a key's row and back), its
// labels, the fields of its rows table and its cells' colors.

/** A hit the layer writes in place (`RowsLay.scratch`): every field set each time, so its shape never changes. */
type ScratchHit = { -readonly [K in keyof Hit]: Hit[K] };

/** No residue number (a row of gaps). */
export const NoNumber = -0x80000000;

/** A rows layer's order: data, tree, or a field of the rows table. */
export type RowsOrder = 'data' | 'tree' | { readonly field: string; readonly direction: 'ascending' | 'descending' };

export interface RowsPrep {
  readonly src?: LayerSource;
  /** The alignment's address: the owner of its clades' and rows' elements. */
  readonly owner: string;
  readonly rowsSystem: string;
  readonly alphabet: Alphabet;
  readonly channels: LayerChannels;
  /** The tree as read and its address, when the layer names one. */
  readonly tree?: { readonly address: string; readonly read: ReturnType<typeof readNewick> };
  readonly order: RowsOrder;
  readonly collapsed: readonly string[];
  readonly pinned: readonly string[];
  readonly hidden: readonly string[];
  readonly summary: 'consensus' | 'bar';
  readonly numbers: 'gutters' | 'none';
  readonly treeWidth: number;
  readonly reference?: string;
  /** The layer's params as written. */
  readonly params: Readonly<Record<string, Value>>;
  /** The layer's theme slots (its channels are made again with the data's measures). */
  readonly channelsForms: ResolvedLayer['theme'];
  /** The row columns (Row columns), as read. */
  readonly columns: readonly RowColumn[];
  /** The sources that may hold the rows' joined fields (a system named as the rows system): every other source. */
  readonly joinSources: readonly string[];
  /** Of those, the ones that could join before they attach (their `system` the rows system, a runtime bound). */
  readonly joinCandidates: readonly string[];
  /** The fields a column or the order names (the only ones joined). */
  readonly joinFields: readonly string[];
}

/** What a rows layer keeps frame to frame: the held cells, the order and folds, the marks in view. */
export class RowsLay {
  key = '';
  section = '';
  system = '';
  rowsSystem = '';
  owner = '';
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  axis!: AxisTransform;
  tokens!: SeqTokens;
  config!: RenderConfig;
  /** Data rows (undefined: not known yet), columns. */
  total: number | undefined = undefined;
  columns = 0;
  held: HeldCells | undefined;
  clades: Clades | undefined;
  cladesOf: unknown;
  readonly display = new DisplayRows();
  displayKey = '';
  /** A band that does not scroll (a text figure's): every row, unfolded, no tree. */
  figure = false;
  /** The layer's params as written (what its commands edit). */
  params: Readonly<Record<string, Value>> = {};
  /** The layer's `collapsed`: the folds a band with none of its own shows. */
  collapsed: readonly string[] = [];
  /** Whether folds apply (tree order with a tree). */
  foldsApply = false;
  /** The folds applied (clade ids) and their names, as the band or the layer's `collapsed` gave them. */
  folds: readonly string[] = [];
  foldIds: number[] = [];
  measures: ColumnMeasures | undefined;
  measuresOf: HeldCells | undefined;
  channel: Channel | undefined;
  channels: LayerChannels | undefined;
  /** Each summary row's consensus (3 bytes a column), by display row; kept by the display's version. */
  readonly summaries = new Map<number, Uint8Array>();
  summariesOf = -1;
  /**
   * The rows' transform: px a row; the band's top (`top`, display rows below the pinned band, fractional: the band's
   * `bandTop`) and its scroll in px; the pinned band's height; the first and last display rows drawn.
   */
  rowScale = 20;
  top = 0;
  bandTop = 0;
  pinnedHeight = 0;
  first = 0;
  last = 0;
  /** Cells drawn as rects: x0, x1, y, h, fill, letter code, ink; `n` of them. */
  n = 0;
  x0 = new Float32Array(256);
  x1 = new Float32Array(256);
  y = new Float32Array(256);
  h = new Float32Array(256);
  fill = new Uint32Array(256);
  code = new Uint8Array(256);
  ink = new Uint32Array(256);
  letters = false;
  font = '10px mono';
  fontSize = 0;
  /** Rects from here on are the pinned rows' (drawn over a band that does not scroll). */
  pinnedFrom = 0;
  /** The band drawn as a neutral fill: lazy rows below a pixel (no row bins served). */
  neutral = false;
  /** Image tiles drawn: their boxes and images. */
  images: ImageData[] = [];
  ib: number[] = [];
  /** The last residue numbers in the right gutter: text, y, measured width; and their x. */
  lastText: string[] = [];
  lastY: number[] = [];
  lastW: number[] = [];
  numberRightX = 0;
  /** Each held row's last residue number (`NoNumber` for none), kept by the held data. */
  lasts = new Int32Array(0);
  lastOf: HeldCells | undefined;
  /** The columns in view (`[viewC0, viewC1)`, set by the cells' layout): the end numbers count the residues in them. */
  viewC0 = 0;
  viewC1 = 0;
  /** The widest last number this layout (the clip of the right gutter). */
  lastMaxW = 0;
  /** The residue numbers in view of the row last asked of `endNumbers`. */
  numFirst = 0;
  numLast = 0;
  /** First residue numbers of lazy rows read from their tiles this frame. */
  readonly lazyFirst = new Map<number, number>();
  /** Lazy rows' last residue numbers, counted once a row (`NoNumber` for none), kept by the source read. */
  readonly lazyLasts = new Map<number, number>();
  lazyLastsOf: unknown;
  /** The node each tree segment hits (an inner node; -1 a leaf's branch). */
  treeNode: number[] = [];
  /** Each segment's clade for the hover: the node it belongs to (a branch: the child it leads to; a vertical: its node). */
  treeClade: number[] = [];
  /** The tree's geometry by node, kept by the display's version (`treeGeometry`). */
  coord = new Float64Array(0);
  nodeX = new Float64Array(0);
  /** The x of a node's deepest leaf (by branch length, as its leaves would sit): a folded clade's wedge ends there. */
  deepX = new Float64Array(0);
  /** The display rows a node's clade takes (its shown rows; a folded clade's one): what the tree's resolution reads. */
  spanRows = new Float64Array(0);
  inFold = new Uint8Array(0);
  folded = new Uint8Array(0);
  treeOf = -1;
  treeClades: Clades | undefined;
  /** The hit a hover reads (one object, written each hit: a hover allocates nothing). */
  readonly scratch: ScratchHit = {
    layer: '',
    section: '',
    system: '',
    position: undefined,
    ySystem: undefined,
    yPosition: undefined,
    target: undefined,
    element: undefined,
    extent: undefined,
    state: undefined,
  };
  /** A clade's element and extent for its hits (a tree node's; a summary row's by display row), kept by the display. */
  readonly cladeHits = new Map<number, { element: Key; extent: Location }>();
  readonly summaryHits = new Map<number, { element: Key; extent: Location }>();
  cladeHitsOf = -1;
  /** The measures' registry entries over the theme's (kept by the held data and the theme's). */
  measureEnv: { measures: ReadonlyMap<string, MeasureEntry>; alphabet: Alphabet } | undefined;
  measureEnvOf: [HeldCells | undefined, unknown, Alphabet] | undefined;
  /** Header: labels (text, y, measured width), numbers, tree segments `[x0, y0, x1, y1]`, wedges `[x0, y0, x1, y1]`. */
  labelText: string[] = [];
  labelY: number[] = [];
  labelW: number[] = [];
  numberText: string[] = [];
  numberY: number[] = [];
  numberW: number[] = [];
  labelFont = '10px sans';
  labelX = 0;
  labelMax = 0;
  numberX = 0;
  tree: number[] = [];
  wedges: number[] = [];
  /** Tree nodes in view for hits: node, x, y. */
  nodes: number[] = [];
  treeX0 = 0;
  treeX1 = 0;
  /** Rows pending (y, h pairs): a placeholder each. */
  pending: number[] = [];
  /** The hidden-rows chip: its text and box (none when no row is hidden). */
  chip = '';
  chipX = 0;
  /** The tree's toggle in a narrow frame (its glyph's left and top; -1: none): its hit reopens the tree. */
  toggleX = -1;
  /** The toggle is drawn pressed: the reader opened the tree, and the toggle folds it back. */
  toggleOpen = false;
  toggleY = 0;
  chipY = 0;
  chipW = 0;
  readonly datum = new PlainDatum();
  /** A cell's fill by (column, letter): computed when first read since the channel or the data changed. */
  colors = new Uint32Array(0);
  known = new Uint8Array(0);
  colorsOf: unknown;
  /** Moves with each new color channel: what the reduction's classes are keyed by. */
  channelVersion = 0;
  /** The reduction over the display rows' theme classes, its class colors, and what it was built from. */
  pyramid: Pyramid | undefined;
  classColors = new Uint32Array(256);
  pyramidOf = '';
  baseBytes: Uint8Array | undefined;
  /** A class by (column, letter): 254 unknown yet. */
  classTable: Uint8Array | undefined;
  /** The shared image cache (the frame's, else the layer's own). */
  cache: ImageCache | undefined;
  /** The prep the display was built from: an edit of the order, the pins or the hidden rows builds it again. */
  prepOf: RowsPrep | undefined;
  /** The joined fields' version an order by a field was built with (-1: not by a field). */
  displayByField = -1;
  /**
   * Row columns (Row columns, `rows-columns.ts`); the fields joined from a records source by row key (kept by its
   * entries' versions, the held rows and the prep), `fieldsVersion` moving with them (an order by a joined field, the
   * columns' values, made again); `joinPending` while a source that could join has not answered.
   */
  readonly cols = new ColumnsLay();
  joined: Map<string, (string | number | null)[]> | undefined;
  readonly joinedOf: unknown[] = [];
  /** The layer sources the joined sources are read through, by source. */
  readonly joinSrc = new Map<string, LayerSource>();
  joinPending = false;
  fieldsVersion = 0;
  /** What joining said (a field shadowed, a field of two sources), made when the join is: said every frame. */
  joinIssue: Issue | undefined;
  /** The rows' band below the layer's title row (its own box, reused). */
  readonly ownBox = { x: 0, y: 0, width: 0, height: 0 };

  /**
   * A field of the rows table at a data row (_spec/README.md "Render", Row columns): `label`, `first`, the rows' own
   * fields, the measures (`length`, `gaps`, `identity`), then the joined ones; undefined when missing.
   */
  fieldValue(row: number, field: string): string | number | undefined {
    const held = this.held;
    if (field === 'label') return this.label(row);
    if (field === 'first') return held ? held.first[row] : this.lazyFirst.get(row);
    const own = held?.fields?.get(field);
    if (own) return own[row] ?? undefined;
    const m = this.measures;
    if (m && held && RowMeasureFields.includes(field)) {
      if (row >= m.rowLength.length) return undefined;
      if (field === 'length') return m.rowLength[row];
      if (field === 'gaps') return m.rowGaps[row];
      return Number.isNaN(m.rowIdentity[row]) ? undefined : m.rowIdentity[row];
    }
    const j = this.joined?.get(field);
    return j ? (j[row] ?? undefined) : undefined;
  }

  /** y of a display row's top in the band (the pinned band fixed, the rest scrolled). */
  rowY(i: number): number {
    const b = this.box;
    if (i < this.display.pinned) return b.y + i * this.rowScale;
    return b.y + this.pinnedHeight + (i - this.display.pinned) * this.rowScale - this.bandTop;
  }

  /** The display row under y in the band, or -1. */
  rowAtY(y: number): number {
    const b = this.box;
    if (y < b.y || y >= b.y + b.height) return -1;
    const pinned = this.display.pinned;
    let i: number;
    if (y < b.y + this.pinnedHeight) i = Math.floor((y - b.y) / this.rowScale);
    else i = pinned + Math.floor((y - b.y - this.pinnedHeight + this.bandTop) / this.rowScale);
    return i >= 0 && i < this.display.count ? i : -1;
  }

  /** The band's top at most (display rows below the pinned band): the last rows fill the band. */
  maxTop(rs: number, count: number, pinned = this.display.pinned): number {
    const area = Math.max(0, this.box.height - pinned * rs);
    return Math.max(0, count - pinned - area / rs);
  }

  /** A data row's key: a held row's label, else `#n` (a lazy row by its index). */
  keyOf(row: number): string {
    return this.held ? this.held.labels[row] : `#${row}`;
  }

  /** The data row of a key (`keyOf`), undefined when none has it. */
  rowOfKey(key: string): number | undefined {
    if (this.held) return this.byLabel.get(key);
    const n = key.startsWith('#') ? Number(key.slice(1)) : Number.NaN;
    return Number.isInteger(n) && n >= 0 && (this.total === undefined || n < this.total) ? n : undefined;
  }

  /**
   * The row at a band's top `top` (display rows below the pinned band) by its key, and the fraction of it scrolled
   * past (_spec/README.md "Render", The rows layer: a kept window): a summary row's first shown row. Undefined past
   * the rows.
   */
  keyAt(top: number): { key: string; offset: number } | undefined {
    const d = this.display;
    const at = Math.max(0, top);
    const i = d.pinned + Math.floor(at);
    if (i >= d.count) return undefined;
    let e = d.rowAt(i);
    if (e < 0) {
      const m = d.members.get(i);
      if (!m) return undefined;
      e = m[0];
    }
    return { key: this.keyOf(e), offset: at - Math.floor(at) };
  }

  /**
   * The band's top putting a key's row at the band's top (`keyAt`'s inverse), under `folds` (the ones shown when
   * absent) and `rowScale`, clamped: a folded row at its summary row, a hidden one at the nearest row shown, a top past
   * the last rows at the last. Undefined when no row has the key.
   */
  topOf(key: string, offset: number, folds?: readonly string[], rowScale?: number): number | undefined {
    const row = this.rowOfKey(key);
    if (row === undefined) return undefined;
    const d = folds !== undefined && folds !== this.folds && this.foldsApply ? this.displayFor(folds) : this.display;
    const i = d.nearestRow(row);
    if (i < 0) return undefined;
    const top = i < d.pinned ? 0 : i - d.pinned + offset;
    return Math.max(0, Math.min(top, this.maxTop(rowScale ?? this.rowScale, d.count, d.pinned)));
  }

  /**
   * The rows system as an axis in display order (_spec/README.md "Render", The rows layer: Shift+click extends in
   * display order): the data rows as the order shows them (a pinned row in its place, a folded clade's rows at its
   * summary row, hidden rows out), consecutive rows one range. Made once a display.
   */
  rowsAxis(): Axis {
    const d = this.display;
    if (this.axisOf === d.version && this.axisKept) return this.axisKept;
    const ranges: number[] = [];
    const push = (r: number) => {
      const n = ranges.length;
      if (n && ranges[n - 1] === r) ranges[n - 1] = r + 1;
      else ranges.push(r, r + 1);
    };
    if (d.identity) {
      if (d.count) ranges.push(0, d.count);
    } else {
      // the pinned band's rows are in their place too (`inPlace`), else only in the band
      const from = this.setArgs?.inPlace === false ? 0 : d.pinned;
      for (let i = from; i < d.count; i++) {
        const e = d.rowAt(i);
        if (e >= 0) push(e);
        else for (const r of d.members.get(i) ?? []) push(r);
      }
    }
    this.axisOf = d.version;
    this.axisKept = {
      system: this.rowsSystem,
      ...(this.total !== undefined ? { length: this.total } : {}),
      cyclic: false,
      ranges,
      dir: 1,
      pieces: 'order',
    };
    return this.axisKept;
  }
  axisOf = -1;
  axisKept: Axis | undefined;

  /** The display rows the layer would show under other folds (its order, pins and hidden rows kept). */
  displayFor(folds: readonly string[]): DisplayRows {
    const a = this.setArgs;
    const out = new DisplayRows();
    if (!a || !this.clades) return this.display;
    const ids = folds.map((f) => this.clades!.byName.get(f) ?? -1).filter((c) => c >= 0);
    out.set(a.n, a.order, this.clades, ids, a.pinned, a.hidden, a.inPlace);
    return out;
  }

  /** What the display was last built from (`displayFor` builds another under other folds). */
  setArgs:
    | {
        n: number;
        order: Int32Array | undefined;
        pinned: number[];
        hidden: ReadonlySet<number>;
        inPlace: boolean;
      }
    | undefined;
  /** The held rows by label (the first of a label), kept by the held data. */
  byLabel = new Map<string, number>();
  byLabelOf: HeldCells | undefined;

  /** The label of a data row: held, else from a loaded tile, else `row n`. */
  label(row: number): string {
    return this.held?.labels[row] ?? this.lazyLabels.get(row) ?? `row ${(row + 1).toLocaleString('en-US')}`;
  }

  /** A summary row's count as its label ends (` +n`: its shown rows past the first); '' for a row of its own. */
  summaryMore(i: number): string {
    const m = this.display.members.get(i);
    return m ? ` +${m.length - 1}` : '';
  }

  /** A display row's label whole: a summary row's first row's and its count (`HBA +1`), else its row's. */
  displayLabel(i: number): string {
    const e = this.display.rowAt(i);
    if (e >= 0) return this.label(e);
    const m = this.display.members.get(i);
    return m ? `${this.label(m[0])}${this.summaryMore(i)}` : '';
  }

  /**
   * The display row the keyboard cursor stands on for a data row: its own (a folded row's summary row; a pinned row's
   * in the band), else, hidden or filtered, the nearest row shown; -1 when none is.
   */
  cursorRow(row: number): number {
    const i = this.display.displayRow(row);
    return i >= 0 ? i : this.display.nearestRow(row);
  }

  /** The display row of the last hit on the rows (a cell's, a row's): which copy of a pinned row a click was on. */
  hitRow = -1;

  /** Labels and first numbers of lazy rows read from their tiles this frame (cleared each layout). */
  readonly lazyLabels = new Map<number, string>();
  readonly lazyCodes = new Map<number, Uint8Array>();

  /** Rebuilds what depends on the held data, the order and the folds when they changed. */
  prepareRows(prep: RowsPrep, ctx: LayoutContext): void {
    const held = this.held;
    // the clades: the tree's leaves matched to the held rows' labels (EC-43)
    if (prep.tree && held && this.cladesOf !== held) {
      this.cladesOf = held;
      const rows = held.labels.map((label) => ({ label }));
      this.clades = buildClades(prep.tree.read, matchLeaves(prep.tree.read, rows), held.rows);
    }
    const folds = this.figure ? NoFolds : (ctx.folds ?? prep.collapsed);
    // an order by a field is built again when the joined fields come; an edit of the layer (another prep) too
    const byField = typeof prep.order === 'object' ? this.fieldsVersion : -1;
    const key = `${held?.version ?? -1}|${this.total}|${this.figure}`;
    if (key === this.displayKey && folds === this.folds && prep === this.prepOf && byField === this.displayByField)
      return;
    this.displayKey = key;
    this.displayByField = byField;
    this.folds = folds;
    this.prepOf = prep;
    if (!held) {
      this.display.setIdentity(this.total ?? 0);
      this.foldIds = [];
      return;
    }
    if (this.byLabelOf !== held) {
      this.byLabelOf = held;
      this.byLabel = new Map();
      held.labels.forEach((l, i) => {
        if (!this.byLabel.has(l)) this.byLabel.set(l, i);
      });
    }
    const byLabel = this.byLabel;
    const rowsOf = (labels: readonly string[]) => labels.map((l) => byLabel.get(l) ?? -1).filter((r) => r >= 0);
    const tree = prep.order === 'tree' ? this.clades : undefined;
    this.foldsApply = tree !== undefined;
    this.foldIds = tree ? folds.map((f) => tree.byName.get(f) ?? -1).filter((c) => c >= 0) : [];
    let order: Int32Array | undefined;
    if (typeof prep.order === 'object') order = fieldOrder(this, held, prep.order);
    const plain = !tree && !order && !prep.pinned.length && !prep.hidden.length;
    this.setArgs = undefined;
    if (plain) this.display.setIdentity(held.rows);
    else {
      this.setArgs = {
        n: held.rows,
        order,
        pinned: rowsOf(prep.pinned),
        hidden: new Set(rowsOf(prep.hidden)),
        inPlace: ctx.config.rows.pinInPlace,
      };
      const a = this.setArgs;
      this.display.set(a.n, a.order, tree, this.foldIds, a.pinned, a.hidden, a.inPlace);
    }
  }

  /** The measures of the held data (the columns' and the rows'), when the data changed. */
  prepareMeasures(prep: RowsPrep): void {
    const held = this.held;
    if (held && this.measuresOf !== held) {
      this.measuresOf = held;
      const ref = prep.reference === undefined ? -1 : held.labels.indexOf(prep.reference);
      this.measures = columnMeasures(held, ref);
      // a new channel set: its measures are this data's
      this.channels = channelForms(prep.channelsForms);
    }
  }

  /** The measures and the color channel (the measures registered for its themes), when the data changed. */
  prepareColors(prep: RowsPrep, ctx: LayoutContext): Channel {
    const held = this.held;
    this.prepareMeasures(prep);
    const channels = this.channels ?? (this.channels = channelForms(prep.channelsForms));
    // the column measures registered over the theme's own (merged, never replacing them), kept while both are
    const theme = ctx.theme.measures;
    const of = this.measureEnvOf;
    if (!of || of[0] !== held || of[1] !== theme || of[2] !== prep.alphabet) {
      this.measureEnvOf = [held, theme, prep.alphabet];
      const measures = new Map(theme ?? []);
      if (held && this.measures) for (const [k, v] of measureEntries(held, this.measures)) measures.set(k, v);
      this.measureEnv = { alphabet: prep.alphabet, measures };
    }
    const channel = channels.get('color', LettersColorFallback[prep.alphabet], ctx, this.measureEnv);
    if (this.colorsOf !== channel || this.colors.length < this.columns * 128) {
      this.colorsOf = channel;
      this.channelVersion++;
      this.colors = new Uint32Array(Math.max(this.columns, 1) * 128);
      this.known = new Uint8Array(Math.max(this.columns, 1) * 128);
    }
    this.channel = channel;
    return channel;
  }

  /** A cell's fill: its letter's (uppercase) at its column, by the color channel; none for a gap. */
  colorOf(code: number, column: number, missing: number): number {
    if (isGap(code) || code > 127) return NoColor;
    const k = column * 128 + upper(code);
    if (this.known[k]) return this.colors[k];
    const ch = this.channel!;
    const d = this.datum;
    d.letter = String.fromCharCode(upper(code));
    d.start = column;
    d.end = column + 1;
    const c = ch.constant ?? ch.at(d);
    const fill = typeof c === 'number' ? c : missing;
    this.colors[k] = fill;
    this.known[k] = 1;
    return fill;
  }
}

/**
 * A permutation of the held rows by a field of the rows table (`label`, `first`, an own field, a joined one, a
 * measure): numbers by value, text by code points, missing last in either direction, ties in data order.
 */
function fieldOrder(
  lay: RowsLay,
  held: HeldCells,
  o: { readonly field: string; readonly direction: string },
): Int32Array {
  const n = held.rows;
  const values: (string | number | undefined)[] = new Array(n);
  for (let i = 0; i < n; i++) values[i] = lay.fieldValue(i, o.field);
  // a field whose values all read as numbers orders by value (a `#=GS` TAXID, an untyped table's field)
  const nums = numericValues(values.map((v) => v ?? null));
  if (nums) for (let i = 0; i < n; i++) values[i] = nums[i] ?? undefined;
  const idx = Array.from({ length: n }, (_, i) => i);
  const sign = o.direction === 'descending' ? -1 : 1;
  idx.sort((a, b) => {
    const x = values[a];
    const y = values[b];
    if (x === undefined || y === undefined) return x === y ? a - b : x === undefined ? 1 : -1;
    const c =
      typeof x === 'number' && typeof y === 'number'
        ? x - y
        : String(x) < String(y)
          ? -1
          : String(x) > String(y)
            ? 1
            : 0;
    return sign * c || a - b;
  });
  return Int32Array.from(idx);
}

export const NoFolds: readonly string[] = [];
