/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  CategoryColumn,
  type Column,
  exactAddress,
  type Issue,
  issueOf,
  type Node,
  type Value,
} from '@molstar/core/repr/tree';
import { type CellsCodes, DefaultMaxBlockValues, type RangesList } from '@molstar/seqx-core/data';
import { matchLeaves, readNewick } from '@molstar/seqx-core/data/alignment-readers';
import { type ElementKey, type Location, rangesLocation } from '@molstar/seqx-core/query';
import { rowsSystemOf, type StaticModel } from '@molstar/seqx-core/resolve';
import type { Channel, MeasureEntry } from '@molstar/seqx-core/theme/index';
import type { Axis, ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, channelLegends, type LayerChannels, PlainDatum } from './channels.js';
import { fontSize, type RenderConfig } from './config.js';
import {
  Anchor,
  drawPass,
  grown,
  type ImageData,
  NoColor,
  type PassName,
  type PassWriter,
  passDraws,
  Round,
} from './draw-lists.js';
import { ImageCache, ImageTile } from './images.js';
import { rowsLevel } from './rows-levels.js';
import { type Alphabet, cellEdge, LettersColorFallback } from './letters.js';
import { RenderIssues } from './issues.js';
import { type LayerSource, layerSource, readTiles } from './layer-data.js';
import {
  type Box,
  type Emphasis,
  type Hit,
  type HitTarget,
  type LayerType,
  type LayoutContext,
  type Overlay,
  sectionId,
  type StoreReader,
  type ItemDescription,
} from './layer-type.js';
import { NoClass, Pyramid, pyramidLevelFor } from './pyramid.js';
import {
  ColumnsLay,
  chipAt,
  drawColumnCells,
  drawTitles,
  elide,
  hoverColumnCells,
  layoutRowColumns,
  setColumns,
  titleAt,
} from './rows-columns.js';
import {
  buildClades,
  type Clades,
  type ColumnMeasures,
  columnMeasures,
  columnsRoom,
  consensus,
  DisplayRows,
  type HeldCells,
  isGap,
  measureEntries,
  RowHeight,
  type RowColumn,
  rowColumnsOf,
  RowMeasureFields,
  numericValues,
  TreeHitPx,
  upper,
} from './rows-model.js';
import { squircleRadius } from './squircle.js';
import { inkOn, type SeqTokens } from './tokens.js';

// The `rows` layer type (_spec/README.md "Render", The rows layer): an alignment's rows over its columns, in data,
// tree or a field's order, clades folded to a summary row, rows pinned in a band at the top, hidden rows out; its own
// y inside its band (`bandTop` the rows' scroll, `rowScale` px a row); its header (the tree, the labels, the left row
// columns, the first residue numbers, the hidden-rows chip) in the view's label column, its right columns after the
// last numbers, the columns' titles in the band above the rows. Held alignments (rows × columns within
// `maxBlockValues`) are read whole once; larger ones lazily, the rows in view only. Below 3 device px a cell on either
// axis the rows draw from the shared reduction as image tiles in the shared image cache.

const AnchorStart: number = Anchor.start;
const AnchorEnd: number = Anchor.end;

/** A hit the layer writes in place (`RowsLay.scratch`): every field set each time, so its shape never changes. */
type ScratchHit = { -readonly [K in keyof Hit]: Hit[K] };

/** No residue number (a row of gaps). */
const NoNumber = -0x80000000;

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
  readonly cladeHits = new Map<number, { element: ElementKey; extent: Location }>();
  readonly summaryHits = new Map<number, { element: ElementKey; extent: Location }>();
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

/** A column's value at `i` (`null` missing: a missing category, a NaN). */
function cellOf(c: Column, i: number): string | number | null {
  if (c instanceof CategoryColumn) return c.get(i);
  const v = (c as ArrayLike<number>)[i];
  return v === undefined || Number.isNaN(v) ? null : v;
}

/**
 * The fields of a records source whose system is the rows system, joined by row key (_spec/README.md "Render", Row
 * columns): a record's `id` names its row (a held row's label), else its position is the row. Only the fields a column
 * or the order names, and neither an own field nor a measure (those shadow it: said once, as a field two sources
 * have). Kept by the held rows, the prep and the entries' versions (compared in place: nothing made once warm),
 * `fieldsVersion` moved when they change. `joinPending` while a source that could join (its `system` the rows system,
 * a runtime bound) has not attached and the view has not settled: one that failed, or never attaches, holds nothing.
 */
function readJoined(lay: RowsLay, prep: Prep, ctx: LayoutContext, data: StoreReader) {
  const held = lay.held;
  lay.joinPending = false;
  if (!held || !prep.joinSources.length || !prep.joinFields.length) return;
  const of = lay.joinedOf;
  let same = of[0] === held && of[1] === prep;
  let at = 2;
  let nt = 0;
  for (let s = 0; s < prep.joinSources.length; s++) {
    const source = prep.joinSources[s];
    const coverage = data.coverage(source);
    if (!coverage) {
      if (!ctx.settled && prep.joinCandidates.includes(source)) lay.joinPending = true;
      continue;
    }
    let length = -1;
    for (let k = 0; k < coverage.systems.length; k++)
      if (coverage.systems[k].name === prep.rowsSystem) length = coverage.systems[k].length ?? held.rows;
    if (length < 0) continue;
    JoinRegion[1] = length;
    const read = readTiles(ctx, data, joinSource(lay, source, prep.rowsSystem, length), JoinRegion, 1);
    if (read.pending.length) lay.joinPending = true;
    for (let k = 0; k < read.entries.length; k++) {
      const e = read.entries[k];
      const t = e.tables.ranges as RangesList | undefined;
      if (t?.form !== 'list') continue;
      if (of[at] !== source || of[at + 1] !== t || of[at + 2] !== e.version) same = false;
      of[at] = source;
      of[at + 1] = t;
      of[at + 2] = e.version;
      at += 3;
      nt++;
    }
  }
  if (same && of.length === at) return;
  of[0] = held;
  of[1] = prep;
  of.length = at;
  lay.fieldsVersion++;
  lay.joinIssue = undefined;
  if (!nt) {
    lay.joined = undefined;
    return;
  }
  const joined = new Map<string, (string | number | null)[]>();
  const from = new Map<string, string>();
  const said = new Set<string>();
  const byLabel = new Map<string, number>();
  held.labels.forEach((l, i) => {
    if (!byLabel.has(l)) byLabel.set(l, i);
  });
  for (let q = 2; q < at; q += 3) {
    const source = of[q] as string;
    const t = of[q + 1] as RangesList;
    for (const f of t.fields ?? []) {
      if (!prep.joinFields.includes(f.name)) continue;
      if (held.fields?.has(f.name) || RowMeasureFields.includes(f.name)) {
        said.add(`${f.name} of ${source} is the alignment’s own ${held.fields?.has(f.name) ? 'field' : 'measure'}`);
        continue;
      }
      const first = from.get(f.name);
      if (first !== undefined && first !== source) {
        said.add(`${f.name} is a field of ${first} and of ${source}: ${first}’s read`);
        continue;
      }
      from.set(f.name, source);
      let a = joined.get(f.name);
      if (!a) joined.set(f.name, (a = new Array(held.rows).fill(null)));
      for (let k = 0; k < t.start.length; k++) {
        const id = t.id?.get(k);
        const row = id != null ? byLabel.get(id) : t.start[k];
        if (row === undefined || row < 0 || row >= held.rows) continue;
        a[row] = cellOf(f.column, k);
      }
    }
  }
  lay.joined = joined;
  if (said.size)
    lay.joinIssue = issueOf(
      RenderIssues,
      'seq.render.rows.field',
      `A joined field is not read: ${[...said].join('; ')}`,
      {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.rows.field'].hint,
      },
    );
}

/** The region a join reads: the whole rows system (its end set before each read: a frame makes no array). */
const JoinRegion = [0, 0];

/** The layer source a joined source is read through, kept by the layer (one a source and length). */
function joinSource(lay: RowsLay, source: string, system: string, length: number): LayerSource {
  let s = lay.joinSrc.get(source);
  if (!s || s.length !== length) {
    s = { source, system, dataSystem: system, shape: 'ranges', length, ratio: 1 };
    lay.joinSrc.set(source, s);
  }
  return s;
}

const NoFolds: readonly string[] = [];

const p = (layer: ResolvedLayer, key: string): Value | undefined => layer.params[key];
const strings = (v: Value | undefined): string[] =>
  Array.isArray(v) ? (v.filter((x) => typeof x === 'string') as string[]) : [];

/** A node of the spec by its exact address (`/data/tree`), for the layer's `tree` param. */
function nodeAt(spec: Node, address: string): Node | undefined {
  const ids = address.split('/').filter(Boolean);
  let n: Node | undefined = spec;
  for (const id of ids) n = n?.children?.find((c) => c.id === id);
  return n;
}

/** The layer's prepared form: its source, its params, its tree read once. */
function prepareRows(layer: ResolvedLayer, model: StaticModel): RowsPrep {
  const src = layerSource(model, layer);
  const owner = layer.data?.address ?? '';
  const a = model.annotations.get(owner);
  const system = model.systems.get(layer.data?.system ?? '');
  const order = p(layer, 'order');
  const treeAt = exactAddress(p(layer, 'tree'));
  const treeNode = treeAt === undefined ? undefined : nodeAt(model.spec, treeAt);
  const newick = treeNode?.params?.newick;
  const summary = p(layer, 'summary');
  const width = p(layer, 'tree_width');
  const columns = rowColumnsOf(p(layer, 'columns'));
  // a records source joined by row key: one whose coverage names the rows system (read at layout); none is looked for
  // unless a column or an order may read one
  const byField = !!order && typeof order === 'object' && !Array.isArray(order);
  const joinSources =
    columns.length || byField
      ? [...model.sources.keys()].filter((s) => s !== src?.source && s !== system?.source)
      : NoSources;
  const rowsSystem = a?.rows ?? rowsSystemOf(owner);
  // the sources that could join (pending until they attach): their own system the rows system, a runtime bound
  const joinCandidates = joinSources.filter((s) => {
    const b = model.sources.get(s);
    return b?.runtime !== undefined && b.entry.system === rowsSystem;
  });
  const joinFields = [
    ...new Set([...columns.map((c) => c.field), ...(byField ? [(order as { field: string }).field] : [])]),
  ];
  return {
    ...(src ? { src } : {}),
    owner,
    params: layer.params,
    rowsSystem: a?.rows ?? rowsSystemOf(owner),
    alphabet: system?.alphabet ?? 'protein',
    channels: channelForms(layer.theme),
    channelsForms: layer.theme,
    ...(treeAt !== undefined && typeof newick === 'string'
      ? { tree: { address: treeAt, read: readNewick(newick) } }
      : {}),
    order:
      order === 'tree' && typeof newick === 'string'
        ? 'tree'
        : order &&
            typeof order === 'object' &&
            !Array.isArray(order) &&
            typeof (order as { field?: unknown }).field === 'string'
          ? {
              field: (order as { field: string }).field,
              direction: (order as { direction?: string }).direction === 'descending' ? 'descending' : 'ascending',
            }
          : 'data',
    collapsed: strings(p(layer, 'collapsed')),
    pinned: strings(p(layer, 'pinned')),
    hidden: strings(p(layer, 'hidden_rows')),
    summary: summary === 'bar' ? 'bar' : 'consensus',
    numbers: p(layer, 'numbers') === 'none' ? 'none' : 'gutters',
    treeWidth: typeof width === 'number' ? width : 120,
    ...(a?.reference !== undefined ? { reference: a.reference } : {}),
    columns,
    joinSources,
    joinCandidates,
    joinFields,
  };
}

const NoSources: readonly string[] = [];

type Prep = RowsPrep;

/** The layer's own image cache when the frame gives none (a test's lone layer). */
const OwnImages = new WeakMap<RenderConfig, ImageCache>();

export const RowsType: LayerType<Prep, RowsLay> = {
  name: 'rows',
  axes: 2,

  prepare: prepareRows,

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: RowsLay): RowsLay {
    const { axis, config, tokens } = ctx;
    const lay = prev ?? new RowsLay();
    lay.key = ctx.key;
    lay.section = sectionId(ctx.section);
    lay.system = axis.system;
    lay.rowsSystem = prep.rowsSystem;
    lay.owner = prep.owner;
    lay.params = prep.params;
    lay.collapsed = prep.collapsed;
    lay.box = ctx.box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.config = config;
    lay.n = 0;
    lay.pinnedFrom = 0;
    lay.neutral = false;
    lay.letters = false;
    lay.images.length = 0;
    lay.ib.length = 0;
    lay.labelText.length = 0;
    lay.labelY.length = 0;
    lay.labelW.length = 0;
    lay.numberText.length = 0;
    lay.numberY.length = 0;
    lay.numberW.length = 0;
    lay.lastText.length = 0;
    lay.lastY.length = 0;
    lay.lastW.length = 0;
    lay.lastMaxW = 0;
    lay.viewC0 = 0;
    lay.viewC1 = 0;
    lay.tree.length = 0;
    lay.treeNode.length = 0;
    lay.treeClade.length = 0;
    lay.wedges.length = 0;
    lay.nodes.length = 0;
    lay.pending.length = 0;
    lay.lazyLabels.clear();
    lay.lazyCodes.clear();
    lay.lazyFirst.clear();
    lay.chip = '';
    // the row columns as written and granted; a layer with columns keeps a title row above its rows
    setColumns(lay.cols, prep, ctx);
    if (lay.cols.titleRow) {
      const b = lay.ownBox;
      b.x = ctx.box.x;
      b.y = ctx.box.y + lay.cols.titleRow;
      b.width = ctx.box.width;
      b.height = Math.max(0, ctx.box.height - lay.cols.titleRow);
      lay.box = b;
    }
    const box = lay.box;
    if (ctx.images) lay.cache = ctx.images;
    else if (!lay.cache) {
      let own = OwnImages.get(config);
      if (!own) OwnImages.set(config, (own = new ImageCache(config.imageBytes)));
      lay.cache = own;
    }
    const src = prep.src;
    if (!src) return lay;
    const at = rowsSource(src, data);
    const rowsSrc = at.src;
    if (lay.lazyLastsOf !== src.source) {
      lay.lazyLastsOf = src.source;
      lay.lazyLasts.clear();
    }
    lay.total = at.rows;
    lay.columns = at.columns ?? axis.length ?? 0;
    // a band that does not scroll (text: a figure's lines repeat the rows): every row in the band, unfolded, no tree;
    // a layout of neither (rings, blocks) draws no rows
    const layout = ctx.section.layout;
    if (layout !== 'tracks' && layout !== 'text') {
      ctx.report(
        issueOf(RenderIssues, 'seq.render.rows.layout', 'A rows layer draws in tracks and text sections only', {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.rows.layout'].hint,
        }),
      );
      return lay;
    }
    lay.figure = layout === 'text';
    lay.rowScale = ctx.rowScale ?? RowHeight;
    const N = lay.total;
    if (N === undefined) {
      // not attached yet: the band is pending
      readTiles(ctx, data, rowsSrc, [0, 1], 1);
      lay.pending.push(box.y, box.height);
      return lay;
    }
    // held: read whole once (within the block bound), joined when every tile is in
    const heldable = N * lay.columns <= DefaultMaxBlockValues;
    if (heldable) {
      const held = readHeld(ctx, data, rowsSrc, N, lay.columns);
      if (held) lay.held = held;
      else if (!lay.held) lay.pending.push(box.y, box.height);
    } else {
      lay.held = undefined;
      // ordering needs every row's field: a lazy alignment draws in data order
      if (prep.order !== 'data')
        ctx.report(
          issueOf(RenderIssues, 'seq.render.rows.order', 'A lazy alignment is drawn in data order', {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.rows.order'].hint,
          }),
        );
      // pins and hidden rows are by key, and a lazy row's key is read with its tile: neither applies
      if (prep.pinned.length || prep.hidden.length)
        ctx.report(
          issueOf(RenderIssues, 'seq.render.rows.lazy', 'A lazy alignment’s pinned and hidden rows are not applied', {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.rows.lazy'].hint,
          }),
        );
    }
    // the fields an order or a column reads: the measures and the joined records, before the order is built
    lay.prepareMeasures(prep);
    readJoined(lay, prep, ctx, data);
    if (lay.joinIssue) ctx.report(lay.joinIssue);
    lay.prepareRows(prep, ctx);
    // a spec's fold naming no clade of the tree is said, never dropped silently (the session's folds are the tree's)
    const clades = lay.foldsApply ? lay.clades : undefined;
    const unknown = clades && prep.collapsed.length ? prep.collapsed.filter((n) => !clades.byName.has(n)) : NoFolds;
    if (unknown.length)
      ctx.report(
        issueOf(RenderIssues, 'seq.render.rows.clade', `No clade of the tree is named ${unknown.join(', ')}`, {
          address: ctx.key,
          hint: RenderIssues.codes['seq.render.rows.clade'].hint,
        }),
      );
    const display = lay.display;
    if (lay.figure) {
      lay.rowScale = box.height / Math.max(1, display.count);
      if ((ctx.folds ?? prep.collapsed).length)
        ctx.report(
          issueOf(RenderIssues, 'seq.render.rows.text.unfolded', 'Every clade is drawn unfolded in a text figure', {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.rows.text.unfolded'].hint,
          }),
        );
    }
    // the rows' transform: the pinned band fixed, the rest scrolled within the band; the band's top in display rows
    // below the pinned band (fractional), clamped, turned into px here
    const rs = lay.rowScale;
    lay.pinnedHeight = display.pinned * rs;
    const area = Math.max(0, box.height - lay.pinnedHeight);
    lay.top = lay.figure ? 0 : Math.max(0, Math.min(ctx.bandTop ?? 0, lay.maxTop(rs, display.count)));
    lay.bandTop = lay.top * rs;
    lay.first = display.pinned + Math.floor(lay.bandTop / rs);
    lay.last = Math.min(display.count, display.pinned + Math.ceil((lay.bandTop + area) / rs));
    if (!heldable) readLazy(lay, ctx, data, rowsSrc);
    lay.prepareColors(prep, ctx);
    if (!lay.held && heldable) return lay;
    layoutCells(lay, prep, ctx);
    layoutHeader(lay, prep, ctx);
    if (prep.columns.length) layoutRowColumns(lay, prep, ctx);
    if (display.hidden) {
      lay.chip = `${display.hidden} hidden`;
      lay.chipW = ctx.text.measure(config.itemFont, lay.chip);
      lay.chipX = config.inset;
      lay.chipY = box.y + box.height - 3;
    }
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(RowsDraws, lay, out, emphasis);
  },

  legend(lay, options) {
    if (!lay.channel) return [];
    // the letters drawn as cells; zoomed out (image tiles) the alphabet's own (an RNA's U, never its T)
    let painted = '';
    const seen = new Set<number>();
    for (let k = 0; k < lay.n; k++) if (lay.code[k] > 32) seen.add(lay.code[k]);
    for (const c of seen) painted += String.fromCharCode(c);
    const alphabet = lay.measureEnv?.alphabet;
    if (lay.images.length || !painted) {
      if (alphabet === 'rna') painted += 'ACGU';
      else if (alphabet === 'dna') painted += 'ACGT';
      else if (lay.images.length) painted = '';
    }
    return channelLegends([lay.channel], {
      ...options,
      layer: lay.key,
      ...(painted ? { painted: painted.toUpperCase() } : {}),
    });
  },

  hitOutside(lay, at): Hit | undefined {
    // a column's title in the layer's title row: `state` its field (a click orders the rows by it); the hidden-columns
    // chip names what it hides (`state` its text: the hover card lists them)
    const t = titleAt(lay, at[0], at[1]);
    if (t < 0) {
      if (!chipAt(lay, at[0], at[1])) return undefined;
      const h = writeHit(lay, lay.rowsSystem, undefined, undefined, 'columns-hidden');
      lay.scratch.state = lay.cols.chip;
      return h;
    }
    const h = writeHit(lay, lay.rowsSystem, undefined, undefined, 'column-title');
    lay.scratch.state = lay.cols.columns[lay.cols.titleCol[t]].field;
    return h;
  },

  items: rowItems,

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const i = lay.rowAtY(y);
    if (x < box.x) {
      // the header: a tree node within reach of its branch, else the row's label
      const node = nearestNode(lay, x, y);
      if (node >= 0 && lay.clades)
        return writeHit(lay, lay.rowsSystem, undefined, undefined, 'tree-node', cladeHit(lay, node));
      if (i < 0 || x < lay.treeX1) return undefined;
      return named(lay, rowHit(lay, i, undefined), lay.cols.cellAt(i, x));
    }
    if (i < 0) return undefined;
    // a right column's cell hits as its row's header, naming its value
    if (x >= box.x + box.width) {
      const q = lay.cols.cellAt(i, x);
      return q < 0 ? undefined : named(lay, rowHit(lay, i, undefined), q);
    }
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return rowHit(lay, i, position);
  },
};

/** A hit on display row `i`: a summary row's (its clade), else its data row's, a cell where `position` is given. */
function rowHit(lay: RowsLay, i: number, position: number | undefined): Hit {
  lay.hitRow = i;
  const e = lay.display.rowAt(i);
  if (e < 0 && lay.clades) {
    const c = -1 - e;
    return writeHit(
      lay,
      position === undefined ? lay.rowsSystem : lay.system,
      position,
      undefined,
      'summary-row',
      cladeHit(lay, c, i),
    );
  }
  if (position === undefined) return writeHit(lay, lay.rowsSystem, e, e, 'row');
  return writeHit(lay, lay.system, position, e, 'cell');
}

/**
 * The inner tree node whose branch is nearest a point within `TreeHitPx` (by the distance to the segment), a folded
 * clade's wedge included; -1 for none.
 */
function nearestNode(lay: RowsLay, x: number, y: number): number {
  const reach = TreeHitPx;
  // the target stops at the label column: a point over a row's label is the row's (a cherry at the tree's edge would
  // take it, the reach being more than the 4 px between the tree and the labels)
  if (x < lay.treeX0 - reach || x > lay.treeX1 + reach || x >= lay.labelX) return -1;
  const w = lay.wedges;
  for (let k = 0; k + 4 < w.length; k += 5)
    if (x >= w[k] - reach && x <= w[k + 2] && y >= w[k + 1] && y <= w[k + 3]) return w[k + 4];
  let best = -1;
  let bestD = reach;
  const t = lay.tree;
  const owner = lay.treeNode;
  for (let k = 0, s = 0; k + 3 < t.length; k += 4, s++) {
    const node = owner[s];
    if (node < 0) continue;
    const xa = t[k] < t[k + 2] ? t[k] : t[k + 2];
    const xb = t[k] < t[k + 2] ? t[k + 2] : t[k];
    const ya = t[k + 1] < t[k + 3] ? t[k + 1] : t[k + 3];
    const yb = t[k + 1] < t[k + 3] ? t[k + 3] : t[k + 1];
    const dx = x < xa ? xa - x : x > xb ? x - xb : 0;
    const dy = y < ya ? ya - y : y > yb ? y - yb : 0;
    const d = Math.sqrt(dx * dx + dy * dy);
    if (d <= bestD) {
      bestD = d;
      best = node;
    }
  }
  return best;
}

/** Joined held cells by store and source, shared by the layers reading one alignment (the rows, its logo). */
const Joined = new WeakMap<StoreReader, Map<string, { key: string; held: HeldCells }>>();

/**
 * An alignment held whole (_spec/README.md "Render", The rows layer): every tile of its rows wanted at level 1, joined
 * once they are all in, kept by the entries' versions and shared by the layers reading it on one store; undefined while
 * some are pending or failed.
 */
export function readHeld(ctx: LayoutContext, data: StoreReader, src: LayerSource, rows: number, columns: number) {
  const read = readTiles(ctx, data, src, [0, rows], 1);
  if (read.pending.length || read.failed.length) return undefined;
  let key = `${rows}|${columns}`;
  for (const e of read.entries) key += `,${e.version}`;
  let bySource = Joined.get(data);
  if (!bySource) Joined.set(data, (bySource = new Map()));
  const kept = bySource.get(src.source);
  if (kept?.key === key) return kept.held;
  const C = columns;
  const codes = new Uint8Array(rows * C);
  const labels: string[] = new Array(rows).fill('');
  const first = new Int32Array(rows).fill(1);
  let fields: Map<string, (string | number | null)[]> | undefined;
  for (const e of read.entries) {
    const t = e.tables.cells;
    if (!t || t.form !== 'codes' || t.columns !== C) continue;
    const n = Math.floor(t.codes.length / Math.max(C, 1));
    codes.set(t.codes.subarray(0, Math.min(n, rows - t.start) * C), t.start * C);
    for (let i = 0; i < n && t.start + i < rows; i++) {
      labels[t.start + i] = t.label?.get(i) ?? `row ${t.start + i + 1}`;
      if (t.first) first[t.start + i] = t.first[i];
    }
    // the rows' own fields (Row columns), a value a row
    for (const f of t.fields ?? []) {
      fields ??= new Map();
      let a = fields.get(f.name);
      if (!a) fields.set(f.name, (a = new Array(rows).fill(null)));
      for (let i = 0; i < n && t.start + i < rows; i++) a[t.start + i] = cellOf(f.column, i);
    }
  }
  const held: HeldCells = {
    rows,
    columns: C,
    codes,
    labels,
    first,
    version: (kept?.held.version ?? 0) + 1,
    ...(fields ? { fields } : {}),
  };
  bySource.set(src.source, { key, held });
  return held;
}

/** The rows system's source of a layer (its data's columns path dropped) and the alignment's size, when attached. */
export function rowsSource(src: LayerSource, data: StoreReader) {
  const coverage = data.coverage(src.source);
  const rows = coverage?.systems.find((s) => s.name === src.system)?.length;
  const columns = coverage?.systems.find((s) => s.name !== src.system)?.length;
  const { path: _p, back: _b, ...rest } = src;
  return { src: { ...rest, ...(rows !== undefined ? { length: rows } : {}) } as LayerSource, rows, columns };
}

/** A lazy alignment: the tiles of the rows in view at level 1 (their labels and cells), the rest pending. */
function readLazy(lay: RowsLay, ctx: LayoutContext, data: StoreReader, src: LayerSource) {
  const { first, last } = lay;
  if (last <= first) return;
  // bounded by the screen: a band of more rows than device pixels reads none; the band a neutral fill (a source's
  // coarse row levels, later)
  if ((last - first) * lay.rowScale * ctx.dpr < (last - first) * 0.99) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.rows.reduced', 'The rows are below a pixel and the source serves no row bins', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.rows.reduced'].hint,
      }),
    );
    lay.neutral = true;
    return;
  }
  const read = readTiles(ctx, data, src, [first, last], 1);
  for (const e of read.entries) {
    const t = e.tables.cells as CellsCodes | undefined;
    if (!t || t.form !== 'codes') continue;
    const C = t.columns;
    const rows = Math.floor(t.codes.length / Math.max(C, 1));
    const a = Math.max(first, t.start);
    const b = Math.min(last, t.start + rows);
    for (let r = a; r < b; r++) {
      lay.lazyLabels.set(r, t.label?.get(r - t.start) ?? `row ${r + 1}`);
      lay.lazyCodes.set(r, t.codes.subarray((r - t.start) * C, (r - t.start + 1) * C));
      if (t.first) lay.lazyFirst.set(r, t.first[r - t.start]);
    }
  }
  for (let i = 0; i + 1 < read.pending.length; i += 2)
    for (let r = read.pending[i]; r < read.pending[i + 1]; r++) {
      lay.pending.push(lay.rowY(r), lay.rowScale);
    }
}

/**
 * The cells in view: rects (letters at the letters step), or image tiles from the reduction below 3 device px; the
 * pinned rows after the scrolled ones (`pinnedFrom`), drawn on their own over a band that does not scroll.
 */
function layoutCells(lay: RowsLay, prep: Prep, ctx: LayoutContext) {
  const { axis, config, tokens } = ctx;
  const dpr = axis.dpr;
  const colPx = axis.devicePx;
  const rowPx = lay.rowScale * dpr;
  const vis = axis.visible([]);
  if (!vis.length) return;
  const c0 = Math.max(0, vis[0]);
  const c1 = Math.min(lay.columns, vis[vis.length - 1]);
  lay.viewC0 = c0;
  lay.viewC1 = c1;
  const missing = tokens.colors['seq-missing'];
  const held = lay.held;
  const C = lay.columns;
  if ((colPx < config.density.cells || rowPx < config.density.cells) && held) layoutImages(lay, ctx, c0, c1);
  else {
    const lettersStep = axis.scale >= config.density.letters && lay.rowScale >= config.density.letters;
    lay.letters = lettersStep;
    if (lettersStep) {
      let size = 0;
      for (const s of config.letterSizes)
        if (s * config.monoAdvance <= axis.scale - 2 && s <= lay.rowScale * 0.75) size = s;
      lay.fontSize = size || (config.letterSizes[0] ?? 9);
      lay.font = `${lay.fontSize}px mono`;
    }
    for (let i = lay.first; i < lay.last; i++) addRow(lay, i, c0, c1, held, C, missing, tokens, prep);
  }
  lay.pinnedFrom = lay.n;
  for (let i = 0; i < lay.display.pinned; i++) addRow(lay, i, c0, c1, held, C, missing, tokens, prep);
}

function addRow(
  lay: RowsLay,
  i: number,
  c0: number,
  c1: number,
  held: HeldCells | undefined,
  C: number,
  missing: number,
  tokens: SeqTokens,
  prep: Prep,
) {
  const axis = lay.axis;
  const dpr = axis.dpr;
  const e = lay.display.rowAt(i);
  const y = lay.rowY(i);
  const h = lay.rowScale;
  let codes: Uint8Array | undefined;
  let base = 0;
  let summary: Uint8Array | undefined;
  if (e < 0) summary = summaryOf(lay, i);
  else if (held) {
    codes = held.codes;
    base = e * C;
  } else codes = lay.lazyCodes.get(e);
  if (!codes && !summary) return;
  let run = -1;
  for (let j = c0; j < c1; j++) {
    let code: number;
    let alpha = 255;
    if (summary) {
      code = summary[3 * j];
      if (prep.summary === 'bar') {
        alpha = summary[3 * j + 2];
        code = alpha ? 65 : 45;
      } else alpha = summary[3 * j + 1];
    } else code = codes![base + j];
    let fill = lay.colorOf(code, j, missing);
    if (fill === NoColor) {
      run = -1;
      continue;
    }
    if (alpha < 255) fill = ((fill & 0xffffff00) | Math.round(((fill & 0xff) * alpha) / 255)) >>> 0;
    if (!axis.cell(j)) continue;
    const x0 = Math.round(axis.sx0 * dpr) / dpr;
    const x1 = Math.round(axis.sx1 * dpr) / dpr;
    // a run of one fill on a row at the cells step is one rect (bounded by the screen)
    if (!lay.letters && run >= 0 && lay.fill[run] === fill && lay.x1[run] >= x0 - 0.01) {
      lay.x1[run] = x1;
      continue;
    }
    const k = lay.n++;
    if (k >= lay.x0.length) {
      lay.x0 = grown(lay.x0, k + 1);
      lay.x1 = grown(lay.x1, k + 1);
      lay.y = grown(lay.y, k + 1);
      lay.h = grown(lay.h, k + 1);
      lay.fill = grown(lay.fill, k + 1);
      lay.code = grown(lay.code, k + 1);
      lay.ink = grown(lay.ink, k + 1);
    }
    lay.x0[k] = x0;
    lay.x1[k] = x1;
    lay.y[k] = y;
    lay.h[k] = h;
    lay.fill[k] = fill;
    lay.code[k] = summary && prep.summary === 'bar' ? 0 : code;
    lay.ink[k] = lay.letters ? inkOn(fill, tokens) : NoColor;
    run = k;
  }
}

/** A folded clade's consensus at display row `i`, computed exactly from its cells when first drawn and kept. */
function summaryOf(lay: RowsLay, i: number): Uint8Array | undefined {
  const held = lay.held;
  const members = lay.display.members.get(i);
  if (!held || !members) return undefined;
  if (lay.display.version !== lay.summariesOf) {
    lay.summaries.clear();
    lay.summariesOf = lay.display.version;
  }
  let s = lay.summaries.get(i);
  if (!s) lay.summaries.set(i, (s = consensus(held, members, new Uint8Array(3 * held.columns))));
  return s;
}

/**
 * Below 3 device px a cell: the shared reduction over the scrolled display rows' theme classes (a class a distinct
 * fill), at the level pair whose sample is at least a device pixel, drawn as image tiles (256 × 256 samples) kept in
 * the shared image cache by (layer, level, tile, base, theme), offset from display row `pinned` (the pinned rows are
 * drawn on their own).
 */
function layoutImages(lay: RowsLay, ctx: LayoutContext, c0: number, c1: number) {
  const { axis, tokens } = ctx;
  const held = lay.held!;
  const display = lay.display;
  const C = held.columns;
  const P = display.pinned;
  const pyr = (lay.pyramid ??= new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' }));
  const baseKey = `${held.version}|${display.version}|${tokens.version}|${lay.channelVersion}`;
  const R = display.count - P;
  if (lay.pyramidOf !== baseKey) {
    lay.pyramidOf = baseKey;
    const need = R * C;
    if (!lay.baseBytes || lay.baseBytes.length < need) lay.baseBytes = new Uint8Array(need);
    const base = lay.baseBytes;
    // one class a distinct fill (at most 254): the theme's classes, so a reduced sample reads as the theme; a table by
    // (column, letter) made once, then one lookup a cell
    const classOf = new Map<number, number>();
    const missing = tokens.colors['seq-missing'];
    const table = lay.classTable && lay.classTable.length >= C * 128 ? lay.classTable : new Uint8Array(C * 128);
    lay.classTable = table;
    table.fill(254, 0, C * 128);
    const classAt = (code: number, j: number): number => {
      const u = upper(code);
      if (u > 127) return NoClass;
      const k = j * 128 + u;
      let c = table[k];
      if (c !== 254) return c;
      const fill = lay.colorOf(code, j, missing);
      if (fill === NoColor) c = NoClass;
      else {
        let id = classOf.get(fill);
        if (id === undefined) {
          id = classOf.size < 253 ? classOf.size : 252;
          if (classOf.size < 253) classOf.set(fill, id);
          lay.classColors[id] = fill;
        }
        c = id;
      }
      table[k] = c;
      return c;
    };
    for (let i = 0; i < R; i++) {
      const e = display.rowAt(P + i);
      const sum = e < 0 ? summaryOf(lay, P + i) : undefined;
      if (sum) {
        for (let j = 0; j < C; j++) base[i * C + j] = classAt(sum[3 * j], j);
        continue;
      }
      const row = e * C;
      for (let j = 0; j < C; j++) {
        const code = held.codes[row + j];
        const c = table[j * 128 + (code & 127)];
        base[i * C + j] = c !== 254 && code < 97 ? c : classAt(code, j);
      }
    }
    pyr.setClasses(R, C, base.subarray(0, need));
  }
  if (R <= 0) return;
  const dpr = axis.dpr;
  const rowLevel = pyramidLevelFor(1 / (lay.rowScale * dpr), 4);
  const colLevel = pyramidLevelFor(1 / axis.devicePx, 4);
  const level = rowsLevel(pyr, rowLevel, colLevel, ctx.settled);
  const fr = 4 ** level.rowLevel;
  const fc = 4 ** level.columnLevel;
  const s0 = Math.floor((lay.first - P) / fr);
  const s1 = Math.ceil((lay.last - P) / fr);
  const t0 = Math.floor(c0 / fc);
  const t1 = Math.ceil(c1 / fc);
  const cache = lay.cache!;
  const whole = level.shares.length === 0;
  for (let ti = Math.floor(s0 / ImageTile); ti * ImageTile < s1; ti++)
    for (let tj = Math.floor(t0 / ImageTile); tj * ImageTile < t1; tj++) {
      const key = `${lay.key}|${level.rowLevel}.${level.columnLevel}|${ti}.${tj}|${lay.pyramidOf}`;
      let img = cache.get(key);
      const rows = Math.min(ImageTile, level.rows - ti * ImageTile);
      const cols = Math.min(ImageTile, level.columns - tj * ImageTile);
      if (rows <= 0 || cols <= 0) continue;
      if (!img) {
        img = cache.slot(cols, rows);
        const d = img.data;
        for (let r = 0; r < rows; r++)
          for (let c = 0; c < cols; c++) {
            const o = (ti * ImageTile + r) * level.columns + tj * ImageTile + c;
            const cl = level.classes[o];
            const q = 4 * (r * cols + c);
            if (cl === NoClass) {
              d[q] = d[q + 1] = d[q + 2] = d[q + 3] = 0;
              continue;
            }
            const fill = lay.classColors[cl];
            // opaque by the majority's share: a mixed sample reads lighter (the base level holds no shares: whole)
            const share = whole ? 255 : level.shares[o];
            const a = ((fill & 0xff) * (0.35 + (0.65 * share) / 255)) | 0;
            d[q] = fill >>> 24;
            d[q + 1] = (fill >>> 16) & 0xff;
            d[q + 2] = (fill >>> 8) & 0xff;
            d[q + 3] = a;
          }
        cache.put(key, img);
      }
      const rowStart = ti * ImageTile * fr;
      const colStart = tj * ImageTile * fc;
      const x0 = axis.x(Math.min(colStart, lay.columns));
      const x1 = axis.x(Math.min(colStart + cols * fc, lay.columns));
      const y0 = lay.rowY(P + rowStart);
      const y1 = lay.rowY(P + Math.min(rowStart + rows * fr, R));
      lay.images.push(img);
      lay.ib.push(Math.min(x0, x1), y0, Math.abs(x1 - x0), y1 - y0);
    }
}

/**
 * The header in the label column: the tree (its own width, the label column widened by it), the labels (cut with an
 * ellipsis past their room), the first residue numbers before the data and the last after it (a right gutter).
 */
function layoutHeader(lay: RowsLay, prep: Prep, ctx: LayoutContext) {
  const { config, box, text } = ctx;
  const rs = lay.rowScale;
  const left = config.inset;
  const right = box.x - config.labelGap;
  const treeW =
    lay.clades && prep.order === 'tree' && !lay.figure ? Math.max(0, Math.min(prep.treeWidth, right - left - 24)) : 0;
  lay.treeX0 = left;
  lay.treeX1 = left + treeW;
  if (treeW > 0) layoutTree(lay, lay.clades!, left, treeW);
  if (rs < config.rows.labelMinPx) return;
  // a row's label is data, as an item's: the items' font (a Pfam label, NAME/1-147, reads whole in the column)
  lay.labelFont = config.itemFont;
  lay.labelX = left + treeW + (treeW ? 4 : 0);
  const numbers = prep.numbers === 'gutters';
  const numberW = numbers ? text.measure(config.itemFont, '00000') : 0;
  lay.numberX = right;
  lay.numberRightX = box.x + box.width + config.labelGap;
  // the left row columns between the label and the first numbers (the label column was widened by them)
  const leftRoom = columnsRoom(prep.columns, lay.cols.width, 'left', config.rows.columns.gapPx);
  lay.labelMax = Math.max(0, right - numberW - 4 - lay.labelX - leftRoom);
  const size = fontSize(config.itemFont);
  const add = (i: number) => {
    const e = lay.display.rowAt(i);
    const y = lay.rowY(i) + rs / 2 + size * 0.35;
    let t: string;
    if (e < 0) {
      // a summary's count is what it says (EC-100): the first label is cut, never ` +n`
      const m = lay.display.members.get(i);
      const more = lay.summaryMore(i);
      const head = m
        ? elide(text, lay.labelFont, lay.label(m[0]), lay.labelMax - text.measure(lay.labelFont, more))
        : '';
      t = m ? (head ? `${head}${more}` : elide(text, lay.labelFont, more.trim(), lay.labelMax)) : '';
    } else {
      t = lay.label(e);
      if (lay.display.pinned && i >= lay.display.pinned && config.rows.pinInPlace && isPinned(lay, e)) t = `· ${t}`;
      t = elide(text, lay.labelFont, t, lay.labelMax);
    }
    lay.labelText.push(t);
    lay.labelY.push(y);
    lay.labelW.push(text.measure(lay.labelFont, t));
    if (!numbers || e < 0 || !endNumbers(lay, e)) return;
    const n = String(lay.numFirst);
    lay.numberText.push(n);
    lay.numberY.push(y);
    lay.numberW.push(text.measure(config.itemFont, n));
    // the last residue number in the right gutter (a figure has none: its lines repeat the rows)
    if (!lay.figure) {
      const s = String(lay.numLast);
      const w = text.measure(config.itemFont, s);
      lay.lastText.push(s);
      lay.lastY.push(y);
      lay.lastW.push(w);
      if (w > lay.lastMaxW) lay.lastMaxW = w;
    }
  };
  for (let i = 0; i < lay.display.pinned; i++) add(i);
  for (let i = lay.first; i < lay.last; i++) add(i);
}

function isPinned(lay: RowsLay, row: number): boolean {
  for (let i = 0; i < lay.display.pinned; i++) if (lay.display.rowAt(i) === row) return true;
  return false;
}

/** A wedge's least length, px: a clade of zero branch lengths still shows one. */
const WedgeMinPx = 6;

/**
 * The least height a clade's rows take for its branch to be drawn, px (_spec/README.md "Render", The rows layer, The
 * tree's look): below it a clade's branches would sit under 2 px apart, and 10,000 leaves in a band read as a solid
 * slab; the tree is drawn down to the clades that resolve.
 */
const TreeMinPx = 4;

/**
 * The tree beside the rows from its kept geometry: only what crosses the band is kept (segments, each with the node a
 * hit on it toggles: an inner node's own vertical and the branch to it; a leaf's branch none), the wedges of folded
 * clades and the inner nodes a hit reads.
 */
function layoutTree(lay: RowsLay, clades: Clades, x0: number, w: number) {
  treeGeometry(lay, clades);
  const n = clades.nodes;
  const { coord, nodeX, inFold, folded, box } = lay;
  const rs = lay.rowScale;
  const P = lay.display.pinned;
  // a display row's y: the scrolled band's (the tree's rows are in the order, never the pinned band's copies)
  const y0 = box.y + lay.pinnedHeight - P * rs - lay.bandTop;
  const top = box.y + lay.pinnedHeight;
  const bottom = box.y + box.height;
  const inner = (k: number) => clades.hi[k] - clades.lo[k] > 1;
  // a clade drawn only when its rows take `TreeMinPx` (a node's segments, the branch to it, a folded clade's wedge)
  const least = TreeMinPx / rs;
  const span = lay.spanRows;
  for (let k = 0; k < n; k++) {
    if (inFold[k] || Number.isNaN(coord[k]) || span[k] < least) continue;
    const x = x0 + nodeX[k] * w;
    const yk = y0 + coord[k] * rs;
    if (folded[k]) {
      if (yk + rs / 2 >= top && yk - rs / 2 <= bottom) {
        // to its deepest leaf's x (a few px at least, so a clade of zero lengths still shows one); the hit's box the
        // row's middle 0.8, the drawn wedge's far side half the row (`drawWedges`)
        const far = Math.min(x0 + w, Math.max(x0 + lay.deepX[k] * w, x + WedgeMinPx));
        // a node within `WedgeMinPx` of the tree's edge: the wedge starts that far before it (its hit box too)
        lay.wedges.push(Math.min(x, far - WedgeMinPx), yk - rs * 0.4, far, yk + rs * 0.4, k);
        lay.nodes.push(k, x, yk);
      }
      continue;
    }
    const kids = clades.children[k];
    // the vertical spans the node and its drawn children only (none drawn: none; zoomed out, a child too small to
    // resolve has no branch to meet)
    let lo = yk;
    let hi = yk;
    let drawn = false;
    for (const c of kids) {
      if (Number.isNaN(coord[c]) || span[c] < least) continue;
      const y = y0 + coord[c] * rs;
      drawn = true;
      if (y < lo) lo = y;
      if (y > hi) hi = y;
      // the child's branch: from this node's x to the child's, at the child's y
      if (y >= top - 1 && y <= bottom + 1) {
        lay.tree.push(x, y, x0 + nodeX[c] * w, y);
        lay.treeNode.push(inner(c) ? c : -1);
        lay.treeClade.push(c);
      }
    }
    if (drawn && hi >= top && lo <= bottom) {
      lay.tree.push(x, Math.max(lo, top), x, Math.min(hi, bottom));
      lay.treeNode.push(inner(k) ? k : -1);
      lay.treeClade.push(k);
      if (yk >= top && yk <= bottom) lay.nodes.push(k, x, yk);
    }
  }
}

/** A clade's first shown row (its summary row's), -1 when none is shown. */
function firstShown(lay: RowsLay, clades: Clades, k: number): number {
  for (let l = clades.lo[k]; l < clades.hi[k]; l++) {
    const r = clades.leafRow[l];
    if (r >= 0 && lay.display.placeRow(r) >= 0) return r;
  }
  return -1;
}

/**
 * The residue numbers at the ends of a data row's part in view (_spec/README.md "Render", The rows layer: the header):
 * the numbers of its first and last residue in the columns in view, into `numFirst` and `numLast`; a gap is never
 * counted; false for a row with no residue in view. At the fit view they are the row's own (kept: counted once a data
 * or a source read), a pan or a zoom counts the residues before and in the columns in view.
 */
function endNumbers(lay: RowsLay, row: number): boolean {
  const held = lay.held;
  let codes: Uint8Array | undefined;
  let base = 0;
  let first: number;
  let C: number;
  if (held) {
    codes = held.codes;
    C = held.columns;
    base = row * C;
    first = held.first[row];
  } else {
    codes = lay.lazyCodes.get(row);
    if (!codes) return false;
    C = codes.length;
    first = lay.lazyFirst.get(row) ?? 1;
  }
  const c0 = lay.figure ? 0 : lay.viewC0;
  const c1 = lay.figure ? C : Math.min(lay.viewC1, C);
  if (c1 <= c0) return false;
  if (c0 === 0 && c1 === C) {
    const last = ownLast(lay, row, codes, base, C, first);
    if (last === NoNumber) return false;
    lay.numFirst = first;
    lay.numLast = last;
    return true;
  }
  let before = 0;
  for (let j = 0; j < c0; j++) if (!isGap(codes[base + j])) before++;
  let inView = 0;
  for (let j = c0; j < c1; j++) if (!isGap(codes[base + j])) inView++;
  if (!inView) return false;
  lay.numFirst = first + before;
  lay.numLast = first + before + inView - 1;
  return true;
}

/** A row's own last residue number: its first plus its residues less one (`NoNumber` for a row of gaps), kept. */
function ownLast(lay: RowsLay, row: number, codes: Uint8Array, base: number, C: number, first: number): number {
  const held = lay.held;
  if (held) {
    if (lay.lastOf !== held) {
      lay.lastOf = held;
      const out = new Int32Array(held.rows);
      for (let r = 0; r < held.rows; r++) {
        let n = 0;
        for (let j = 0; j < C; j++) if (!isGap(held.codes[r * C + j])) n++;
        out[r] = n ? held.first[r] + n - 1 : NoNumber;
      }
      lay.lasts = out;
    }
    return lay.lasts[row];
  }
  // a lazy row's: counted once when its tile is read, kept (bounded: let go past 4,096 rows)
  const kept = lay.lazyLasts.get(row);
  if (kept !== undefined) return kept;
  let n = 0;
  for (let j = 0; j < C; j++) if (!isGap(codes[base + j])) n++;
  const last = n ? first + n - 1 : NoNumber;
  if (lay.lazyLasts.size >= 4096) lay.lazyLasts.clear();
  lay.lazyLasts.set(row, last);
  return last;
}

/**
 * The tree's geometry in display rows, kept while the display is (_spec/README.md "Render", The rows layer: nothing
 * sized by the nodes made a pan): each node's x as a share of the tree's width (branch lengths; depth when none), its
 * y as a display row (a leaf at its row's place in the order, a pinned one in place; a folded clade at its summary
 * row; an inner node between its first and last child), whether it is folded or inside a fold.
 */
function treeGeometry(lay: RowsLay, clades: Clades) {
  const d = lay.display;
  if (lay.treeOf === d.version && lay.treeClades === clades) return;
  lay.treeOf = d.version;
  lay.treeClades = clades;
  const n = clades.nodes;
  if (lay.coord.length < n) {
    lay.coord = new Float64Array(n);
    lay.nodeX = new Float64Array(n);
    lay.deepX = new Float64Array(n);
    lay.spanRows = new Float64Array(n);
    lay.inFold = new Uint8Array(n);
    lay.folded = new Uint8Array(n);
  }
  const { coord, nodeX, deepX, inFold, folded } = lay;
  let maxD = 0;
  for (let k = 0; k < n; k++) if (clades.distance[k] > maxD) maxD = clades.distance[k];
  let maxDepth = 0;
  for (let k = 0; k < n; k++) if (clades.depth[k] > maxDepth) maxDepth = clades.depth[k];
  for (let k = 0; k < n; k++) nodeX[k] = maxD > 0 ? clades.distance[k] / maxD : clades.depth[k] / Math.max(1, maxDepth);
  // a node's children come after it: its deepest leaf's x from theirs
  for (let k = n - 1; k >= 0; k--) {
    const kids = clades.children[k];
    let d = nodeX[k];
    for (const c of kids) if (deepX[c] > d) d = deepX[c];
    deepX[k] = d;
  }
  // folded: a clade shown as a summary row (one folded to its one shown row is that row)
  folded.fill(0, 0, n);
  for (const i of d.members.keys()) {
    const e = d.rowAt(i);
    if (e < 0) folded[-1 - e] = 1;
  }
  for (let k = 0; k < n; k++) {
    const p = clades.parent[k];
    inFold[k] = p >= 0 && (inFold[p] || folded[p]) ? 1 : 0;
  }
  const { spanRows } = lay;
  for (let k = n - 1; k >= 0; k--) {
    coord[k] = Number.NaN;
    spanRows[k] = 0;
    if (inFold[k]) continue;
    const kids = clades.children[k];
    if (!kids.length || folded[k]) {
      const leafRow = kids.length ? -1 : clades.leafRow[clades.lo[k]];
      const row = leafRow >= 0 ? leafRow : firstShown(lay, clades, k);
      const at = row >= 0 ? d.placeRow(row) : -1;
      if (at >= 0) {
        coord[k] = at + 0.5;
        spanRows[k] = 1;
      }
      continue;
    }
    let lo = Number.NaN;
    let hi = Number.NaN;
    for (const c of kids) {
      const y = coord[c];
      if (Number.isNaN(y)) continue;
      if (Number.isNaN(lo) || y < lo) lo = y;
      if (Number.isNaN(hi) || y > hi) hi = y;
      spanRows[k] += spanRows[c];
    }
    if (!Number.isNaN(lo)) coord[k] = (lo + hi) / 2;
  }
}

/**
 * The hit written into the layer's scratch (a hover allocates nothing): every field set. A summary row's or a tree
 * node's element is its clade's, its extent the clade's rows (kept by the display).
 */
function writeHit(
  lay: RowsLay,
  system: string,
  position: number | undefined,
  yPosition: number | undefined,
  target: HitTarget,
  clade?: { element: ElementKey; extent: Location },
): Hit {
  const h = lay.scratch;
  h.layer = lay.key;
  h.section = lay.section;
  h.system = system;
  h.position = position;
  h.ySystem = yPosition === undefined ? undefined : lay.rowsSystem;
  h.yPosition = yPosition;
  h.target = target;
  h.element = clade?.element;
  h.extent = clade?.extent;
  h.state = undefined;
  return h;
}

/** A header's hit naming the column cell it is on (`title: value`, made at layout: a hover allocates nothing). */
function named(lay: RowsLay, h: Hit, q: number): Hit {
  if (q >= 0 && h === lay.scratch) {
    const c = lay.cols;
    lay.scratch.state = c.name[c.col[q]][c.row[q]];
  }
  return h;
}

/**
 * A clade's element and its rows as a Location on the rows system (its shown rows when folded at display row `i`),
 * made once a display.
 */
function cladeHit(lay: RowsLay, c: number, i?: number): { element: ElementKey; extent: Location } {
  if (lay.cladeHitsOf !== lay.display.version) {
    lay.cladeHitsOf = lay.display.version;
    lay.cladeHits.clear();
    lay.summaryHits.clear();
  }
  const kept = i === undefined ? lay.cladeHits.get(c) : lay.summaryHits.get(i);
  if (kept) return kept;
  const members = i === undefined ? undefined : lay.display.members.get(i);
  const ranges: number[] = [];
  if (members) for (const r of members) ranges.push(r, r + 1);
  else if (lay.clades)
    for (let k = lay.clades.lo[c]; k < lay.clades.hi[c]; k++) {
      const r = lay.clades.leafRow[k];
      if (r >= 0) ranges.push(r, r + 1);
    }
  const made = {
    element: { owner: lay.owner, key: lay.clades?.name[c] ?? '' },
    extent: rangesLocation(lay.rowsSystem, ranges),
  };
  if (i === undefined) lay.cladeHits.set(c, made);
  else lay.summaryHits.set(i, made);
  return made;
}

// ---- drawing (nothing allocated once warm: numbers through the writer's `args`) ----

/** Clips to a box from its top `y` down to the band's foot (the scrolled rows below the pinned band). */
function clipBelow(lay: RowsLay, out: PassWriter, y: number) {
  const g = out.args;
  g[0] = lay.box.x;
  g[1] = y;
  g[2] = lay.box.width;
  g[3] = Math.max(0, lay.box.y + lay.box.height - y);
  out.clipArgs();
}

function drawData(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  // the scrolled rows below the pinned band, clipped to it: the image tiles, the rects, their letters
  clipBelow(lay, out, lay.box.y + lay.pinnedHeight);
  if (lay.neutral) drawNeutral(lay, out);
  for (let k = 0; k < lay.images.length; k++) {
    g[0] = lay.ib[4 * k];
    g[1] = lay.ib[4 * k + 1];
    g[2] = lay.ib[4 * k + 2];
    g[3] = lay.ib[4 * k + 3];
    out.imageArgs(lay.images[k]);
  }
  drawRects(lay, out, 0, lay.pinnedFrom);
  if (lay.pending.length) drawPending(lay, out);
  if (lay.letters) drawLetters(lay, out, 0, lay.pinnedFrom);
  out.unclip();
  if (!lay.display.pinned) return;
  // the pinned rows on their own, over a band that does not scroll
  clipBelow(lay, out, lay.box.y);
  drawRects(lay, out, lay.pinnedFrom, lay.n);
  if (lay.letters) drawLetters(lay, out, lay.pinnedFrom, lay.n);
  drawPinnedRule(lay, out);
  out.unclip();
}

function drawRects(lay: RowsLay, out: PassWriter, from: number, to: number) {
  const g = out.args;
  // a letter cell's hairline (`seq-cell-edge`: dark only, none in light), as the letters layer's
  const edge = lay.letters ? lay.tokens.colors['seq-cell-edge'] : NoColor;
  const px = edge === NoColor ? 0 : 1 / lay.axis.dpr;
  for (let k = from; k < to; k++) {
    g[0] = lay.x0[k];
    g[1] = lay.y[k];
    g[2] = lay.x1[k] - lay.x0[k];
    g[3] = lay.h[k];
    g[4] = lay.fill[k];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
    if (edge !== NoColor) {
      g[5] = edge;
      g[6] = px;
      cellEdge(out);
    }
  }
}

function drawLetters(lay: RowsLay, out: PassWriter, from: number, to: number) {
  const g = out.args;
  let y = Number.NaN;
  let open = false;
  for (let k = from; k < to; k++) {
    const code = lay.code[k];
    if (!code) continue;
    const base = lay.y[k] + lay.h[k] / 2 + lay.fontSize * 0.35;
    if (base !== y) {
      if (open) out.endGlyphs();
      g[1] = base;
      out.beginGlyphsArgs(lay.font, NoColor);
      open = true;
      y = base;
    }
    g[0] = (lay.x0[k] + lay.x1[k]) / 2;
    g[4] = lay.ink[k];
    out.letterArgs(code);
  }
  if (open) out.endGlyphs();
}

function drawPending(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  for (let k = 0; k + 1 < lay.pending.length; k += 2) {
    g[0] = lay.box.x;
    g[1] = lay.pending[k];
    g[2] = lay.box.width;
    g[3] = lay.pending[k + 1];
    g[4] = lay.tokens.colors['seq-pending'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

/** Lazy rows below a pixel: the band a neutral fill (never blank, never a pending stripe that never resolves). */
function drawNeutral(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  g[0] = lay.box.x;
  g[1] = lay.box.y + lay.pinnedHeight;
  g[2] = lay.box.width;
  g[3] = Math.max(0, lay.box.height - lay.pinnedHeight);
  g[4] = lay.tokens.colors['seq-missing'];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** A hairline under the pinned band. */
function drawPinnedRule(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  g[0] = lay.box.x;
  g[1] = lay.box.y + lay.pinnedHeight - 0.5;
  g[2] = lay.box.width;
  g[3] = 1;
  g[4] = lay.tokens.colors['msui-line'];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** The header: the tree's branches and wedges, the labels, the first numbers, the chip; the last numbers after the data. */
function drawHeader(lay: RowsLay, out: PassWriter) {
  const { box, config, tokens } = lay;
  const g = out.args;
  g[0] = config.inset;
  g[1] = box.y;
  g[2] = Math.max(0, box.x - config.inset);
  g[3] = box.height;
  out.clipArgs();
  const ink = tokens.colors['msui-text-muted'];
  const t = lay.tree;
  branchWidth(lay, 1);
  for (let k = 0; k + 3 < t.length; k += 4) branchRect(lay, out, k, ink);
  if (lay.wedges.length) drawWedges(lay, out, ink, -1, false);
  const text = tokens.colors['msui-text'];
  for (let k = 0; k < lay.labelText.length; k++) {
    g[0] = lay.labelX;
    g[1] = lay.labelY[k];
    g[2] = lay.labelW[k];
    g[3] = lay.labelMax;
    g[4] = text;
    out.textArgs(lay.labelText[k], AnchorStart, lay.labelFont);
  }
  for (let k = 0; k < lay.numberText.length; k++) {
    g[0] = lay.numberX;
    g[1] = lay.numberY[k];
    g[2] = lay.numberW[k];
    g[3] = lay.numberW[k];
    g[4] = ink;
    out.textArgs(lay.numberText[k], AnchorEnd, config.itemFont);
  }
  if (lay.chip) {
    g[0] = lay.chipX;
    g[1] = lay.chipY;
    g[2] = lay.chipW;
    g[3] = lay.chipW;
    g[4] = tokens.colors['msui-accent'];
    out.textArgs(lay.chip, AnchorStart, config.itemFont);
  }
  if (lay.cols.n) drawColumnCells(lay, out, 0);
  out.unclip();
  if (lay.cols.n) {
    // the right columns, clipped to the band as the rows are
    g[0] = box.x + box.width;
    g[1] = box.y;
    g[2] = Math.max(0, lay.cols.rightEnd - g[0]);
    g[3] = box.height;
    out.clipArgs();
    drawColumnCells(lay, out, 1);
    out.unclip();
  }
  if (lay.cols.titleN || lay.cols.chip) drawTitles(lay, out);
  // the last numbers in the right gutter, clipped to the band as its rows are (a partly shown row ends at the foot)
  if (lay.lastText.length) {
    g[0] = box.x + box.width;
    g[1] = box.y;
    g[2] = lay.numberRightX - g[0] + lay.lastMaxW;
    g[3] = box.height;
    out.clipArgs();
  }
  for (let k = 0; k < lay.lastText.length; k++) {
    g[0] = lay.numberRightX;
    g[1] = lay.lastY[k];
    g[2] = lay.lastW[k];
    g[3] = lay.lastW[k];
    g[4] = ink;
    out.textArgs(lay.lastText[k], AnchorStart, config.itemFont);
  }
  if (lay.lastText.length) out.unclip();
}

/**
 * A row column's cell cut short, read whole in the `tip` pass (_spec/README.md "Render", The rows layer, Emphasis along
 * y): the one cell under the pointer (a header hit naming it), never from a hover on the row's data or label; the
 * keyboard cursor on a row's header reads all of its row's. A function of its own (measured: no boxing).
 */
function hoverColumns(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  if (!lay.cols.n) return;
  const h = emphasis?.hover;
  if (h && h.layer === lay.key && h.target === 'row' && h.state !== undefined && h.position !== undefined) {
    // the copy the pointer is on (a pinned row kept in place is shown twice): the last hit's, when it shows that row
    const r = lay.hitRow;
    const i =
      r >= 0 && r < lay.display.count && lay.display.rowAt(r) === h.position ? r : lay.display.displayRow(h.position);
    if (i >= 0) hoverColumnCells(lay, out, i, h.state);
  }
  const c = emphasis?.cursor;
  if (!c || c.part !== 'header' || c.ySystem !== lay.rowsSystem || c.yPosition === undefined) return;
  // the layer's address is its section's and more: `startsWith` and a code read allocate nothing
  const s = c.section;
  if (s !== undefined && !(lay.key.startsWith(s) && lay.key.charCodeAt(s.length) === 47)) return;
  // on every copy the cursor is drawn on (as `drawCursor` outlines both of a pinned row's)
  const i = lay.cursorRow(c.yPosition);
  if (i >= 0) hoverColumnCells(lay, out, i, undefined);
  const own = lay.display.placeRow(c.yPosition);
  if (own >= 0 && own !== i) hoverColumnCells(lay, out, own, undefined);
}

/** Sets a branch's width `t` (px) in whole device pixels into `BranchW[0]`, one at least (_spec/README.md "Render", The rows layer, The tree's look). */
function branchWidth(lay: RowsLay, t: number) {
  const dpr = lay.axis.dpr;
  BranchW[0] = Math.max(1, Math.round(t * dpr)) / dpr;
}

/**
 * The branches' width `branchWidth` set, read by `branchRect` and `drawWedges`; a wedge's corner radius (`softCorners`):
 * fractional numbers through an array (one passed to a call that is not inlined is boxed: the hover pass allocates
 * nothing).
 */
const BranchW = new Float64Array(2);

/**
 * Tree segment `k` (of `lay.tree`) as a rect of thickness `t` in `ink`, its edges on device pixels: across the segment
 * its centre's pixel run of the branch's width, along it from each end less half the width to the end plus half (so a
 * branch and its node's vertical meet square, and every branch draws one width).
 */
function branchRect(lay: RowsLay, out: PassWriter, k: number, ink: number) {
  const tr = lay.tree;
  const g = out.args;
  const dpr = lay.axis.dpr;
  const w = BranchW[0];
  const vertical = tr[k] === tr[k + 2];
  const across = vertical ? tr[k] : tr[k + 1];
  const a = vertical ? Math.min(tr[k + 1], tr[k + 3]) : Math.min(tr[k], tr[k + 2]);
  const b = vertical ? Math.max(tr[k + 1], tr[k + 3]) : Math.max(tr[k], tr[k + 2]);
  const c0 = Math.round((across - w / 2) * dpr) / dpr;
  const a0 = Math.round((a - w / 2) * dpr) / dpr;
  const b1 = Math.round((b - w / 2) * dpr) / dpr + w;
  g[0] = vertical ? c0 : a0;
  g[1] = vertical ? a0 : c0;
  g[2] = vertical ? w : b1 - a0;
  g[3] = vertical ? b1 - a0 : w;
  g[4] = ink;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** A wedge's corners and its outline with each corner softened (3 corners, 4 points each: two cuts, two between). */
const WedgeCorners = new Float64Array(6);
const Wedge = new Float64Array(24);

/**
 * The folded clades' wedges (_spec/README.md "Render", The rows layer, The tree's look): a triangle from the clade's
 * node to its deepest leaf's x, its far side half the row high, its corners softened by the shape rule, filled with a
 * light tint of `ink` (a step darker `hovered`) and outlined in `ink` at the branches' width `t`; all, or (`clade` ≥ 0)
 * those inside that clade.
 */
function drawWedges(lay: RowsLay, out: PassWriter, ink: number, clade: number, hovered: boolean) {
  const w = lay.wedges;
  const g = out.args;
  const cl = lay.clades;
  const { config } = lay;
  const half = lay.rowScale / 4;
  const stroke = BranchW[0];
  // a tint of `ink` at a share of its own opacity
  const fill = ((ink & 0xffffff00) | Math.round((ink & 0xff) * (hovered ? 0.36 : 0.18))) >>> 0;
  for (let k = 0; k + 4 < w.length; k += 5) {
    if (clade >= 0 && cl && !within(cl, w[k + 4], clade)) continue;
    const x0 = w[k];
    const yc = (w[k + 1] + w[k + 3]) / 2;
    const x1 = w[k + 2];
    WedgeCorners[0] = x0;
    WedgeCorners[1] = yc;
    WedgeCorners[2] = x1;
    WedgeCorners[3] = yc - half;
    WedgeCorners[4] = x1;
    WedgeCorners[5] = yc + half;
    // the squircle rule's radius, inline (`squircleRadius` given fractional numbers boxes them)
    const len = x1 - x0;
    const r =
      len < config.squircleMin || 2 * half < config.squircleMin
        ? 0
        : Math.min(2 * half * config.squircleRatio, config.squircleMax, len / 2, half);
    BranchW[1] = r;
    const n = softCorners();
    g[4] = fill;
    g[5] = ink;
    g[6] = stroke;
    g[7] = 0;
    out.pathArgs(Wedge, n, true);
  }
}

/**
 * The triangle in `WedgeCorners` into `Wedge`, each corner cut `r` along its two sides (at most a third of the shorter)
 * and joined by a quadratic curve through the corner, sampled: its point count (3, the corners, when `r` is 0).
 */
function softCorners(): number {
  const r = BranchW[1];
  const p = WedgeCorners;
  if (!(r > 0)) {
    for (let i = 0; i < 6; i++) Wedge[i] = p[i];
    return 3;
  }
  let o = 0;
  for (let c = 0; c < 3; c++) {
    const px = p[2 * c];
    const py = p[2 * c + 1];
    const ax = p[2 * ((c + 2) % 3)];
    const ay = p[2 * ((c + 2) % 3) + 1];
    const bx = p[2 * ((c + 1) % 3)];
    const by = p[2 * ((c + 1) % 3) + 1];
    const la = Math.sqrt((ax - px) * (ax - px) + (ay - py) * (ay - py));
    const lb = Math.sqrt((bx - px) * (bx - px) + (by - py) * (by - py));
    const d = Math.min(r, la / 3, lb / 3);
    const sx = px + ((ax - px) * d) / la;
    const sy = py + ((ay - py) * d) / la;
    const ex = px + ((bx - px) * d) / lb;
    const ey = py + ((by - py) * d) / lb;
    for (let s = 0; s < 4; s++) {
      const u = s / 3;
      const v = 1 - u;
      Wedge[o++] = v * v * sx + 2 * u * v * px + u * u * ex;
      Wedge[o++] = v * v * sy + 2 * u * v * py + u * u * ey;
    }
  }
  return 12;
}

/** Node `k` is `c` or below it (its leaves are within `c`'s; a chain of one child shares them: by depth too). */
function within(cl: Clades, k: number, c: number): boolean {
  return cl.lo[k] >= cl.lo[c] && cl.hi[k] <= cl.hi[c] && cl.depth[k] >= cl.depth[c];
}

/**
 * A light band over `rows` display rows from `i` (one by default) as far as it is in the band: the scrolled rows
 * (`i` past the pinned band) never cover the pinned band.
 */
function hoverBand(lay: RowsLay, out: PassWriter, i: number, ink: number, rows = 1) {
  const { box } = lay;
  const y0 = Math.max(lay.rowY(i), i >= lay.display.pinned ? box.y + lay.pinnedHeight : box.y);
  const y1 = Math.min(lay.rowY(i) + rows * lay.rowScale, box.y + box.height);
  if (y1 <= y0) return;
  const g = out.args;
  g[0] = box.x;
  g[1] = y0;
  g[2] = box.width;
  g[3] = y1 - y0;
  g[4] = ((ink & 0xffffff00) | 0x1a) >>> 0;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** The hovered cell's x extents (`hoverCell` reads them: set before the call). */
const CellX = new Float64Array(2);

/**
 * The hovered cell's outline at display row `i`, between `CellX[0]` and `CellX[1]`: a squircle by the cell's size (every emphasis a
 * squircle), its radius plus the offset, square on a side the view or the band's edge cuts and below the minimum.
 */
function hoverCell(lay: RowsLay, out: PassWriter, i: number, ink: number) {
  const { box, config } = lay;
  // the cell's x through an array: a fractional number passed to a call that is not inlined is boxed
  const a = CellX[0];
  const b = CellX[1];
  const y = lay.rowY(i);
  const rh = lay.rowScale;
  if (y + rh <= box.y || y >= box.y + box.height) return;
  const g = out.args;
  const d = config.outline.hover / 2 + config.emphasis.hoverOffset;
  const w = b - a;
  const r =
    w < config.squircleMin || rh < config.squircleMin
      ? 0
      : Math.min(rh * config.squircleRatio, config.squircleMax, w / 2, rh / 2);
  // cut: past the data box's sides (the view), above or below the band (its rows scroll inside it)
  let round = (a >= box.x ? 1 : 0) | (b <= box.x + box.width ? 2 : 0);
  if (y < box.y) round |= 4;
  if (y + rh > box.y + box.height) round |= 8;
  g[0] = a - d;
  g[1] = y - d;
  g[2] = w + 2 * d;
  g[3] = rh + 2 * d;
  g[4] = NoColor;
  g[5] = ink;
  g[6] = config.outline.hover;
  g[7] = r > 0 && (round & 3) !== 0 ? r + d : 0;
  g[8] = round;
  out.rectArgs();
}

/**
 * A hovered clade (_spec/README.md "Render", The rows layer, Hover): its branches in the hover ink, one band over each
 * run of its rows in the order (a folded clade's summary row alone), a pinned row's band too.
 */
function hoverClade(lay: RowsLay, out: PassWriter, key: string, ink: number) {
  const cl = lay.clades;
  const c = cl?.byName.get(key);
  if (!cl || c === undefined) return;
  const g = out.args;
  const { box, config } = lay;
  const P = lay.display.pinned;
  let start = -1;
  let end = -1;
  for (let l = cl.lo[c]; l < cl.hi[c]; l++) {
    const r = cl.leafRow[l];
    const p = r >= 0 ? lay.display.placeRow(r) : -1;
    if (p < 0) continue;
    // a pinned row in the band is lit there (a run is the scrolled rows and its place in the order)
    if (P) {
      const d = lay.display.displayRow(r);
      if (d >= 0 && d < P) hoverBand(lay, out, d, ink);
    }
    if (p < P) continue;
    if (start < 0) start = end = p;
    else if (p === end)
      continue; // a folded clade's rows share its summary row
    else if (p === end + 1) end = p;
    else {
      hoverBand(lay, out, start, ink, end - start + 1);
      start = end = p;
    }
  }
  if (start >= 0) hoverBand(lay, out, start, ink, end - start + 1);
  g[0] = config.inset;
  g[1] = box.y;
  g[2] = Math.max(0, box.x - config.inset);
  g[3] = box.height;
  out.clipArgs();
  const t = lay.tree;
  const owner = lay.treeClade;
  branchWidth(lay, config.outline.hover);
  for (let k = 0, s = 0; k + 3 < t.length; k += 4, s++) if (within(cl, owner[s], c)) branchRect(lay, out, k, ink);
  if (lay.wedges.length) drawWedges(lay, out, ink, c, true);
  out.unclip();
}

/**
 * The hover: its row's band tinted lightly (a pinned row at both its places: the pinned band and the order), its
 * column within the rows, the cell outlined.
 */
function drawHover(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key) return;
  const { box, tokens } = lay;
  const g = out.args;
  const ink = tokens.colors['seq-hover'];
  if ((h.target === 'tree-node' || h.target === 'summary-row') && h.element?.owner === lay.owner)
    hoverClade(lay, out, String(h.element.key), ink);
  const row = h.yPosition ?? (h.system === lay.rowsSystem ? h.position : undefined);
  const i = row === undefined ? -1 : lay.display.displayRow(row);
  const j = row === undefined ? -1 : lay.display.placeRow(row);
  if (i >= 0) hoverBand(lay, out, i, ink);
  if (j >= 0 && j !== i) hoverBand(lay, out, j, ink);
  if (h.position === undefined || h.system !== lay.system || !lay.axis.cell(h.position)) return;
  const dpr = lay.axis.dpr;
  const a = Math.round(lay.axis.sx0 * dpr) / dpr;
  const b = Math.max(Math.round(lay.axis.sx1 * dpr) / dpr, a + 1 / dpr);
  g[0] = a;
  g[1] = box.y;
  g[2] = b - a;
  g[3] = box.height;
  g[4] = ((ink & 0xffffff00) | 0x14) >>> 0;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
  CellX[0] = a;
  CellX[1] = b;
  if (i >= 0) hoverCell(lay, out, i, ink);
  if (j >= 0 && j !== i) hoverCell(lay, out, j, ink);
}

/** The selection's tint and outline colors and the x span (`drawSelection` sets them): read by `paintRows`. */
const SelPaint = new Float64Array(4);

/**
 * The selection along y: a rows part (data rows) as row bands over the data box, tinted with an outline; a rect of
 * rows × columns as its cells' boxes. A summary row draws the partial mark: the tint over the share of its rows
 * selected, the outline when all are. O(rows in view) whatever the selection's size.
 */
function drawSelection(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const sel = emphasis?.selection;
  if (!sel) return;
  const { box, tokens, config } = lay;
  const accent = tokens.colors['msui-accent'];
  // the colors through a typed array: a number past a small integer passed to a call is boxed
  SelPaint[0] = ((accent & 0xffffff00) | Math.round(255 * config.emphasis.selectionTint)) >>> 0;
  SelPaint[1] = ((accent & 0xffffff00) | Math.round(255 * config.emphasis.selectionOutline)) >>> 0;
  for (let p = 0; p < sel.parts.length; p++) {
    const part = sel.parts[p];
    if (part.system !== lay.rowsSystem) continue;
    SelPaint[2] = box.x;
    SelPaint[3] = box.x + box.width;
    paintRows(lay, out, part.ranges);
  }
  for (let r = 0; r < sel.rects.length; r++) {
    const rect = sel.rects[r];
    const rows = rect.a.system === lay.rowsSystem ? rect.a : rect.b.system === lay.rowsSystem ? rect.b : undefined;
    const cols = rows === rect.a ? rect.b : rect.a;
    if (!rows || cols.system !== lay.system) continue;
    for (let c = 0; c + 1 < cols.ranges.length; c += 2) {
      const xa = lay.axis.x(cols.ranges[c]);
      const xb = lay.axis.x(cols.ranges[c + 1]);
      if (Number.isNaN(xa) || Number.isNaN(xb)) continue;
      SelPaint[2] = Math.min(xa, xb);
      SelPaint[3] = Math.max(xa, xb);
      paintRows(lay, out, rows.ranges);
    }
  }
}

/** Whether a data row is in sorted half-open `ranges` (a binary search). */
function inRanges(ranges: ArrayLike<number>, row: number): boolean {
  let lo = 0;
  let hi = ranges.length >> 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (ranges[2 * mid + 1] <= row) lo = mid + 1;
    else hi = mid;
  }
  return lo < ranges.length >> 1 && ranges[2 * lo] <= row;
}

/** The display rows in view (the pinned band, then the scrolled rows) whose data rows `ranges` holds, each once. */
function paintRows(lay: RowsLay, out: PassWriter, ranges: ArrayLike<number>) {
  for (let i = 0; i < lay.display.pinned; i++) paintRow(lay, out, i, ranges);
  for (let i = lay.first; i < lay.last; i++) paintRow(lay, out, i, ranges);
}

function paintRow(lay: RowsLay, out: PassWriter, i: number, ranges: ArrayLike<number>) {
  const e = lay.display.rowAt(i);
  let share = 0;
  if (e >= 0) share = inRanges(ranges, e) ? 1 : 0;
  else {
    const m = lay.display.members.get(i);
    if (!m) return;
    let k = 0;
    for (let q = 0; q < m.length; q++) if (inRanges(ranges, m[q])) k++;
    share = k / m.length;
  }
  if (share <= 0) return;
  const { box } = lay;
  const y = lay.rowY(i);
  if (y + lay.rowScale <= box.y || y >= box.y + box.height) return;
  const g = out.args;
  g[0] = SelPaint[2];
  g[1] = y;
  g[2] = (SelPaint[3] - SelPaint[2]) * share;
  g[3] = lay.rowScale;
  g[4] = SelPaint[0];
  g[5] = share === 1 ? SelPaint[1] : NoColor;
  g[6] = share === 1 ? lay.config.outline.selection : 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The keyboard cursor on a cell of two axes (_spec/README.md "Render", The rows layer), on its own section: its cell, or
 * its row's header (the label column, `labelX` to `numberX`), outlined in the neutral ink, a squircle by the cell (every
 * emphasis a squircle) inside it by half its width, square on a side the view or the band cuts and below the minimum. A
 * folded row's at its summary row, a hidden one's at the nearest row shown (`cursorRow`); a pinned row's in the band and
 * in its place too (by data row, as the selection).
 */
function drawCursor(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const c = emphasis?.cursor;
  if (!c || c.ySystem !== lay.rowsSystem || c.yPosition === undefined) return;
  // the layer's address is its section's and more: `startsWith` and a code read allocate nothing
  const s = c.section;
  if (s !== undefined && !(lay.key.startsWith(s) && lay.key.charCodeAt(s.length) === 47)) return;
  const i = lay.cursorRow(c.yPosition);
  if (i < 0) return;
  const header = c.part === 'header';
  if (header) {
    CellX[0] = lay.labelX;
    CellX[1] = lay.numberX;
  } else {
    if (c.system !== lay.system || !lay.axis.cell(c.position)) return;
    CellX[0] = Math.min(lay.axis.sx0, lay.axis.sx1);
    CellX[1] = Math.max(lay.axis.sx0, lay.axis.sx1);
  }
  cursorCell(lay, out, i, header);
  const own = lay.display.placeRow(c.yPosition);
  if (own >= 0 && own !== i) cursorCell(lay, out, own, header);
}

/** The cursor's outline at display row `i` between `CellX[0]` and `CellX[1]` (`drawCursor`). */
function cursorCell(lay: RowsLay, out: PassWriter, i: number, header: boolean) {
  const { box, tokens, config } = lay;
  const a = CellX[0];
  const b = CellX[1];
  const y = lay.rowY(i);
  const rh = lay.rowScale;
  if (y + rh <= box.y || y >= box.y + box.height || !(b > a)) return;
  const w = config.outline.cursor;
  const cw = b - a;
  const r = squircleRadius(cw, rh, config.squircleRatio, config.squircleMax, config.squircleMin);
  // cut: past the data box's sides (the view; the header is whole), above or below the band (its rows scroll in it)
  let round = header ? Round.both : (a >= box.x ? Round.left : 0) | (b <= box.x + box.width ? Round.right : 0);
  if (y < box.y) round |= Round.cutTop;
  if (y + rh > box.y + box.height) round |= Round.cutBottom;
  const g = out.args;
  g[0] = a + w / 2;
  g[1] = y + w / 2;
  g[2] = Math.max(cw - w, 0.5);
  g[3] = Math.max(rh - w, 0.5);
  g[4] = NoColor;
  g[5] = tokens.colors[config.emphasis.cursorInk];
  g[6] = w;
  g[7] = r > 0 && (round & Round.both) !== 0 ? Math.max(r - w / 2, 0) : 0;
  g[8] = round;
  out.rectArgs();
}

const RowsDraws = passDraws<RowsLay, Emphasis | undefined>({
  data: drawData,
  labels: drawHeader,
  hover: drawHover,
  tip: hoverColumns,
  selection: drawSelection,
  cursor: drawCursor,
});

/**
 * The rows in view as items (a row an item: its label, its place), a folded clade's summary row by its summary label
 * (`HBA +12`, its extent its first row's: EC-100); none until the rows are placed.
 */
function rowItems(lay: RowsLay): ItemDescription[] {
  const out: ItemDescription[] = [];
  const d = lay.display;
  if (!lay.axis || !(lay.rowScale > 0)) return out;
  const first = d.pinned + Math.max(0, Math.floor(lay.top));
  for (let i = 0; i < d.count; i++) {
    // the pinned band, then the scrolled rows from the band's top
    if (i === d.pinned) i = first;
    if (i >= d.count) break;
    const y = lay.rowY(i);
    if (y >= lay.box.y + lay.box.height) break;
    const e = d.rowAt(i);
    const row = e >= 0 ? e : (d.members.get(i)?.[0] ?? -1);
    if (row < 0 || y + lay.rowScale <= lay.box.y) continue;
    out.push({
      key: e >= 0 ? e : lay.displayLabel(i),
      label: lay.displayLabel(i),
      start: row,
      end: row + 1,
      x0: lay.box.x,
      x1: lay.box.x + lay.box.width,
      y: Math.max(y, lay.box.y),
      height: Math.min(y + lay.rowScale, lay.box.y + lay.box.height) - Math.max(y, lay.box.y),
      system: lay.rowsSystem,
    });
  }
  return out;
}
