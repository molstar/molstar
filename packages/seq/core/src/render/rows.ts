/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, exactAddress, issueOf, type Value } from '@molstar/core/repr/tree';
import { type CellsCodes, DefaultMaxBlockValues } from '@molstar/seqx-core/data';
import { matchLeaves, readNewick } from '@molstar/seqx-core/data/alignment-readers';
import { rangesLocation } from '@molstar/seqx-core/query';
import { rowsSystemOf, type StaticModel } from '@molstar/seqx-core/resolve';
import { type Channel, lettersDefault, type ThemeBody } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, channelLegends, type LayerChannels, PlainDatum } from './channels.js';
import type { RenderConfig } from './config.js';
import {
  Anchor,
  drawPass,
  grown,
  type ImageData,
  NoColor,
  type PassName,
  type PassWriter,
  passDraws,
} from './draw-lists.js';
import { ImageCache, ImageTile } from './images.js';
import { RenderIssues } from './issues.js';
import { type LayerSource, layerSource, readTiles } from './layer-data.js';
import {
  type Box,
  type Emphasis,
  type Hit,
  type LayerType,
  type LayoutContext,
  type Overlay,
  sectionId,
  type StoreReader,
} from './layer-type.js';
import { NoClass, Pyramid, pyramidLevelFor } from './pyramid.js';
import {
  buildClades,
  type Clades,
  type ColumnMeasures,
  columnMeasures,
  consensus,
  DisplayRows,
  type HeldCells,
  isGap,
  measureEntries,
  upper,
} from './rows-model.js';
import { inkOn, type SeqTokens } from './tokens.js';

// The `rows` layer type (_spec/README.md "Render", The rows layer): an alignment's rows over its columns, in data,
// tree or a field's order, clades folded to a summary row, rows pinned in a band at the top, hidden rows out; its own
// y inside its band (`bandTop` the rows' scroll, `rowScale` px a row); its header (the tree, the labels, the first
// residue numbers, the hidden-rows chip) in the view's label column. Held alignments (rows × columns within
// `maxBlockValues`) are read whole once; larger ones lazily, the rows in view only. Below 3 device px a cell on either
// axis the rows draw from the shared reduction as image tiles in the shared image cache.

const AnchorStart: number = Anchor.start;
const AnchorEnd: number = Anchor.end;
const AnchorMiddle: number = Anchor.middle;

type Alphabet = 'protein' | 'dna' | 'rna';

const ColorFallback: Record<Alphabet, ThemeBody> = {
  protein: { kind: 'scheme', params: { name: lettersDefault('protein') } },
  dna: { kind: 'scheme', params: { name: lettersDefault('dna') } },
  rna: { kind: 'scheme', params: { name: lettersDefault('rna') } },
};

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
  readonly summary: 'consensus' | 'logo' | 'bar';
  readonly numbers: 'gutters' | 'none';
  readonly treeWidth: number;
  readonly reference?: string;
  /** The layer's theme slots (its channels are made again with the data's measures). */
  readonly channelsForms: ResolvedLayer['theme'];
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
  heldKey = '';
  clades: Clades | undefined;
  cladesOf: unknown;
  readonly display = new DisplayRows();
  displayKey = '';
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
  /** The rows' transform: px a row, the scroll, the pinned band's height, the first and last display rows drawn. */
  rowScale = 20;
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
  /** Image tiles drawn: their boxes and images. */
  images: ImageData[] = [];
  ib: number[] = [];
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
  /** The shared image cache (the frame's, else the layer's own). */
  cache: ImageCache | undefined;

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

  /** The label of a data row: held, else from a loaded tile, else `row n`. */
  label(row: number): string {
    return this.held?.labels[row] ?? this.lazyLabels.get(row) ?? `row ${(row + 1).toLocaleString('en-US')}`;
  }

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
    const folds = ctx.folds ?? prep.collapsed;
    const key = `${held?.version ?? -1}|${this.total}`;
    if (key === this.displayKey && folds === this.folds) return;
    this.displayKey = key;
    this.folds = folds;
    if (!held) {
      this.display.setIdentity(this.total ?? 0);
      this.foldIds = [];
      return;
    }
    const byLabel = new Map<string, number>();
    held.labels.forEach((l, i) => {
      if (!byLabel.has(l)) byLabel.set(l, i);
    });
    const rowsOf = (labels: readonly string[]) => labels.map((l) => byLabel.get(l) ?? -1).filter((r) => r >= 0);
    const tree = prep.order === 'tree' ? this.clades : undefined;
    this.foldIds = tree ? folds.map((f) => tree.name.indexOf(f)).filter((c) => c >= 0) : [];
    let order: Int32Array | undefined;
    if (typeof prep.order === 'object') order = fieldOrder(held, prep.order);
    const plain = !tree && !order && !prep.pinned.length && !prep.hidden.length;
    if (plain) this.display.setIdentity(held.rows);
    else
      this.display.set(
        held.rows,
        order,
        tree,
        this.foldIds,
        rowsOf(prep.pinned),
        new Set(rowsOf(prep.hidden)),
        ctx.config.rows.pinInPlace,
      );
  }

  /** The measures and the color channel (the measures registered for its themes), when the data changed. */
  prepareColors(prep: RowsPrep, ctx: LayoutContext): Channel {
    const held = this.held;
    if (held && this.measuresOf !== held) {
      this.measuresOf = held;
      const ref = prep.reference === undefined ? -1 : held.labels.indexOf(prep.reference);
      this.measures = columnMeasures(held, ref);
      // a new channel set: its measures are this data's
      this.channels = channelForms(prep.channelsForms);
    }
    const channels = this.channels ?? (this.channels = channelForms(prep.channelsForms));
    const env = {
      alphabet: prep.alphabet,
      ...(held && this.measures ? { measures: measureEntries(held, this.measures) } : {}),
    };
    const channel = channels.get('color', ColorFallback[prep.alphabet], ctx, env);
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

/** A permutation of the held rows by a field of the rows table (`label`, `first`): missing last, ties in data order. */
function fieldOrder(held: HeldCells, o: { readonly field: string; readonly direction: string }): Int32Array {
  const idx = Array.from({ length: held.rows }, (_, i) => i);
  const sign = o.direction === 'descending' ? -1 : 1;
  if (o.field === 'first') idx.sort((a, b) => sign * (held.first[a] - held.first[b]) || a - b);
  else
    idx.sort(
      (a, b) => sign * (held.labels[a] < held.labels[b] ? -1 : held.labels[a] > held.labels[b] ? 1 : 0) || a - b,
    );
  return Int32Array.from(idx);
}

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
  return {
    ...(src ? { src } : {}),
    owner,
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
    summary: summary === 'logo' || summary === 'bar' ? summary : 'consensus',
    numbers: p(layer, 'numbers') === 'none' ? 'none' : 'gutters',
    treeWidth: typeof width === 'number' ? width : 120,
    ...(a?.reference !== undefined ? { reference: a.reference } : {}),
  };
}

type Prep = RowsPrep;

/** The layer's own image cache when the frame gives none (a test's lone layer). */
const OwnImages = new WeakMap<RenderConfig, ImageCache>();

export const RowsType: LayerType<Prep, RowsLay> = {
  name: 'rows',
  axes: 2,

  prepare: prepareRows,

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: RowsLay): RowsLay {
    const { axis, box, config, tokens } = ctx;
    const lay = prev ?? new RowsLay();
    lay.key = ctx.key;
    lay.section = sectionId(ctx.section);
    lay.system = axis.system;
    lay.rowsSystem = prep.rowsSystem;
    lay.owner = prep.owner;
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.config = config;
    lay.n = 0;
    lay.images.length = 0;
    lay.ib.length = 0;
    lay.labelText.length = 0;
    lay.labelY.length = 0;
    lay.labelW.length = 0;
    lay.numberText.length = 0;
    lay.numberY.length = 0;
    lay.numberW.length = 0;
    lay.tree.length = 0;
    lay.wedges.length = 0;
    lay.nodes.length = 0;
    lay.pending.length = 0;
    lay.lazyLabels.clear();
    lay.lazyCodes.clear();
    lay.chip = '';
    if (!lay.cache) {
      let own = OwnImages.get(config);
      if (!own) OwnImages.set(config, (own = new ImageCache(config.imageBytes)));
      lay.cache = own;
    }
    const src = prep.src;
    if (!src) return lay;
    const coverage = data.coverage(src.source);
    const rowsInfo = coverage?.systems.find((s) => s.name === src.system);
    const colsInfo = coverage?.systems.find((s) => s.name !== src.system);
    lay.total = rowsInfo?.length;
    lay.columns = colsInfo?.length ?? axis.length ?? 0;
    const rowsSrc: LayerSource = { ...src, ...(lay.total !== undefined ? { length: lay.total } : {}) };
    delete (rowsSrc as { path?: unknown }).path;
    delete (rowsSrc as { back?: unknown }).back;
    lay.rowScale = ctx.rowScale ?? config.rows.rowHeight;
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
      const read = readTiles(ctx, data, rowsSrc, [0, N], 1);
      if (!read.pending.length && !read.failed.length) joinHeld(lay, read.entries, N);
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
    }
    lay.prepareRows(prep, ctx);
    const display = lay.display;
    // the rows' transform: the pinned band fixed, the rest scrolled within the band, clamped
    const rs = lay.rowScale;
    lay.pinnedHeight = display.pinned * rs;
    const area = Math.max(0, box.height - lay.pinnedHeight);
    const content = (display.count - display.pinned) * rs;
    lay.bandTop = Math.max(0, Math.min(ctx.bandTop ?? 0, Math.max(0, content - area)));
    lay.first = display.pinned + Math.floor(lay.bandTop / rs);
    lay.last = Math.min(display.count, display.pinned + Math.ceil((lay.bandTop + area) / rs));
    if (!heldable) readLazy(lay, ctx, data, rowsSrc);
    const channel = lay.prepareColors(prep, ctx);
    if (!lay.held && heldable) return lay;
    layoutCells(lay, prep, ctx, channel);
    layoutHeader(lay, prep, ctx);
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
    return lay.channel ? channelLegends([lay.channel], { ...options, layer: lay.key }) : [];
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const i = lay.rowAtY(y);
    if (x < box.x) {
      // the header: a tree node within reach of its branch, else the row's label
      const node = nearestNode(lay, x, y);
      if (node >= 0 && lay.clades)
        return {
          layer: lay.key,
          section: lay.section,
          system: lay.rowsSystem,
          target: 'tree-node',
          element: { owner: lay.owner, key: lay.clades.name[node] },
          extent: cladeExtent(lay, node),
        };
      if (i < 0 || x < lay.treeX1) return undefined;
      return rowHit(lay, i, undefined);
    }
    if (i < 0) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return rowHit(lay, i, position);
  },
};

/** A hit on display row `i`: a summary row's (its clade), else its data row's, a cell where `position` is given. */
function rowHit(lay: RowsLay, i: number, position: number | undefined): Hit {
  const e = lay.display.rowAt(i);
  if (e < 0 && lay.clades) {
    const c = -1 - e;
    return {
      layer: lay.key,
      section: lay.section,
      system: position === undefined ? lay.rowsSystem : lay.system,
      ...(position !== undefined ? { position } : {}),
      target: 'summary-row',
      element: { owner: lay.owner, key: lay.clades.name[c] },
      extent: cladeExtent(lay, c, i),
    };
  }
  if (position === undefined)
    return {
      layer: lay.key,
      section: lay.section,
      system: lay.rowsSystem,
      position: e,
      ySystem: lay.rowsSystem,
      yPosition: e,
      target: 'row',
    };
  return {
    layer: lay.key,
    section: lay.section,
    system: lay.system,
    position,
    ySystem: lay.rowsSystem,
    yPosition: e,
    target: 'cell',
  };
}

/** A clade's rows as a Location on the rows system (its shown rows when folded at display row `i`). */
function cladeExtent(lay: RowsLay, c: number, i?: number) {
  const members = i === undefined ? undefined : lay.display.members.get(i);
  const ranges: number[] = [];
  if (members) for (const r of members) ranges.push(r, r + 1);
  else if (lay.clades)
    for (let k = lay.clades.lo[c]; k < lay.clades.hi[c]; k++) {
      const r = lay.clades.leafRow[k];
      if (r >= 0) ranges.push(r, r + 1);
    }
  return rangesLocation(lay.rowsSystem, ranges);
}

/** The tree node nearest a point within `rows.treeHitPx` of its branch (a wedge's triangle included), or -1. */
function nearestNode(lay: RowsLay, x: number, y: number): number {
  if (x < lay.treeX0 || x > lay.treeX1 + 2) return -1;
  const reach = lay.config.rows.treeHitPx;
  let best = -1;
  let bestD = reach;
  const w = lay.wedges;
  for (let k = 0; k + 4 < w.length; k += 5)
    if (x >= w[k] - reach && x <= w[k + 2] && y >= w[k + 1] && y <= w[k + 3]) return w[k + 4];
  const n = lay.nodes;
  for (let k = 0; k + 2 < n.length; k += 3) {
    const d = Math.max(Math.abs(x - n[k + 1]), Math.abs(y - n[k + 2]));
    if (d <= bestD) {
      bestD = d;
      best = n[k];
    }
  }
  return best;
}

/** The held cells joined from the tiles in order, made again only when an entry changed. */
function joinHeld(
  lay: RowsLay,
  entries: readonly { readonly version: number; readonly tables: { cells?: unknown } }[],
  N: number,
) {
  let key = `${N}`;
  for (const e of entries) key += `,${e.version}`;
  if (lay.held && key === lay.heldKey) return;
  const C = lay.columns;
  const codes = new Uint8Array(N * C);
  const labels: string[] = new Array(N).fill('');
  const first = new Int32Array(N).fill(1);
  for (const e of entries) {
    const t = e.tables.cells as CellsCodes | undefined;
    if (!t || t.form !== 'codes' || t.columns !== C) continue;
    const rows = Math.floor(t.codes.length / Math.max(C, 1));
    codes.set(t.codes.subarray(0, Math.min(rows, N - t.start) * C), t.start * C);
    for (let i = 0; i < rows && t.start + i < N; i++) {
      labels[t.start + i] = t.label?.get(i) ?? `row ${t.start + i + 1}`;
      if (t.first) first[t.start + i] = t.first[i];
    }
  }
  lay.heldKey = key;
  lay.held = { rows: N, columns: C, codes, labels, first, version: (lay.held?.version ?? 0) + 1 };
}

/** A lazy alignment: the tiles of the rows in view at level 1 (their labels and cells), the rest pending. */
function readLazy(lay: RowsLay, ctx: LayoutContext, data: StoreReader, src: LayerSource) {
  const { first, last } = lay;
  if (last <= first) return;
  // bounded by the screen: a band of more rows than device pixels reads none (a navigator's coarse rows, later)
  if ((last - first) * lay.rowScale * ctx.dpr < (last - first) * 0.99) {
    ctx.report(
      issueOf(RenderIssues, 'seq.render.rows.reduced', 'The rows are below a pixel and the source serves no row bins', {
        address: ctx.key,
        hint: RenderIssues.codes['seq.render.rows.reduced'].hint,
      }),
    );
    lay.pending.push(lay.box.y, lay.box.height);
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
    }
  }
  for (let i = 0; i + 1 < read.pending.length; i += 2)
    for (let r = read.pending[i]; r < read.pending[i + 1]; r++) {
      lay.pending.push(lay.rowY(r), lay.rowScale);
    }
}

/** The cells in view: rects (letters at the letters step), or image tiles from the reduction below 3 device px. */
function layoutCells(lay: RowsLay, prep: Prep, ctx: LayoutContext, channel: Channel) {
  const { axis, config, tokens, box } = ctx;
  const dpr = axis.dpr;
  const colPx = axis.devicePx;
  const rowPx = lay.rowScale * dpr;
  const vis = axis.visible([]);
  if (!vis.length) return;
  const c0 = Math.max(0, vis[0]);
  const c1 = Math.min(lay.columns, vis[vis.length - 1]);
  const missing = tokens.colors['seq-missing'];
  if ((colPx < config.density.cells || rowPx < config.density.cells) && lay.held) {
    layoutImages(lay, ctx, c0, c1);
    return;
  }
  const lettersStep = axis.scale >= config.density.letters && lay.rowScale >= config.density.letters;
  lay.letters = lettersStep;
  if (lettersStep) {
    let size = 0;
    for (const s of config.letterSizes)
      if (s * config.monoAdvance <= axis.scale - 2 && s <= lay.rowScale * 0.75) size = s;
    lay.fontSize = size || (config.letterSizes[0] ?? 9);
    lay.font = `${lay.fontSize}px mono`;
  }
  const held = lay.held;
  const C = lay.columns;
  for (let i = lay.first - (lay.first > lay.display.pinned ? 0 : 0); i < lay.last; i++)
    addRow(lay, i, c0, c1, held, C, missing, tokens, prep);
  for (let i = 0; i < lay.display.pinned; i++) addRow(lay, i, c0, c1, held, C, missing, tokens, prep);
  void channel;
  void box;
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
 * Below 3 device px a cell: the shared reduction over the display rows' theme classes (`rows.reduce`: a class a
 * distinct fill), at the level pair whose sample is at least a device pixel, drawn as image tiles (256 × 256 samples)
 * kept in the shared image cache by (layer, level, tile, base, theme).
 */
function layoutImages(lay: RowsLay, ctx: LayoutContext, c0: number, c1: number) {
  const { axis, tokens } = ctx;
  const held = lay.held!;
  const display = lay.display;
  const C = held.columns;
  const pyr = (lay.pyramid ??= new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' }));
  const baseKey = `${held.version}|${display.version}|${tokens.version}|${lay.channelVersion}`;
  const own = lay;
  if (own.pyramidOf !== baseKey) {
    own.pyramidOf = baseKey;
    const R = display.count;
    const need = R * C;
    if (!own.baseBytes || own.baseBytes.length < need) own.baseBytes = new Uint8Array(need);
    const base = own.baseBytes;
    // one class a distinct fill (at most 254): the theme's classes, so a reduced sample reads as the theme
    const classOf = new Map<number, number>();
    const missing = tokens.colors['seq-missing'];
    for (let i = 0; i < R; i++) {
      const e = display.rowAt(i);
      const sum = e < 0 ? summaryOf(lay, i) : undefined;
      for (let j = 0; j < C; j++) {
        const code = sum ? sum[3 * j] : held.codes[e * C + j];
        const fill = lay.colorOf(code, j, missing);
        let k: number;
        if (fill === NoColor) k = NoClass;
        else {
          let c = classOf.get(fill);
          if (c === undefined) {
            c = classOf.size < 254 ? classOf.size : 253;
            if (classOf.size < 254) classOf.set(fill, c);
            lay.classColors[c] = fill;
          }
          k = c;
        }
        base[i * C + j] = k;
      }
    }
    pyr.setClasses(R, C, base.subarray(0, need));
  }
  const dpr = axis.dpr;
  const rowLevel = pyramidLevelFor(1 / (lay.rowScale * dpr), 4);
  const colLevel = pyramidLevelFor(1 / axis.devicePx, 4);
  const level = pyr.level(rowLevel, colLevel);
  const fr = 4 ** level.rowLevel;
  const fc = 4 ** level.columnLevel;
  const s0 = Math.floor(lay.first / fr);
  const s1 = Math.ceil(lay.last / fr);
  const t0 = Math.floor(c0 / fc);
  const t1 = Math.ceil(c1 / fc);
  const cache = lay.cache!;
  for (let ti = Math.floor(s0 / ImageTile); ti * ImageTile < s1; ti++)
    for (let tj = Math.floor(t0 / ImageTile); tj * ImageTile < t1; tj++) {
      const key = `${lay.key}|${level.rowLevel}.${level.columnLevel}|${ti}.${tj}|${own.pyramidOf}`;
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
            // opaque by the majority's share: a mixed sample reads lighter
            const a = ((fill & 0xff) * (0.35 + (0.65 * level.shares[o]) / 255)) | 0;
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
      const y0 = lay.rowY(rowStart);
      const y1 = lay.rowY(Math.min(rowStart + rows * fr, display.count));
      lay.images.push(img);
      lay.ib.push(Math.min(x0, x1), y0, Math.abs(x1 - x0), y1 - y0);
    }
}

/** The header in the label column: the tree, the labels, the first residue numbers. */
function layoutHeader(lay: RowsLay, prep: Prep, ctx: LayoutContext) {
  const { config, box, text } = ctx;
  const rs = lay.rowScale;
  const left = config.inset;
  const right = box.x - config.labelGap;
  const treeW = lay.clades && prep.order === 'tree' ? Math.min(prep.treeWidth, (right - left) * 0.45) : 0;
  lay.treeX0 = left;
  lay.treeX1 = left + treeW;
  if (treeW > 0) layoutTree(lay, lay.clades!, left, treeW);
  if (rs < config.rows.labelMinPx) return;
  lay.labelFont = rs >= 16 ? config.labelFont : config.itemFont;
  lay.labelX = left + treeW + (treeW ? 4 : 0);
  const numberW = prep.numbers === 'gutters' ? text.measure(config.itemFont, '00000') : 0;
  lay.numberX = right;
  lay.labelMax = Math.max(0, right - numberW - 4 - lay.labelX);
  const size = rs >= 16 ? 12 : 10;
  const add = (i: number) => {
    const e = lay.display.rowAt(i);
    const y = lay.rowY(i) + rs / 2 + size * 0.35;
    let t: string;
    if (e < 0) {
      const m = lay.display.members.get(i);
      t = m ? `${lay.label(m[0])} +${m.length - 1}` : '';
    } else t = lay.label(e);
    if (lay.display.pinned && i >= lay.display.pinned && e >= 0 && config.rows.pinInPlace && isPinned(lay, e))
      t = `· ${t}`;
    lay.labelText.push(t);
    lay.labelY.push(y);
    lay.labelW.push(text.measure(lay.labelFont, t));
    if (prep.numbers === 'gutters' && e >= 0) {
      const n = lay.held ? String(lay.held.first[e]) : '';
      if (n) {
        lay.numberText.push(n);
        lay.numberY.push(y);
        lay.numberW.push(text.measure(config.itemFont, n));
      }
    }
  };
  for (let i = 0; i < lay.display.pinned; i++) add(i);
  for (let i = lay.first; i < lay.last; i++) add(i);
}

function isPinned(lay: RowsLay, row: number): boolean {
  for (let i = 0; i < lay.display.pinned; i++) if (lay.display.rowAt(i) === row) return true;
  return false;
}

/**
 * The tree beside the rows: a node at its distance from the root (branch lengths; depth when none), a leaf at its row's
 * middle, an inner node between its first and last child; a folded clade a wedge at its summary row. Only what crosses
 * the band is kept (segments and the nodes a hit reads).
 */
function layoutTree(lay: RowsLay, clades: Clades, x0: number, w: number) {
  const n = clades.nodes;
  const ys = new Float64Array(n).fill(Number.NaN);
  let maxD = 0;
  for (let k = 0; k < n; k++) if (clades.distance[k] > maxD) maxD = clades.distance[k];
  const byDepth = maxD <= 0;
  let maxDepth = 0;
  if (byDepth) for (let k = 0; k < n; k++) if (clades.depth[k] > maxDepth) maxDepth = clades.depth[k];
  const xOf = (k: number) => x0 + (byDepth ? clades.depth[k] / Math.max(1, maxDepth) : clades.distance[k] / maxD) * w;
  const rs = lay.rowScale;
  const folded = new Set(lay.display.members.size ? lay.foldIds : []);
  const box = lay.box;
  // a folded clade's rows show at its summary row: its node's y is that row's middle, nothing below it is drawn
  const inFold = new Uint8Array(n);
  for (let k = 0; k < n; k++) {
    const p = clades.parent[k];
    if (p >= 0 && (inFold[p] || folded.has(p))) inFold[k] = 1;
  }
  for (let k = n - 1; k >= 0; k--) {
    if (inFold[k]) continue;
    const kids = clades.children[k];
    if (!kids.length || folded.has(k)) {
      const leafRow = kids.length ? -1 : clades.leafRow[clades.lo[k]];
      const row = leafRow >= 0 ? leafRow : firstShown(lay, clades, k);
      const d = row >= 0 ? lay.display.displayRow(row) : -1;
      if (d >= 0) ys[k] = lay.rowY(d) + rs / 2;
      continue;
    }
    let lo = Number.NaN;
    let hi = Number.NaN;
    for (const c of kids) {
      const y = ys[c];
      if (Number.isNaN(y)) continue;
      if (Number.isNaN(lo) || y < lo) lo = y;
      if (Number.isNaN(hi) || y > hi) hi = y;
    }
    if (!Number.isNaN(lo)) ys[k] = (lo + hi) / 2;
  }
  const top = box.y + lay.pinnedHeight;
  const bottom = box.y + box.height;
  for (let k = 0; k < n; k++) {
    if (inFold[k] || Number.isNaN(ys[k])) continue;
    const x = xOf(k);
    const kids = folded.has(k) ? [] : clades.children[k];
    if (folded.has(k)) {
      const y = ys[k];
      if (y + rs / 2 >= top && y - rs / 2 <= bottom) {
        lay.wedges.push(x, y - rs * 0.4, x0 + w, y + rs * 0.4, k);
        lay.nodes.push(k, x, y);
      }
      continue;
    }
    let lo = ys[k];
    let hi = ys[k];
    for (const c of kids) {
      const y = ys[c];
      if (Number.isNaN(y)) continue;
      if (y < lo) lo = y;
      if (y > hi) hi = y;
      // the child's branch: from this node's x to the child's, at the child's y
      if (y >= top - 1 && y <= bottom + 1) lay.tree.push(x, y, xOf(c), y);
    }
    if (kids.length && hi >= top && lo <= bottom) {
      lay.tree.push(x, Math.max(lo, top), x, Math.min(hi, bottom));
      if (ys[k] >= top && ys[k] <= bottom) lay.nodes.push(k, x, ys[k]);
    }
  }
}

/** A clade's first shown row (its summary row's), -1 when none is shown. */
function firstShown(lay: RowsLay, clades: Clades, k: number): number {
  for (let l = clades.lo[k]; l < clades.hi[k]; l++) {
    const r = clades.leafRow[l];
    if (r >= 0 && lay.display.displayRow(r) >= 0) return r;
  }
  return -1;
}

// ---- drawing (nothing allocated once warm: numbers through the writer's `args`) ----

function drawData(lay: RowsLay, out: PassWriter) {
  const { box } = lay;
  const g = out.args;
  g[0] = box.x;
  g[1] = box.y;
  g[2] = box.width;
  g[3] = box.height;
  out.clipArgs();
  for (let k = 0; k < lay.images.length; k++) {
    g[0] = lay.ib[4 * k];
    g[1] = lay.ib[4 * k + 1];
    g[2] = lay.ib[4 * k + 2];
    g[3] = lay.ib[4 * k + 3];
    out.imageArgs(lay.images[k]);
  }
  for (let k = 0; k < lay.n; k++) {
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
  }
  if (lay.pending.length) drawPending(lay, out);
  if (lay.letters && lay.n) drawLetters(lay, out);
  if (lay.display.pinned) drawPinnedRule(lay, out);
  out.unclip();
}

function drawLetters(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  let y = Number.NaN;
  let open = false;
  for (let k = 0; k < lay.n; k++) {
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

/** The header: the tree's branches and wedges, the labels faded past their width, the first numbers, the chip. */
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
  for (let k = 0; k + 3 < t.length; k += 4) {
    const vertical = t[k] === t[k + 2];
    g[0] = vertical ? t[k] - 0.5 : Math.min(t[k], t[k + 2]);
    g[1] = vertical ? Math.min(t[k + 1], t[k + 3]) : t[k + 1] - 0.5;
    g[2] = vertical ? 1 : Math.abs(t[k + 2] - t[k]);
    g[3] = vertical ? Math.abs(t[k + 3] - t[k + 1]) : 1;
    g[4] = ink;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  if (lay.wedges.length) drawWedges(lay, out);
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
  out.unclip();
}

const Wedge = new Float64Array(6);
function drawWedges(lay: RowsLay, out: PassWriter) {
  const w = lay.wedges;
  const g = out.args;
  for (let k = 0; k + 4 < w.length; k += 5) {
    Wedge[0] = w[k];
    Wedge[1] = (w[k + 1] + w[k + 3]) / 2;
    Wedge[2] = w[k + 2];
    Wedge[3] = w[k + 1];
    Wedge[4] = w[k + 2];
    Wedge[5] = w[k + 3];
    g[4] = lay.tokens.colors['msui-text-muted'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    out.pathArgs(Wedge, 3, true);
  }
}

/** The hover: its row's band tinted lightly, its column within the rows, the cell outlined. */
function drawHover(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key) return;
  const { box, tokens, config } = lay;
  const g = out.args;
  const ink = tokens.colors['seq-hover'];
  const row = h.yPosition ?? (h.system === lay.rowsSystem ? h.position : undefined);
  const i = row === undefined ? -1 : lay.display.displayRow(row);
  const y = i >= 0 ? lay.rowY(i) : Number.NaN;
  if (i >= 0 && y + lay.rowScale > box.y && y < box.y + box.height) {
    g[0] = box.x;
    g[1] = y;
    g[2] = box.width;
    g[3] = lay.rowScale;
    g[4] = ((ink & 0xffffff00) | 0x1a) >>> 0;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
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
  if (Number.isNaN(y)) return;
  const d = config.outline.hover / 2 + config.emphasis.hoverOffset;
  g[0] = a - d;
  g[1] = y - d;
  g[2] = b - a + 2 * d;
  g[3] = lay.rowScale + 2 * d;
  g[4] = NoColor;
  g[5] = ink;
  g[6] = config.outline.hover;
  out.rectArgs();
}

/**
 * The selection along y: a rows part (data rows) as row bands over the data box, tinted with an outline; a rect of
 * rows × columns as its cells' boxes. A row inside a folded clade lights its summary row.
 */
function drawSelection(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const sel = emphasis?.selection;
  if (!sel) return;
  const { box, tokens, config } = lay;
  const accent = tokens.colors['msui-accent'];
  const tint = ((accent & 0xffffff00) | Math.round(255 * config.emphasis.selectionTint)) >>> 0;
  const line = ((accent & 0xffffff00) | Math.round(255 * config.emphasis.selectionOutline)) >>> 0;
  const g = out.args;
  for (let p = 0; p < sel.parts.length; p++) {
    const part = sel.parts[p];
    if (part.system !== lay.rowsSystem) continue;
    for (let k = 0; k + 1 < part.ranges.length; k += 2)
      paintRows(lay, out, part.ranges[k], part.ranges[k + 1], box.x, box.x + box.width, tint, line);
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
      for (let k = 0; k + 1 < rows.ranges.length; k += 2)
        paintRows(lay, out, rows.ranges[k], rows.ranges[k + 1], Math.min(xa, xb), Math.max(xa, xb), tint, line);
    }
  }
  void g;
}

/** Data rows `[a, b)` as the display rows they are shown at, between x0 and x1: a tint and an outline each. */
function paintRows(
  lay: RowsLay,
  out: PassWriter,
  a: number,
  b: number,
  x0: number,
  x1: number,
  tint: number,
  line: number,
) {
  const g = out.args;
  const { box } = lay;
  const w = lay.config.outline.selection;
  let last = -1;
  for (let r = a; r < b; r++) {
    const i = lay.display.displayRow(r);
    if (i < 0 || i === last) continue;
    last = i;
    const y = lay.rowY(i);
    if (y + lay.rowScale <= box.y || y >= box.y + box.height) continue;
    g[0] = x0;
    g[1] = y;
    g[2] = x1 - x0;
    g[3] = lay.rowScale;
    g[4] = tint;
    g[5] = line;
    g[6] = w;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
    // a row in view past many: rows of a huge selection are bounded by the band (each row once)
    if (b - a > 4 * (lay.last - lay.first + lay.display.pinned) && lay.display.identity) {
      r = Math.max(r, lay.first - 1);
      if (r >= lay.last) break;
    }
  }
}

/** The keyboard cursor on a cell of two axes: its cell outlined in the neutral ink. */
function drawCursor(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const c = emphasis?.cursor;
  if (!c || c.ySystem !== lay.rowsSystem || c.yPosition === undefined) return;
  const i = lay.display.displayRow(c.yPosition);
  if (i < 0) return;
  const { box, tokens, config } = lay;
  const y = lay.rowY(i);
  const g = out.args;
  const w = config.outline.cursor;
  let x0 = box.x - (lay.treeX1 > 0 ? box.x - lay.labelX : 0);
  let x1 = box.x;
  if (c.system === lay.system && lay.axis.cell(c.position)) {
    x0 = lay.axis.sx0;
    x1 = lay.axis.sx1;
  }
  g[0] = x0 + w / 2;
  g[1] = y + w / 2;
  g[2] = Math.max(x1 - x0 - w, 0.5);
  g[3] = lay.rowScale - w;
  g[4] = NoColor;
  g[5] = tokens.colors[config.emphasis.cursorInk];
  g[6] = w;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

const RowsDraws = passDraws<RowsLay, Emphasis | undefined>({
  data: drawData,
  labels: drawHeader,
  hover: drawHover,
  selection: drawSelection,
  cursor: drawCursor,
});

void AnchorMiddle;
