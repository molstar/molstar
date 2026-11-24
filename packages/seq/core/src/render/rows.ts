/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { exactAddress, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import { DefaultMaxBlockValues } from '@molstar/seqx-core/data';
import { readNewick } from '@molstar/seqx-core/data/alignment-readers';
import { rowsSystemOf, type StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { channelForms, channelLegends } from './channels.js';
import type { RenderConfig } from './config.js';
import { drawPass, passDraws, type PassName, type PassWriter } from './draw-lists.js';
import { ImageCache } from './images.js';
import { RenderIssues } from './issues.js';
import { layerSource, readTiles } from './layer-data.js';
import {
  type Emphasis,
  type Hit,
  type ItemDescription,
  type LayerType,
  type LayoutContext,
  type Overlay,
  sectionId,
  type StoreReader,
} from './layer-type.js';
import { chipAt, layoutRowColumns, setColumns, titleAt } from './rows-columns.js';
import { rowColumnsOf, RowHeight } from './rows-model.js';
import { NoFolds, RowsLay, type RowsPrep } from './rows-lay.js';
import { readHeld, readJoined, readLazy, rowsSource } from './rows-read.js';
import { layoutCells } from './rows-cells.js';
import { inToggle, layoutHeader } from './rows-header.js';
import { cladeHit, named, nearestNode, rowHit, writeHit } from './rows-hit.js';
import { DefaultRowsZoom } from './rail.js';
import { drawData, drawHeader } from './rows-draw.js';
import { drawCursor, drawHover, drawSelection, hoverColumns } from './rows-emphasis.js';

export { RowsLay, type RowsOrder, type RowsPrep } from './rows-lay.js';
export { readHeld, rowsSource } from './rows-read.js';

// The `rows` layer type (_spec/README.md "Render", The rows layer): an alignment's rows over its columns, in data,
// tree or a field's order, clades folded to a summary row, rows pinned in a band at the top, hidden rows out; its own
// y inside its band (`bandTop` the rows' scroll, `rowScale` px a row); its header (the tree, the labels, the left row
// columns, the first residue numbers, the hidden-rows chip) in the view's label column, its right columns after the
// last numbers, the columns' titles in the band above the rows. Held alignments (rows × columns within
// `maxBlockValues`) are read whole once; larger ones lazily, the rows in view only. Below 3 device px a cell on either
// axis the rows draw from the shared reduction as image tiles in the shared image cache. This module is the layer type
// (its prepared form, its layout's steps, its passes, hits and items); what the layer keeps, reads, lays out (cells,
// header), hits, draws and emphasizes are the `rows-*.ts` modules beside it (`rows-model.ts`: what it computes without
// drawing).

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

/** The layer's own image cache when the frame gives none (a test's lone layer). */
const OwnImages = new WeakMap<RenderConfig, ImageCache>();

export const RowsType: LayerType<RowsPrep, RowsLay> = {
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
    lay.toggleX = -1;
    lay.toggleOpen = false;
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
    lay.rowScale = RowHeight;
    lay.rowFloor = 0;
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
      // no row levels are served: below a device pixel a row the band is blank, so the rows' zoom stops there
      lay.rowFloor = 1 / ctx.dpr;
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
    // a stored scale (a saved view, a story, a motion) is clamped where it is read: the host's bounds and the source's
    // floor (`rail.ts`), and with the zoom off ignored (the fit's)
    const zoom = ctx.rowsZoom ?? DefaultRowsZoom;
    if (ctx.rowScale !== undefined && zoom.enabled)
      lay.rowScale = Math.max(zoom.minPx, lay.rowFloor, Math.min(zoom.maxPx, ctx.rowScale));
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
    // the tree's toggle in a narrow frame: a click reopens the tree (session state)
    if (inToggle(lay, at[0], at[1])) return writeHit(lay, lay.rowsSystem, undefined, undefined, 'tree-toggle');
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
