/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArrayRange } from '@molstar/core/data/util/sort';
import { issueOf } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { type AxisTransform, Cut } from './axis.js';
import { fontSize } from './config.js';
import { drawPass, grown, type PassName, type PassWriter, Round } from './draw-lists.js';
import {
  beginItems,
  describeItems,
  itemExtents,
  itemHit,
  itemsUsed,
  markEnds,
  readItems,
  RowImages,
  tickHit,
  tickOn,
} from './items.js';
import { RenderIssues } from './issues.js';
import { channelForms, channelLegends } from './channels.js';
import { layerSource, toSection } from './layer-data.js';
import type { Emphasis, Hit, LayerType, LayoutContext, Overlay, StoreReader } from './layer-type.js';
import { squircleRadius } from './squircle.js';
import { geneGroups, geneLabelAt, geneMarks, placeGeneLabels } from './ranges-gene.js';
import {
  Gene,
  itemPiece,
  laneId,
  LineStyles,
  mappedPieces,
  newLay,
  type RangesLay,
  type RangesPrep,
  Solid,
  ThinStripShare,
} from './ranges-lay.js';
import { placeLabels, sinkLabels } from './ranges-beside.js';
import { knownLanesVersion, noteLanes } from './known.js';
import { noteReachable } from './ranges-reach.js';
import { rangesContent } from './ranges-content.js';
import { connect, connectorAt, startConnectors } from './ranges-connector.js';
import { RangesDraws, growWave, WaveStep } from './ranges-draw.js';
import { densityStrip } from './ranges-density.js';
import { themeItems } from './ranges-theme.js';

export { addToStrip, paintStrip, resetStrip, type StripBuffers } from './ranges-density.js';
export { LineStyles, type RangesLay, type RangesPrep, ThinStripShare } from './ranges-lay.js';

// The `ranges` layer type (_spec/README.md "Render", Layer types): items in sticky lanes. Each item is carried onto
// the section along the layer's path (a domain on its gene: a piece an exon) and drawn as a squircle where it is big
// enough, square below, with square ends where the view or a piece's edge cuts it (the frame rule). Below a device
// pixel a position, the marks of a lane merge into one bar a run; at the density step the layer is one heat strip, a
// device column a pixel. Items past the lanes collapse into a density row at the bottom (`ranges.overflow`:
// `extra-lanes`), or the whole layer does (`layer`). An item's key is its `id`, else its `row`, else its tile and
// index; a key and start seen in an earlier tile is not drawn again (a range repeated in the tiles it touches), and a
// key repeated with another start is another piece of one item (a wrapped row, the exons of a transcript). This module
// is the layer type (its prepared form, its layout's lanes and marks, its hit and legend); what the layer keeps, its
// content, labels beside, draw, heat strip and channels are the `ranges-*.ts` modules beside it.

/** A mark this wide (CSS px) or wider draws its item's arrowhead, at any step: narrower, the head would be all of it. */
const ArrowMinPx = 4;

/**
 * The next mark's index, its arrays grown to hold it and its label, radius and ink cleared; the caller writes its extent,
 * lane, item and ends (a number passed to a call that is not inlined is boxed: an allocation a mark).
 */
export function nextMark(lay: RangesLay): number {
  const j = lay.m++;
  if (j >= lay.mx0.length) {
    lay.mx0 = grown(lay.mx0, j + 1);
    lay.mx1 = grown(lay.mx1, j + 1);
    lay.mlane = grown(lay.mlane, j + 1);
    lay.mitem = grown(lay.mitem, j + 1);
    lay.mround = grown(lay.mround, j + 1);
    lay.mlabel = grown(lay.mlabel, j + 1);
    lay.mcut = grown(lay.mcut, j + 1);
    lay.mradius = grown(lay.mradius, j + 1);
    lay.mink = grown(lay.mink, j + 1);
    lay.mtip = grown(lay.mtip, j + 1);
    lay.manchor = grown(lay.manchor, j + 1);
    lay.mthin = grown(lay.mthin, j + 1);
  }
  lay.mthin[j] = 0;
  lay.mlabel[j] = 0;
  lay.mradius[j] = 0;
  lay.mink[j] = 0;
  lay.mtip[j] = 0;
  return j;
}

// made at first use: `items.js` and this module import each other, and a composition may load either first (a registry
// that lists bands before ranges), so nothing of `items.js` is read while the modules load
let rowImages: RowImages | undefined;

/** The style `line` names (an index of `LineStyles`; 0: blocks); another value is solid, said once. */
function lineStyle(line: string | undefined, ctx: LayoutContext): number {
  if (line === undefined) return 0;
  const i = LineStyles.indexOf(line as (typeof LineStyles)[number]);
  if (i > 0) return i;
  ctx.report(
    issueOf(RenderIssues, 'seq.render.line-style', `The line style ${JSON.stringify(line)} is not drawn`, {
      address: ctx.key,
      hint: RenderIssues.codes['seq.render.line-style'].hint,
    }),
  );
  return Solid;
}

/**
 * Where the content of piece `a` of the axis starts on the screen, unclipped: the left edge of its low boundary on a
 * forward axis, the right edge on a reversed one (the dash and the wave run from it).
 */
function anchorOf(axis: AxisTransform, a: number, s: number, e: number): number {
  const ua = axis.xOfUnit(axis.unitOf(a, Math.max(s, axis.pieces[2 * a])));
  const ub = axis.xOfUnit(axis.unitOf(a, Math.min(e, axis.pieces[2 * a + 1])));
  return axis.dir === 1 ? Math.min(ua, ub) : Math.max(ua, ub);
}

/**
 * The band's lanes (_spec/README.md "Render", Track heights): in a band tall enough, lanes of `laneHeight` below its
 * padding, the track's label lane above them when its content says labels go beside, the block of them centred in the
 * band by the lanes the whole data needs (so a pan never moves them), the marks in those lanes only (a band taller than
 * its content draws none below them). A band too short for its lanes (a track's own height, a text line's underline, a
 * navigator's miniature; not one at the most lanes, whose rest go to its density row) drops the label lane, then
 * compacts its lanes over the whole band, as many as a mark's height (`laneHeight − laneGap`) fits, at least one. A
 * tiled source's lanes (not known whole) are centred by the most the layer has shown (`seen`: it only grows, so its
 * marks never move back on a pan). Returns the lanes the marks may take.
 */
function laneGeometry(lay: RangesLay, ctx: LayoutContext): number {
  const { box, config } = ctx;
  const lh = config.laneHeight;
  const pad = config.lanePadding;
  const content = ctx.content;
  lay.labelY = -1;
  lay.scrolls = false;
  lay.top = 0;
  lay.scrollMax = 0;
  const room = box.height >= lh + 2 * pad ? Math.floor((box.height - 2 * pad) / lh) : 0;
  const need = content ? Math.max(1, content.lanes) : 1;
  // lanes past the band scroll inside it (tracks; a ring's annulus and a text's lines compact, as `overflow: compact`)
  if (
    (room < need || content?.partial === true) &&
    config.track.overflow === 'scroll' &&
    ctx.section.layout === 'tracks' &&
    (ctx.labels === undefined || ctx.labels === 'inside')
  ) {
    lay.scrolls = true;
    lay.room = Math.max(1, room);
    lay.laneHeight = lh;
    // a band under a lane and its padding centres its one lane, one under a lane shows the lane's top (EC-249)
    lay.laneY0 = box.y + Math.max(0, Math.min(pad, (box.height - lh) / 2));
    lay.laneY = lay.laneY0;
    return config.ranges.maxLanes;
  }
  // `density`: the band at its cap, the rest of the lanes in a density row at its bottom
  if (room >= need || (config.track.overflow === 'density' && room >= 1)) {
    // a label lane only beside marks the band holds too, and only where labels are placed in the band (tracks)
    const label = content?.labels && ctx.labels === 'inside' && room > need ? 1 : 0;
    // a source's lanes known so far are not its whole: its marks may take every lane the band holds
    const lanes = content && !content.partial ? Math.max(1, Math.min(need, room - label)) : room;
    const block = (content && !content.partial ? lanes : Math.max(1, Math.min(lay.seen, lanes))) + label;
    lay.laneHeight = lh;
    // centred in a band with a little to spare (one lane's marks in a thin band); a tall band's lanes start at its top,
    // evenly (not a blank stretch above and below them: a transcripts track of 47 lanes)
    const slack = box.height - block * lh;
    const top = box.y + (slack > 2 * lh ? pad : slack / 2);
    if (label) lay.labelY = top;
    lay.laneY = top + label * lh;
    lay.laneY0 = lay.laneY;
    // a source's lanes are packed past the band's (a ring's, a text's): the lanes it needs are counted, those past the band not drawn
    return content?.partial && room >= need ? config.ranges.maxLanes : lanes;
  }
  const fit = Math.floor(box.height / Math.max(1, lh - config.laneGap));
  const lanes = Math.max(1, Math.min(content ? need : fit, fit));
  lay.laneHeight = box.height / lanes;
  lay.laneY = box.y;
  lay.laneY0 = box.y;
  return lanes;
}

export const RangesType: LayerType<RangesPrep, RangesLay> = {
  name: 'ranges',

  prepare(layer: ResolvedLayer, model: StaticModel): RangesPrep {
    const src = layerSource(model, layer);
    const line = layer.params.line;
    const arrow = layer.params.arrow;
    return {
      ...(typeof line === 'string' ? { line } : {}),
      arrow: arrow === 'end' ? 1 : arrow === 'start' ? 2 : arrow === 'both' ? 3 : 0,
      connector: layer.params.connector !== 'none',
      cyclic: model.systems.get(layer.data?.system ?? '')?.cyclic ?? false,
      ...(itemExtents(model, layer) ? { extents: itemExtents(model, layer) } : {}),
      ...(src ? { src } : {}),
      owner: layer.data?.address ?? '',
      channels: channelForms(layer.theme),
    };
  },

  content: rangesContent,
  known: knownLanesVersion,

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: RangesLay): RangesLay {
    const lay = prev ?? newLay(ctx, prep);
    const { axis, box, config, tokens } = ctx;
    beginItems(lay, ctx, prep);
    lay.font = config.itemFont;
    lay.fontSize = fontSize(config.itemFont);
    lay.items = 0;
    lay.m = 0;
    lay.hn = 0;
    lay.overflow = 0;
    lay.laneCount = 0;
    startConnectors(lay, ctx, prep.connector);
    const src = prep.src;
    if (!src) return lay;

    lay.style = lineStyle(prep.line, ctx);
    lay.arrow = prep.arrow;
    // the items in view and in the stretches its breaks hide (an item there is its break's tick)
    const densities = readItems(lay, ctx, data, src, ['category']);
    const step = lay.step;
    const rows = lay.rows;
    lay.items = rows.n;
    themeItems(lay, prep, ctx);
    for (let j = 0; j < rows.n; j++) lay.ids[j] = laneId(lay.idCache, rows.keys[j], rows.occ[j]);
    const n = lay.items;
    if (lay.idCache.size > 4 * n + 4096) lay.idCache.clear();
    if (lay.arrow) {
      if (lay.first.length < n) {
        lay.first = new Uint8Array(Math.max(n, 2 * lay.first.length));
        lay.last = new Uint8Array(lay.first.length);
      }
      markEnds(rows, prep.extents, prep.cyclic, src.length ?? Number.POSITIVE_INFINITY, lay.first, lay.last);
    }
    const maxLanes = laneGeometry(lay, ctx);
    let row = false;
    // at the density step a list of few rows (at most one a CSS pixel of the band) stays items: each its own mark, so
    // wide ones keep their labels, hits and arrows (a chromosome's bands at its fit); more rows are the density strip
    const dense = step === 'density' && n > box.width;
    // the gene line lays an item's rows in one lane, by its whole extent (`ranges-gene.ts`)
    const gene = lay.style === Gene;
    if (!dense && !densities.length) {
      lay.order = grown(lay.order, n);
      let count = n;
      if (gene) {
        // a label's room before each item: positions a CSS px on the data's system (the path's step, `ratio`)
        const perPx = ((axis.u1 - axis.u0) * src.ratio) / Math.max(1, box.width);
        const text = ctx.text;
        const font = config.itemFont;
        count = geneGroups(
          lay,
          n,
          lay.order,
          perPx,
          lay.laneHeight >= 10 && (ctx.labels === undefined || ctx.labels === 'inside')
            ? (l) => text.measure(font, l)
            : undefined,
        );
      } else for (let i = 0; i < n; i++) lay.order[i] = i;
      const order = lay.order.subarray(0, count);
      const starts: ArrayLike<number> = gene ? lay.gene.packStart : lay.starts;
      const ends: ArrayLike<number> = gene ? lay.gene.gEnd : lay.ends;
      sortArrayRange(order, 0, count, (o, a, b) => starts[o[a]] - starts[o[b]] || o[a] - o[b]);
      lay.lanes = grown(lay.lanes, n);
      const lanes = lay.lanes.subarray(0, n);
      const repOf = lay.gene.repOf;
      const assign = (cap: number, fresh: boolean) => {
        lay.laneCount = lay.sticky.assign(lay.ids, starts, ends, order, cap, fresh, lanes);
        // an item's other rows take its lane
        if (gene) for (let i = 0; i < n; i++) if (repOf[i] !== i) lanes[i] = lanes[repOf[i]];
        lay.overflow = 0;
        for (let k = 0; k < count; k++) if (lanes[order[k]] < 0) lay.overflow++;
      };
      assign(maxLanes, ctx.settled);
      // the whole data packs in its lanes: a sticky lane that leaves an item out is packed fresh before any overflow
      if (lay.overflow && ctx.content && !ctx.content.partial && !ctx.settled) assign(maxLanes, true);
      // past the lanes: the last lane becomes a density row of what does not fit; the rest stay items
      if (lay.overflow && config.ranges.overflow === 'extra-lanes' && maxLanes > 1) {
        assign(maxLanes - 1, ctx.settled);
        row = true;
      }
      // the whole data's lanes not known (a tiled source's): centred by the most lanes the layer has shown, every lane
      // when the last is the density row (else the lanes stay centred by the first frame's and run past the band)
      const shown = row ? maxLanes : lay.laneCount;
      if ((!ctx.content || ctx.content.partial) && lay.laneHeight === config.laneHeight && shown > lay.seen) {
        lay.seen = shown;
        lay.laneY = box.y + (box.height - Math.min(lay.seen, maxLanes) * lay.laneHeight) / 2;
      }
    }
    if (lay.scrolls) {
      // the lanes packed (as known, when the data is not whole): the band's top lane, kept inside them
      const total = Math.max(ctx.content?.lanes ?? 1, lay.laneCount);
      lay.scrollMax = Math.max(0, total - lay.room);
      lay.top = Math.max(0, Math.min(ctx.bandTop ?? 0, lay.scrollMax));
      lay.laneY = lay.laneY0 - lay.top * lay.laneHeight;
    }
    // the lanes out of the band are no marks (the draw lists, the hit index and the labels hold the lanes in view)
    if (lay.scrolls || ctx.content?.partial) {
      if (!dense && !densities.length) {
        // the lanes whose rows meet the box
        const first = Math.floor((box.y - lay.laneY) / lay.laneHeight);
        const past = Math.ceil((box.y + box.height - lay.laneY) / lay.laneHeight);
        const lanes = lay.lanes;
        for (let i = 0; i < n; i++) if (lanes[i] >= 0 && (lanes[i] < first || lanes[i] >= past)) lanes[i] = -1;
      }
    }
    // a source's lanes as known grow once the view settles (its track grows with them, `ranges-content.ts`)
    // its reachable region's, once read (the track keeps the height a pan will, `ranges-reach.ts`), and the view's
    if (!prep.extents) {
      noteReachable(prep, src, ctx, data);
      if (ctx.settled) noteLanes(prep, lay.laneCount);
    }
    lay.mode = dense || densities.length || (lay.overflow && !row) ? 'density' : step === 'bars' ? 'bars' : 'items';
    // below a position a pixel a run has no ends, no dash and no wave: a plain line or block
    if (lay.mode !== 'items' && lay.style >= 2 && lay.style <= 5) lay.style = Solid;
    // an arrowhead is an item's, by its own width (below): a run of bars keeps the arrows of the items wide enough
    if (lay.mode === 'density') lay.arrow = 0;
    if (lay.style === 4) growWave(2 * Math.ceil(box.width / WaveStep + 8));

    if (lay.mode === 'density') {
      lay.densityY = box.y;
      lay.densityH = box.height;
      densityStrip(lay, src, axis, densities, tokens, false, box.height < config.laneHeight ? ThinStripShare : 1);
      return lay;
    }
    if (row) {
      lay.densityY = lay.laneY + (maxLanes - 1) * lay.laneHeight;
      lay.densityH = lay.laneHeight;
      densityStrip(lay, src, axis, [], tokens, true, 1);
    }

    // marks: each item's pieces on the section, square where cut, round where the item really ends (on its own
    // system an item is its one piece: no mapping and no array an item)
    const px = 1 / axis.dpr;
    if (gene) geneMarks(lay, axis, src, n, box.x, box.x + box.width, nextMark);
    else
      for (let i = 0; i < n; i++) {
        const lane = lay.lanes[i];
        if (lane < 0) continue;
        // a span row (an item with no piece in its tile) is no mark: its item is drawn by its other rows
        if (lay.ends[i] <= lay.starts[i] && rows.itemStarts[i] >= 0 && rows.itemStarts[i] < rows.itemEnds[i]) continue;
        let pieces = itemPiece;
        if (src.path) pieces = toSection(src, lay.starts[i], lay.ends[i], mappedPieces);
        else {
          itemPiece[0] = lay.starts[i];
          itemPiece[1] = lay.ends[i];
        }
        let lo = Number.POSITIVE_INFINITY;
        let hi = Number.NEGATIVE_INFINITY;
        for (let k = 0; k + 1 < pieces.length; k += 2) {
          lo = Math.min(lo, pieces[k]);
          hi = Math.max(hi, pieces[k + 1]);
        }
        const marks = lay.m;
        // the arrowheads go on the piece holding the end (the start) the item really has, in reading order
        const atEnd = (lay.arrow & 1) !== 0 && lay.last[i] === 1;
        const atStart = (lay.arrow & 2) !== 0 && lay.first[i] === 1;
        const images = (rowImages ??= new RowImages());
        if (atEnd || atStart) images.set(src, lay.starts[i], lay.ends[i]);
        for (let k = 0; k + 1 < pieces.length; k += 2) {
          axis.piecesIn(pieces[k], pieces[k + 1]);
          const ka = axis.ka;
          const kb = axis.kb;
          for (let a = ka; a <= kb; a++) {
            if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
            // an end the item really has is round; on a reversed axis its start is on the right
            const atLo = pieces[k] === lo;
            const atHi = pieces[k + 1] === hi;
            let round = 0;
            if (!(axis.scut & Cut.left) && (axis.dir === -1 ? atHi : atLo)) round |= Round.left;
            if (!(axis.scut & Cut.right) && (axis.dir === -1 ? atLo : atHi)) round |= Round.right;
            const x0 = axis.snap(axis.sx0);
            const j = nextMark(lay);
            lay.mx0[j] = x0;
            lay.mx1[j] = Math.max(axis.snap(axis.sx1), x0 + px);
            lay.mlane[j] = lane;
            lay.mitem[j] = i;
            lay.mround[j] = round;
            if (lay.style > 0) lay.manchor[j] = anchorOf(axis, a, pieces[k], pieces[k + 1]);
            if ((atEnd || atStart) && lay.mx1[j] - x0 >= ArrowMinPx) {
              let tips = 0;
              const pa = Math.max(pieces[k], axis.pieces[2 * a]);
              const pb = Math.min(pieces[k + 1], axis.pieces[2 * a + 1]);
              const t = atEnd ? images.tip(axis, box, pa, pb, true) : 0;
              const u = atStart ? images.tip(axis, box, pa, pb, false) : 0;
              if (t) tips |= t === 1 ? 2 : 1;
              if (u) tips |= u === 1 ? 2 : 1;
              lay.mtip[j] = tips;
            }
          }
        }
        if (lay.mode === 'items' && pieces.length > 2) connect(lay, axis, pieces, lane, i, box.x, box.x + box.width);
        // an item wholly in a hidden stretch: its break's tick in its lane
        if (lay.m === marks && pieces.length >= 2) tickOn(lay, axis, pieces[0], lane);
      }
    // the gene line's boxes are never merged: a box an exon's part, thick or thin, at every step but the density's
    const merged = lay.mode === 'bars' && !gene;
    if (merged) mergeRuns(lay, px);
    // a layout's label sink (rings: placed outside with leaders) takes every shown item's label, an item once, with its
    // hit; else each label inside its marks where it fits whole (a modified residue's code over its one position too),
    // else beside in the label lane (drawn in the labels pass)
    const sink = ctx.labels !== undefined && ctx.labels !== 'inside' ? ctx.labels : undefined;
    lay.ln = 0;
    if (sink) sinkLabels(lay, src, sink, box.y);
    // a mark cut at both ends is square
    lay.gap = lay.laneHeight >= 6 ? 2 : 0;
    const h = lay.laneHeight - lay.gap;
    for (let j = 0; j < lay.m; j++)
      if (lay.mround[j])
        lay.mradius[j] = squircleRadius(lay.mx1[j] - lay.mx0[j], h, lay.radius.ratio, lay.radius.max, lay.radius.min);
    // the hit marks: an item's own marks, also inside a run merged with others (_spec/README.md "Render", Lanes and Hits)
    if (merged) {
      const s = lay.merge;
      for (let k = 0; k < lay.hn; k++)
        s.radius[k] = s.round[k]
          ? squircleRadius(s.x1[k] - s.x0[k], h, lay.radius.ratio, lay.radius.max, lay.radius.min)
          : 0;
      lay.hx0 = s.x0;
      lay.hx1 = s.x1;
      lay.hlane = s.lane;
      lay.hitem = s.item;
      lay.hround = s.round;
      lay.hradius = s.radius;
      lay.hslack = px;
    } else {
      lay.hn = lay.m;
      lay.hx0 = lay.mx0;
      lay.hx1 = lay.mx1;
      lay.hlane = lay.mlane;
      lay.hitem = lay.mitem;
      lay.hround = lay.mround;
      lay.hradius = lay.mradius;
      lay.hslack = lay.mode === 'bars' ? px : 0;
    }
    lay.index.build(lay.hn, lay.hlane, lay.hx0, lay.hx1, lay.laneCount);
    // below a position a pixel too: a run's items keep their labels beside (their own marks are the hit marks)
    // the gene line's labels are text beside its models (`ranges-gene.ts`), never inside a box
    if (gene) {
      lay.gene.ln = 0;
      if (!sink) placeGeneLabels(lay, n, box.x, lay.fontSize);
    } else if (!sink && lay.laneHeight >= 10) placeLabels(lay, ctx);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(RangesDraws, lay, out, emphasis);
  },

  items: describeItems,

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, section: lay.section, system: lay.system, ...(position >= 0 ? { position } : {}) };
    // a label beside its mark is its item's
    if (lay.labelY >= 0 && y >= lay.labelY && y < lay.labelY + lay.laneHeight) {
      for (let k = 0; k < lay.ln; k++)
        if (x >= lay.lx[k] && x < lay.lx[k] + lay.lw[k]) {
          const i = lay.litem[k];
          return itemHit(base, lay, lay.keys[i], lay.starts[i], lay.ends[i], i);
        }
      return base;
    }
    // a gene model's label is its item's
    if (lay.gene.ln) {
      const i = geneLabelAt(lay, x, y);
      if (i >= 0) return itemHit(base, lay, lay.keys[i], lay.starts[i], lay.ends[i], i);
    }
    const lane = Math.floor((y - lay.laneY) / lay.laneHeight);
    // a hidden items' tick is hit as its break (within the break's slop), the count of what it hides with it
    const tick = tickHit(lay, x, lane);
    if (tick) return tick;
    if (lay.mode === 'density' || !lay.hn) return base;
    const j = lay.index.near(lane, x, lay.hslack);
    if (j < 0) {
      // between an item's pieces: its connector is the item's
      const c = connectorAt(lay, lane, x);
      if (c < 0) return base;
      const i = lay.citem[c];
      return itemHit(base, lay, lay.keys[i], lay.starts[i], lay.ends[i], i);
    }
    const item = lay.hitem[j];
    if (item < 0) return base;
    const hit = itemHit(base, lay, lay.keys[item], lay.starts[item], lay.ends[item], item);
    const mixed = lay.hn > lay.m ? categoriesAt(lay, lane, x) : undefined;
    return mixed ? { ...hit, categories: mixed } : hit;
  },

  filtered: (lay) => lay.rows.filtered,
  used: itemsUsed,

  legend(lay, options) {
    // the density strip is a count, not a theme: its key is the strip itself
    if (lay.mode === 'density') return [];
    const painted = new Set<string>();
    if (lay.color?.legend?.kind === 'palette')
      for (let j = 0; j < lay.items; j++) {
        const v = lay.painted[j];
        if (v !== undefined) painted.add(v);
      }
    return channelLegends([lay.color, lay.size], { ...options, layer: lay.key, painted });
  },
};

/**
 * The categories of the items under the pointer within the merge slack, in the legend's order, when they are several
 * (a merged run of mixed categories): their names joined by a comma, else undefined. A hover's work, not a frame's.
 */
function categoriesAt(lay: RangesLay, lane: number, x: number): string | undefined {
  const { index, hx0, hx1, hitem } = lay;
  if (lane < 0 || lane + 1 >= index.laneStart.length || !lay.groupNames.length) return undefined;
  const slack = lay.hslack;
  const from = index.laneStart[lane];
  let lo = from;
  let hi = index.laneStart[lane + 1];
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (hx0[index.order[mid]] <= x + slack) lo = mid + 1;
    else hi = mid;
  }
  const groups: number[] = [];
  for (let k = lo - 1; k >= from && index.runMax[k] >= x - slack; k--) {
    const i = index.order[k];
    if (hx1[i] < x - slack || hitem[i] < 0) continue;
    const v = lay.painted[hitem[i]];
    const g = v === undefined ? -1 : lay.group[hitem[i]];
    if (g >= 0 && !groups.includes(g)) groups.push(g);
  }
  if (groups.length < 2) return undefined;
  groups.sort((a, b) => a - b);
  return groups.map((g) => lay.groupNames[g]).join(', ');
}

/**
 * The item of the run `[from, to)` (marks in the merge scratch) whose color group covers most of the run's width (each mark weighs its width), the lower group
 * (the legend's order) on a tie, the first such item. Counts in a kept scratch, zeroed after: allocates nothing warm.
 */
function majority(lay: RangesLay, item: Int32Array, from: number, to: number): number {
  const count = lay.merge.count;
  const { x0, x1 } = lay.merge;
  const group = lay.group;
  for (let q = from; q < to; q++) count[group[item[q]]] += x1[q] - x0[q];
  let best = item[from];
  let bestCount = -1;
  let bestGroup = 0;
  for (let q = from; q < to; q++) {
    const g = group[item[q]];
    const c = count[g];
    // widths are sums of floats: equal within a thousandth of a pixel is a tie
    if (c > bestCount + 1e-3 || (c >= bestCount - 1e-3 && g < bestGroup)) {
      best = item[q];
      bestCount = c;
      bestGroup = g;
    }
  }
  for (let q = from; q < to; q++) count[group[item[q]]] = 0;
  return best;
}

/**
 * Merges a lane's marks that touch within a device pixel into one bar a run: at most one mark a device column. The marks
 * are sorted by lane and left edge (ties in the order they were made) into the layout's kept scratch, then written back.
 */
function mergeRuns(lay: RangesLay, px: number) {
  const m = lay.m;
  const s = lay.merge;
  if (s.order.length < m) {
    const n = Math.max(m, 2 * s.order.length);
    s.order = new Int32Array(n);
    s.x0 = new Float32Array(n);
    s.x1 = new Float32Array(n);
    s.lane = new Int32Array(n);
    s.item = new Int32Array(n);
    s.round = new Uint8Array(n);
    s.radius = new Float64Array(n);
    s.tip = new Uint8Array(n);
  }
  if (s.count.length < lay.groupNames.length + lay.groupIds.size + 1)
    s.count = new Float64Array(2 * (lay.groupNames.length + lay.groupIds.size + 1));
  const order = s.order;
  for (let j = 0; j < m; j++) order[j] = j;
  const { mlane, mx0 } = lay;
  sortArrayRange(order, 0, m, (o, i, j) => mlane[o[i]] - mlane[o[j]] || mx0[o[i]] - mx0[o[j]] || o[i] - o[j]);
  const { x0, x1, lane, item, round, tip } = s;
  for (let k = 0; k < m; k++) {
    const j = order[k];
    x0[k] = lay.mx0[j];
    x1[k] = lay.mx1[j];
    lane[k] = lay.mlane[j];
    item[k] = lay.mitem[j];
    round[k] = lay.mround[j];
    tip[k] = lay.mtip[j];
  }
  // the item's own marks stay in the scratch: the hit marks (`hn`)
  lay.hn = m;
  lay.m = 0;
  for (let k = 0; k < m; ) {
    let end = x1[k];
    let it = item[k];
    // a run of several items whose fills are all equal keeps that fill (`−2 − item`, its first item's: drawn so, never
    // hit, sized or labelled as one item); a run of mixed fills takes the majority group's (ties by the legend's order)
    let same = it >= 0;
    let r = round[k];
    let e = k + 1;
    // an arrow is its item's alone: never merged into a run
    while (!tip[k] && e < m && !tip[e] && lane[e] === lane[k] && x0[e] <= end + px) {
      end = Math.max(end, x1[e]);
      if (item[e] !== it) it = -1;
      if (same && (item[e] < 0 || lay.fills[item[e]] !== lay.fills[item[k]])) same = false;
      r = (r & Round.left) | (round[e] & Round.right);
      e++;
    }
    if (it < 0 && same) it = -2 - item[k];
    else if (it < 0) it = -2 - majority(lay, item, k, e);
    const j = nextMark(lay);
    lay.mx0[j] = x0[k];
    lay.mx1[j] = end;
    lay.mlane[j] = lane[k];
    lay.mitem[j] = it;
    lay.mround[j] = r;
    lay.mtip[j] = tip[k];
    k = e;
  }
}
