/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArrayRange } from '@molstar/core/data/util/sort';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { type AxisTransform, Cut } from './axis.js';
import { drawPass, grown, NoColor, type PassName, type PassWriter, passDraws, Round } from './draw-lists.js';
import {
  beginItems,
  countStrip,
  drawTicks,
  type ItemExtents,
  type ItemsBase,
  type ItemsPrep,
  itemExtents,
  itemHit,
  markEnds,
  newItemsBase,
  readItems,
  RowImages,
  tickHit,
  tickOn,
  itemsUsed,
} from './items.js';
import { LaneIndex } from './lanes.js';
import { type LayerSource, layerSource, toSection } from './layer-data.js';
import type { Emphasis, Hit, LayerType, LayoutContext, Overlay, StoreReader } from './layer-type.js';
import { ArrowBody, arrowArgs, drawFailed, drawSpans } from './marks.js';
import { squircleRadius } from './squircle.js';
import type { SeqTokenName } from './tokens.js';
import type { Channel, ThemeNode } from '@molstar/seqx-core/theme/index';
import { channelForms, channelLegends, type LayerChannels, RowDatum } from './channels.js';

// The `secondary` layer type (_spec/README.md "Render", Layer types): helices and strands as blocks in one band, a turn
// a thin bar. The strand's arrowhead is on the piece holding its C-terminal end (the end of the item's last row in
// reading order: across the origin of a cyclic system the row that starts at 0), pointing the way the residues run on
// screen (the axis's direction, the layer's path's), and never where the view clips that end. An end cut by the view, a
// piece of the axis or a piece of the item's image is square. Items equal in extent and kind (the copies of a chain)
// are drawn once; a hit is the first of them.

/** A kind's class: the `kind` field's value (DSSP codes or words) as helix, strand or turn. */
export const Helix = 0;
export const Strand = 1;
export const Turn = 2;
type Kind = typeof Helix | typeof Strand | typeof Turn;

/** The colors of the classes: the `secondary` scheme's tokens, drawn directly when the layer has no color of its own. */
const KindTokens: readonly SeqTokenName[] = ['seq-ss-helix', 'seq-ss-strand', 'seq-ss-turn'];
/** The classes as the `kind` field reads to a theme (_spec/README.md "Themes"): whatever the data's spelling. */
const KindWords = ['helix', 'strand', 'turn'] as const;
/** The color's default: the `secondary` scheme, so its legend keys the classes. */
const SecondaryScheme: ThemeNode = { kind: 'scheme', params: { name: 'secondary' } };

/** A row's datum: its table's fields, `kind` read as its class's word (a DSSP code, an mmCIF one or a word alike). */
class KindDatum extends RowDatum {
  kind = 'turn';
  override field(name: string): string | number | undefined {
    return name === 'kind' ? this.kind : super.field(name);
  }
}

// whole codes and words: DSSP (H G I, E B, T S), mmCIF (HELX_*, STRN, TURN_*), helix, strand, turn
const HelixKind = /^(h|g|i|helix|helx(_\w+)?)$/i;
const StrandKind = /^(e|b|strand|strn|sheet)$/i;

/** The class of a `kind` value; anything else (a turn, a coil, no value) is a turn. */
export function secondaryKind(v: string | undefined): Kind {
  if (v === undefined) return Turn;
  return HelixKind.test(v) ? Helix : StrandKind.test(v) ? Strand : Turn;
}

export interface SecondaryPrep extends ItemsPrep {
  readonly cyclic: boolean;
  /** The rows of every item when the data is inline (an item's end is known), else undefined: no head. */
  readonly extents?: ItemExtents;
  /** Its `color` channel; `themed`: the layer has a color of its own (else the class tokens, directly). */
  readonly channels: LayerChannels;
  readonly themed: boolean;
}

export interface SecondaryLay extends ItemsBase {
  /** Per row: its class, whether it is the last in reading order of its item, whether it is drawn (not a copy). */
  kind: Uint8Array;
  first: Uint8Array;
  last: Uint8Array;
  shown: Uint8Array;
  fills: Uint32Array;
  /** The color channel read, and the datum it reads each row through. */
  color?: Channel;
  datum: KindDatum;
  /** Marks: extents, row, round ends, the arrowhead's side (0 none, 1 right, 2 left). */
  m: number;
  mx0: Float32Array;
  mx1: Float32Array;
  mitem: Int32Array;
  mround: Uint8Array;
  mtip: Uint8Array;
  mradius: Float64Array;
  mhead: Float32Array;
  /** The block's height and top, px. */
  bandY: number;
  bandH: number;
  index: LaneIndex;
  lane: Int32Array;
  order: Int32Array;
}

function nextMark(lay: SecondaryLay): number {
  const j = lay.m++;
  if (j >= lay.mx0.length) {
    lay.mx0 = grown(lay.mx0, j + 1);
    lay.mx1 = grown(lay.mx1, j + 1);
    lay.mitem = grown(lay.mitem, j + 1);
    lay.mround = grown(lay.mround, j + 1);
    lay.mtip = grown(lay.mtip, j + 1);
    lay.mradius = grown(lay.mradius, j + 1);
    lay.mhead = grown(lay.mhead, j + 1);
  }
  lay.mtip[j] = 0;
  lay.mradius[j] = 0;
  lay.mhead[j] = 0;
  return j;
}

/** A turn is this share of the band's height. */
const TurnShare = 0.25;
/** A head is this share of the band's height long. */
const HeadShare = 0.9;

const itemPiece: number[] = [0, 0];
const mapped: number[] = [];
const images = new RowImages();

export const SecondaryType: LayerType<SecondaryPrep, SecondaryLay> = {
  name: 'secondary',
  fallback: 'ranges',

  prepare(layer: ResolvedLayer, model: StaticModel): SecondaryPrep {
    const src = layerSource(model, layer);
    return {
      ...(src ? { src } : {}),
      owner: layer.data?.address ?? '',
      cyclic: model.systems.get(layer.data?.system ?? '')?.cyclic ?? false,
      ...(itemExtents(model, layer) ? { extents: itemExtents(model, layer) } : {}),
      channels: channelForms(layer.theme),
      themed: layer.theme.color !== undefined,
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: SecondaryLay): SecondaryLay {
    const { axis, box, config, tokens } = ctx;
    const lay: SecondaryLay = prev ?? {
      ...newItemsBase(ctx, prep),
      kind: new Uint8Array(64),
      first: new Uint8Array(64),
      last: new Uint8Array(64),
      shown: new Uint8Array(64),
      fills: new Uint32Array(64),
      datum: new KindDatum(),
      m: 0,
      mx0: new Float32Array(64),
      mx1: new Float32Array(64),
      mitem: new Int32Array(64),
      mround: new Uint8Array(64),
      mtip: new Uint8Array(64),
      mradius: new Float64Array(64),
      mhead: new Float32Array(64),
      bandY: 0,
      bandH: 0,
      index: new LaneIndex(),
      lane: new Int32Array(64),
      order: new Int32Array(64),
    };
    beginItems(lay, ctx, prep);
    lay.m = 0;
    const src = prep.src;
    if (!src) return lay;
    const densities = readItems(lay, ctx, data, src, ['kind', 'category']);
    const rows = lay.rows;
    const n = rows.n;
    if (lay.kind.length < n) {
      const size = Math.max(n, 2 * lay.kind.length);
      lay.kind = new Uint8Array(size);
      lay.first = new Uint8Array(size);
      lay.last = new Uint8Array(size);
      lay.shown = new Uint8Array(size);
      lay.fills = new Uint32Array(size);
    }
    if (densities.length) {
      lay.density = countStrip(lay.strip, src, axis, rows, densities, tokens);
      lay.bandY = box.y;
      lay.bandH = box.height;
      lay.index.build(0, lay.lane, lay.mx0, lay.mx1, 1);
      return lay;
    }
    classify(lay, prep, ctx);
    marks(lay, src, axis, config.squircleRatio, config.squircleMax, config.squircleMin);
    lay.lane = grown(lay.lane, lay.m);
    lay.lane.fill(0, 0, lay.m);
    lay.index.build(lay.m, lay.lane, lay.mx0, lay.mx1, 1);
    return lay;
  },

  draw(lay, _pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(SecondaryDraws, lay, out, emphasis);
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, section: lay.section, system: lay.system, ...(position >= 0 ? { position } : {}) };
    const tick = tickHit(lay, x, 0);
    if (tick) return tick;
    if (lay.density || !lay.m) return base;
    const j = lay.index.near(0, x, 0);
    if (j < 0) return base;
    const i = lay.mitem[j];
    return itemHit(base, lay, lay.rows.keys[i], lay.rows.starts[i], lay.rows.ends[i]);
  },

  filtered: (lay) => lay.rows.filtered,
  used: itemsUsed,

  legend(lay, options) {
    if (lay.density) return [];
    const painted = new Set<string>();
    for (let i = 0; i < lay.rows.n; i++) painted.add(KindWords[lay.kind[i]]);
    return channelLegends([lay.color], { ...options, layer: lay.key, painted });
  },
};

/**
 * Each row's class, fill, whether it is drawn once (the copies of a chain are rows of one extent and class: only the
 * first is) and whether it is the last of its item in reading order (the one holding the C-terminal end).
 */
function classify(lay: SecondaryLay, prep: SecondaryPrep, ctx: LayoutContext) {
  const { rows, tokens } = lay;
  const n = rows.n;
  const datum = lay.datum;
  datum.bind(n ? rows.fields[0] : undefined);
  const color = prep.channels.get('color', SecondaryScheme, ctx);
  lay.color = color;
  const fixed = typeof color.constant === 'number' ? color.constant : undefined;
  const missing = tokens.colors['seq-missing'];
  for (let i = 0; i < n; i++) {
    const k = (lay.kind[i] = secondaryKind(rows.values[i]));
    if (fixed !== undefined) lay.fills[i] = fixed;
    else if (!prep.themed) lay.fills[i] = tokens.colors[KindTokens[k]];
    else {
      datum.bind(rows.fields[i]);
      datum.row = rows.row[i];
      datum.start = rows.starts[i];
      datum.end = rows.ends[i];
      datum.kind = KindWords[k];
      const c = color.at(datum);
      lay.fills[i] = typeof c === 'number' ? c : missing;
    }
  }
  markEnds(rows, prep.extents, prep.cyclic, lay.src?.length ?? Number.POSITIVE_INFINITY, lay.first, lay.last);
  // copies: one row of each extent and class is drawn (the first), the order by extent then class
  lay.order = grown(lay.order, n);
  const order = lay.order.subarray(0, n);
  for (let i = 0; i < n; i++) order[i] = i;
  const { starts, ends } = rows;
  const kinds = lay.kind;
  sortArrayRange(
    order,
    0,
    n,
    (o, a, b) => starts[o[a]] - starts[o[b]] || ends[o[a]] - ends[o[b]] || kinds[o[a]] - kinds[o[b]] || o[a] - o[b],
  );
  for (let k = 0; k < n; k++) {
    const i = order[k];
    const p = k ? order[k - 1] : -1;
    lay.shown[i] = p >= 0 && starts[i] === starts[p] && ends[i] === ends[p] && kinds[i] === kinds[p] ? 0 : 1;
  }
}

/** The marks: each shown row's pieces in view, with their round ends and the arrowhead of the strand holding the C end. */
function marks(lay: SecondaryLay, src: LayerSource, axis: AxisTransform, ratio: number, max: number, min: number) {
  const { rows, box } = lay;
  const px = 1 / axis.dpr;
  const gap = box.height >= 6 ? 1 : 0;
  lay.bandY = box.y + gap;
  lay.bandH = box.height - 2 * gap;
  const head = lay.bandH * HeadShare;
  for (let i = 0; i < rows.n; i++) {
    if (!lay.shown[i]) continue;
    const s = rows.starts[i];
    const e = rows.ends[i];
    let pieces = itemPiece;
    if (src.path) pieces = toSection(src, s, e, mapped);
    else {
      itemPiece[0] = s;
      itemPiece[1] = e;
    }
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    for (let k = 0; k + 1 < pieces.length; k += 2) {
      lo = Math.min(lo, pieces[k]);
      hi = Math.max(hi, pieces[k + 1]);
    }
    const wantsHead = lay.last[i] === 1 && lay.kind[i] === Strand;
    if (wantsHead) images.set(src, s, e);
    const before = lay.m;
    for (let k = 0; k + 1 < pieces.length; k += 2) {
      axis.piecesIn(pieces[k], pieces[k + 1]);
      const ka = axis.ka;
      const kb = axis.kb;
      for (let a = ka; a <= kb; a++) {
        if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
        const atLo = pieces[k] === lo;
        const atHi = pieces[k + 1] === hi;
        let round = 0;
        if (!(axis.scut & Cut.left) && (axis.dir === -1 ? atHi : atLo)) round |= Round.left;
        if (!(axis.scut & Cut.right) && (axis.dir === -1 ? atLo : atHi)) round |= Round.right;
        const x0 = axis.snap(axis.sx0);
        const j = nextMark(lay);
        lay.mx0[j] = x0;
        lay.mx1[j] = Math.max(axis.snap(axis.sx1), x0 + px);
        lay.mitem[j] = i;
        lay.mround[j] = round;
        // the piece holding the C end: the head's tip at its boundary, unless the view clips it
        if (wantsHead) {
          const tip = images.tip(
            axis,
            box,
            Math.max(pieces[k], axis.pieces[2 * a]),
            Math.min(pieces[k + 1], axis.pieces[2 * a + 1]),
            true,
          );
          if (tip && lay.mx1[j] - lay.mx0[j] > px) {
            lay.mtip[j] = tip;
            lay.mhead[j] = Math.min(head, lay.mx1[j] - lay.mx0[j]);
            round &= tip === 1 ? Round.left : Round.right;
            lay.mround[j] = round;
          }
        }
        // only a block with a round end has a radius; a mark cut at both ends is square
        if (lay.kind[i] === Helix && round)
          lay.mradius[j] = squircleRadius(lay.mx1[j] - lay.mx0[j], lay.bandH, ratio, max, min);
      }
    }
    // an item wholly in a hidden stretch: its break's tick
    if (lay.m === before && pieces.length >= 2) tickOn(lay, axis, pieces[0], 0);
  }
}

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: placeholders, the strip, the blocks (a helix full height, a strand's body and head, a turn thin). */
function drawBlocks(lay: SecondaryLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  if (lay.tn) drawTicks(lay, out);
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  const a = out.args;
  if (lay.density) {
    a[0] = box.x;
    a[1] = box.y;
    a[2] = box.width;
    a[3] = box.height;
    out.imageArgs(lay.density);
    return;
  }
  const turnH = lay.bandH * TurnShare;
  for (let j = 0; j < lay.m; j++) {
    const i = lay.mitem[j];
    const kind = lay.kind[i];
    a[0] = lay.mx0[j];
    a[2] = lay.mx1[j] - lay.mx0[j];
    a[4] = lay.fills[i];
    a[5] = NoColor;
    a[6] = 0;
    if (kind === Turn) {
      a[1] = lay.bandY + (lay.bandH - turnH) / 2;
      a[3] = turnH;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
    } else if (kind === Helix) {
      a[1] = lay.bandY;
      a[3] = lay.bandH;
      a[7] = lay.mradius[j];
      a[8] = lay.mround[j];
      out.rectArgs();
    } else if (lay.mtip[j]) {
      a[1] = lay.bandY;
      a[3] = lay.bandH;
      a[7] = lay.mhead[j];
      a[8] = lay.mtip[j] === 1 ? 1 : -1;
      arrowArgs(out);
    } else {
      // a strand's body where its head is not here (its N end, a piece, or its C end cut by the view)
      const inset = (lay.bandH * (1 - ArrowBody)) / 2;
      a[1] = lay.bandY + inset;
      a[3] = lay.bandH - 2 * inset;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
    }
  }
}

/** The hover pass: an outline around each of the hovered item's marks, outside it by the offset. */
function drawHover(lay: SecondaryLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || !h.element) return;
  const half = lay.hoverWidth / 2 + lay.hoverOffset;
  const a = out.args;
  for (let j = 0; j < lay.m; j++) {
    const i = lay.mitem[j];
    if (lay.rows.keys[i] !== h.element.key) continue;
    a[0] = lay.mx0[j] - half;
    a[1] = lay.bandY - half;
    a[2] = lay.mx1[j] - lay.mx0[j] + 2 * half;
    a[3] = lay.bandH + 2 * half;
    a[4] = NoColor;
    a[5] = lay.tokens.colors['seq-hover'];
    a[6] = lay.hoverWidth;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  }
}

const SecondaryDraws = passDraws<SecondaryLay, Emphasis | undefined>({ data: drawBlocks, hover: drawHover });
