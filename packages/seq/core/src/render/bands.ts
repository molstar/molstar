/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ThemeNode } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { Channel } from '@molstar/seqx-core/theme/index';
import { channelForms, channelLegends, type LayerChannels, RowDatum } from './channels.js';
import { drawPass, grown, NoColor, type PassName, type PassWriter, passDraws } from './draw-lists.js';
import {
  beginItems,
  countStrip,
  drawTicks,
  type ItemsBase,
  type ItemsPrep,
  itemHit,
  itemsUsed,
  newItemsBase,
  readItems,
  tickHit,
  tickOn,
} from './items.js';
import { layerSource, toSection } from './layer-data.js';
import type { Emphasis, Hit, LayerType, LayoutContext, Overlay, StoreReader } from './layer-type.js';
import { drawFailed, drawSpans } from './marks.js';

// The `bands` layer type (_spec/README.md "Render", Layer types): a chromosome's stained bands in one band, an
// ideogram. Each item a block colored by its `stain` field through the `giemsa` scheme (grays derived from the text
// ink, so dark flips them), the centromere's `acen` bands as two triangles meeting at it, a hairline where no band is.
// A reader without the type draws the bands as ranges.

/** The color's default: the `giemsa` scheme, so its legend keys the stains. */
const GiemsaScheme: ThemeNode = { kind: 'scheme', params: { name: 'giemsa' } };

/** A row's datum: its table's fields, `stain` read as written. */
class StainDatum extends RowDatum {
  stain: string | undefined = undefined;
  override field(name: string): string | number | undefined {
    return name === 'stain' ? this.stain : super.field(name);
  }
}

export interface BandsPrep extends ItemsPrep {
  readonly channels: LayerChannels;
}

export interface BandsLay extends ItemsBase {
  /** Per row: its fill, and whether it is a centromere band (`acen`) and which way its tip points (1 right, 2 left). */
  fills: Uint32Array;
  tip: Uint8Array;
  color?: Channel;
  datum: StainDatum;
  /** Marks: extents and row. */
  m: number;
  mx0: Float32Array;
  mx1: Float32Array;
  mitem: Int32Array;
}

/** A gneg band's outline: a band of no stain still shows its edges against the surface. */
const Outline = 0.5;

const mapped: number[] = [];
const tri = new Float64Array(6);

export const BandsType: LayerType<BandsPrep, BandsLay> = {
  name: 'bands',
  fallback: 'ranges',

  prepare(layer: ResolvedLayer, model: StaticModel): BandsPrep {
    const src = layerSource(model, layer);
    return { ...(src ? { src } : {}), owner: layer.data?.address ?? '', channels: channelForms(layer.theme) };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: BandsLay): BandsLay {
    const { axis, tokens } = ctx;
    const lay: BandsLay = prev ?? {
      ...newItemsBase(ctx, prep),
      fills: new Uint32Array(64),
      tip: new Uint8Array(64),
      datum: new StainDatum(),
      m: 0,
      mx0: new Float32Array(64),
      mx1: new Float32Array(64),
      mitem: new Int32Array(64),
    };
    beginItems(lay, ctx, prep);
    lay.m = 0;
    const src = prep.src;
    if (!src) return lay;
    const densities = readItems(lay, ctx, data, src, ['stain']);
    const rows = lay.rows;
    const n = rows.n;
    if (lay.fills.length < n) {
      lay.fills = new Uint32Array(Math.max(n, 2 * lay.fills.length));
      lay.tip = new Uint8Array(lay.fills.length);
    }
    if (densities.length) {
      lay.density = countStrip(lay.strip, src, axis, rows, densities, tokens);
      return lay;
    }
    const color = prep.channels.get('color', GiemsaScheme, ctx);
    lay.color = color;
    const fixed = typeof color.constant === 'number' ? color.constant : undefined;
    const missing = tokens.colors['seq-missing'];
    const datum = lay.datum;
    for (let i = 0; i < n; i++) {
      const stain = rows.values[i];
      datum.bind(rows.fields[i]);
      datum.row = rows.row[i];
      datum.start = rows.starts[i];
      datum.end = rows.ends[i];
      datum.stain = stain;
      const c = fixed ?? color.at(datum);
      lay.fills[i] = typeof c === 'number' ? c : missing;
      lay.tip[i] = 0;
    }
    // the centromere: the p arm's `acen` band (the one an `acen` band follows in positions) points to its end, the
    // q arm's to its start; on screen by the axis's direction
    for (let i = 0; i < n; i++) {
      if (rows.values[i] !== 'acen') continue;
      let next = false;
      for (let j = 0; j < n && !next; j++)
        next = j !== i && rows.values[j] === 'acen' && rows.starts[j] >= rows.ends[i];
      const towardEnd = next;
      lay.tip[i] = (towardEnd ? 1 : 2) ^ (axis.dir === -1 ? 3 : 0);
    }
    for (let i = 0; i < n; i++) {
      const s = rows.starts[i];
      const e = rows.ends[i];
      let pieces: readonly number[] = mapped;
      if (src.path) pieces = toSection(src, s, e, mapped);
      else {
        mapped.length = 2;
        mapped[0] = s;
        mapped[1] = e;
      }
      const before = lay.m;
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        axis.piecesIn(pieces[k], pieces[k + 1]);
        for (let a = axis.ka; a <= axis.kb; a++) {
          if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
          const j = lay.m++;
          if (j >= lay.mx0.length) {
            lay.mx0 = grown(lay.mx0, j + 1);
            lay.mx1 = grown(lay.mx1, j + 1);
            lay.mitem = grown(lay.mitem, j + 1);
          }
          const x0 = axis.snap(axis.sx0);
          lay.mx0[j] = x0;
          lay.mx1[j] = Math.max(axis.snap(axis.sx1), x0 + 1 / axis.dpr);
          lay.mitem[j] = i;
        }
      }
      if (lay.m === before && pieces.length >= 2) tickOn(lay, axis, pieces[0], 0);
    }
    return lay;
  },

  draw(lay, _pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(BandsDraws, lay, out, emphasis);
  },

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, section: lay.section, system: lay.system, ...(position >= 0 ? { position } : {}) };
    const tick = tickHit(lay, x, 0);
    if (tick) return tick;
    if (lay.density) return base;
    for (let j = 0; j < lay.m; j++) {
      if (x < lay.mx0[j] || x >= lay.mx1[j]) continue;
      const i = lay.mitem[j];
      const label = lay.rows.labels[i];
      const h = itemHit(base, lay, lay.rows.keys[i], lay.rows.starts[i], lay.rows.ends[i]);
      // the band's name, what a navigator's readout says
      return label !== undefined ? { ...h, state: label } : h;
    }
    return base;
  },

  filtered: (lay) => lay.rows.filtered,
  used: itemsUsed,

  legend(lay, options) {
    if (lay.density) return [];
    const painted = new Set<string>();
    // the bands with a mark in view, never the rows only read
    for (let j = 0; j < lay.m; j++) {
      const v = lay.rows.values[lay.mitem[j]];
      if (v !== undefined) painted.add(v);
    }
    return channelLegends([lay.color], { ...options, layer: lay.key, painted });
  },
};

/** The data pass: the hairline where no band is, then each band (a block, or a triangle for `acen`). */
function drawBands(lay: BandsLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  const a = out.args;
  const ink = tokens.colors['msui-text-muted'];
  const mid = box.y + box.height / 2;
  a[0] = box.x;
  a[1] = mid - 0.5;
  a[2] = box.width;
  a[3] = 1;
  a[4] = ink;
  a[5] = NoColor;
  a[6] = 0;
  a[7] = 0;
  a[8] = 0;
  out.rectArgs();
  if (lay.tn) drawTicks(lay, out);
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  if (lay.density) {
    a[0] = box.x;
    a[1] = box.y;
    a[2] = box.width;
    a[3] = box.height;
    out.imageArgs(lay.density);
    return;
  }
  for (let j = 0; j < lay.m; j++) {
    const i = lay.mitem[j];
    const x0 = lay.mx0[j];
    const x1 = lay.mx1[j];
    const tip = lay.tip[i];
    if (tip) {
      // a triangle: its base at the arm's side, its tip at the centromere
      const wide = tip === 1 ? x0 : x1;
      const narrow = tip === 1 ? x1 : x0;
      tri[0] = wide;
      tri[1] = box.y;
      tri[2] = narrow;
      tri[3] = mid;
      tri[4] = wide;
      tri[5] = box.y + box.height;
      a[4] = lay.fills[i];
      a[5] = NoColor;
      a[6] = 0;
      a[7] = 0;
      out.pathArgs(tri, 3, true);
      continue;
    }
    a[0] = x0;
    a[1] = box.y;
    a[2] = x1 - x0;
    a[3] = box.height;
    a[4] = lay.fills[i];
    a[5] = ink;
    a[6] = Outline;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  }
}

/**
 * The hover pass: an outline around the hovered band, a squircle by the band's size (every emphasis a squircle), its
 * radius plus the offset; square on a side at the view's edge (a band there may go on past it) and below the minimum.
 */
function drawHover(lay: BandsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || !h.element) return;
  const half = lay.hoverWidth / 2 + lay.hoverOffset;
  const a = out.args;
  const c = lay.config;
  const { box } = lay;
  const bh = box.height;
  for (let j = 0; j < lay.m; j++) {
    if (lay.rows.keys[lay.mitem[j]] !== h.element.key) continue;
    const x0 = lay.mx0[j];
    const x1 = lay.mx1[j];
    const w = x1 - x0;
    // `squircleRadius` inline (a number returned from a call is boxed)
    const r =
      w < c.squircleMin || bh < c.squircleMin ? 0 : Math.min(bh * c.squircleRatio, c.squircleMax, w / 2, bh / 2);
    const round = (x0 > box.x + 0.5 ? 1 : 0) | (x1 < box.x + box.width - 0.5 ? 2 : 0);
    a[0] = x0 - half;
    a[1] = box.y - half;
    a[2] = w + 2 * half;
    a[3] = bh + 2 * half;
    a[4] = NoColor;
    a[5] = lay.tokens.colors['seq-hover'];
    a[6] = lay.hoverWidth;
    a[7] = r > 0 && round !== 0 ? r + half : 0;
    a[8] = round;
    out.rectArgs();
  }
}

const BandsDraws = passDraws<BandsLay, Emphasis | undefined>({ data: drawBands, hover: drawHover });
