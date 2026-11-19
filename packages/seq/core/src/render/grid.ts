/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The grid layer (_spec/README.md "Render", The grid layer): a value per pair of positions drawn over a grid section's
// plot, x across (`ctx.axis`) and y down (`ctx.axisY`), its `color` channel at the value (by value over the grid's
// whole range, held). A cell of `grid.rectMinPx` and more is a rect with an outline; below, image tiles (256 × 256
// samples, never smoothed) from the shared reduction (`Pyramid`: square levels, factor 2, the layer's `aggregate`),
// built once from the held grid and kept by its version, at the level whose sample is at least a device pixel, in the
// shared image cache by (layer, level, tile, theme). While the view moves a frame builds at most `grid.buildBudget`
// samples and draws a coarser tile meanwhile, never empty; a settled frame builds the rest. A missing value is
// transparent: later layers over earlier merge. `triangle` keeps to the cells above (`upper`: y before x) or below the
// diagonal. The y path (`axes: 2`): the layer paints the emphasis on its y system (a row at each lit y, the crosshair,
// a selected rect), the strip paints x. A grid draws only on its own systems: on others it says so and draws nothing.

import { issueOf } from '@molstar/core/repr/tree';
import { theme as themeForm } from '@molstar/seqx-core/lang';
import type { Location, LocationPart } from '@molstar/seqx-core/query';
import type { GridData, StaticModel } from '@molstar/seqx-core/resolve';
import type { ChannelName, Channel } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, channelLegends, Extents, type LayerChannels, PlainDatum } from './channels.js';
import type { RenderConfig } from './config.js';
import { grown, type ImageData, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { walkStripes } from './highlight-stripes.js';
import { ImageCache, ImageTile } from './images.js';
import { RenderIssues } from './issues.js';
import {
  type Box,
  type Emphasis,
  type EmphasisHighlight,
  type Hit,
  type LayerType,
  type LayoutContext,
  sectionId,
} from './layer-type.js';
import { Pyramid, pyramidLevelFor } from './pyramid.js';
import type { SeqTokens } from './tokens.js';

/** The color channel's environment: its one field is the value. */
const ValueEnv = { numeric: (n: string) => n === 'value' };

/** The default color: by value, a continuous palette over the grid's range. */
const DefaultColor = themeForm.field('value', { list: false });

/** The crosshair's tint (the hover's ink at this opacity) over a lit row and column. */
const CrossTint = 0.14;
/** The crosshair's edge: the hover's ink at this opacity, a hairline each side of its column and row. */
const CrossEdge = 0.85;

type Triangle = 'both' | 'upper' | 'lower';

export interface GridPrep {
  readonly owner: string;
  readonly grid?: GridData;
  readonly xSystem: string;
  readonly ySystem: string;
  readonly triangle: Triangle;
  readonly channels: LayerChannels;
  /** The values' range (NaN with none): the color's domain is held to it. */
  readonly lo: number;
  readonly hi: number;
  /** The shared reduction over the values, built once (`rows` the x positions, `columns` the y). */
  readonly pyramid?: Pyramid;
}

type ScratchHit = { -readonly [K in keyof Hit]: Hit[K] };

export class GridLay {
  key = '';
  section = '';
  sectionAddress = '';
  system = '';
  ySystem = '';
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  axis: AxisTransform | undefined;
  axisY: AxisTransform | undefined;
  tokens: SeqTokens | undefined;
  config: RenderConfig | undefined;
  triangle: Triangle = 'both';
  /**
   * A strip of a grid over blocks on another pair of systems: nothing of the layer is drawn there, its emphasis neither
   * (the strip of its own pair draws it once).
   */
  off = false;
  grid: GridData | undefined;
  color: Channel | undefined;
  /** Moves when the color changes (its channel or its held domain): part of an image tile's key. */
  colorVersion = 0;
  readonly datum = new PlainDatum();
  readonly extents = new Extents();
  /** Cells as rects (`rectMinPx` and more): x, y, w, h each, and the fill. */
  n = 0;
  rx = new Float64Array(64);
  rc = new Uint32Array(16);
  /** Image tiles: their boxes (x, y, w, h each) and images. */
  readonly images: ImageData[] = [];
  readonly ib: number[] = [];
  cache: ImageCache | undefined;
  /** The level the tiles were drawn at (−1: rects). */
  level = -1;
  readonly scratch: ScratchHit = { layer: '', section: '', system: '' };
}

const OwnImages = new WeakMap<RenderConfig, ImageCache>();

function prepare(layer: ResolvedLayer, model: StaticModel): GridPrep {
  const owner = layer.data?.address ?? '';
  const a = model.annotations.get(owner);
  const t = layer.params.triangle;
  const triangle: Triangle = t === 'upper' || t === 'lower' ? t : 'both';
  const agg = layer.params.aggregate;
  const reducer = agg === 'max' || agg === 'min' ? agg : 'mean';
  const grid = a?.shape === 'grid' ? a.grid : undefined;
  let lo = Number.NaN;
  let hi = Number.NaN;
  let pyramid: Pyramid | undefined;
  if (grid) {
    const v = grid.values;
    for (let i = 0; i < v.length; i++) {
      const x = v[i];
      if (Number.isNaN(x)) continue;
      if (Number.isNaN(lo) || x < lo) lo = x;
      if (Number.isNaN(hi) || x > hi) hi = x;
    }
    // a triangle reduces its own cells alone: the others missing, so a coarse sample on the diagonal keeps its half
    let kept = v;
    if (triangle !== 'both') {
      kept = Float32Array.from(v);
      for (let i = 0; i < grid.nx; i++)
        for (let j = 0; j < grid.ny; j++) if (!keeps(triangle, i, j)) kept[i * grid.ny + j] = Number.NaN;
    }
    pyramid = new Pyramid({ factor: 2, pairing: 'square', reducer });
    pyramid.setValues(grid.nx, grid.ny, kept);
  }
  return {
    owner,
    ...(grid ? { grid } : {}),
    xSystem: a?.target ?? '',
    ySystem: a?.ySystem ?? '',
    triangle,
    channels: channelForms(layer.theme),
    lo,
    hi,
    ...(pyramid ? { pyramid } : {}),
  };
}

function say(ctx: LayoutContext, code: 'seq.render.layer.axes' | 'seq.render.grid.systems', message: string) {
  ctx.report(issueOf(RenderIssues, code, message, { address: ctx.key, hint: RenderIssues.codes[code].hint }));
}

/** Whether cell (i, j) is kept by the triangle (`upper`: y before x, `lower`: y after x). */
function keeps(t: Triangle, i: number, j: number): boolean {
  return t === 'both' || (t === 'upper' ? j < i : j > i);
}

export const GridType: LayerType<GridPrep, GridLay> = {
  name: 'grid',
  axes: 2,
  prepare,

  layout(prep, ctx: LayoutContext, _data, _overlay, prev?: GridLay): GridLay {
    const { axis, box, config, tokens } = ctx;
    const lay = prev ?? new GridLay();
    lay.key = ctx.key;
    lay.section = sectionId(ctx.section);
    lay.sectionAddress = ctx.section.address;
    lay.system = axis.system;
    lay.ySystem = ctx.axisY?.system ?? prep.ySystem;
    lay.box = box;
    lay.axis = axis;
    lay.axisY = ctx.axisY;
    lay.tokens = tokens;
    lay.config = config;
    lay.triangle = prep.triangle;
    lay.n = 0;
    lay.off = false;
    lay.images.length = 0;
    lay.ib.length = 0;
    lay.level = -1;
    lay.grid = undefined;
    if (ctx.images) lay.cache = ctx.images;
    else if (!lay.cache) {
      let own = OwnImages.get(config);
      if (!own) OwnImages.set(config, (own = new ImageCache(config.imageBytes)));
      lay.cache = own;
    }
    const ay = ctx.axisY;
    if (!ay) {
      say(ctx, 'seq.render.layer.axes', 'A grid layer draws in a grid section only');
      return lay;
    }
    if (!prep.grid) return lay;
    if (prep.xSystem !== axis.system || prep.ySystem !== ay.system) {
      // a grid over blocks is a strip a pair of blocks: a grid draws on the pairs of its own systems only, silently
      if (!ctx.section.blocks)
        say(ctx, 'seq.render.grid.systems', 'The grid is not on its section’s two systems: it is not drawn');
      else lay.off = true;
      return lay;
    }
    lay.grid = prep.grid;
    const color = prep.channels.get('color' as ChannelName, DefaultColor, ctx, ValueEnv);
    if (color !== lay.color) {
      lay.color = color;
      lay.colorVersion++;
    }
    // the color's open domain is the grid's whole range, held: a tile's colors never shimmer with the view
    lay.extents.reset([color]);
    if (lay.extents.any && !Number.isNaN(prep.lo)) {
      lay.datum.value = prep.lo;
      lay.extents.add(lay.datum);
      lay.datum.value = prep.hi;
      lay.extents.add(lay.datum);
      lay.extents.flush(ctx.theme);
    }
    if (Math.min(axis.scale, ay.scale) >= config.grid.rectMinPx) layoutRects(lay, prep, axis, ay);
    else layoutImages(lay, prep, ctx, axis, ay);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    if (!lay.axis || !lay.axisY || !lay.tokens || !lay.config || lay.off) return;
    switch (pass) {
      case 'data':
        drawData(lay, out);
        break;
      case 'highlight':
      case 'selection':
      case 'hover':
      case 'cursor':
        if (emphasis) drawEmphasis(lay, pass, out, emphasis);
        break;
    }
  },

  hit(lay, at) {
    const a = lay.axis;
    const ay = lay.axisY;
    if (!a || !ay || !lay.grid) return undefined;
    const i = a.positionOf(at, 0);
    const j = ay.positionOf(at, 1);
    // the cell under the point, drawn or not (a triangle's other half, a missing value): a plot's every point is a cell
    if (i < 0 || j < 0) return undefined;
    const h = lay.scratch;
    h.layer = lay.key;
    h.section = lay.section;
    h.system = lay.system;
    h.position = i;
    h.ySystem = lay.ySystem;
    h.yPosition = j;
    h.target = 'cell';
    return h;
  },

  legend(lay, options) {
    return channelLegends([lay.color], { ...options, layer: lay.key });
  },
};

/** The color of a value (missing: none). */
function colorOf(lay: GridLay, v: number): number {
  if (Number.isNaN(v) || !lay.color) return NoColor;
  lay.datum.value = v;
  const c = lay.color.at(lay.datum);
  return typeof c === 'number' ? c : NoColor;
}

/** Every cell in view as a rect (few: each at least `rectMinPx` on both axes). */
function layoutRects(lay: GridLay, prep: GridPrep, a: AxisTransform, ay: AxisTransform) {
  const g = prep.grid!;
  const i0 = Math.floor(a.u0);
  const i1 = Math.ceil(a.u1);
  const j0 = Math.floor(ay.u0);
  const j1 = Math.ceil(ay.u1);
  for (let ui = i0; ui < i1; ui++) {
    const i = a.elementAtUnit(ui + 0.5);
    if (i < 0 || i >= g.nx || !a.cell(i)) continue;
    const x0 = Math.min(a.sx0, a.sx1);
    const w = Math.abs(a.sx1 - a.sx0);
    for (let uj = j0; uj < j1; uj++) {
      const j = ay.elementAtUnit(uj + 0.5);
      if (j < 0 || j >= g.ny || !keeps(prep.triangle, i, j) || !ay.cell(j)) continue;
      const c = colorOf(lay, g.values[i * g.ny + j]);
      if (c === NoColor) continue;
      const k = lay.n++;
      lay.rx = grown(lay.rx, 4 * k + 4);
      lay.rc = grown(lay.rc, k + 1);
      lay.rx[4 * k] = x0;
      lay.rx[4 * k + 1] = Math.min(ay.sx0, ay.sx1);
      lay.rx[4 * k + 2] = w;
      lay.rx[4 * k + 3] = Math.abs(ay.sx1 - ay.sx0);
      lay.rc[k] = c;
    }
  }
}

/** A position's sample range on an axis's view: the first and past the last position in view. */
function viewRange(a: AxisTransform, n: number, out: Int32Array, o: number) {
  const p0 = a.elementAtUnit(Math.max(a.u0, 0));
  const p1 = a.elementAtUnit(Math.max(a.u0, a.u1 - 1e-6));
  const lo = Math.max(0, Math.min(p0 < 0 ? 0 : p0, p1 < 0 ? n - 1 : p1));
  const hi = Math.min(n, Math.max(p0 < 0 ? 0 : p0, p1 < 0 ? n - 1 : p1) + 1);
  out[o] = lo;
  out[o + 1] = hi;
}

const Range = new Int32Array(4);

/**
 * The image tiles in view at the level whose sample is at least a device pixel on both axes, each from the cache or
 * built into it; past the frame's budget while the view moves, the coarser level's tile there (built if none is held).
 */
function layoutImages(lay: GridLay, prep: GridPrep, ctx: LayoutContext, a: AxisTransform, ay: AxisTransform) {
  const g = prep.grid!;
  const pyr = prep.pyramid!;
  const cache = lay.cache!;
  const tokens = ctx.tokens;
  const want = pyramidLevelFor(Math.max(1 / a.devicePx, 1 / ay.devicePx), 2);
  const top = pyr.levels('rows') > pyr.levels('columns') ? pyr.levels('rows') : pyr.levels('columns');
  const L = Math.min(want, top - 1);
  lay.level = L;
  viewRange(a, g.nx, Range, 0);
  viewRange(ay, g.ny, Range, 2);
  let budget = ctx.settled ? Number.POSITIVE_INFINITY : ctx.config.grid.buildBudget;
  let lastCoarse = '';
  const f = 2 ** L;
  const side = ImageTile * f;
  for (let ti = Math.floor(Range[0] / side); ti * side < Range[1]; ti++)
    for (let tj = Math.floor(Range[2] / side); tj * side < Range[3]; tj++) {
      let level = L;
      let key = tileKey(lay, prep, tokens.version, level, ti, tj);
      let img = cache.get(key);
      let tx = ti;
      let ty = tj;
      if (!img) {
        const cost = tileCost(pyr, level, ti, tj);
        if (cost > budget) {
          // a coarser level's tile over this one meanwhile: the first held, else the first the budget builds, else the
          // coarsest (one sample's worth: never empty); drawn once
          for (let c = L + 1; c < top; c++) {
            const k = 2 ** (c - L);
            const ci = Math.floor(ti / k);
            const cj = Math.floor(tj / k);
            const ck = tileKey(lay, prep, tokens.version, c, ci, cj);
            const held = cache.get(ck);
            const fits = !held && tileCost(pyr, c, ci, cj) <= budget;
            if (held || fits || c === top - 1) {
              if (fits) budget -= tileCost(pyr, c, ci, cj);
              level = c;
              tx = ci;
              ty = cj;
              key = ck;
              img = held;
              break;
            }
          }
          if (key === lastCoarse) continue;
          lastCoarse = key;
        } else budget -= cost;
      }
      if (!img) {
        img = buildTile(lay, prep, cache, level, tx, ty);
        if (!img) continue;
        cache.put(key, img);
      }
      placeTile(lay, g, a, ay, level, tx, ty, img);
    }
}

/** The samples a tile of a level holds. */
function tileCost(pyr: Pyramid, level: number, ti: number, tj: number): number {
  const lv = pyr.level(level);
  return (
    Math.max(0, Math.min(ImageTile, lv.rows - ti * ImageTile)) *
    Math.max(0, Math.min(ImageTile, lv.columns - tj * ImageTile))
  );
}

function tileKey(lay: GridLay, prep: GridPrep, tokens: number, level: number, ti: number, tj: number): string {
  return `${lay.key}|${level}|${ti}.${tj}|${prep.pyramid!.version}|${lay.colorVersion}.${tokens}|${prep.triangle}`;
}

/**
 * A tile's image: x samples across, y samples down, in display order (a reversed axis's samples right to left); a
 * sample outside the triangle or missing is transparent.
 */
function buildTile(
  lay: GridLay,
  prep: GridPrep,
  cache: ImageCache,
  level: number,
  ti: number,
  tj: number,
): ImageData | undefined {
  const lv = prep.pyramid!.level(level);
  const w = Math.min(ImageTile, lv.rows - ti * ImageTile);
  const h = Math.min(ImageTile, lv.columns - tj * ImageTile);
  if (w <= 0 || h <= 0) return undefined;
  const img = cache.slot(w, h);
  const d = img.data;
  const flipX = lay.axis!.dir === -1;
  const flipY = lay.axisY!.dir === -1;
  for (let r = 0; r < h; r++) {
    const sy = tj * ImageTile + r;
    const row = flipY ? h - 1 - r : r;
    for (let c = 0; c < w; c++) {
      const sx = ti * ImageTile + c;
      const q = 4 * (row * w + (flipX ? w - 1 - c : c));
      // the pyramid of a triangle holds its own cells alone (a sample of none of them is missing)
      const fill = colorOf(lay, lv.values[sx * lv.columns + sy]);
      if (fill === NoColor) {
        d[q] = d[q + 1] = d[q + 2] = d[q + 3] = 0;
        continue;
      }
      d[q] = fill >>> 24;
      d[q + 1] = (fill >>> 16) & 0xff;
      d[q + 2] = (fill >>> 8) & 0xff;
      d[q + 3] = fill & 0xff;
    }
  }
  return img;
}

/** A tile's box on the plot: the positions its samples cover on each axis, through the transforms. */
function placeTile(
  lay: GridLay,
  g: GridData,
  a: AxisTransform,
  ay: AxisTransform,
  level: number,
  ti: number,
  tj: number,
  img: ImageData,
) {
  const f = 2 ** level;
  const xs = Math.min(ti * ImageTile * f, g.nx);
  const xe = Math.min(xs + img.width * f, g.nx);
  const ys = Math.min(tj * ImageTile * f, g.ny);
  const ye = Math.min(ys + img.height * f, g.ny);
  // the edges on the device grid: two tiles meeting at a fractional pixel would each antialias it, a seam the page
  // shows through (a thin line where a tile ends)
  const d = a.dpr;
  const x0 = Math.round(a.x(xs) * d) / d;
  const x1 = Math.round(a.x(xe) * d) / d;
  const y0 = Math.round(ay.x(ys) * d) / d;
  const y1 = Math.round(ay.x(ye) * d) / d;
  if (Number.isNaN(x0) || Number.isNaN(x1) || Number.isNaN(y0) || Number.isNaN(y1)) return;
  lay.images.push(img);
  lay.ib.push(Math.min(x0, x1), Math.min(y0, y1), Math.abs(x1 - x0), Math.abs(y1 - y0));
}

function drawData(lay: GridLay, out: PassWriter) {
  const g = out.args;
  for (let k = 0; k < lay.images.length; k++) {
    g[0] = lay.ib[4 * k];
    g[1] = lay.ib[4 * k + 1];
    g[2] = lay.ib[4 * k + 2];
    g[3] = lay.ib[4 * k + 3];
    out.imageArgs(lay.images[k]);
  }
  // rects with a thin outline in the page's ink: cells read apart
  const line = lay.tokens!.colors['msui-bg'];
  const dpr = lay.axis!.dpr;
  Ink[1] = line;
  Ink[2] = 1 / dpr;
  for (let k = 0; k < lay.n; k++) {
    R[0] = lay.rx[4 * k];
    R[1] = lay.rx[4 * k + 1];
    R[2] = lay.rx[4 * k + 2];
    R[3] = lay.rx[4 * k + 3];
    Ink[0] = lay.rc[k];
    rect(out);
  }
}

/**
 * Whether a selection holds a rect over a section's two systems: such a section draws the rect itself, never its sides
 * projected (they light the sections of one of its systems).
 */
export function rectOver(loc: Location | undefined, x: string, y: string | undefined): boolean {
  if (!loc || y === undefined) return false;
  for (let i = 0; i < loc.rects.length; i++)
    if (loc.rects[i].a.system === x && loc.rects[i].b.system === y) return true;
  return false;
}

/** A Location's part on a system. */
function partOn(loc: Location | undefined, system: string): LocationPart | undefined {
  // an indexed loop: `partOf`'s `for...of` makes an iterator on this redraw path (nothing allocated once warm)
  if (!loc) return undefined;
  for (let i = 0; i < loc.parts.length; i++) if (loc.parts[i].system === system) return loc.parts[i];
  return undefined;
}

/**
 * What the emphasis paints with: the fill, the stroke and its width (`Ink[0..2]`), set before each paint. The colors
 * (above 2^30, no small integers) and the widths go through a typed array: a number passed to a call that is not
 * inlined is boxed, and an emphasis redraw allocates nothing once warm.
 */
const Ink = new Float64Array(3);

/**
 * Sets `Ink` from `Src` (a fill color and its opacity, a stroke color and its opacity, a width): a color at an opacity
 * (0 to 1), none at 0. The colors are read from an array, never passed (a color above 2^30 passed to a call is boxed).
 */
const Src = new Float64Array(5);

function ink() {
  const f = Src[0];
  const s = Src[2];
  Ink[0] = Src[1] > 0 ? ((f & 0xffffff00) | Math.round(Math.min(1, Src[1]) * (f & 0xff))) >>> 0 : NoColor;
  Ink[1] = Src[3] > 0 ? ((s & 0xffffff00) | Math.round(Math.min(1, Src[3]) * (s & 0xff))) >>> 0 : NoColor;
  Ink[2] = Src[4];
}

/** A range's extent on an axis, its first cell to its last, into `Box[o]` (low) and `Box[o + 1]` (high); false: hidden. */
const Box4 = new Float64Array(4);

function span(a: AxisTransform, s: number, e: number, o: number): boolean {
  if (!a.cell(s)) return false;
  let lo = a.sx0 < a.sx1 ? a.sx0 : a.sx1;
  let hi = a.sx0 < a.sx1 ? a.sx1 : a.sx0;
  if (!a.cell(e - 1)) return false;
  if (a.sx0 < lo) lo = a.sx0;
  if (a.sx1 < lo) lo = a.sx1;
  if (a.sx0 > hi) hi = a.sx0;
  if (a.sx1 > hi) hi = a.sx1;
  Box4[o] = lo;
  Box4[o + 1] = hi;
  return true;
}

/** A rect's x, y, width and height, set before `rect`. */
const R = new Float64Array(4);

/** The rect `R` holds in `Ink`, square. */
function rect(out: PassWriter) {
  const g = out.args;
  g[0] = R[0];
  g[1] = R[1];
  g[2] = R[2];
  g[3] = R[3];
  g[4] = Ink[0];
  g[5] = Ink[1];
  g[6] = Ink[2];
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * Rows across the plot at a part's ranges on y (each range's span down the plot), in `Ink`: the y half of the
 * crosshair, a selection's, a highlight's (_spec/README.md "Render", The grid layer).
 */
function paintRows(lay: GridLay, out: PassWriter, ranges: Int32Array) {
  const ay = lay.axisY!;
  const a = lay.axis!;
  const b = lay.box;
  // a stripe across the plot (a highlight's, where channels overlap): its edges on device pixels, shared with its
  // neighbors; a band is a bar along its leading edge
  stripeEdges(a);
  const w = RowX[1] - RowX[0];
  if (w <= 0) return;
  for (let r = 0; r + 1 < ranges.length; r += 2) {
    if (!span(ay, ranges[r], ranges[r + 1], 2)) continue;
    const top = Math.max(b.y, Box4[2]);
    const bottom = Math.min(b.y + b.height, Box4[3]);
    if (bottom <= top) continue;
    // across the plot: x's extent (the plot may be narrower than its band)
    R[0] = RowX[0];
    R[1] = top;
    R[2] = RowX[1] - RowX[0];
    R[3] = bottom - top;
    rect(out);
  }
}

/**
 * The x edges `[RowX[0], RowX[1]]` of stripe `slot` of `count` across the plot, on device pixels (shared with its
 * neighbors); a bar (`RowLook[2]` px) is along the leading edge. Through an array: a double returned or kept is boxed.
 */
const RowX = new Float64Array(2);

function stripeEdges(a: AxisTransform) {
  const slot = RowLook[0];
  const count = RowLook[1];
  const dpr = a.dpr;
  if (count > 1) {
    RowX[0] = Math.round((a.x0 + (slot * a.width) / count) * dpr) / dpr;
    RowX[1] = Math.round((a.x0 + ((slot + 1) * a.width) / count) * dpr) / dpr;
  } else {
    RowX[0] = a.x0;
    RowX[1] = a.x0 + a.width;
  }
  if (RowLook[2] > 0)
    RowX[1] = RowX[0] + Math.min(Math.max(Math.round(RowLook[2] * dpr) / dpr, 1 / dpr), RowX[1] - RowX[0]);
}

/** A stripe's look for `paintRows`: its slot, the stripe count, and a bar's width (0: the whole stripe). */
const RowLook = new Float64Array([0, 1, 0]);

/** The rows' highlight channels (`walkStripes`' context: set before it, since it is given no closure). */
const GParts: (LocationPart | undefined)[] = [];
const GIndex: number[] = [];
const GOutlines: number[] = [];
const GRun = new Int32Array(2);
let GLay!: GridLay;
let GOut!: PassWriter;
let GList!: readonly EmphasisHighlight[];

/**
 * The highlight channels' rows (the plot's y part; the x part is the strip's, drawn with the tracks'): where channels
 * overlap the row is split into equal stripes across the plot, in channel order (as in tracks, a band's stripes
 * perpendicular to its axis); `band` is a bar along a stripe's left edge; `outline` takes no stripe, one outline over
 * the whole row.
 */
function drawRowHighlights(lay: GridLay, out: PassWriter, hs: readonly EmphasisHighlight[]) {
  const config = lay.config!;
  let n = 0;
  let o = 0;
  for (let i = 0; i < hs.length; i++) {
    const part = partOn(hs[i].location, lay.ySystem);
    if (!part) continue;
    if (hs[i].style === 'outline') GOutlines[o++] = i;
    else {
      GParts[n] = part;
      GIndex[n] = i;
      n++;
    }
  }
  GLay = lay;
  GOut = out;
  GList = hs;
  if (n === 1) {
    rowsOf(lay, out, hs[GIndex[0]], GParts[0]!.ranges, 0, 1);
  } else if (n > 1) walkStripes(GParts, n, Math.max(1, config.emphasis.highlightStripes | 0), emitRows);
  for (let k = 0; k < o; k++) {
    const h = hs[GOutlines[k]];
    rowsOf(lay, out, h, partOn(h.location, lay.ySystem)!.ranges, 0, 1);
  }
  for (let c = 0; c < n; c++) GParts[c] = undefined;
}

function emitRows(channel: number, start: number, end: number, slot: number, count: number, _first: number) {
  GRun[0] = start;
  GRun[1] = end;
  rowsOf(GLay, GOut, GList[GIndex[channel]], GRun, slot, count);
}

/** A channel's ranges as rows in its style, in stripe `slot` of `count`. */
function rowsOf(lay: GridLay, out: PassWriter, h: EmphasisHighlight, ranges: Int32Array, slot: number, count: number) {
  const em = lay.config!.emphasis;
  const style = h.style ?? 'tint';
  Src[0] = h.color;
  Src[1] = style === 'tint' ? em.highlightTint : style === 'band' ? 1 : 0;
  Src[2] = h.color;
  Src[3] = style === 'outline' ? 1 : 0;
  Src[4] = style === 'outline' ? lay.config!.outline.selection : 0;
  ink();
  RowLook[0] = slot;
  RowLook[1] = count;
  RowLook[2] = style === 'band' ? em.highlightBand : 0;
  paintRows(lay, out, ranges);
  RowLook[0] = 0;
  RowLook[1] = 1;
  RowLook[2] = 0;
}

/**
 * A cell's box (x position `Pos[0]`, y position `Pos[1]`) in `Ink`; nothing when it is not shown. The positions through
 * an array: a hover's may be held as doubles, and a double passed to a call is boxed.
 */
const Pos = new Float64Array(2);

function paintCell(lay: GridLay, out: PassWriter) {
  const a = lay.axis!;
  const ay = lay.axisY!;
  const i = Pos[0];
  const j = Pos[1];
  if (!a.cell(i)) return;
  R[0] = Math.min(a.sx0, a.sx1);
  R[2] = Math.abs(a.sx1 - a.sx0);
  if (!ay.cell(j)) return;
  R[1] = Math.min(ay.sx0, ay.sx1);
  R[3] = Math.abs(ay.sx1 - ay.sx0);
  rect(out);
}

/** A column down the plot at x position `Pos[0]` in `Ink` (the crosshair's x half, on the hover's own section). */
function paintColumn(lay: GridLay, out: PassWriter) {
  const a = lay.axis!;
  const i = Pos[0];
  if (!a.cell(i)) return;
  const b = lay.box;
  R[0] = Math.min(a.sx0, a.sx1);
  R[1] = b.y;
  R[2] = Math.abs(a.sx1 - a.sx0);
  R[3] = b.height;
  rect(out);
}

/** Each rect over the plot's two systems in `Ink` (a selected rectangle). */
function paintRects(lay: GridLay, out: PassWriter, loc: Location) {
  const a = lay.axis!;
  const ay = lay.axisY!;
  for (let k = 0; k < loc.rects.length; k++) {
    const r = loc.rects[k];
    if (r.a.system !== lay.system || r.b.system !== lay.ySystem) continue;
    for (let p = 0; p + 1 < r.a.ranges.length; p += 2)
      for (let q = 0; q + 1 < r.b.ranges.length; q += 2) {
        if (!span(a, r.a.ranges[p], r.a.ranges[p + 1], 0) || !span(ay, r.b.ranges[q], r.b.ranges[q + 1], 2)) continue;
        R[0] = Box4[0];
        R[1] = Box4[2];
        R[2] = Box4[1] - Box4[0];
        R[3] = Box4[3] - Box4[2];
        rect(out);
      }
  }
}

/**
 * The emphasis on y (the strip paints x): highlights' rows; the selection's rows and rects (its own, else projected,
 * lighter); the hover's crosshair on this section (a column and a row and the cell's outline) or a projected hover's
 * rows; the cursor's cell when it is on a cell of the plot.
 */
function drawEmphasis(lay: GridLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const tokens = lay.tokens!;
  const config = lay.config!;
  const em = config.emphasis;
  switch (pass) {
    case 'highlight': {
      if (e.highlights && em.selectionPlace === 'over') drawRowHighlights(lay, out, e.highlights);
      return;
    }
    case 'selection': {
      const accent = tokens.colors['msui-accent'];
      const own = partOn(e.selection, lay.ySystem);
      const part =
        own ??
        (rectOver(e.selection, lay.system, lay.ySystem) ? undefined : partOn(e.projected?.selection, lay.ySystem));
      const k = own ? 1 : em.projected;
      if (part) {
        Src[0] = accent;
        Src[1] = em.selectionPlace === 'over' ? em.selectionTint * k : 0;
        Src[2] = accent;
        Src[3] = em.selectionOutline * k;
        Src[4] = config.outline.selection;
        ink();
        paintRows(lay, out, part.ranges);
      }
      if (e.selection?.rects.length) {
        Src[0] = accent;
        Src[1] = em.selectionTint;
        Src[2] = accent;
        Src[3] = em.selectionOutline;
        Src[4] = config.outline.selection;
        ink();
        paintRects(lay, out, e.selection);
      }
      return;
    }
    case 'hover': {
      const hover = tokens.colors['seq-hover'];
      const h = e.hover;
      if (
        h &&
        h.section === lay.section &&
        h.yPosition !== undefined &&
        h.position !== undefined &&
        h.element === undefined &&
        (h.ySystem === lay.ySystem || h.system === lay.system)
      ) {
        // on a grid over blocks the crosshair crosses every strip of the hovered x block (its column) and of its y
        // block (its row), the cell outlined on the one strip of both; the column and row tinted and edged in the
        // hover's ink (its contrast token: dark in light, light in dark), so they read on either end of a palette
        Src[0] = hover;
        Src[1] = CrossTint;
        Src[2] = hover;
        Src[3] = CrossEdge;
        Src[4] = config.outline.selection;
        ink();
        Pos[0] = h.position;
        Pos[1] = h.yPosition;
        if (h.system === lay.system) paintColumn(lay, out);
        Cell[0] = h.yPosition;
        Cell[1] = h.yPosition + 1;
        if (h.ySystem === lay.ySystem) paintRows(lay, out, Cell);
        if (h.ySystem !== lay.ySystem || h.system !== lay.system) return;
        Src[0] = 0;
        Src[1] = 0;
        Src[2] = hover;
        Src[3] = 1;
        Src[4] = config.outline.hover;
        ink();
        paintCell(lay, out);
        return;
      }
      const p = e.projected;
      if (p?.hover && p.hoverSection !== lay.sectionAddress) {
        const part = partOn(p.hover, lay.ySystem);
        if (!part) return;
        Src[0] = 0;
        Src[1] = 0;
        Src[2] = hover;
        Src[3] = em.projected;
        Src[4] = config.outline.hover;
        ink();
        paintRows(lay, out, part.ranges);
      }
      return;
    }
    case 'cursor': {
      const c = e.cursor;
      if (!c || c.ySystem !== lay.ySystem || c.yPosition === undefined || c.system !== lay.system) return;
      Src[0] = 0;
      Src[1] = 0;
      Src[2] = tokens.colors[em.cursorInk];
      Src[3] = 1;
      Src[4] = config.outline.cursor;
      ink();
      Pos[0] = c.position;
      Pos[1] = c.yPosition;
      paintCell(lay, out);
      return;
    }
  }
}

/** One cell's y range, reused (a hover's row). */
const Cell = new Int32Array(2);

/** What `gridSummary` gives: the mean, least and greatest value of the cells counted, and their count (NaN when none). */
export interface GridSummary {
  readonly mean: number;
  readonly min: number;
  readonly max: number;
  readonly count: number;
}

/**
 * A summary of a view's grids over a rect (_spec/README.md "Render", The grid layer, A summary over a rect): the cells
 * of every grid a grid layer of the view draws (each grid once) inside a rect of `rect` on the grid's two systems, each
 * cell once however many rects hold it, within the layer's `triangle`, missing values left out. A host's readout reads it (a selected rect's mean PAE); a library function, no state.
 */
export function gridSummary(
  view: { readonly model: StaticModel; readonly view: ResolvedView },
  rect: Location,
): GridSummary {
  const seen = new Set<string>();
  let sum = 0;
  let count = 0;
  let min = Number.NaN;
  let max = Number.NaN;
  const visit = (layer: ResolvedLayer) => {
    const owner = layer.data?.address;
    if (layer.type !== 'grid' || !owner || seen.has(owner)) return;
    seen.add(owner);
    const a = view.model.annotations.get(owner);
    const grid = a?.shape === 'grid' ? a.grid : undefined;
    if (!a || !grid) return;
    const t = layer.params.triangle;
    const triangle: Triangle = t === 'upper' || t === 'lower' ? t : 'both';
    // a cell once however many rects hold it (rects added by a drag overlap): a mark a cell of this grid
    let counted: Uint8Array | undefined;
    for (const r of rect.rects) {
      if (r.a.system !== a.target || r.b.system !== a.ySystem) continue;
      counted ??= new Uint8Array(grid.nx * grid.ny);
      const xs = r.a.ranges;
      const ys = r.b.ranges;
      for (let p = 0; p + 1 < xs.length; p += 2)
        for (let i = Math.max(0, xs[p]); i < Math.min(grid.nx, xs[p + 1]); i++)
          for (let q = 0; q + 1 < ys.length; q += 2)
            for (let j = Math.max(0, ys[q]); j < Math.min(grid.ny, ys[q + 1]); j++) {
              const k = i * grid.ny + j;
              if (counted[k] || !keeps(triangle, i, j)) continue;
              counted[k] = 1;
              const v = grid.values[k];
              if (Number.isNaN(v)) continue;
              sum += v;
              count++;
              if (!(v >= min)) min = v;
              if (!(v <= max)) max = v;
            }
    }
  };
  for (const s of view.view.sections) {
    if (s.layout !== 'grid') continue;
    for (const t of s.tracks) for (const l of t.layers) visit(l);
    for (const b of s.blocks ?? []) for (const ls of b.layers) for (const l of ls) visit(l);
  }
  return { mean: count ? sum / count : Number.NaN, min, max, count };
}
