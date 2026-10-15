/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the built-in renderers share, part of the renderer API (draw/api.ts; seqspec:viewer#extension-points):
// `kernelEntry` (a kernel type's entry from its prepare in `PreparedParts`), legends over channels, item hits, the lane
// geometry, `residue_form`, and the pairs layer's geometry, one function for the draw, the
// overlay and the hit (an arc lit by a hover is the arc drawn). The entries themselves are one module each under
// `renderers/`. Headless and allocation-free.

import { CoordinateSystemKind } from '@molstar/seq-core/model';
import type { AxisView } from './axis.js';
import { DOT_R } from './draw/sizes.js';
import { type Hit, HoverKind } from './hit.js';
import type { RenderProps } from './layout.js';
import { pairsInView } from './pairs-index.js';
import { type PairsState, type PreparedPart, PreparedParts } from './prepare.js';
import type { LayerRenderer } from './registry.js';
import type { FontId, TextManagerApi } from './text.js';
import { type LegendSink, NO_SIZE, SIZE_ONE, type Theme } from './themes.js';
import type { TreeRows } from './tree.js';

type Params = Readonly<Record<string, unknown>>;

/** A layer's legend: each channel's theme's, in its type's channel order. */
export function legendOf(_p: unknown, themes: readonly (Theme | null)[], out: LegendSink) {
  for (const t of themes) t?.legend?.(out);
}

/** Lanes in a band: the padding and each lane's height (the draw and the hit share them). */
export function laneHeight(height: number, lanes: number): number {
  return Math.max(1, (height - 2 * lanePad(height)) / Math.max(1, lanes));
}
export const lanePad = (height: number) => (height > 10 ? 2 : 0);
/** The gap above and below a feature in its lane (px). */
export const laneInset = (laneH: number) => (laneH >= 6 ? 1 : 0);

/**
 * The line a letter is centered on, by code (seqspec:viewer#the-render-backend, text placement): 1, the x-height, for
 * lower case and the gaps (`-` and `.`, drawn for lower case); 0, a capital's height, for every other.
 */
export const LETTER_LINE = Uint8Array.from({ length: 128 }, (_, c) =>
  (c >= 97 && c <= 122) || c === 45 || c === 46 ? 1 : 0,
);

/**
 * The alphabetic baselines that center a letter on the band [y0, y1) of device px (its cell): `out[0]` a capital's ink,
 * `out[1]` a lower-case letter's x-height (`LETTER_LINE`), from the font's ink heights at the view's DPR
 * (`TextManager.ink`). The em box's middle sits 0.1 em under a capital's (ui-monospace), so a `middle` baseline drew
 * capitals high. Whole device px: Canvas2D snaps an upright baseline to one (Chromium), and the glyph atlas then draws
 * each letter at one phase.
 */
export function letterLines(
  v: { readonly text: TextManagerApi; readonly dpr: number },
  font: FontId,
  y0: number,
  y1: number,
  out: Float64Array,
): void {
  const ink = v.text.ink(font, Math.round(v.dpr * 64));
  // Half way between two px, the higher: the optical center is a little above the middle.
  out[0] = Math.ceil((y0 + y1 + ink[0]) / 2 - 0.5);
  out[1] = Math.ceil((y0 + y1 + ink[1]) / 2 - 0.5);
}

/** An item hit: the annotation's row as its item (kind item), with its location and coordinate system. */
export function hitItem(out: Hit, set: PairsState['set'], annotation: number, row: number): boolean {
  const a = set.annotations[annotation];
  if (!a) return false;
  out.kind = HoverKind.item;
  out.annotation = annotation;
  out.item = a.item[row];
  out.itemLoc.set = set;
  out.itemLoc.annotation = annotation as never;
  out.itemLoc.item = out.item as never;
  out.system.kind = CoordinateSystemKind.item;
  out.system.id = annotation;
  out.system.row = out.item;
  return true;
}

/** A pair end's line from the band's edge, CSS px: the end's dot and a margin (read once: see `glyphFrame`). */
const ENDS_INSET = Math.ceil(DOT_R + 2.5);
let inView = new Int32Array(64);
/** How far past the data area's edge an arc's piece runs (CSS px), so its stroke's cap lies outside the clip. */
export const ARC_PAST = 4;
/** How far past the window's edge an arc's end is taken at most (`Arcs.virtual`), as a share of the window's width. */
export const ARC_MARGIN = 0.25;
/** The shown ranges the window and its margins cover, per block (`Arcs.frame`'s query). */
let marginLo = new Int32Array(0);
let marginHi = new Int32Array(0);
/** The pairs layer's arcs (seqspec:format#layer-types `pairs`): one geometry in any view, for every pass. */
export const Arcs = {
  /**
   * Into `s.frame`: the ends' line (`side`, else toward the residues: `v.facing`), the direction arcs rise in, the
   * band they fill and the span that fills it (heights are proportional to span, so nested arcs nest): the widest
   * virtual span (`virtual`: each end taken at most `ARC_MARGIN` of the window past its edge) over the pairs within the
   * window and its margins, from the interval tree. Every arc's apex is then within the band, and as a pan moves the
   * ends, the virtual ends and this scale move continuously. Recomputed only when the view moved.
   */
  frame(v: AxisView, s: PairsState, p: Params) {
    const f = s.frame;
    const below = p.side === 'below' || (p.side !== 'above' && v.facing > 0);
    const ends = p.display !== 'arcs';
    const inset = ends ? ENDS_INSET : 3;
    const base = below ? v.top + inset : v.top + v.height - inset;
    const left = v.blockX.length > 0 ? v.blockX[0] : 0;
    const at = s.at;
    if (
      at.axis !== v.axis ||
      f[4] !== left ||
      f[5] !== v.cell ||
      f[6] !== v.left ||
      f[7] !== v.right ||
      f[8] !== v.onlyBlock
    ) {
      at.axis = v.axis;
      f[4] = left;
      f[5] = v.cell;
      f[6] = v.left;
      f[7] = v.right;
      f[8] = v.onlyBlock;
      if (inView.length < s.count) inView = new Int32Array(s.count);
      // The pairs within the window and its margins, by x: the whole view's in a pan's uncovered strip too (the frame
      // is kept, and the layer's redraw after the copy reuses it; over the strip's few pairs it was another scale). A
      // pair with an end within a margin counts before it shows, so none enters the scale at once.
      const lo = v.visLo;
      const hi = v.visHi;
      const n0 = lo.length;
      if (marginLo.length < n0) {
        marginLo = new Int32Array(n0);
        marginHi = new Int32Array(n0);
      }
      const margin = ARC_MARGIN * (v.right - v.left);
      for (let b = 0; b < n0; b++) {
        const w = v.cell * v.span[b];
        const max = v.axis && b < v.axis.count ? v.axis.shown[b] : hi[b];
        const a = w > 0 ? Math.floor((v.left - margin - v.blockX[b]) / w) - 1 : 0;
        const z = w > 0 ? Math.ceil((v.right + margin - v.blockX[b]) / w) + 1 : max;
        marginLo[b] = a < 0 ? 0 : a > max ? max : a;
        marginHi[b] = z < marginLo[b] ? marginLo[b] : z > max ? max : z;
      }
      v.visLo = marginLo;
      v.visHi = marginHi;
      const n = pairsInView(v, s, inView);
      v.visLo = lo;
      v.visHi = hi;
      let widest = 0;
      for (let j = 0; j < n; j++) {
        if (!Arcs.ends(v, s, inView[j])) continue;
        Arcs.virtual(s);
        if (s.xs[3] - s.xs[2] > widest) widest = s.xs[3] - s.xs[2];
      }
      f[3] = widest;
    }
    f[0] = Math.round(base * v.dpr) / v.dpr;
    f[1] = below ? 1 : -1;
    f[2] = Math.max(0, v.height - 9);
  },
  /**
   * One end's center x (side 0: a, 1: b) into `s.xs[side]`; NaN when filtered out, not placed (block -1) or on another
   * row (stacked).
   */
  end(v: AxisView, s: PairsState, k: number, side: number) {
    const e = 2 * k + side;
    const b = s.block[e];
    const k0 = b < 0 || (v.onlyBlock >= 0 && v.onlyBlock !== b) ? -1 : v.shownOf(b, s.start[e]);
    const k1 = k0 < 0 ? -1 : v.shownOf(b, s.stop[e] - 1);
    s.xs[side] = k0 < 0 || k1 < 0 ? Number.NaN : (v.x(b, k0 as never) + v.x(b, (k1 + 1) as never)) / 2;
  },
  /** Both ends into `s.xs`; true when both are drawn here. */
  ends(v: AxisView, s: PairsState, k: number): boolean {
    Arcs.end(v, s, k, 0);
    Arcs.end(v, s, k, 1);
    return !Number.isNaN(s.xs[0]) && !Number.isNaN(s.xs[1]);
  },
  /**
   * The virtual ends of the pair in `s.xs` (after `ends`) into `s.xs[2]` and `s.xs[3]`, the lower first: each end
   * taken at most `ARC_MARGIN` of the window's width (the frame's, `s.frame[6..7]`) past its edge. An arc's curve is
   * the parabola over them: a partner far beyond the window rises from its end in view toward the edge as to a
   * partner just past the margin, and moves on continuously as it comes into view.
   */
  virtual(s: PairsState) {
    const xs = s.xs;
    const f = s.frame;
    const margin = ARC_MARGIN * (f[7] - f[6]);
    const lo = f[6] - margin;
    const hi = f[7] + margin;
    const a = xs[0] < xs[1] ? xs[0] : xs[1];
    const b = xs[0] < xs[1] ? xs[1] : xs[0];
    xs[2] = a < lo ? lo : a > hi ? hi : a;
    xs[3] = b < lo ? lo : b > hi ? hi : b;
  },
  /**
   * Into `s.frame[9]`: the apex height of the arc between the ends in `s.xs` (after `frame`), at most the band
   * (`s.frame[2]`): by its virtual span (`virtual`, into `s.xs[2..3]`) against the widest in view, or with a `height`
   * (a size theme's int, `SIZE_ONE` the band; `NO_SIZE`: none) that share of the band.
   */
  apex(s: PairsState, height = NO_SIZE) {
    const f = s.frame;
    Arcs.virtual(s);
    if (height >= 0) f[9] = (f[2] * Math.min(height, SIZE_ONE)) / SIZE_ONE;
    else f[9] = f[3] > 0 ? Math.min(f[2], (f[2] * (s.xs[3] - s.xs[2])) / f[3]) : 0;
  },
  /**
   * The height over the ends' line at x (`s.frame[10]`) of the arc in `s.xs` (after `apex`) into `s.frame[11]`; -1:
   * none there: the parabola over its virtual ends, `4 · apex · t(1 − t)`. In and out through the frame (no double
   * crosses a call). The draw (`pushArc`), its labels and the hit share it.
   */
  heightAt(s: PairsState) {
    const f = s.frame;
    const x = f[10];
    const a = s.xs[2];
    const b = s.xs[3];
    const t = a === b ? (x === a ? 0 : -1) : (x - a) / (b - a);
    f[11] = t < 0 || t > 1 ? -1 : 4 * f[9] * t * (1 - t);
  },
};

/** A collapsed clade's row label (laid out with the tree, "SYN001 +19"); "" when none. */
export const wedgeLabel = (s: TreeRows, node: number) => s.tree?.wedgeText.get(node) ?? '';

/** Kernel types that draw images (matrix tiles, dense rows' tiles, the map's miniature): no ring form. */
const IMAGES = new Set(['matrix', 'rows', 'seqv:map']);

/**
 * A kernel type's entry: its prepare and what goes with it (`PreparedParts`), its members over the defaults (linear,
 * unstable, 6 px, the legend over its channels). Stable (a copied strip equals a full redraw): what draws only
 * rectangles on device pixels and letters within their cells. Paths (glyphs, curves, a values line), labels wider than
 * a cell and what stays in view are redrawn whole on a copy. Its linear draw holds on rings, mapped by the engine
 *, so it declares `circular` too, unless it draws images (`IMAGES`), which rings do not show.
 */
export function kernelEntry<S>(
  name: string,
  height: number | ((s: S, p: unknown, props: RenderProps, rowScale?: number, width?: number, room?: number) => number),
  o: Omit<Partial<LayerRenderer<S>>, 'name' | 'height'> = {},
): LayerRenderer<S> {
  const part = PreparedParts[name] as PreparedPart;
  const geometries = o.geometries ?? ['linear'];
  return {
    name,
    stable: false,
    minHeight: 6,
    legend: legendOf,
    ...part,
    ...o,
    geometries: IMAGES.has(name) ? geometries : [...geometries, 'circular'],
    height: typeof height === 'number' ? () => height : height,
  } as LayerRenderer<S>;
}
