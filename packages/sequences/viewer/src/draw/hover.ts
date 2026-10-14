/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover emphasis (seqspec:interaction#hover; HV1), part of the renderer API: the element under the pointer (an item of
// the layer, `marks.item`, or its own hit, `marks.element`) outlined on the overlay in an ink derived from the color it
// is drawn in. `hoverInk` moves that color in OKLab lightness away from the skin's surface (darker on a light skin,
// lighter on a dark one; toward it for a fill already at the far end) until it holds 3:1 against the surface and 1.5:1
// against the fill as seen, cached per color and skin in a two-way table; `HoverOutline` batches outlines `HOVER_WIDTH`
// wide inside a shape's own box (David, 2026-10-12: an outline over the edge makes the item look bigger), a small box
// filled in the ink, a glyph ringed only inside its cell, else outlined inset; `outlineItem` and `ringItem` are one
// call each for an extension's overlay. Allocation-free once warm.

import type { Axis } from "../axis.js";
import type { Skin } from "../graph.js";
import type { Features } from "../prepare.js";
import { GlyphBatch, NO_COLOR, PathBatch, PathMode, PathVerb, RectBatch, SkinToken } from "./backend.js";
import type { DrawView } from "./view.js";

/** The hover outline's width, CSS px: the pairs' hovered arc with `arcs: "hover"` too. */
export const HOVER_WIDTH = 1.5;
/** The contrast (WCAG 2 ratios) the hover ink keeps against the surface, and against the fill where the scale allows. */
export const HoverContrast = { background: 3, fill: 1.5 } as const;
/** The outline's width in device px: `HOVER_WIDTH` at the DPR, a device px at least. */
export const hoverWidth = (dpr: number): number => Math.max(1, Math.round(HOVER_WIDTH * dpr)) | 0;
/** Whether the section outlines the hovered element (`RenderProps.hoverEmphasis`, the style key `hover_emphasis`). */
export const hoverEmphasis = (v: DrawView): boolean => v.props.hoverEmphasis !== false;

// ---- The ink ---------------------------------------------------------------------------------------------------------

/** One skin's inks: its surface and tokens as packed colors, and a two-way table of key → ink and fill. */
interface InkTable {
  readonly skin: Skin | null;
  readonly dark: boolean;
  readonly bg: number;
  readonly bgY: number;
  /** By `-token - 2`: a token's color and alpha. */
  readonly tokens: Int32Array;
  readonly tokenAlpha: Uint8Array;
  readonly keys: Int32Array;
  readonly inks: Int32Array;
  /** The fill as seen (its alpha over the surface), opaque: what an inset glyph outline fills inside. */
  readonly fills: Int32Array;
  readonly filled: Uint8Array;
}
/** Two ways of 256 sets: a pair of a hover's colors that hash to one set keeps both. */
const TABLE = 512;
/** Skin keys by `-token - 2` (`SkinToken`). */
const NAMES: (keyof Skin)[] = [];
for (const [name, id] of Object.entries(SkinToken)) NAMES[-id - 2] = name as keyof Skin;
const held: (InkTable | null)[] = [null, null, null, null];
let next = 0;
const rgba = new Float64Array(4);

/**
 * A CSS color as the skin reader writes it (`rgb(r g b)`, `rgb(r g b / a)`, `rgba(r, g, b, a)`) or `#rgb`, `#rrggbb`,
 * `#rrggbbaa`, into `rgba` (0..255, alpha 0..1); false: another form.
 */
function parse(css: unknown): boolean {
  if (typeof css !== "string") return false;
  const t = css.trim().toLowerCase();
  if (t[0] === "#") {
    const h = t.length === 4 ? t[1] + t[1] + t[2] + t[2] + t[3] + t[3] : t.slice(1);
    if (!/^[0-9a-f]{6}([0-9a-f]{2})?$/.test(h)) return false;
    for (let i = 0; i < 3; i++) rgba[i] = Number.parseInt(h.slice(2 * i, 2 * i + 2), 16);
    rgba[3] = h.length === 8 ? Number.parseInt(h.slice(6), 16) / 255 : 1;
    return true;
  }
  const m = /^rgba?\(([^)]*)\)$/.exec(t);
  const p = m ? m[1].split(/[\s,/]+/).filter(Boolean) : [];
  if (p.length < 3) return false;
  for (let i = 0; i < 3; i++) rgba[i] = Math.max(0, Math.min(255, Number.parseFloat(p[i]) || 0));
  rgba[3] = p.length > 3 ? Math.max(0, Math.min(1, Number.parseFloat(p[3]) / (p[3].endsWith("%") ? 100 : 1))) : 1;
  return true;
}
const packed = () => (Math.round(rgba[0]) << 16) | (Math.round(rgba[1]) << 8) | Math.round(rgba[2]);

/** Edge rate: a skin's first hover (a skin change, or a second viewer's skin), apart from the hover's path. */
function newTable(skin: Skin | null): InkTable {
  const dark = skin?.dark === true;
  const bg = parse(skin?.surface) ? packed() : dark ? 0x181a1e : 0xffffff;
  const tokens = new Int32Array(NAMES.length);
  const tokenAlpha = new Uint8Array(NAMES.length);
  for (let i = 0; i < NAMES.length; i++) {
    const ok = parse(skin?.[NAMES[i]]);
    tokens[i] = ok ? packed() : dark ? 0xe6e8eb : 0x1b1f24;
    tokenAlpha[i] = ok ? Math.round(rgba[3] * 255) : 255;
  }
  const t: InkTable = {
    skin,
    dark,
    bg,
    bgY: luminance(bg),
    tokens,
    tokenAlpha,
    keys: new Int32Array(TABLE),
    inks: new Int32Array(TABLE),
    fills: new Int32Array(TABLE),
    filled: new Uint8Array(TABLE),
  };
  held[next] = t;
  next = (next + 1) % held.length;
  return t;
}
function tableOf(skin: Skin | null): InkTable {
  for (let i = 0; i < held.length; i++) {
    const t = held[i];
    if (t !== null && t.skin === skin) return t;
  }
  return newTable(skin);
}

/**
 * The hover outline's ink for what is drawn in `key` (0xRRGGBB, a skin token, a color key with its alpha folded in as
 * `featureKey` returns; `NO_COLOR` counts as the surface) on `skin`: a packed color, cached per key and skin.
 */
export function hoverInk(skin: Skin | null, key: number): number {
  return entry(skin, key).inks[slot];
}
/** What `key` looks like on `skin`: its color with its alpha over the surface, opaque (cached with `hoverInk`). */
export function hoverFill(skin: Skin | null, key: number): number {
  return entry(skin, key).fills[slot];
}
let slot = 0;
function entry(skin: Skin | null, key: number): InkTable {
  const t = tableOf(skin);
  const set = (Math.imul(key, 0x9e3779b1) >>> 24) << 1;
  if ((t.filled[set] & 1) === 1 && t.keys[set] === key) {
    slot = set;
    return t;
  }
  if ((t.filled[set + 1] & 1) === 1 && t.keys[set + 1] === key) {
    slot = set + 1;
    return t;
  }
  // A miss: into the way not filled last (bit 1 of the set's first entry marks its second as the newer).
  const h = (t.filled[set] & 2) === 2 ? set : set + 1;
  t.filled[set] = (t.filled[set] & 1) | (h === set + 1 ? 2 : 0);
  slot = h;
  t.inks[h] = derive(t, key);
  t.fills[h] = seen;
  t.keys[h] = key;
  t.filled[h] |= 1;
  return t;
}

const lin = (c: number): number => (c <= 10.31475 ? c / 3294.6 : ((c / 255 + 0.055) / 1.055) ** 2.4);
const enc = (x: number): number => {
  const v = x <= 0.0031308 ? 12.92 * x : 1.055 * x ** (1 / 2.4) - 0.055;
  return Math.max(0, Math.min(255, Math.round(v * 255)));
};
function luminance(c: number): number {
  return 0.2126 * lin((c >> 16) & 255) + 0.7152 * lin((c >> 8) & 255) + 0.0722 * lin(c & 255);
}
/** The WCAG 2 contrast ratio of two relative luminances. */
export const contrastRatio = (y0: number, y1: number): number =>
  y0 > y1 ? (y0 + 0.05) / (y1 + 0.05) : (y1 + 0.05) / (y0 + 0.05);

/** OKLab (L, a, b) and linear sRGB, through scratch. */
const lab = new Float64Array(3);
const rgb = new Float64Array(3);
function toLab(c: number): void {
  const r = lin((c >> 16) & 255);
  const g = lin((c >> 8) & 255);
  const b = lin(c & 255);
  const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b);
  const m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b);
  const s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
  lab[0] = 0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s;
  lab[1] = 1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s;
  lab[2] = 0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s;
}
/** Linear sRGB of OKLab (L, a, b) into `rgb`; whether it is in gamut. */
function fromLab(L: number, a: number, b: number): boolean {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3;
  const m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3;
  const s = (L - 0.0894841775 * a - 1.291485548 * b) ** 3;
  rgb[0] = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
  rgb[1] = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
  rgb[2] = -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s;
  const e = 1e-7;
  return rgb[0] >= -e && rgb[0] <= 1 + e && rgb[1] >= -e && rgb[1] <= 1 + e && rgb[2] >= -e && rgb[2] <= 1 + e;
}
/** OKLab lightness L at hue (ca, cb) and the largest chroma up to `c` in sRGB (12 bisection steps), 8 bits. */
function inGamut(L: number, c: number, ca: number, cb: number): number {
  if (!fromLab(L, c * ca, c * cb)) {
    let lo = 0;
    let hi = c;
    for (let i = 0; i < 12; i++) {
      const mid = (lo + hi) / 2;
      if (fromLab(L, mid * ca, mid * cb)) lo = mid;
      else hi = mid;
    }
    fromLab(L, lo * ca, lo * cb);
  }
  return (enc(rgb[0]) << 16) | (enc(rgb[1]) << 8) | enc(rgb[2]);
}
/** The last `derive`'s fill as seen. */
let seen = 0;
/** The fill as seen (over the surface), then lightness steps of 0.01 away from the surface. */
function derive(t: InkTable, key: number): number {
  let c = t.bg;
  let a = 1;
  if (key !== NO_COLOR && key < 0 && key >= -64) {
    const i = -key - 2;
    c = i < t.tokens.length ? t.tokens[i] : t.tokens[-SkinToken.ink - 2];
    a = i < t.tokens.length ? t.tokenAlpha[i] / 255 : 1;
  } else if (key !== NO_COLOR) {
    c = key & 0xffffff;
    a = (255 - (key >>> 24)) / 255;
  }
  const bg = t.bg;
  const fill =
    a >= 1
      ? c
      : (Math.round(((c >> 16) & 255) * a + ((bg >> 16) & 255) * (1 - a)) << 16) |
        (Math.round(((c >> 8) & 255) * a + ((bg >> 8) & 255) * (1 - a)) << 8) |
        Math.round((c & 255) * a + (bg & 255) * (1 - a));
  seen = fill;
  const yf = luminance(fill);
  toLab(fill);
  const L = lab[0];
  const chroma = Math.hypot(lab[1], lab[2]);
  const ca = chroma > 1e-9 ? lab[1] / chroma : 0;
  const cb = chroma > 1e-9 ? lab[2] / chroma : 0;
  // Away from the surface first; a fill already at that end (a near-black bar on a light skin) has no 1.5:1 there, and
  // its outline, inside it, goes the other way.
  const away = t.dark ? 1 : -1;
  const ink = walk(t, L, chroma, ca, cb, yf, away);
  if (ink >= 0) return ink;
  const back = walk(t, L, chroma, ca, cb, yf, -away);
  return back >= 0 ? back : away < 0 ? 0 : 0xffffff;
}
/** The last `walk`'s color. */
let walked = 0;
/** Lightness steps of 0.01 in direction `dir` until 3:1 against the surface and 1.5:1 against the fill: the color; -1 at the end. */
function walk(t: InkTable, L: number, chroma: number, ca: number, cb: number, yf: number, dir: number): number {
  for (let i = 1; i <= 100; i++) {
    const L1 = Math.min(1, Math.max(0, L + dir * 0.01 * i));
    walked = inGamut(L1, chroma, ca, cb);
    const y = luminance(walked);
    if (contrastRatio(y, t.bgY) >= HoverContrast.background && contrastRatio(y, yf) >= HoverContrast.fill)
      return walked;
    if (L1 <= 0 || L1 >= 1) return -1;
  }
  return -1;
}

// ---- The outline ------------------------------------------------------------------------------------------------------

const boxes = new RectBatch(16);
const rings = new PathBatch(8);
const inkGlyphs = new GlyphBatch(8);
const fillGlyphs = new GlyphBatch(8);

/**
 * What the built-ins' overlays batch through: `reset`, then boxes and glyphs (device px; `key` the color drawn, as
 * `hoverInk` takes it), then `hand`, which draws nothing where the section turns hover emphasis off. Everything stays
 * inside the item's own box.
 */
export const HoverOutline = {
  reset(): void {
    boxes.reset();
    rings.reset();
    inkGlyphs.reset();
    fillGlyphs.reset();
  },
  /**
   * An outline `hoverWidth` wide inside the rectangle (x, y, w, h): its outer edge the item's drawn edge. A box
   * narrower than twice the width, or lower than twice the width or 3 device px, is filled in the ink instead.
   */
  box(v: DrawView, key: number, x: number, y: number, w: number, h: number): void {
    if (w <= 0 || h <= 0) return;
    const d = hoverWidth(v.dpr);
    const fill = w < 2 * d || h < 2 * d || h < 3;
    boxes.begin(hoverInk(v.skin, key), 255, fill ? 0 : d);
    boxes.push(x | 0, y | 0, w | 0, h | 0);
  },
  /**
   * Glyph `id` of radius `r` at (cx, cy) in a cell reaching `half` device px from its center each way: a ring around it
   * where the ring stays inside the cell, else the glyph's own outline inset (the glyph in the ink, then in its fill as
   * seen, `hoverWidth` smaller).
   */
  glyph(v: DrawView, key: number, id: number, cx: number, cy: number, r: number, half: number): void {
    const d = hoverWidth(v.dpr);
    if (r + d + 1 <= half) {
      const rr = r + 0.5 + d / 2;
      rings.begin(hoverInk(v.skin, key), 255, PathMode.stroke, 2 * d);
      let o = rings.verb(PathVerb.move);
      rings.coords[o] = cx + rr;
      rings.coords[o + 1] = cy;
      o = rings.verb(PathVerb.arc);
      rings.coords[o] = cx;
      rings.coords[o + 1] = cy;
      rings.coords[o + 2] = rr;
      rings.coords[o + 3] = 0;
      rings.coords[o + 4] = 2 * Math.PI;
      rings.coords[o + 5] = 0;
      rings.verb(PathVerb.close);
      return;
    }
    inkGlyphs.begin(hoverInk(v.skin, key));
    let i = inkGlyphs.add(id);
    inkGlyphs.x[i] = cx;
    inkGlyphs.y[i] = cy;
    inkGlyphs.radius[i] = r;
    if (r <= 2 * d) return;
    fillGlyphs.begin(hoverFill(v.skin, key));
    i = fillGlyphs.add(id);
    fillGlyphs.x[i] = cx;
    fillGlyphs.y[i] = cy;
    fillGlyphs.radius[i] = r - d;
  },
  hand(v: DrawView): void {
    if (!hoverEmphasis(v)) return;
    if (boxes.n > 0) v.backend.rects(boxes);
    if (rings.n > 0) v.backend.paths(rings);
    if (inkGlyphs.n === 0) return;
    inkGlyphs.font = fillGlyphs.font = v.text.font("ui", 12, 600);
    v.backend.glyphs(inkGlyphs);
    if (fillGlyphs.n > 0) v.backend.glyphs(fillGlyphs);
  },
};

/** One call for an extension's `overlay`: the hovered item's box (device px) outlined inside, in the hover ink of `key`. */
export function outlineItem(v: DrawView, key: number, x: number, y: number, w: number, h: number): void {
  HoverOutline.reset();
  HoverOutline.box(v, key, x, y, w, h);
  HoverOutline.hand(v);
}
/**
 * One call for an extension's `overlay`: the hovered glyph `id` (radius `r` at (cx, cy), in a cell `half` device px
 * each way) ringed inside its cell, else outlined inset (`HoverOutline.glyph`).
 */
export function ringItem(v: DrawView, key: number, id: number, cx: number, cy: number, r: number, half: number): void {
  HoverOutline.reset();
  HoverOutline.glyph(v, key, id, cx, cy, r, half);
  HoverOutline.hand(v);
}

// ---- The hovered item's features ------------------------------------------------------------------------------------

/** Features by item: per feature in (item, annotation) order, its item, annotation and index. */
export interface ItemOrder {
  readonly item: Int32Array;
  readonly annotation: Int32Array;
  readonly order: Int32Array;
}
const orders = new WeakMap<Features, ItemOrder>();
/** Edge rate: on the layer's first item hover. */
function newOrder(f: Features, set: Axis["set"]): ItemOrder {
  const id = (k: number) => set.annotations[f.annotation[k]]?.item[f.row[k]] ?? -1;
  const ks = Array.from({ length: f.count }, (_, k) => k).sort(
    (a, b) => id(a) - id(b) || f.annotation[a] - f.annotation[b] || a - b,
  );
  const x = {
    item: Int32Array.from(ks, id),
    annotation: Int32Array.from(ks, (k) => f.annotation[k]),
    order: Int32Array.from(ks),
  };
  orders.set(f, x);
  return x;
}
/**
 * The hovered item of the layer drawn (`marks.item`, `itemAnnotation`) among `f`: into `span` the range [lo, hi) of its
 * features in the returned order (`order[j]`, those of another annotation skipped by the caller through
 * `annotation[j]`); null: no item of this layer is hovered.
 */
export function hoveredItem(v: DrawView, f: Features, span: Int32Array): ItemOrder | null {
  const item = v.marks.item(v.layer);
  if (item < 0 || f.count === 0 || v.axis === null) return null;
  const x = orders.get(f) ?? newOrder(f, v.axis.set);
  let lo = 0;
  let hi = f.count;
  while (lo < hi) {
    const m = (lo + hi) >> 1;
    if (x.item[m] < item) lo = m + 1;
    else hi = m;
  }
  let end = lo;
  while (end < f.count && x.item[end] === item) end++;
  span[0] = lo;
  span[1] = end;
  return end > lo ? x : null;
}
