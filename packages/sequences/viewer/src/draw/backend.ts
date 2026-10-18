/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The backend contract and the batches every renderer draws through (seqspec:viewer#the-render-backend), a leaf module
// the renderer API (api.ts) re-exports. Batches are struct-of-arrays and grow only, one of each kind per renderer
// module, reused by every draw, so a frame allocates nothing: a draw calls `begin(color, alpha)`, pushes, and hands the
// batch to the backend once. Positions are device px: rectangle edges integers (Smis), the rest doubles in typed
// arrays.
// V8 boxes a double that crosses a call it does not inline, so a draw with a large loop writes doubles through the
// slot `add()` returns rather than calling `push` with them (`push` is for small draws and tests).

/** The composer's passes, one target each. */
export const Pass = { main: 0, overlay: 1, nav: 2, navOverlay: 3, vnav: 4, vnavOverlay: 5 } as const;
/** `RenderBackend.caps` bits: without `Copy` the composer never copies, without `Back` it draws whole frames. */
export const BackendCaps = { Copy: 1, Progressive: 2, RotatedText: 4, Back: 8 } as const;
/** A color that draws nothing (a theme's "no color", seqspec:viewer#extension-points): a backend skips its items. */
export const NO_COLOR = -1;
/**
 * Skin tokens as colors (seqspec:interaction#the-skin): negative ids a backend resolves through `DrawResources.token`.
 */
export const SkinToken = {
  surface: -2,
  surfaceAlt: -3,
  ink: -4,
  inkQuiet: -5,
  inkFaint: -6,
  grid: -7,
  accent: -8,
  hover: -9,
  highlight: -10,
  selection: -11,
  selectionBorder: -12,
  unobserved: -13,
  conflict: -14,
  pending: -15,
  grip: -16,
  gripHover: -17,
  focusRing: -18,
  window: -19,
} as const;
/** The shape GlyphIds, the same in every glyph table; ids from `SHAPE_COUNT` on are characters and paths. */
export const GlyphShape = {
  circle: 0,
  diamond: 1,
  triangle: 2,
  square: 3,
  pin: 4,
  star: 5,
  cross: 6,
  tick: 7,
} as const;
export const SHAPE_COUNT = 8;
/** Path verbs and the coordinates each takes: move, line 2; quad 4; cubic 6; arc 6 (cx cy r a0 a1 ccw); round-rect 5. */
export const PathVerb = { move: 0, line: 1, quad: 2, cubic: 3, arc: 4, roundRect: 5, close: 6 } as const;
/** The coordinates each `PathVerb` takes, by verb. */
export const PathVerbCoords: Readonly<Int32Array> = new Int32Array([2, 2, 4, 6, 6, 5, 0]);
export const PathMode = { fill: 0, stroke: 1 } as const;
export const TextAlign = { left: 0, center: 1, right: 2 } as const;
export const TextBaseline = { alphabetic: 0, middle: 1, top: 2, bottom: 3 } as const;

export interface RenderBackend {
  readonly caps: number;
  readonly images: BackendImages;
  // The composer's: passes, clears, clips, copies and presentation; device px, integers.
  beginPass(pass: number, width: number, height: number, dpr: number): void;
  endPass(): void;
  clear(x: number, y: number, w: number, h: number): void;
  pushClip(x: number, y: number, w: number, h: number): void;
  popClip(): void;
  /** Only with Copy: the region (x, y, w, h) set to the pass's image as it was at its first copy, moved by (dx, dy). */
  copy(dx: number, dy: number, x: number, y: number, w: number, h: number): void;
  /** Only with Back: a finished region of the back target, shown at `endPass`. */
  present(x: number, y: number, w: number, h: number): void;
  /** Only with Back: the last presented image's (sx, sy, sw, sh) shown at (dx, dy, dw, dh), at `endPass`, under the presents. */
  stale(sx: number, sy: number, sw: number, sh: number, dx: number, dy: number, dw: number, dh: number): void;
  // The renderers': one call per batch, never per element.
  rects(b: RectBatch): void;
  sectors(b: SectorBatch): void;
  paths(b: PathBatch): void;
  glyphs(b: GlyphBatch): void;
  text(b: TextBatch): void;
  /** `r`: sx sy sw sh dx dy dw dh, device px (integers draw without boxing). */
  image(img: number, r: Float64Array, smooth: boolean): void;
}
export interface BackendImages {
  /** An image id from pixels in ImageData's byte order (R, G, B, A in memory). */
  fromPixels(w: number, h: number, rgba: Uint32Array): number;
  /** A backend drawing a w × h image (its pass begun), until `endLayer`; one at a time. */
  beginLayer(w: number, h: number): RenderBackend;
  endLayer(): number;
  dispose(img: number): void;
}
/**
 * A character or path glyph (seqspec:format#themes): `path` is SVG path data in a `box`, scaled so its larger side is
 * 2r.
 */
export interface GlyphRecord {
  readonly char: string | null;
  readonly path: string | null;
  readonly box: readonly [number, number];
  /** The shape drawn where the record has neither. */
  readonly fallback: number;
}
/** What a backend resolves ids through: the context's text manager, skin and glyph table. Called per batch or group. */
export interface DrawResources {
  /** A FontId's CSS font, its size in CSS px. */
  font(id: number): string;
  fontPx(id: number): number;
  /** A skin token's CSS color. */
  token(token: number): string;
  /** A glyph from SHAPE_COUNT on; null: draw the circle. */
  glyph(id: number): GlyphRecord | null;
}

type Typed = Int32Array | Uint8Array | Uint16Array | Float64Array;
/** `a` copied into a new array of length `n` (growth only: never on a warm frame). */
function grown<T extends Typed>(a: T, n: number): T {
  const b = new (a.constructor as new (n: number) => T)(n);
  b.set(a);
  return b;
}

/** Rectangles, integer device px; `stroke` 0 fills, else it is the outline's width in device px, inside the rect. */
export class RectBatch {
  n = 0;
  x: Int32Array;
  y: Int32Array;
  w: Int32Array;
  h: Int32Array;
  color: Int32Array;
  alpha: Uint8Array;
  stroke: Uint8Array;
  private c = 0;
  private a = 255;
  private s = 0;
  constructor(capacity = 256) {
    this.x = new Int32Array(capacity);
    this.y = new Int32Array(capacity);
    this.w = new Int32Array(capacity);
    this.h = new Int32Array(capacity);
    this.color = new Int32Array(capacity);
    this.alpha = new Uint8Array(capacity);
    this.stroke = new Uint8Array(capacity);
  }
  /** The color, alpha (0..255) and stroke of the pushes that follow. */
  begin(color: number, alpha = 255, stroke = 0): void {
    this.c = color;
    this.a = alpha;
    this.s = stroke;
  }
  push(x: number, y: number, w: number, h: number): void {
    if (this.n === this.x.length) this.grow(this.n + 1);
    const i = this.n++;
    this.x[i] = x;
    this.y[i] = y;
    this.w[i] = w;
    this.h[i] = h;
    this.color[i] = this.c;
    this.alpha[i] = this.a;
    this.stroke[i] = this.s;
  }
  reset(): void {
    this.n = 0;
  }
  /** Items [start, start + count) of `b` added, their colors with them. */
  append(b: RectBatch, start: number, count: number): void {
    if (this.n + count > this.x.length) this.grow(this.n + count);
    for (let k = 0; k < count; k++) {
      const i = this.n++;
      const j = start + k;
      this.x[i] = b.x[j];
      this.y[i] = b.y[j];
      this.w[i] = b.w[j];
      this.h[i] = b.h[j];
      this.color[i] = b.color[j];
      this.alpha[i] = b.alpha[j];
      this.stroke[i] = b.stroke[j];
    }
  }
  private grow(min: number): void {
    const m = Math.max(min, this.x.length * 2);
    this.x = grown(this.x, m);
    this.y = grown(this.y, m);
    this.w = grown(this.w, m);
    this.h = grown(this.h, m);
    this.color = grown(this.color, m);
    this.alpha = grown(this.alpha, m);
    this.stroke = grown(this.stroke, m);
  }
}

/** Annular sectors (the circular geometry): center, radii r0 < r1, angles a0 to a1 clockwise, radians. */
export class SectorBatch {
  n = 0;
  cx: Float64Array;
  cy: Float64Array;
  r0: Float64Array;
  r1: Float64Array;
  a0: Float64Array;
  a1: Float64Array;
  color: Int32Array;
  alpha: Uint8Array;
  private c = 0;
  private a = 255;
  constructor(capacity = 64) {
    this.cx = new Float64Array(capacity);
    this.cy = new Float64Array(capacity);
    this.r0 = new Float64Array(capacity);
    this.r1 = new Float64Array(capacity);
    this.a0 = new Float64Array(capacity);
    this.a1 = new Float64Array(capacity);
    this.color = new Int32Array(capacity);
    this.alpha = new Uint8Array(capacity);
  }
  begin(color: number, alpha = 255): void {
    this.c = color;
    this.a = alpha;
  }
  /** A new item in the current color: its slot, to write the doubles into. */
  add(): number {
    if (this.n === this.cx.length) this.grow(this.n + 1);
    const i = this.n++;
    this.color[i] = this.c;
    this.alpha[i] = this.a;
    return i;
  }
  push(cx: number, cy: number, r0: number, r1: number, a0: number, a1: number): void {
    const i = this.add();
    this.cx[i] = cx;
    this.cy[i] = cy;
    this.r0[i] = r0;
    this.r1[i] = r1;
    this.a0[i] = a0;
    this.a1[i] = a1;
  }
  reset(): void {
    this.n = 0;
  }
  append(b: SectorBatch, start: number, count: number): void {
    if (this.n + count > this.cx.length) this.grow(this.n + count);
    for (let k = 0; k < count; k++) {
      const i = this.n++;
      const j = start + k;
      this.cx[i] = b.cx[j];
      this.cy[i] = b.cy[j];
      this.r0[i] = b.r0[j];
      this.r1[i] = b.r1[j];
      this.a0[i] = b.a0[j];
      this.a1[i] = b.a1[j];
      this.color[i] = b.color[j];
      this.alpha[i] = b.alpha[j];
    }
  }
  private grow(min: number): void {
    const m = Math.max(min, this.cx.length * 2);
    this.cx = grown(this.cx, m);
    this.cy = grown(this.cy, m);
    this.r0 = grown(this.r0, m);
    this.r1 = grown(this.r1, m);
    this.a0 = grown(this.a0, m);
    this.a1 = grown(this.a1, m);
    this.color = grown(this.color, m);
    this.alpha = grown(this.alpha, m);
  }
}

/**
 * Subpaths of verbs over coordinates (device px). `begin` starts a subpath with its style; consecutive subpaths of one
 * style may be drawn as one path, so overlapping translucent fills of one style do not darken where they meet. A
 * stroke may be dashed (`dashed`): `dash` px on, `gap` px off (device px; 0: solid), the pattern starting `phase` px
 * into it at the subpath's first point.
 */
export class PathBatch {
  /** Subpaths, and per subpath its first verb and coordinate, color, alpha, mode and width (half device px). */
  n = 0;
  start: Int32Array;
  cstart: Int32Array;
  color: Int32Array;
  alpha: Uint8Array;
  mode: Uint8Array;
  width: Uint8Array;
  dash: Uint16Array;
  gap: Uint16Array;
  phase: Float64Array;
  nv = 0;
  verbs: Uint8Array;
  nc = 0;
  coords: Float64Array;
  constructor(capacity = 64) {
    this.start = new Int32Array(capacity);
    this.cstart = new Int32Array(capacity);
    this.color = new Int32Array(capacity);
    this.alpha = new Uint8Array(capacity);
    this.mode = new Uint8Array(capacity);
    this.width = new Uint8Array(capacity);
    this.dash = new Uint16Array(capacity);
    this.gap = new Uint16Array(capacity);
    this.phase = new Float64Array(capacity);
    this.verbs = new Uint8Array(capacity * 4);
    this.coords = new Float64Array(capacity * 8);
  }
  /** A new subpath; `width` in half device px (2: a 1 px stroke). */
  begin(color: number, alpha = 255, mode: number = PathMode.fill, width = 2): void {
    if (this.n === this.start.length) this.growPaths(this.n + 1);
    const p = this.n++;
    this.start[p] = this.nv;
    this.cstart[p] = this.nc;
    this.color[p] = color;
    this.alpha[p] = alpha;
    this.mode[p] = mode;
    this.width[p] = width;
    this.dash[p] = this.gap[p] = 0;
    this.phase[p] = 0;
  }
  /** The current subpath's stroke dashed: `dash` px on, `gap` off, from `phase` px into the pattern (device px). */
  dashed(dash: number, gap: number, phase: number): void {
    const p = this.n - 1;
    this.dash[p] = dash;
    this.gap[p] = gap;
    this.phase[p] = phase;
  }
  /** A verb added: the offset of its coordinates in `coords`, to write them into. */
  verb(v: number): number {
    const k = PathVerbCoords[v];
    // At least one more: a batch made with no room grows too (VL-20).
    if (this.nv === this.verbs.length) this.verbs = grown(this.verbs, Math.max(this.nv + 1, this.nv * 2));
    if (this.nc + k > this.coords.length) this.coords = grown(this.coords, Math.max(this.nc + k, this.nc * 2));
    this.verbs[this.nv++] = v;
    const o = this.nc;
    this.nc += k;
    return o;
  }
  moveTo(x: number, y: number): void {
    const o = this.verb(PathVerb.move);
    this.coords[o] = x;
    this.coords[o + 1] = y;
  }
  lineTo(x: number, y: number): void {
    const o = this.verb(PathVerb.line);
    this.coords[o] = x;
    this.coords[o + 1] = y;
  }
  quadTo(cx: number, cy: number, x: number, y: number): void {
    const o = this.verb(PathVerb.quad);
    this.coords[o] = cx;
    this.coords[o + 1] = cy;
    this.coords[o + 2] = x;
    this.coords[o + 3] = y;
  }
  close(): void {
    this.verb(PathVerb.close);
  }
  reset(): void {
    this.n = this.nv = this.nc = 0;
  }
  /** Subpaths [start, start + count) of `b` added with their verbs and coordinates. */
  append(b: PathBatch, start: number, count: number): void {
    if (count === 0) return;
    if (this.n + count > this.start.length) this.growPaths(this.n + count);
    const v0 = b.start[start];
    const c0 = b.cstart[start];
    const v1 = start + count < b.n ? b.start[start + count] : b.nv;
    const c1 = start + count < b.n ? b.cstart[start + count] : b.nc;
    for (let k = 0; k < count; k++) {
      const p = this.n++;
      const j = start + k;
      this.start[p] = b.start[j] - v0 + this.nv;
      this.cstart[p] = b.cstart[j] - c0 + this.nc;
      this.color[p] = b.color[j];
      this.alpha[p] = b.alpha[j];
      this.mode[p] = b.mode[j];
      this.width[p] = b.width[j];
      this.dash[p] = b.dash[j];
      this.gap[p] = b.gap[j];
      this.phase[p] = b.phase[j];
    }
    if (this.nv + v1 - v0 > this.verbs.length) this.verbs = grown(this.verbs, (this.nv + v1 - v0) * 2);
    if (this.nc + c1 - c0 > this.coords.length) this.coords = grown(this.coords, (this.nc + c1 - c0) * 2);
    // Copied by loops: a `subarray` view would be an allocation per call.
    for (let v = v0; v < v1; v++) this.verbs[this.nv++] = b.verbs[v];
    for (let c = c0; c < c1; c++) this.coords[this.nc++] = b.coords[c];
  }
  private growPaths(min: number): void {
    const m = Math.max(min, this.start.length * 2);
    this.start = grown(this.start, m);
    this.cstart = grown(this.cstart, m);
    this.color = grown(this.color, m);
    this.alpha = grown(this.alpha, m);
    this.mode = grown(this.mode, m);
    this.width = grown(this.width, m);
    this.dash = grown(this.dash, m);
    this.gap = grown(this.gap, m);
    this.phase = grown(this.phase, m);
  }
}

/** Glyph marks centered at (x, y), `radius` device px; a character glyph in the batch's `font`, 2r high. */
export class GlyphBatch {
  n = 0;
  /** The FontId character glyphs are drawn in, scaled to their radius (set once per draw). */
  font = 0;
  x: Float64Array;
  y: Float64Array;
  radius: Float64Array;
  glyph: Int32Array;
  color: Int32Array;
  alpha: Uint8Array;
  private c = 0;
  private a = 255;
  constructor(capacity = 64) {
    this.x = new Float64Array(capacity);
    this.y = new Float64Array(capacity);
    this.radius = new Float64Array(capacity);
    this.glyph = new Int32Array(capacity);
    this.color = new Int32Array(capacity);
    this.alpha = new Uint8Array(capacity);
  }
  begin(color: number, alpha = 255): void {
    this.c = color;
    this.a = alpha;
  }
  /** A new item of glyph `id` in the current color: its slot, to write x, y and radius into. */
  add(id: number): number {
    if (this.n === this.x.length) this.grow(this.n + 1);
    const i = this.n++;
    this.glyph[i] = id;
    this.color[i] = this.c;
    this.alpha[i] = this.a;
    return i;
  }
  push(x: number, y: number, radius: number, id: number): void {
    const i = this.add(id);
    this.x[i] = x;
    this.y[i] = y;
    this.radius[i] = radius;
  }
  reset(): void {
    this.n = 0;
  }
  append(b: GlyphBatch, start: number, count: number): void {
    if (this.n + count > this.x.length) this.grow(this.n + count);
    this.font = b.font;
    for (let k = 0; k < count; k++) {
      const i = this.n++;
      const j = start + k;
      this.x[i] = b.x[j];
      this.y[i] = b.y[j];
      this.radius[i] = b.radius[j];
      this.glyph[i] = b.glyph[j];
      this.color[i] = b.color[j];
      this.alpha[i] = b.alpha[j];
    }
  }
  private grow(min: number): void {
    const m = Math.max(min, this.x.length * 2);
    this.x = grown(this.x, m);
    this.y = grown(this.y, m);
    this.radius = grown(this.radius, m);
    this.glyph = grown(this.glyph, m);
    this.color = grown(this.color, m);
    this.alpha = grown(this.alpha, m);
  }
}

/** A one-character string per 7-bit code (`CHARS[c & 127]`): a letter's text item without a string made per frame. */
export const CHARS: readonly string[] = Array.from({ length: 128 }, (_, i) => String.fromCharCode(i));

/**
 * Text runs in one font, alignment, baseline and angle (`style`); per item its anchor, string and color, and its scale
 * about the anchor (`sx` along the run, `sy` across it, before the angle; 1 from `add`): a logo's letters stretched.
 */
export class TextBatch {
  n = 0;
  font = 0;
  align: number = TextAlign.left;
  baseline: number = TextBaseline.alphabetic;
  angle = 0;
  x: Float64Array;
  y: Float64Array;
  sx: Float64Array;
  sy: Float64Array;
  /** References to strings the draw holds (prepared labels, a table of characters): never built per draw. */
  text: string[];
  color: Int32Array;
  alpha: Uint8Array;
  private c = 0;
  private a = 255;
  constructor(capacity = 256) {
    this.x = new Float64Array(capacity);
    this.y = new Float64Array(capacity);
    this.sx = new Float64Array(capacity);
    this.sy = new Float64Array(capacity);
    this.text = new Array<string>(capacity).fill('');
    this.color = new Int32Array(capacity);
    this.alpha = new Uint8Array(capacity);
  }
  /** The batch's font (FontId), alignment, baseline and angle (radians, clockwise): one per batch. */
  style(font: number, align: number, baseline: number, angle = 0): void {
    this.font = font;
    this.align = align;
    this.baseline = baseline;
    this.angle = angle;
  }
  begin(color: number, alpha = 255): void {
    this.c = color;
    this.a = alpha;
  }
  /** A new run of `text` in the current color, unscaled: its slot, to write x and y (and a scale) into. */
  add(text: string): number {
    if (this.n === this.x.length) this.grow(this.n + 1);
    const i = this.n++;
    this.sx[i] = this.sy[i] = 1;
    this.text[i] = text;
    this.color[i] = this.c;
    this.alpha[i] = this.a;
    return i;
  }
  push(x: number, y: number, text: string): void {
    const i = this.add(text);
    this.x[i] = x;
    this.y[i] = y;
  }
  reset(): void {
    this.n = 0;
  }
  append(b: TextBatch, start: number, count: number): void {
    if (this.n + count > this.x.length) this.grow(this.n + count);
    for (let k = 0; k < count; k++) {
      const i = this.n++;
      const j = start + k;
      this.x[i] = b.x[j];
      this.y[i] = b.y[j];
      this.sx[i] = b.sx[j];
      this.sy[i] = b.sy[j];
      this.text[i] = b.text[j];
      this.color[i] = b.color[j];
      this.alpha[i] = b.alpha[j];
    }
  }
  private grow(min: number): void {
    const m = Math.max(min, this.x.length * 2);
    this.x = grown(this.x, m);
    this.y = grown(this.y, m);
    this.sx = grown(this.sx, m);
    this.sy = grown(this.sy, m);
    while (this.text.length < m) this.text.push('');
    this.color = grown(this.color, m);
    this.alpha = grown(this.alpha, m);
  }
}
