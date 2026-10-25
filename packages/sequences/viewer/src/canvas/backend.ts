/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Canvas2D backend (seqspec:viewer#the-render-backend): the live path, over `Paint2D`, the subset of Canvas2D that
// CanvasRenderingContext2D and OffscreenCanvasRenderingContext2D share. Device px under the identity transform. Batches
// are grouped by (color, alpha) with a counting sort over their distinct keys (`ColorGroups`, as the SVG backend groups
// them), one path per group, so style switches are at most the distinct keys. Text is placed by `setTransform` and
// drawn at (0, 0) (V8 boxes fractional arguments to fillText, drawImage and roundRect, not to fillRect, rect, arc,
// setTransform, clip or clearRect), with the batch's font, alignment and baseline set here, so no draw inherits
// another's. A back pass draws into its own offscreen canvas and shows only what is presented, at `endPass`.
// What a context's backends share, `Canvas2DShared`: one scratch canvas for copies and stale images, the images, the
// path glyphs, each a Path2D built once, and the glyph atlas: a one-character text item drawn upright at unit scale (a
// letter) is a `drawImage` of a cell rasterized once, not a `setTransform` and a `fillText`.

import {
  BackendCaps,
  type BackendImages,
  type DrawResources,
  type GlyphBatch,
  type GlyphRecord,
  NO_COLOR,
  Pass,
  type PathBatch,
  PathMode,
  PathVerb,
  PathVerbCoords,
  type RectBatch,
  type RenderBackend,
  type SectorBatch,
  SHAPE_COUNT,
  type TextBatch,
} from '../draw/backend.js';
import { ColorGroups } from '../draw/group.js';

/** An opaque stroke's items per path (`paths`): Chrome's closePath slows as a path grows. */
const STROKE_CHUNK = 256;

export type Paint2D = Pick<
  CanvasRenderingContext2D,
  | 'fillStyle'
  | 'strokeStyle'
  | 'globalAlpha'
  | 'lineWidth'
  | 'setLineDash'
  | 'lineDashOffset'
  | 'font'
  | 'textAlign'
  | 'textBaseline'
  | 'imageSmoothingEnabled'
  | 'setTransform'
  | 'clearRect'
  | 'beginPath'
  | 'closePath'
  | 'rect'
  | 'moveTo'
  | 'lineTo'
  | 'quadraticCurveTo'
  | 'bezierCurveTo'
  | 'arc'
  | 'arcTo'
  | 'fill'
  | 'stroke'
  | 'clip'
  | 'save'
  | 'restore'
  | 'fillText'
  | 'drawImage'
> & { readonly canvas: HTMLCanvasElement | OffscreenCanvas };

const ALIGN: CanvasTextAlign[] = ['left', 'center', 'right'];
const BASELINE: CanvasTextBaseline[] = ['alphabetic', 'middle', 'top', 'bottom'];
const ALPHA = Float64Array.from({ length: 256 }, (_, i) => i / 255);

/**
 * A packed color's CSS string, made once per color and kept, never per draw: a frame under a continuous palette holds
 * thousands of colors, and a direct-mapped table evicted colliding ones every frame. A lookup by an int key allocates
 * nothing; past `CSS_MAX` colors the table starts again (edge rate).
 */
const CSS_MAX = 1 << 16;
const cssText = new Map<number, string>();
function css(c: number): string {
  let s = cssText.get(c);
  if (s === undefined) {
    if (cssText.size >= CSS_MAX) cssText.clear();
    s = `#${c.toString(16).padStart(6, '0')}`;
    cssText.set(c, s);
  }
  return s;
}

const surface = (w: number, h: number) => {
  const c = new OffscreenCanvas(Math.max(1, w), Math.max(1, h));
  return { c, g: c.getContext('2d') as OffscreenCanvasRenderingContext2D };
};

export interface SharedOptions {
  /** False: every text item through `fillText` (the atlas's pixel tests draw their reference so). */
  atlas?: boolean;
  /** The context's text version (bumped by a font load): the glyph atlas is dropped when it moves. */
  fonts?: () => number;
}

/** One per context: its draw resources, the scratch canvas, the images, the path glyphs and the glyph atlas. */
export class Canvas2DShared implements BackendImages {
  readonly scratch = surface(1, 1);
  dpr = 1;
  /** Letters' cells (null: disabled). */
  readonly atlas: GlyphAtlas | null;
  private readonly bitmaps: (ImageBitmap | null)[] = [];
  private readonly free: number[] = [];
  private pixels: { c: OffscreenCanvas; g: OffscreenCanvasRenderingContext2D } | null = null;
  private layer: { c: OffscreenCanvas; g: OffscreenCanvasRenderingContext2D; b: Canvas2DBackend } | null = null;
  private readonly glyphPaths: (Path2D | null)[] = [];
  private readonly glyphRecs: (GlyphRecord | null)[] = [];
  readonly resources: DrawResources;
  constructor(resources: DrawResources, opts: SharedOptions = {}) {
    this.resources = resources;
    this.atlas = opts.atlas === false ? null : new GlyphAtlas(resources, opts.fonts ?? null);
  }

  bitmap(img: number): ImageBitmap | null {
    return this.bitmaps[img] ?? null;
  }
  fromPixels(w: number, h: number, rgba: Uint32Array): number {
    if (!this.pixels) this.pixels = surface(w, h);
    const p = this.pixels;
    if (p.c.width !== w || p.c.height !== h) {
      p.c.width = w;
      p.c.height = h;
    }
    p.g.putImageData(
      new ImageData(new Uint8ClampedArray(rgba.buffer as ArrayBuffer, rgba.byteOffset, w * h * 4), w, h),
      0,
      0,
    );
    return this.keep(p.c.transferToImageBitmap());
  }
  beginLayer(w: number, h: number): RenderBackend {
    if (!this.layer) {
      const s = surface(w, h);
      this.layer = { ...s, b: new Canvas2DBackend([s.g], this, 0) };
    }
    const l = this.layer;
    l.c.width = Math.max(1, w);
    l.c.height = Math.max(1, h);
    l.b.beginPass(Pass.main, w, h, this.dpr);
    return l.b;
  }
  endLayer(): number {
    const l = this.layer as NonNullable<Canvas2DShared['layer']>;
    l.b.endPass();
    return this.keep(l.c.transferToImageBitmap());
  }
  dispose(img: number): void {
    const b = this.bitmaps[img];
    if (!b) return;
    b.close();
    this.bitmaps[img] = null;
    this.free.push(img);
  }
  /** A path glyph's Path2D, built once per glyph record. */
  glyphPath(id: number, rec: GlyphRecord): Path2D {
    let p = this.glyphPaths[id];
    if (!p || this.glyphRecs[id] !== rec) {
      while (this.glyphPaths.length <= id) {
        this.glyphPaths.push(null);
        this.glyphRecs.push(null);
      }
      this.glyphPaths[id] = p = new Path2D(rec.path ?? '');
      this.glyphRecs[id] = rec;
    }
    return p;
  }
  private keep(b: ImageBitmap): number {
    const id = this.free.length > 0 ? (this.free.pop() as number) : this.bitmaps.length;
    this.bitmaps[id] = b;
    return id;
  }
}

/** The style a group or subpath sets: cached, so a context switch is made only on a change. */
const UNSET = -0x80000000;
/** A solid stroke's dash pattern. */
const NO_DASH: number[] = [];

export class Canvas2DBackend implements RenderBackend {
  readonly caps = BackendCaps.Copy | BackendCaps.Progressive | BackendCaps.RotatedText | BackendCaps.Back;
  readonly images: Canvas2DShared;
  /** Fill and stroke style switches: at most the distinct colors of each batch. */
  styles = 0;
  private g: Paint2D;
  private view: Paint2D | null = null;
  private dpr = 1;
  private readonly backs: ({ c: OffscreenCanvas; g: OffscreenCanvasRenderingContext2D } | null)[] = [];
  private copied = false;
  private presents = new Int32Array(64);
  private np = 0;
  private stales = new Float64Array(64);
  private ns = 0;
  // the cached context state
  private fillC = UNSET;
  private strokeC = UNSET;
  private alphaV = -1;
  private widthV = -1;
  private fontV = '';
  private alignV = -1;
  private baseV = -1;
  private smoothV = -1;
  private readonly num = new Float64Array(1);
  /** A dashed stroke's pattern, reused (`setLineDash` copies it). */
  private readonly dashArr = [0, 0];
  /** Batches grouped by (color, alpha): `group`'s order and starts are read after each grouping (they can grow). */
  private readonly groups = new ColorGroups();

  /**
   * `targets`: the visible context of each pass (by `Pass`); `back`: a bit per pass (1 << pass) that draws into a
   * back target and shows only what is presented.
   */
  private readonly targets: readonly (Paint2D | null)[];
  private readonly back: number;
  constructor(targets: readonly (Paint2D | null)[], shared: Canvas2DShared, back = 1 << Pass.main) {
    this.targets = targets;
    this.back = back;
    this.images = shared;
    this.g = targets[0] as Paint2D;
  }

  beginPass(pass: number, width: number, height: number, dpr: number): void {
    const view = this.targets[pass];
    if (!view) throw new Error(`Canvas2DBackend: no target for pass ${pass}.`);
    this.dpr = this.images.dpr = dpr;
    if (this.back & (1 << pass)) {
      let b = this.backs[pass];
      if (!b) {
        while (this.backs.length <= pass) this.backs.push(null);
        this.backs[pass] = b = surface(width, height);
      }
      if (b.c.width !== width || b.c.height !== height) {
        b.c.width = width;
        b.c.height = height;
      }
      this.g = b.g;
      this.view = view;
    } else {
      this.g = view;
      this.view = null;
    }
    this.copied = false;
    this.np = this.ns = 0;
    this.unknown();
    this.g.setTransform(1, 0, 0, 1, 0, 0);
  }
  /** A back pass shows its stale regions, then its presented ones over them. */
  endPass(): void {
    const v = this.view;
    if (!v || (this.np === 0 && this.ns === 0)) return;
    if (this.ns > 0) {
      const s = this.snapshot(v.canvas);
      v.imageSmoothingEnabled = true;
      const r = this.stales;
      for (let k = 0; k < this.ns * 8; k += 8) {
        v.clearRect(r[k + 4], r[k + 5], r[k + 6], r[k + 7]);
        v.drawImage(s, r[k], r[k + 1], r[k + 2], r[k + 3], r[k + 4], r[k + 5], r[k + 6], r[k + 7]);
      }
    }
    const src = this.g.canvas;
    const p = this.presents;
    for (let k = 0; k < this.np * 4; k += 4) {
      v.clearRect(p[k], p[k + 1], p[k + 2], p[k + 3]);
      v.drawImage(src, p[k], p[k + 1], p[k + 2], p[k + 3], p[k], p[k + 1], p[k + 2], p[k + 3]);
    }
    this.np = this.ns = 0;
  }
  clear(x: number, y: number, w: number, h: number): void {
    this.g.clearRect(x, y, w, h);
  }
  pushClip(x: number, y: number, w: number, h: number): void {
    const g = this.g;
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
  }
  popClip(): void {
    this.g.restore();
    this.unknown();
  }
  copy(dx: number, dy: number, x: number, y: number, w: number, h: number): void {
    const g = this.g;
    if (!this.copied) {
      this.snapshot(g.canvas);
      this.copied = true;
    }
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    g.clearRect(x, y, w, h);
    g.drawImage(this.images.scratch.c, dx, dy);
    g.restore();
    this.unknown();
  }
  present(x: number, y: number, w: number, h: number): void {
    if (!this.view) return;
    if (this.np * 4 === this.presents.length) this.presents = grow(this.presents);
    const p = this.presents;
    const k = this.np++ * 4;
    p[k] = x;
    p[k + 1] = y;
    p[k + 2] = w;
    p[k + 3] = h;
  }
  stale(sx: number, sy: number, sw: number, sh: number, dx: number, dy: number, dw: number, dh: number): void {
    if (!this.view) return;
    if (this.ns * 8 === this.stales.length) this.stales = grow(this.stales);
    const r = this.stales;
    const k = this.ns++ * 8;
    r[k] = sx;
    r[k + 1] = sy;
    r[k + 2] = sw;
    r[k + 3] = sh;
    r[k + 4] = dx;
    r[k + 5] = dy;
    r[k + 6] = dw;
    r[k + 7] = dh;
  }

  /** Fills, and outlines as four filled edges inside the rect: whole device px at every width, grouped with fills. */
  rects(b: RectBatch): void {
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    const g = this.g;
    const o = G.order;
    const st = G.start;
    for (let q = 0; q < groups; q++) {
      g.beginPath();
      for (let j = st[q]; j < st[q + 1]; j++) {
        const i = o[j];
        const s = b.stroke[i];
        if (s === 0 || 2 * s >= b.w[i] || 2 * s >= b.h[i]) g.rect(b.x[i], b.y[i], b.w[i], b.h[i]);
        else {
          g.rect(b.x[i], b.y[i], b.w[i], s);
          g.rect(b.x[i], b.y[i] + b.h[i] - s, b.w[i], s);
          g.rect(b.x[i], b.y[i] + s, s, b.h[i] - 2 * s);
          g.rect(b.x[i] + b.w[i] - s, b.y[i] + s, s, b.h[i] - 2 * s);
        }
      }
      this.fillStyle(b.color[o[st[q]]], b.alpha[o[st[q]]]);
      g.fill();
    }
  }
  sectors(b: SectorBatch): void {
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    const g = this.g;
    const o = G.order;
    const st = G.start;
    for (let q = 0; q < groups; q++) {
      g.beginPath();
      for (let j = st[q]; j < st[q + 1]; j++) {
        const i = o[j];
        g.moveTo(b.cx[i] + b.r1[i] * Math.cos(b.a0[i]), b.cy[i] + b.r1[i] * Math.sin(b.a0[i]));
        g.arc(b.cx[i], b.cy[i], b.r1[i], b.a0[i], b.a1[i], false);
        g.arc(b.cx[i], b.cy[i], b.r0[i], b.a1[i], b.a0[i], true);
        // No closePath: the fill closes each sector, and Chrome's slows as the path grows (`paths`).
      }
      this.fillStyle(b.color[o[st[q]]], b.alpha[o[st[q]]]);
      g.fill();
    }
  }
  paths(b: PathBatch): void {
    const g = this.g;
    const k = b.coords;
    for (let p = 0; p < b.n; p++) {
      if (b.color[p] === NO_COLOR) continue;
      g.beginPath();
      let c = b.cstart[p];
      // Chrome's closePath costs as much as the path holds so far (16,000 closed subpaths in one path: 574 ms, 0.7 ms
      // unclosed), so a path never closes what it need not: a fill closes every subpath itself, so a close before the
      // next move or at the path's end is left out; an opaque stroke is drawn every `STROKE_CHUNK` items.
      const fill = b.mode[p] === PathMode.fill;
      const chunked = !fill && b.alpha[p] === 255;
      // Subpaths of one style are one path: drawn when the next differs.
      let q = p;
      for (; q < b.n && (q === p || this.sameStyle(b, p, q)); q++) {
        if (chunked && q > p && (q - p) % STROKE_CHUNK === 0) {
          this.drawPath(b, p);
          g.beginPath();
        }
        const end = q + 1 < b.n ? b.start[q + 1] : b.nv;
        for (let v = b.start[q]; v < end; v++) {
          switch (b.verbs[v]) {
            case PathVerb.move:
              g.moveTo(k[c], k[c + 1]);
              break;
            case PathVerb.line:
              g.lineTo(k[c], k[c + 1]);
              break;
            case PathVerb.quad:
              g.quadraticCurveTo(k[c], k[c + 1], k[c + 2], k[c + 3]);
              break;
            case PathVerb.cubic:
              g.bezierCurveTo(k[c], k[c + 1], k[c + 2], k[c + 3], k[c + 4], k[c + 5]);
              break;
            case PathVerb.arc:
              g.arc(k[c], k[c + 1], k[c + 2], k[c + 3], k[c + 4], k[c + 5] !== 0);
              break;
            case PathVerb.roundRect: {
              const x = k[c];
              const y = k[c + 1];
              const r = k[c + 4];
              g.moveTo(x + r, y);
              g.arcTo(x + k[c + 2], y, x + k[c + 2], y + k[c + 3], r);
              g.arcTo(x + k[c + 2], y + k[c + 3], x, y + k[c + 3], r);
              g.arcTo(x, y + k[c + 3], x, y, r);
              g.arcTo(x, y, x + k[c + 2], y, r);
              if (!fill || !this.closedByFill(b, v, end, q, p)) g.closePath();
              break;
            }
            default:
              if (!fill || !this.closedByFill(b, v, end, q, p)) g.closePath();
          }
          c += PathVerbCoords[b.verbs[v]];
        }
      }
      this.drawPath(b, p);
      p = q - 1;
    }
  }
  /**
   * Whether a fill closes the subpath that verb v (of item q, `end` its items' end; p the path's first item) closes: the
   * path ends there, or the next verb starts a subpath (a move, a round rect), as the next item's first verb does.
   */
  private closedByFill(b: PathBatch, v: number, end: number, q: number, p: number): boolean {
    let next = v + 1;
    if (next === end) {
      if (q + 1 >= b.n || !this.sameStyle(b, p, q + 1)) return true;
      next = b.start[q + 1];
    }
    const verb = b.verbs[next];
    return verb === PathVerb.move || verb === PathVerb.roundRect;
  }
  /** The path built so far, in item p's style. */
  private drawPath(b: PathBatch, p: number): void {
    const g = this.g;
    if (b.mode[p] === PathMode.stroke && b.dash[p] > 0) {
      // Dashed: the pattern set for this stroke alone (setLineDash copies the scratch), then solid again.
      this.strokeStyle(b.color[p], b.alpha[p], b.width[p] / 2);
      this.dashArr[0] = b.dash[p];
      this.dashArr[1] = b.gap[p];
      g.setLineDash(this.dashArr);
      g.lineDashOffset = b.phase[p];
      g.stroke();
      g.setLineDash(NO_DASH);
      g.lineDashOffset = 0;
    } else if (b.mode[p] === PathMode.stroke) {
      this.strokeStyle(b.color[p], b.alpha[p], b.width[p] / 2);
      g.stroke();
    } else {
      this.fillStyle(b.color[p], b.alpha[p]);
      g.fill();
    }
  }
  glyphs(b: GlyphBatch): void {
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    const g = this.g;
    const o = G.order;
    const st = G.start;
    const res = this.images.resources;
    // A double from a call V8 does not inline is boxed: read once per batch, into a typed slot.
    this.num[0] = res.fontPx(b.font);
    for (let q = 0; q < groups; q++) {
      g.beginPath();
      this.fillStyle(b.color[o[st[q]]], b.alpha[o[st[q]]]);
      for (let j = st[q]; j < st[q + 1]; j++) {
        const i = o[j];
        const id = b.glyph[i];
        const rec = id < SHAPE_COUNT ? null : res.glyph(id);
        if (rec?.char) this.charGlyph(b, i, rec.char);
        else if (rec?.path) {
          // Scaled so the box's larger side is 2r, centered on the point (seqspec:format#themes).
          const s = (2 * b.radius[i]) / Math.max(rec.box[0], rec.box[1]);
          g.setTransform(s, 0, 0, s, b.x[i] - (rec.box[0] * s) / 2, b.y[i] - (rec.box[1] * s) / 2);
          g.fill(this.images.glyphPath(id, rec));
          g.setTransform(1, 0, 0, 1, 0, 0);
        } else shape(g, id < SHAPE_COUNT ? id : rec ? rec.fallback : 0, b, i);
      }
      g.fill();
    }
  }
  text(b: TextBatch): void {
    if (b.n === 0) return;
    this.textStyle(b.font, b.align, b.baseline);
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    const g = this.g;
    const o = G.order;
    const st = G.start;
    const cos = Math.cos(b.angle) * this.dpr;
    const sin = Math.sin(b.angle) * this.dpr;
    // Upright batches through the atlas: its letters (one character, unit scale) from cells, the rest by fillText.
    const at = b.angle === 0 ? this.images.atlas : null;
    if (at !== null) {
      at.px[0] = this.dpr;
      at.begin(b.font, b.align, b.baseline);
      // Every cell the batch lacks made first, so a page changed is uploaded once a batch, not once a new letter.
      if (at.prepass)
        for (let q = 0; q < groups; q++) {
          at.color(b.color[o[st[q]]]);
          for (let j = st[q]; j < st[q + 1]; j++) {
            const i = o[j];
            if (b.sx[i] === 1 && b.sy[i] === 1 && b.text[i].length === 1) at.cellOf(b, i);
          }
        }
    }
    let moved = false;
    for (let q = 0; q < groups; q++) {
      this.fillStyle(b.color[o[st[q]]], b.alpha[o[st[q]]]);
      if (at !== null) at.color(b.color[o[st[q]]]);
      for (let j = st[q]; j < st[q + 1]; j++) {
        const i = o[j];
        const sx = b.sx[i];
        const sy = b.sy[i];
        if (at !== null && sx === 1 && sy === 1 && b.text[i].length === 1) {
          if (moved) {
            g.setTransform(1, 0, 0, 1, 0, 0);
            moved = false;
          }
          if (at.draw(g, b, i)) continue;
        }
        g.setTransform(cos * sx, sin * sx, -sin * sy, cos * sy, b.x[i], b.y[i]);
        g.fillText(b.text[i], 0, 0);
        moved = true;
      }
    }
    g.setTransform(1, 0, 0, 1, 0, 0);
  }
  image(img: number, r: Float64Array, smooth: boolean): void {
    const bmp = this.images.bitmap(img);
    if (!bmp) return;
    const g = this.g;
    if (+smooth !== this.smoothV) {
      g.imageSmoothingEnabled = smooth;
      this.smoothV = +smooth;
    }
    if (this.alphaV !== 255) {
      g.globalAlpha = 1;
      this.alphaV = 255;
    }
    g.drawImage(bmp, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
  }

  /** A character glyph in the batch's font, 2r high, centered on its point. */
  private charGlyph(b: GlyphBatch, i: number, ch: string): void {
    const g = this.g;
    this.textStyle(b.font, 1, 1);
    const s = (2 * b.radius[i]) / this.num[0];
    g.setTransform(s, 0, 0, s, b.x[i], b.y[i]);
    g.fillText(ch, 0, 0);
    g.setTransform(1, 0, 0, 1, 0, 0);
  }
  private sameStyle(b: PathBatch, p: number, q: number): boolean {
    return (
      b.color[q] === b.color[p] &&
      b.alpha[q] === b.alpha[p] &&
      b.mode[q] === b.mode[p] &&
      b.width[q] === b.width[p] &&
      b.dash[q] === b.dash[p] &&
      b.gap[q] === b.gap[p] &&
      b.phase[q] === b.phase[p]
    );
  }
  /** The whole of `src` copied into the scratch canvas. */
  private snapshot(src: HTMLCanvasElement | OffscreenCanvas): OffscreenCanvas {
    const s = this.images.scratch;
    if (s.c.width < src.width || s.c.height < src.height) {
      s.c.width = Math.max(s.c.width, src.width);
      s.c.height = Math.max(s.c.height, src.height);
    } else s.g.clearRect(0, 0, src.width, src.height);
    s.g.drawImage(src, 0, 0);
    return s.c;
  }
  private color(c: number): string {
    return c < NO_COLOR ? this.images.resources.token(c) : css(c);
  }
  private alpha(a: number): void {
    if (a === this.alphaV) return;
    this.g.globalAlpha = ALPHA[a];
    this.alphaV = a;
  }
  private fillStyle(c: number, a: number): void {
    this.alpha(a);
    if (c === this.fillC) return;
    this.g.fillStyle = this.color(c);
    this.fillC = c;
    this.styles++;
  }
  private strokeStyle(c: number, a: number, width: number): void {
    this.alpha(a);
    if (width !== this.widthV) {
      this.g.lineWidth = width;
      this.widthV = width;
    }
    if (c === this.strokeC) return;
    this.g.strokeStyle = this.color(c);
    this.strokeC = c;
    this.styles++;
  }
  private textStyle(font: number, align: number, baseline: number): void {
    const f = this.images.resources.font(font);
    if (f !== this.fontV) {
      this.g.font = f;
      this.fontV = f;
    }
    if (align !== this.alignV) {
      this.g.textAlign = ALIGN[align];
      this.alignV = align;
    }
    if (baseline !== this.baseV) {
      this.g.textBaseline = BASELINE[baseline];
      this.baseV = baseline;
    }
  }
  /** After a restore or a new target: every cached style is unknown. */
  private unknown(): void {
    this.fillC = this.strokeC = UNSET;
    this.alphaV = this.widthV = this.alignV = this.baseV = this.smoothV = -1;
    this.fontV = '';
  }
}

/** Atlas pages: device px a side, at most this many, made on first use. */
const PAGE = 1024;
const PAGES = 4;
/** Cell sides (device px), each page's slots one of them; a glyph needing more than the last is drawn by fillText. */
const SIDES = [24, 32, 48, 64, 96, 128];
/** Slots a page holds at most (the smallest side), cells in all, and the hash table (a power of two, under half full). */
const SLOTS = Math.floor(PAGE / SIDES[0]) ** 2;
const CELLS = PAGES * SLOTS;
const TABLE = 1 << 14;
/**
 * A cell's margin (device px), and where the anchor sits in it, a share of its side: a left or top anchor `EDGE` in
 * from its edge (a glyph may cross the em box's), an alphabetic baseline at `ALPHABETIC` down, a center in the middle.
 */
const PAD = 2;
const EDGE = 0.2;
const ALPHABETIC = 0.7;

/**
 * The glyph atlas (Canvas2D): one-character items drawn upright at unit scale come from cells, each rasterized once by
 * `fillText` with the item's font, color, alignment and baseline, and drawn with one `drawImage` of whole device px
 * (CS2 measured 2.7 µs a letter as a `setTransform` and a `fillText`). A cell is keyed by font, DPR
 * (so the device font size), color, alignment, baseline, character and the anchor's quarter-pixel phase on each axis:
 * rasterized at its phase and drawn at the anchor's whole pixel, its pixels are `fillText`'s at the anchor (exactly at
 * full alpha; translucent, a channel at most 1 apart: `_test/canvas/atlas.test.ts`). Anchors off the quarter pixels
 * (none of the letters') go to `fillText`. Cells sit in square slots of a few sizes on at most `PAGES` pages; a
 * full atlas evicts the page drawn least recently, never one the current batch drew (the item is drawn by `fillText`
 * then). A font load (the context's text version) or a skin token's new color drops every cell; a DPR is part of the
 * key. A draw allocates nothing: keys, slots and places in typed arrays, pages made on first use (edge rate).
 *
 * A batch's missing cells are made before any of its items is drawn (`prepass`): a page is rasterized on the CPU and
 * the target is the GPU's, so a draw from a page changed since its last draw uploads the whole page (4 MB, about
 * 1.5 ms). Made as each item was drawn, every letter new to a frame cost an upload: a zoom reaching a new letter size
 * made some 80 cells in one frame, a 50 to 110 ms frame (case 15, 2026-10-13). Now a page is uploaded at most once a
 * batch that changed it (`uploads` counts them), and the pixels are as before.
 */
export class GlyphAtlas {
  /** Cells rasterized and pages evicted, ever; draws from a page changed since its last draw (each an upload). */
  made = 0;
  evicted = 0;
  uploads = 0;
  /**
   * True: a batch's missing cells made before its draws (`cellOf`); false: each as its item is drawn, an upload after
   * each (the old behavior: the upload guard's positive control).
   */
  prepass = true;
  /** The DPR of the batch in progress (a double passed to a call is boxed: written here by the backend). */
  readonly px = new Float64Array(2);
  private readonly res: DrawResources;
  private readonly fonts: (() => number) | null;
  private fontsAt = -1;
  /** Per token (by -token), the color its cells were made in: another drops the atlas. */
  private readonly tokens: string[] = new Array<string>(32).fill('');
  private readonly pages: ({ c: OffscreenCanvas; g: OffscreenCanvasRenderingContext2D } | null)[] = [];
  /** Per page: its slots' side (0: free), slots used and the batch it was last drawn in. */
  private readonly side = new Int32Array(PAGES);
  private readonly used = new Int32Array(PAGES);
  private readonly stamp = new Int32Array(PAGES);
  /** Per page: changed since it was last drawn from (its next draw uploads it). */
  private readonly changed = new Int32Array(PAGES);
  private clock = 0;
  /** Per cell (page · SLOTS + slot): its key and the anchor's place in the cell (device px). */
  private readonly k0 = new Int32Array(CELLS);
  private readonly k1 = new Int32Array(CELLS);
  private readonly k2 = new Int32Array(CELLS);
  private readonly ax = new Int32Array(CELLS);
  private readonly ay = new Int32Array(CELLS);
  /** Open addressing: a cell per entry (-1: empty); `at` the slot the last miss stopped at. */
  private readonly table = new Int32Array(TABLE).fill(-1);
  private at = 0;
  /** The batch in progress: font, its cell side (-1: too large), the key's style bits, the group's color. */
  private font = 0;
  private cell = 0;
  private bits = 0;
  private ink = 0;
  /** Cell sides by (font, DPR), direct-mapped: the font's px read once per pair, not per batch. */
  private readonly sideKey = new Int32Array(64).fill(-1);
  private readonly sideVal = new Int32Array(64);

  constructor(res: DrawResources, fonts: (() => number) | null) {
    this.res = res;
    this.fonts = fonts;
  }

  /** A batch's style; the DPR in `px[0]`. A font load since the last batch drops every cell. */
  begin(font: number, align: number, baseline: number): void {
    if (this.fonts !== null) {
      const v = this.fonts();
      if (v !== this.fontsAt) {
        if (this.fontsAt >= 0) this.drop();
        this.fontsAt = v;
      }
    }
    this.clock++;
    const dpr = Math.round(this.px[0] * 64);
    this.font = font;
    this.bits = (dpr << 8) | (align << 6) | (baseline << 4);
    const h = (font * 31 + dpr) & 63;
    const key = (font << 12) | dpr;
    if (this.sideKey[h] !== key) {
      // Edge rate: a font's px at this DPR. A glyph spans about an em; its cell a half more and a margin.
      this.px[1] = this.res.fontPx(font) * this.px[0];
      const need = Math.ceil(1.5 * this.px[1]) + 2 * PAD;
      let s = -1;
      for (let k = 0; k < SIDES.length && s < 0; k++) if (SIDES[k] >= need) s = SIDES[k];
      this.sideKey[h] = key;
      this.sideVal[h] = s;
    }
    this.cell = this.sideVal[h];
  }
  /** A group's color; a token whose color changed (a skin) drops every cell. */
  color(c: number): void {
    this.ink = c;
    if (c >= NO_COLOR) return;
    const css = this.res.token(c);
    const t = -c & 31;
    if (this.tokens[t] === css) return;
    if (this.tokens[t] !== '') this.drop();
    this.tokens[t] = css;
  }
  /**
   * Item i's cell (one character), made on a miss, its page kept from eviction this batch; -1: none (too large, an
   * anchor off the quarter pixels, or a full atlas whose pages the batch draws).
   */
  cellOf(b: TextBatch, i: number): number {
    const side = this.cell;
    if (side < 0) return -1;
    const code = b.text[i].charCodeAt(0);
    // Anchors on quarter pixels only (letters sit on whole and half ones): a cell's phase is then the anchor's own, and
    // its pixels are fillText's exactly. Another anchor rounds to another glyph position after the baseline's offset.
    const qx = Math.round(b.x[i] * 4) | 0;
    const qy = Math.round(b.y[i] * 4) | 0;
    if (qx !== b.x[i] * 4 || qy !== b.y[i] * 4) return -1;
    const k0 = (this.font << 16) | code;
    const k2 = this.bits | ((qx & 3) << 2) | (qy & 3);
    let e = this.find(k0, this.ink, k2);
    if (e < 0) e = this.make(b.text[i], k0, k2, side, qx & 3, qy & 3);
    if (e >= 0) this.stamp[(e / SLOTS) | 0] = this.clock;
    return e;
  }
  /** Item i of `b` (one character) onto `g` from its cell, made on a miss; false: not drawn (fillText draws it). */
  draw(g: Paint2D, b: TextBatch, i: number): boolean {
    const e = this.cellOf(b, i);
    if (e < 0) return false;
    const qx = Math.round(b.x[i] * 4) | 0;
    const qy = Math.round(b.y[i] * 4) | 0;
    const p = (e / SLOTS) | 0;
    if (this.changed[p] !== 0) {
      this.changed[p] = 0;
      this.uploads++;
    }
    const s = this.side[p];
    const slot = e - p * SLOTS;
    const n = (PAGE / s) | 0;
    const x0 = (slot % n) * s;
    const y0 = ((slot / n) | 0) * s;
    const src = (this.pages[p] as { c: OffscreenCanvas }).c;
    g.drawImage(src, x0, y0, s, s, (qx >> 2) - this.ax[e], (qy >> 2) - this.ay[e], s, s);
    return true;
  }
  /** Every cell dropped (a font load, a skin); the pages are kept for reuse. */
  drop(): void {
    this.side.fill(0);
    this.used.fill(0);
    this.table.fill(-1);
  }

  private hash(k0: number, k1: number, k2: number): number {
    return (Math.imul(k0, 0x9e3779b1) ^ Math.imul(k1, 0x85ebca77) ^ Math.imul(k2, 0xc2b2ae3d)) >>> 18;
  }
  /** The cell of a key, or -1 with `at` the empty entry it would take. */
  private find(k0: number, k1: number, k2: number): number {
    let h = this.hash(k0, k1, k2);
    const t = this.table;
    for (let e = t[h]; e >= 0; e = t[h]) {
      if (this.k0[e] === k0 && this.k1[e] === k1 && this.k2[e] === k2) return e;
      h = (h + 1) & (TABLE - 1);
    }
    this.at = h;
    return -1;
  }
  /** A cell for the key, rasterized: a slot on a page of its side, a free page, or the least recent page evicted. */
  private make(text: string, k0: number, k2: number, side: number, fx: number, fy: number): number {
    let p = -1;
    const n = (PAGE / side) | 0;
    for (let q = 0; q < PAGES && p < 0; q++) if (this.side[q] === side && this.used[q] < n * n) p = q;
    for (let q = 0; q < PAGES && p < 0; q++) if (this.side[q] === 0) p = q;
    if (p < 0) {
      for (let q = 0; q < PAGES; q++) if (this.stamp[q] < this.clock && (p < 0 || this.stamp[q] < this.stamp[p])) p = q;
      if (p < 0) return -1;
      this.evict(p);
      this.find(k0, this.ink, k2);
    }
    if (this.side[p] === 0) {
      this.side[p] = side;
      this.used[p] = 0;
    }
    while (this.pages.length <= p) this.pages.push(null);
    let pg = this.pages[p];
    if (!pg) {
      // Rasterized on the CPU: a GPU page drawn into the back targets left their translucent pixels 1 to 3 apart
      // between two surfaces of one document (the composer's copy test), a GPU or CPU raster chosen by history.
      const c = new OffscreenCanvas(PAGE, PAGE);
      const g = c.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
      this.pages[p] = pg = { c, g };
    }
    const slot = this.used[p]++;
    const e = p * SLOTS + slot;
    this.changed[p] = 1;
    const x0 = (slot % n) * side;
    const y0 = ((slot / n) | 0) * side;
    const align = (k2 >> 6) & 3;
    const base = (k2 >> 4) & 3;
    const edge = Math.round(side * EDGE);
    const ax = align === 0 ? edge : align === 2 ? side - edge : side >> 1;
    const ay = base === 0 ? Math.round(side * ALPHABETIC) : base === 1 ? side >> 1 : base === 2 ? edge : side - edge;
    const g = pg.g;
    g.save();
    g.beginPath();
    g.rect(x0, y0, side, side);
    g.clip();
    g.clearRect(x0, y0, side, side);
    g.font = this.res.font(this.font);
    g.textAlign = ALIGN[align];
    g.textBaseline = BASELINE[base];
    g.fillStyle = this.ink < NO_COLOR ? this.res.token(this.ink) : css(this.ink);
    g.setTransform(this.px[0], 0, 0, this.px[0], x0 + ax + fx / 4, y0 + ay + fy / 4);
    g.fillText(text, 0, 0);
    g.restore();
    this.k0[e] = k0;
    this.k1[e] = this.ink;
    this.k2[e] = k2;
    this.ax[e] = ax;
    this.ay[e] = ay;
    this.table[this.at] = e;
    this.made++;
    return e;
  }
  /** Page p's cells gone: the table rebuilt from the others' (edge rate). */
  private evict(p: number): void {
    this.side[p] = 0;
    this.used[p] = 0;
    this.evicted++;
    this.table.fill(-1);
    for (let q = 0; q < PAGES; q++)
      for (let s = 0; s < this.used[q]; s++) {
        const e = q * SLOTS + s;
        this.find(this.k0[e], this.k1[e], this.k2[e]);
        this.table[this.at] = e;
      }
  }
}

function grow<T extends Int32Array | Float64Array>(a: T): T {
  const b = new (a.constructor as new (n: number) => T)(a.length * 2);
  b.set(a);
  return b;
}

/** A shape glyph added to the current path (`drawGlyph`'s shapes, in device px), item i of the batch. */
function shape(g: Paint2D, kind: number, b: GlyphBatch, i: number): void {
  const x = b.x[i];
  const y = b.y[i];
  const r = b.radius[i];
  switch (kind) {
    case 1: // diamond
      g.moveTo(x, y - r);
      g.lineTo(x + r, y);
      g.lineTo(x, y + r);
      g.lineTo(x - r, y);
      break;
    case 2: // triangle
      g.moveTo(x, y - r);
      g.lineTo(x + r, y + r * 0.8);
      g.lineTo(x - r, y + r * 0.8);
      break;
    case 3: // square
      g.rect(x - r * 0.8, y - r * 0.8, r * 1.6, r * 1.6);
      break;
    case 4: // pin: a head over a point
      g.moveTo(x + r * 0.65 * Math.cos(Math.PI * 0.85), y - r * 0.35 + r * 0.65 * Math.sin(Math.PI * 0.85));
      g.arc(x, y - r * 0.35, r * 0.65, Math.PI * 0.85, Math.PI * 2.15);
      g.lineTo(x, y + r);
      break;
    case 5: // star
      for (let k = 0; k < 10; k++) {
        const a = -Math.PI / 2 + (k * Math.PI) / 5;
        const q = k % 2 === 0 ? r : r * 0.45;
        if (k === 0) g.moveTo(x + q * Math.cos(a), y + q * Math.sin(a));
        else g.lineTo(x + q * Math.cos(a), y + q * Math.sin(a));
      }
      break;
    case 6: // cross
      g.rect(x - r, y - r * 0.25, 2 * r, r * 0.5);
      g.rect(x - r * 0.25, y - r, r * 0.5, 2 * r);
      break;
    case 7: // tick
      g.rect(x - 0.5, y - r, 1, 2 * r);
      break;
    default:
      g.moveTo(x + r * 0.85, y);
      g.arc(x, y, r * 0.85, 0, Math.PI * 2);
  }
  g.closePath();
}
