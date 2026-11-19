/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  Anchor,
  type DrawLists,
  type GlyphList,
  type ImageData as DrawImage,
  NoColor,
  type PassName,
  PathCommand,
  PathList,
  type PassVisitor,
  squirclePath,
  walkPass,
} from '@molstar/seqx-core/render';
import { FontCache } from './font.js';

// The Canvas2D backend (_spec/README.md "Canvas2D backend"): draw lists on three stacked canvases, DPR-aware. The
// background and data passes are drawn on a base canvas; the emphasis passes (window, highlight, selection, hover, cursor) on an
// overlay canvas above it, so a hover redraws the overlay only; the label pass on a labels canvas above both, so the
// labels stay over the emphasis. Squircles are the core's curve (`squirclePath`), so this and
// the SVG backend draw the same corner. Buffers are reused: colors and fonts are strings kept per value, glyphs are
// strings kept per code, a path is replayed from the lists, an image goes through one reusable pixel buffer; once warm,
// a frame allocates nothing (a text faded past its width makes a gradient: rare).

/** The members of a 2D context the backend calls: a recording context in tests implements exactly these. */
export interface Context2D {
  fillStyle: unknown;
  strokeStyle: unknown;
  lineWidth: number;
  font: string;
  textAlign: string;
  textBaseline: string;
  imageSmoothingEnabled: boolean;
  save(): void;
  restore(): void;
  setTransform(a: number, b: number, c: number, d: number, e: number, f: number): void;
  clearRect(x: number, y: number, w: number, h: number): void;
  fillRect(x: number, y: number, w: number, h: number): void;
  strokeRect(x: number, y: number, w: number, h: number): void;
  beginPath(): void;
  closePath(): void;
  moveTo(x: number, y: number): void;
  lineTo(x: number, y: number): void;
  bezierCurveTo(x1: number, y1: number, x2: number, y2: number, x: number, y: number): void;
  rect(x: number, y: number, w: number, h: number): void;
  fill(): void;
  stroke(): void;
  clip(): void;
  setLineDash(segments: number[]): void;
  fillText(text: string, x: number, y: number): void;
  createLinearGradient(
    x0: number,
    y0: number,
    x1: number,
    y1: number,
  ): { addColorStop(offset: number, color: string): void };
  createImageData(w: number, h: number): { readonly data: Uint8ClampedArray };
  putImageData(image: unknown, x: number, y: number): void;
  drawImage(source: unknown, x: number, y: number, w: number, h: number): void;
}

/** What the backend needs of a canvas element (or an off-screen canvas). */
export interface CanvasLike {
  width: number;
  height: number;
  style?: { width: string; height: string };
  getContext(kind: '2d'): Context2D | null;
}

export interface CanvasBackendOptions {
  /** Makes an off-screen canvas for images (dense strips): the DOM's by default. */
  readonly createCanvas?: (width: number, height: number) => CanvasLike;
  readonly fonts?: FontCache;
  /** How far a text past its max width fades out, CSS px (the layout's `fade`). */
  readonly fade?: number;
}

const OverlayPasses: readonly PassName[] = ['window', 'highlight', 'selection', 'hover', 'cursor'];
const BasePasses: readonly PassName[] = ['background', 'data'];
const LabelPasses: readonly PassName[] = ['labels'];

/** The most image strips kept: a strip past it shares the last slot (an image is put and drawn at once). */
const MaxImageSlots = 16;

interface ImageSlot {
  readonly canvas: CanvasLike;
  readonly ctx: Context2D;
  width: number;
  height: number;
  image: { readonly data: Uint8ClampedArray };
}

/** A color `0xrrggbbaa` as CSS, kept per value. */
class ColorCache {
  private readonly _css = new Map<number, string>();

  css(c: number): string {
    let s = this._css.get(c);
    if (s === undefined) {
      if (this._css.size > 8192) this._css.clear();
      const a = (c & 0xff) / 255;
      s = `rgba(${(c >>> 24) & 0xff}, ${(c >>> 16) & 0xff}, ${(c >>> 8) & 0xff}, ${a === 1 ? 1 : Math.round(a * 1000) / 1000})`;
      this._css.set(c, s);
    }
    return s;
  }

  /** The color with no opacity, for a gradient's clear end. */
  clear(c: number): string {
    return `rgba(${(c >>> 24) & 0xff}, ${(c >>> 16) & 0xff}, ${(c >>> 8) & 0xff}, 0)`;
  }
}

const Aligns = ['start', 'center', 'end'] as const;

export class CanvasBackend {
  readonly fonts: FontCache;
  /** The device pixel ratio and the size (CSS px) the canvases were last sized for. */
  dpr = 1;
  width = 0;
  height = 0;
  /** Draws counted: what the tests (and a profile) read. */
  readonly counts = { base: 0, overlay: 0, labels: 0 };
  /** The canvases' scale: the backing size over the CSS size (the rounded backing size, not the ratio asked for). */
  private _sx = 1;
  private _sy = 1;
  private readonly _base: Context2D;
  private readonly _overlay: Context2D;
  private readonly _labels: Context2D;
  private readonly _colors = new ColorCache();
  private readonly _glyphs: string[] = [];
  private readonly _scratch = new PathList();
  /** A slot per image strip of a pass (its index there): the strips a view draws are few and their sizes stay. */
  private readonly _slots: ImageSlot[] = [];
  private readonly _fade: number;
  private _overlayDirty = false;
  private _list: DrawLists | undefined;
  private _visitor: PassVisitor;
  private _ctx: Context2D;

  constructor(
    readonly base: CanvasLike,
    readonly overlay: CanvasLike,
    readonly labels: CanvasLike,
    private readonly options: CanvasBackendOptions = {},
  ) {
    const b = base.getContext('2d');
    const o = overlay.getContext('2d');
    const t = labels.getContext('2d');
    if (!b || !o || !t) throw new Error('A 2D canvas context is not available');
    this._base = this._ctx = b;
    this._overlay = o;
    this._labels = t;
    this.fonts = options.fonts ?? new FontCache();
    this._fade = options.fade ?? 12;
    this._visitor = this._makeVisitor();
  }

  /**
   * Sizes the three canvases for `width` × `height` CSS px at `dpr`; true when anything changed (a redraw follows). The
   * backing size is rounded, so the scale drawn with is the backing size over the CSS size.
   */
  resize(width: number, height: number, dpr: number): boolean {
    if (width === this.width && height === this.height && dpr === this.dpr) return false;
    this.width = width;
    this.height = height;
    this.dpr = dpr;
    for (const c of [this.base, this.overlay, this.labels]) {
      // assigning a canvas's size clears it
      c.width = Math.max(1, Math.round(width * dpr));
      c.height = Math.max(1, Math.round(height * dpr));
      if (c.style) {
        c.style.width = `${width}px`;
        c.style.height = `${height}px`;
      }
    }
    this._sx = this.base.width / width;
    this._sy = this.base.height / height;
    this._overlayDirty = false;
    return true;
  }

  /**
   * Draws the lists: the base canvas (background, data) and the labels canvas above everything when `base` (a layout
   * was rebuilt), the overlay between them always (a hover, a selection or the cursor changed, or the base did). The
   * canvases stack as the passes do, so the picture is the SVG's. An overlay with nothing to clear and nothing to draw is
   * left alone.
   */
  draw(lists: DrawLists, base: boolean): void {
    if (base) {
      this._drawPasses(lists, this._base, BasePasses, true);
      this.counts.labels++;
      this._labels.setTransform(this._sx, 0, 0, this._sy, 0, 0);
      this._labels.clearRect(0, 0, this.width, this.height);
      this._drawPasses(lists, this._labels, LabelPasses, false);
    }
    this._drawOverlay(lists);
  }

  /** Redraws the overlay alone. */
  drawOverlay(lists: DrawLists): void {
    this._drawOverlay(lists);
  }

  private _drawOverlay(lists: DrawLists): void {
    let any = false;
    for (const p of OverlayPasses) if (lists.passes[p].size) any = true;
    if (!any && !this._overlayDirty) return;
    this.counts.overlay++;
    const ctx = this._overlay;
    ctx.setTransform(this._sx, 0, 0, this._sy, 0, 0);
    ctx.clearRect(0, 0, this.width, this.height);
    this._drawPasses(lists, ctx, OverlayPasses, false);
    this._overlayDirty = any;
  }

  private _drawPasses(lists: DrawLists, ctx: Context2D, passes: readonly PassName[], countBase: boolean): void {
    if (countBase) this.counts.base++;
    ctx.setTransform(this._sx, 0, 0, this._sy, 0, 0);
    ctx.textBaseline = 'alphabetic';
    this._list = lists;
    this._ctx = ctx;
    for (const name of passes) {
      const pass = lists.passes[name];
      if (!pass.size) continue;
      walkPass(pass, this._visitor);
      // a pass leaves no clip behind
      for (let d = pass.depth; d > 0; d--) ctx.restore();
    }
  }

  /**
   * A placed run: each letter moved to its place, turned and scaled about it, through the context's transform (the base
   * transform composed with the letter's, then set back); drawn on some views only.
   */
  private _placedGlyphs(l: GlyphList, i: number): void {
    const ctx = this._ctx;
    const sx = this._sx;
    const sy = this._sy;
    const one = l.color[i] !== NoColor;
    if (one) ctx.fillStyle = this._colors.css(l.color[i]);
    for (let j = l.start[i], end = j + l.count[i]; j < end; j++) {
      if (!one) ctx.fillStyle = this._colors.css(l.colors[j]);
      const c = Math.cos(l.angle[j]);
      const s = Math.sin(l.angle[j]);
      const v = l.vs[j];
      ctx.setTransform(sx * c, sy * s, -sx * s * v, sy * c * v, sx * l.x[j], sy * l.ly[j]);
      ctx.fillText(this._glyph(l.codes[j]), 0, 0);
    }
    ctx.setTransform(sx, 0, 0, sy, 0, 0);
  }

  private _glyph(code: number): string {
    return this._glyphs[code] ?? (this._glyphs[code] = String.fromCharCode(code));
  }

  private _tracePath(l: PathList, i: number): void {
    const ctx = this._ctx;
    ctx.beginPath();
    let pt = l.pointStart[i];
    for (let c = l.cmdStart[i], end = c + l.cmdCount[i]; c < end; c++) {
      switch (l.commands[c]) {
        case PathCommand.move:
          ctx.moveTo(l.points[2 * pt], l.points[2 * pt + 1]);
          pt++;
          break;
        case PathCommand.line:
          ctx.lineTo(l.points[2 * pt], l.points[2 * pt + 1]);
          pt++;
          break;
        case PathCommand.cubic:
          ctx.bezierCurveTo(
            l.points[2 * pt],
            l.points[2 * pt + 1],
            l.points[2 * pt + 2],
            l.points[2 * pt + 3],
            l.points[2 * pt + 4],
            l.points[2 * pt + 5],
          );
          pt += 3;
          break;
        case PathCommand.close:
          ctx.closePath();
          break;
      }
    }
  }

  private _paint(fill: number, stroke: number, strokeWidth: number): void {
    const ctx = this._ctx;
    if (fill !== NoColor) {
      ctx.fillStyle = this._colors.css(fill);
      ctx.fill();
    }
    if (stroke !== NoColor && strokeWidth > 0) {
      ctx.strokeStyle = this._colors.css(stroke);
      ctx.lineWidth = strokeWidth;
      ctx.stroke();
    }
  }

  private _image(img: DrawImage, strip: number, x: number, y: number, w: number, h: number): void {
    const at = Math.min(strip, MaxImageSlots - 1);
    let slot = this._slots[at];
    if (!slot) {
      const canvas = this.options.createCanvas
        ? this.options.createCanvas(img.width, img.height)
        : (Object.assign(document.createElement('canvas'), { width: img.width, height: img.height }) as CanvasLike);
      const ctx = canvas.getContext('2d');
      if (!ctx) return;
      this._slots[at] = slot = {
        canvas,
        ctx,
        width: img.width,
        height: img.height,
        image: ctx.createImageData(img.width, img.height),
      };
    } else if (slot.width !== img.width || slot.height !== img.height) {
      // the strip's size changed: its canvas and pixel buffer are made again, in place
      slot.canvas.width = slot.width = img.width;
      slot.canvas.height = slot.height = img.height;
      slot.image = slot.ctx.createImageData(img.width, img.height);
    }
    // the layer owns and reuses its bytes: they are copied into the slot's buffer at each draw
    slot.image.data.set(img.data);
    slot.ctx.putImageData(slot.image, 0, 0);
    const ctx = this._ctx;
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(slot.canvas, x, y, w, h);
  }

  private _makeVisitor(): PassVisitor {
    return {
      rect: (l, i) => {
        const ctx = this._ctx;
        const fill = l.fill[i];
        const stroke = l.stroke[i];
        if (l.radius[i] > 0) {
          const s = this._scratch;
          s.n = 0;
          s.commandCount = 0;
          s.pointCount = 0;
          s.begin();
          squirclePath(s, l.x[i], l.y[i], l.w[i], l.h[i], l.radius[i], l.round[i]);
          s.end(fill);
          this._tracePath(s, 0);
          this._paint(fill, stroke, l.strokeWidth[i]);
          return;
        }
        if (fill !== NoColor) {
          ctx.fillStyle = this._colors.css(fill);
          ctx.fillRect(l.x[i], l.y[i], l.w[i], l.h[i]);
        }
        if (stroke !== NoColor && l.strokeWidth[i] > 0) {
          ctx.strokeStyle = this._colors.css(stroke);
          ctx.lineWidth = l.strokeWidth[i];
          ctx.strokeRect(l.x[i], l.y[i], l.w[i], l.h[i]);
        }
      },
      path: (l, i) => {
        const ctx = this._ctx;
        this._tracePath(l, i);
        // a dash applies to this path alone
        const dashed = l.dash[i] > 0;
        if (dashed) ctx.setLineDash([l.dash[i], l.dash[i]]);
        this._paint(l.fill[i], l.stroke[i], l.strokeWidth[i]);
        if (dashed) ctx.setLineDash(NoDash);
      },
      glyphs: (l, i) => {
        const ctx = this._ctx;
        const lists = this._list!;
        ctx.font = this.fonts.css(lists.string(l.font[i]));
        ctx.textAlign = 'center';
        const start = l.start[i];
        const end = start + l.count[i];
        const y = l.y[i];
        if (l.placed[i] === 1) {
          this._placedGlyphs(l, i);
          return;
        }
        if (l.color[i] !== NoColor) {
          ctx.fillStyle = this._colors.css(l.color[i]);
          for (let j = start; j < end; j++) ctx.fillText(this._glyph(l.codes[j]), l.x[j], y);
          return;
        }
        for (let j = start; j < end; j++) {
          ctx.fillStyle = this._colors.css(l.colors[j]);
          ctx.fillText(this._glyph(l.codes[j]), l.x[j], y);
        }
      },
      text: (l, i) => {
        const ctx = this._ctx;
        const lists = this._list!;
        ctx.font = this.fonts.css(lists.string(l.font[i]));
        ctx.textAlign = Aligns[l.anchor[i]];
        const x = l.x[i];
        const color = l.color[i];
        if (l.width[i] > l.maxWidth[i]) {
          // faded past its max width, never ellipsized: a gradient from opaque to clear over the last `fade` px
          const anchor = l.anchor[i];
          const x0 = anchor === Anchor.start ? x : anchor === Anchor.end ? x - l.maxWidth[i] : x - l.maxWidth[i] / 2;
          const x1 = x0 + l.maxWidth[i];
          const fade = Math.min(this._fade, l.maxWidth[i] / 2);
          const g = ctx.createLinearGradient(x1 - fade, 0, x1, 0);
          g.addColorStop(0, this._colors.css(color));
          g.addColorStop(1, this._colors.clear(color));
          ctx.save();
          ctx.beginPath();
          ctx.rect(x0, l.y[i] - 1000, l.maxWidth[i], 2000);
          ctx.clip();
          ctx.fillStyle = g;
          ctx.fillText(lists.string(l.str[i]), x, l.y[i]);
          ctx.restore();
          return;
        }
        ctx.fillStyle = this._colors.css(color);
        ctx.fillText(lists.string(l.str[i]), x, l.y[i]);
      },
      image: (l, i) => {
        this._image(l.images[i], i, l.x[i], l.y[i], l.w[i], l.h[i]);
      },
      clip: (l, i) => {
        const ctx = this._ctx;
        ctx.save();
        if (l.radius[i] > 0) {
          // a squircle clip (a track's well): the core's curve, as the SVG's clip path
          const s = this._scratch;
          s.n = 0;
          s.commandCount = 0;
          s.pointCount = 0;
          s.begin();
          squirclePath(s, l.x[i], l.y[i], l.w[i], l.h[i], l.radius[i], l.round[i]);
          s.end(NoColor);
          this._tracePath(s, 0);
        } else {
          ctx.beginPath();
          ctx.rect(l.x[i], l.y[i], l.w[i], l.h[i]);
        }
        ctx.clip();
      },
      unclip: () => {
        this._ctx.restore();
      },
    };
  }
}

const NoDash: number[] = [];
