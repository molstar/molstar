/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SVG backend (seqspec:viewer#the-render-backend): each batch written as SVG text, platform-free. Device px are
// divided by the pass's DPR, so the document is in CSS px and an export at DPR 2 lays out as the screen does. Batches
// are grouped by (color, alpha) as Canvas2D groups them (one `<path>` per group, so the file stays small): rectangles
// are subpaths with `shape-rendering="crispEdges"`, sectors and paths use arcs, glyph shapes are subpaths of their
// group, path glyphs are `<path>`s in `<defs>` and a `<use>` per mark, text is `<text>` per item under a `<g>` carrying
// the font, anchor and baseline, clips are `<clipPath>`s, passes are `<g>`s in order. Images are PNG data URIs in
// `<defs>` drawn by `<use>` (a source rectangle smaller than the image is cut by a clip); a layer (`beginLayer`, the
// map miniature) stays vector, a `<g>` in `<defs>`. The writer appends to arrays of parts and joins them once, in
// `toString`. It draws whole frames (`caps`: RotatedText only), so `copy` is an error and `present` and `stale` have
// nothing to do.

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
  SkinToken,
  TextAlign,
  TextBaseline,
  type TextBatch,
} from "../draw/backend.js";
import { ColorGroups } from "../draw/group.js";
import type { TextMeasurer } from "../text.js";
import { cssColor, dec, dec3, esc, hex } from "./format.js";
import { CSS_FONT, type EmBox } from "./metrics.js";
import { base64, encodePng } from "./png.js";

export interface SvgFont {
  family: string;
  weight: number;
  style: "normal" | "italic";
  data: Uint8Array;
  format: "woff2" | "woff" | "truetype";
}
/** What `toString` and `renderSvg` take from the host: the platform's pieces, so this entry imports none. */
export interface SvgPorts {
  /**
   * Text widths for fitting (seqspec:viewer#text-measuring), for a context `renderSvg` makes; absent: the context's
   * own.
   */
  readonly measurer?: TextMeasurer;
  /** PNG compression, a zlib stream out (`zlib.deflateSync`, `CompressionStream("deflate")`); absent: stored blocks. */
  readonly deflate?: (bytes: Uint8Array) => Uint8Array | Promise<Uint8Array>;
  /** Embedded as `@font-face` rules; absent: the skin's families are referenced by name. */
  readonly fonts?: readonly SvgFont[];
}
export interface SvgOptions {
  /** The file's pixels per CSS px (its `width` and `height` over its `viewBox`); default the first pass's DPR. */
  scale?: number;
  /** The overlay's passes (marks, selection, the cursor); default off. */
  overlay?: boolean;
  /** `renderSvg` only: the context's hover lights its element in the overlay; default off. */
  hover?: boolean;
  /** The navigator's passes; default off. */
  nav?: boolean;
  /** A rectangle in the skin's surface color behind everything. */
  background?: boolean;
  title?: string;
  /** Keep a line per backend call in `trace` (the conformance kit compares it with the recording's op log). */
  trace?: boolean;
  /**
   * The em box of a CSS font string (`metricsMeasurer(...).metrics`, a canvas measurer's): with it, `middle`, `top` and
   * `bottom` are written as an explicit `dy` from the alphabetic baseline, as Canvas2D places them and every reader
   * honors; without it they are `dominant-baseline`, which readers place by their own metrics (Inkscape ignores it).
   */
  metrics?: (font: string) => EmBox | null;
}
/** Items written, by kind (NO_COLOR items are not), and numbers that were not finite (written as 0). */
export interface SvgStats {
  rects: number;
  sectors: number;
  paths: number;
  glyphs: number;
  texts: number;
  images: number;
  clips: number;
  nonFinite: number;
}

const PASS_CLASS = ["main", "overlay", "nav", "nav-overlay", "vnav", "vnav-overlay"];
const TAU = 2 * Math.PI;
const MIME = { woff2: "font/woff2", woff: "font/woff", truetype: "font/ttf" } as const;

/** An image: pixels the renderer made, or a vector layer; kept (unused ones dropped) until the file is written. */
export interface SvgImage {
  w: number;
  h: number;
  px: Uint32Array | null;
  layer: string[] | null;
  /** Drawn in this document; disposed by its renderer (kept until written if drawn). */
  used: boolean;
  gone: boolean;
}
interface Font {
  family: string;
  px: number;
  weight: string;
  /** The em box the host's metrics give this font; null: none. */
  box: EmBox | null;
}

/** What the backends of one document share: its images, `<defs>` and clip ids (`images`, as Canvas2DShared is). */
export class SvgShared implements BackendImages {
  readonly imgs: SvgImage[] = [];
  defs: string[] = [];
  clips = new Map<string, number>();
  glyphs = new Set<number>();
  dpr = 1;
  private layer: SvgBackend | null = null;
  private lw = 0;
  private lh = 0;
  readonly res: DrawResources;
  readonly opts: SvgOptions;
  constructor(res: DrawResources, opts: SvgOptions) {
    this.res = res;
    this.opts = opts;
  }
  fromPixels(w: number, h: number, rgba: Uint32Array): number {
    return this.keep({ w, h, px: rgba.slice(0, w * h), layer: null, used: false, gone: false });
  }
  beginLayer(w: number, h: number): RenderBackend {
    this.layer = new SvgBackend(this.res, { ...this.opts, overlay: true, nav: true, trace: false }, this);
    this.lw = w;
    this.lh = h;
    this.layer.beginPass(Pass.main, w, h, this.dpr);
    return this.layer;
  }
  endLayer(): number {
    const l = this.layer as SvgBackend;
    l.endPass();
    // Clipped to the layer's own bounds, as a canvas of its size is.
    const id = this.clip(0, 0, this.lw, this.lh, 100 / this.dpr);
    return this.keep({
      w: this.lw,
      h: this.lh,
      px: null,
      layer: [`<g clip-path="url(#c${id})">`, ...l.out, "</g>"],
      used: false,
      gone: false,
    });
  }
  dispose(img: number): void {
    const m = this.imgs[img];
    if (!m) return;
    m.gone = true;
    // Drawn images stay until the file is written; others are not needed.
    if (!m.used) m.px = m.layer = null;
  }
  /** A clip id for a rectangle in device px at `q` (hundredths per device px), defined once. */
  clip(x: number, y: number, w: number, h: number, q: number): number {
    const x0 = Math.round(x * q);
    const y0 = Math.round(y * q);
    const key = `${x0} ${y0} ${Math.round((x + w) * q)} ${Math.round((y + h) * q)}`;
    let id = this.clips.get(key);
    if (id === undefined) {
      id = this.clips.size;
      this.clips.set(key, id);
      const [, , x1, y1] = key.split(" ").map(Number);
      this.defs.push(
        `<clipPath id="c${id}"><rect x="${dec(x0)}" y="${dec(y0)}" width="${dec(x1 - x0)}" height="${dec(y1 - y0)}"/></clipPath>`,
      );
    }
    return id;
  }
  /** A new document: the defs and clip ids dropped, images kept (renderers hold their ids across frames). */
  reset(): void {
    this.defs = [];
    this.clips.clear();
    this.glyphs.clear();
    for (const m of this.imgs) {
      m.used = false;
      if (m.gone) m.px = m.layer = null;
    }
  }
  private keep(m: SvgImage): number {
    return this.imgs.push(m) - 1;
  }
}

export class SvgBackend implements RenderBackend {
  readonly caps = BackendCaps.RotatedText;
  readonly images: BackendImages;
  readonly stats: SvgStats = { rects: 0, sectors: 0, paths: 0, glyphs: 0, texts: 0, images: 0, clips: 0, nonFinite: 0 };
  /** One line per call received (`SvgOptions.trace`), muted passes included. */
  readonly trace: string[] | null;
  /** The body so far (a nested layer's, for its image). */
  out: string[] = [];
  private readonly sh: SvgShared;
  private readonly res: DrawResources;
  private readonly opts: SvgOptions;
  private readonly groups = new ColorGroups();
  private readonly fonts = new Map<string, Font>();
  /** The pass: hundredths per device px (100 / DPR), whether its output is dropped, its first part, open clips. */
  private q = 100;
  private dpr = 1;
  private mute = false;
  private mark = 0;
  private depth = 0;
  private pw = 0;
  private ph = 0;
  /** The section offset (CSS px), the document's extent (CSS px) and the first pass's DPR. */
  private ox = 0;
  private oy = 0;
  /** The offset `beginOffset` added (CSS px). */
  private bx = 0;
  private by = 0;
  private width = 0;
  private height = 0;
  private first = 0;

  constructor(res: DrawResources, opts: SvgOptions = {}, shared?: SvgShared) {
    this.res = res;
    this.opts = opts;
    this.sh = shared ?? new SvgShared(res, opts);
    this.images = this.sh;
    this.trace = opts.trace ? [] : null;
  }

  /** A new document; the images stay (renderers hold their ids across frames). */
  reset(): void {
    this.out = [];
    this.sh.reset();
    if (this.trace) this.trace.length = 0;
    this.width = this.height = this.first = this.ox = this.oy = this.bx = this.by = this.depth = 0;
    for (const k of Object.keys(this.stats) as (keyof SvgStats)[]) this.stats[k] = 0;
  }
  /** Passes that follow are placed at (x, y) CSS px, until `endSection`: sections stacked. */
  beginSection(x: number, y: number): void {
    this.ox = x;
    this.oy = y;
    this.out.push(
      `<g class="seqv-section" transform="translate(${dec(Math.round(x * 100))} ${dec(Math.round(y * 100))})">`,
    );
  }
  endSection(): void {
    this.out.push("</g>");
  }
  /** Passes that follow are offset (x, y) CSS px more inside the section, until `endOffset`: a box of their own (`VerticalNav`). */
  beginOffset(x: number, y: number): void {
    this.bx = x;
    this.by = y;
    this.ox += x;
    this.oy += y;
    this.out.push(`<g transform="translate(${dec(Math.round(x * 100))} ${dec(Math.round(y * 100))})">`);
  }
  endOffset(): void {
    this.ox -= this.bx;
    this.oy -= this.by;
    this.bx = this.by = 0;
    this.out.push("</g>");
  }

  beginPass(pass: number, width: number, height: number, dpr: number): void {
    this.trace?.push(`beginPass ${pass} ${width} ${height} ${dpr}`);
    this.dpr = this.sh.dpr = dpr;
    this.q = 100 / dpr;
    this.pw = width;
    this.ph = height;
    this.mute = pass === Pass.main ? false : !(pass === Pass.overlay ? this.opts.overlay : this.opts.nav);
    if (this.mute) return;
    if (this.first === 0) this.first = dpr;
    this.width = Math.max(this.width, this.ox + width / dpr);
    this.height = Math.max(this.height, this.oy + height / dpr);
    this.out.push(`<g class="seqv-${PASS_CLASS[pass] ?? pass}">`);
    this.mark = this.out.length;
    this.depth = 0;
  }
  endPass(): void {
    this.trace?.push("endPass");
    if (this.mute) return;
    for (; this.depth > 0; this.depth--) this.out.push("</g>");
    this.out.push("</g>");
  }
  /** Whole frames never clear part of a pass: a clear of all of it drops what was drawn, any other is ignored. */
  clear(x: number, y: number, w: number, h: number): void {
    this.trace?.push(`clear ${x} ${y} ${w} ${h}`);
    if (!this.mute && this.depth === 0 && x <= 0 && y <= 0 && x + w >= this.pw && y + h >= this.ph)
      this.out.length = this.mark;
  }
  pushClip(x: number, y: number, w: number, h: number): void {
    this.trace?.push(`pushClip ${x} ${y} ${w} ${h}`);
    if (this.mute) return;
    this.out.push(`<g clip-path="url(#c${this.sh.clip(x, y, w, h, this.q)})">`);
    this.depth++;
    this.stats.clips++;
  }
  popClip(): void {
    this.trace?.push("popClip");
    if (this.mute || this.depth === 0) return;
    this.out.push("</g>");
    this.depth--;
  }
  copy(dx: number, dy: number, x: number, y: number, w: number, h: number): void {
    this.trace?.push(`copy ${dx} ${dy} ${x} ${y} ${w} ${h}`);
    throw new Error("SvgBackend has no Copy: it draws whole frames (caps).");
  }
  present(x: number, y: number, w: number, h: number): void {
    this.trace?.push(`present ${x} ${y} ${w} ${h}`);
  }
  stale(sx: number, sy: number, sw: number, sh: number, dx: number, dy: number, dw: number, dh: number): void {
    this.trace?.push(`stale ${sx} ${sy} ${sw} ${sh} ${dx} ${dy} ${dw} ${dh}`);
  }

  rects(b: RectBatch): void {
    this.trace?.push(`rects ${b.n}`);
    if (this.mute || b.n === 0) return;
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    for (let g = 0; g < groups; g++) {
      const first = G.order[G.start[g]];
      let d = "";
      for (let j = G.start[g]; j < G.start[g + 1]; j++) {
        const i = G.order[j];
        const s = b.stroke[i];
        // An outline is four filled edges inside the rectangle.
        if (s === 0 || 2 * s >= b.w[i] || 2 * s >= b.h[i]) d += this.box(b.x[i], b.y[i], b.w[i], b.h[i]);
        else {
          d += this.box(b.x[i], b.y[i], b.w[i], s) + this.box(b.x[i], b.y[i] + b.h[i] - s, b.w[i], s);
          d +=
            this.box(b.x[i], b.y[i] + s, s, b.h[i] - 2 * s) +
            this.box(b.x[i] + b.w[i] - s, b.y[i] + s, s, b.h[i] - 2 * s);
        }
      }
      this.out.push(
        `<path${this.paint(b.color[first], b.alpha[first], "fill")} shape-rendering="crispEdges" d="${d}"/>`,
      );
    }
    this.stats.rects += G.start[groups];
  }
  sectors(b: SectorBatch): void {
    this.trace?.push(`sectors ${b.n}`);
    if (this.mute || b.n === 0) return;
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    for (let g = 0; g < groups; g++) {
      const first = G.order[G.start[g]];
      let d = "";
      for (let j = G.start[g]; j < G.start[g + 1]; j++) {
        const i = G.order[j];
        d += this.sector(b.cx[i], b.cy[i], b.r0[i], b.r1[i], b.a0[i], b.a1[i]);
      }
      this.out.push(`<path${this.paint(b.color[first], b.alpha[first], "fill")} d="${d}"/>`);
    }
    this.stats.sectors += G.start[groups];
  }
  paths(b: PathBatch): void {
    this.trace?.push(`paths ${b.n}`);
    if (this.mute) return;
    for (let p = 0; p < b.n; p++) {
      if (b.color[p] === NO_COLOR) continue;
      // Subpaths of one style are one path, as Canvas2D draws them: translucent fills do not darken where they overlap.
      let end = p + 1;
      while (
        end < b.n &&
        b.color[end] === b.color[p] &&
        b.alpha[end] === b.alpha[p] &&
        b.mode[end] === b.mode[p] &&
        b.width[end] === b.width[p]
      )
        end++;
      let d = "";
      let has = false;
      for (let s = p; s < end; s++) {
        let c = b.cstart[s];
        const v1 = s + 1 < b.n ? b.start[s + 1] : b.nv;
        for (let v = b.start[s]; v < v1; v++) {
          const k = b.coords;
          switch (b.verbs[v]) {
            case PathVerb.move:
              d += `M${this.pt(k[c], k[c + 1])}`;
              has = true;
              break;
            case PathVerb.line:
              d += `${has ? "L" : "M"}${this.pt(k[c], k[c + 1])}`;
              has = true;
              break;
            case PathVerb.quad:
              d += `${has ? "" : `M${this.pt(k[c], k[c + 1])}`}Q${this.pt(k[c], k[c + 1])} ${this.pt(k[c + 2], k[c + 3])}`;
              has = true;
              break;
            case PathVerb.cubic:
              d += `${has ? "" : `M${this.pt(k[c], k[c + 1])}`}C${this.pt(k[c], k[c + 1])} ${this.pt(k[c + 2], k[c + 3])} ${this.pt(k[c + 4], k[c + 5])}`;
              has = true;
              break;
            case PathVerb.arc:
              d += this.arc(k[c], k[c + 1], k[c + 2], k[c + 3], k[c + 4], k[c + 5] !== 0, has);
              has = true;
              break;
            case PathVerb.roundRect:
              d += this.roundRect(k[c], k[c + 1], k[c + 2], k[c + 3], k[c + 4]);
              has = true;
              break;
            default:
              d += "Z";
          }
          c += PathVerbCoords[b.verbs[v]];
        }
      }
      const style =
        b.mode[p] === PathMode.stroke
          ? `${this.paint(b.color[p], b.alpha[p], "stroke")} fill="none" stroke-width="${this.n(Math.round((b.width[p] / 2) * this.q))}"`
          : this.paint(b.color[p], b.alpha[p], "fill");
      this.out.push(`<path${style} d="${d}"/>`);
      this.stats.paths += end - p;
      p = end - 1;
    }
  }
  glyphs(b: GlyphBatch): void {
    this.trace?.push(`glyphs ${b.n} ${b.font}`);
    if (this.mute || b.n === 0) return;
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    const res = this.res;
    const font = this.fontOf(res.font(b.font));
    for (let g = 0; g < groups; g++) {
      const first = G.order[G.start[g]];
      let shapes = "";
      let marks = "";
      for (let j = G.start[g]; j < G.start[g + 1]; j++) {
        const i = G.order[j];
        const id = b.glyph[i];
        const rec = id < SHAPE_COUNT ? null : res.glyph(id);
        if (rec?.char) {
          // A character in the batch's font, 2r high, centered on its point.
          marks += `<text x="${this.p(b.x[i])}" y="${this.p(b.y[i])}" font-size="${this.p(2 * b.radius[i])}"${this.dy(font, TextBaseline.middle, (2 * b.radius[i]) / this.dpr) ?? ""}>${esc(rec.char)}</text>`;
        } else if (rec?.path) {
          // Scaled so the box's larger side is 2r, centered on the point.
          const s = (2 * b.radius[i]) / Math.max(rec.box[0], rec.box[1]);
          this.glyphPath(id, rec);
          marks += `<use xlink:href="#p${id}" transform="translate(${this.p(b.x[i] - (rec.box[0] * s) / 2)} ${this.p(b.y[i] - (rec.box[1] * s) / 2)}) scale(${sig((s * this.q) / 100)})"/>`;
        } else shapes += this.shape(id < SHAPE_COUNT ? id : rec ? rec.fallback : 0, b.x[i], b.y[i], b.radius[i]);
      }
      const paint = this.paint(b.color[first], b.alpha[first], "fill");
      if (shapes) this.out.push(`<path${paint} d="${shapes}"/>`);
      if (marks)
        this.out.push(
          `<g${paint} font-family="${esc(font.family)}"${weight(font)} text-anchor="middle"${font.box ? "" : ' dominant-baseline="central"'}>${marks}</g>`,
        );
    }
    this.stats.glyphs += G.start[groups];
  }
  text(b: TextBatch): void {
    this.trace?.push(`text ${b.n} ${b.font} ${b.align} ${b.baseline} ${b.angle}`);
    if (this.mute || b.n === 0) return;
    const G = this.groups;
    const groups = G.group(b.color, b.alpha, b.n);
    const font = this.fontOf(this.res.font(b.font));
    const turned = b.angle !== 0 ? ` rotate(${sig((b.angle * 180) / Math.PI)})` : "";
    const head =
      `<g font-family="${esc(font.family)}" font-size="${this.n(Math.round(font.px * 100))}"${weight(font)}` +
      `${b.align === TextAlign.left ? "" : ` text-anchor="${b.align === TextAlign.center ? "middle" : "end"}"`}` +
      `${b.baseline === TextBaseline.alphabetic || font.box ? "" : ` dominant-baseline="${BASELINE[b.baseline]}"`} xml:space="preserve"`;
    const dy = this.dy(font, b.baseline, font.px) ?? "";
    for (let g = 0; g < groups; g++) {
      const first = G.order[G.start[g]];
      let items = "";
      for (let j = G.start[g]; j < G.start[g + 1]; j++) {
        const i = G.order[j];
        const sx = b.sx[i];
        const sy = b.sy[i];
        // Placed by `translate(x y) rotate(a) scale(sx sy)`: the scale about the anchor, before the angle.
        items +=
          turned === "" && sx === 1 && sy === 1
            ? `<text x="${this.p(b.x[i])}" y="${this.p(b.y[i])}"${dy}>${esc(b.text[i])}</text>`
            : `<text transform="translate(${this.pt(b.x[i], b.y[i])})${turned}${sx === 1 && sy === 1 ? "" : ` scale(${sig(sx)} ${sig(sy)})`}"${dy}>${esc(b.text[i])}</text>`;
      }
      this.out.push(`${head}${this.paint(b.color[first], b.alpha[first], "fill")}>${items}</g>`);
    }
    this.stats.texts += G.start[groups];
  }
  image(img: number, r: Float64Array, smooth: boolean): void {
    this.trace?.push(`image ${img} ${r[0]} ${r[1]} ${r[2]} ${r[3]} ${r[4]} ${r[5]} ${r[6]} ${r[7]} ${+smooth}`);
    const m = this.sh.imgs[img];
    if (this.mute || !m || (!m.px && !m.layer)) return;
    m.used = true;
    const [sx, sy, sw, sh, dx, dy, dw, dh] = r;
    const kx = dw / sw;
    const ky = dh / sh;
    // A pixel image is w × h pixels: one pixel is kx device px, so kx / dpr CSS px. A layer's content is in CSS px
    // already, scaled by kx only.
    const u = m.px ? 1 / this.dpr : 1;
    const use =
      `<use xlink:href="#i${img}" transform="translate(${this.p(dx - sx * kx)} ${this.p(dy - sy * ky)}) scale(${sig(kx * u)} ${sig(ky * u)})"` +
      `${smooth || !m.px ? "" : ' image-rendering="pixelated"'}/>`;
    // The source rectangle is the whole image, or the destination cuts it out.
    const whole = sx <= 0 && sy <= 0 && sx + sw >= m.w && sy + sh >= m.h;
    this.out.push(whole ? use : `<g clip-path="url(#c${this.sh.clip(dx, dy, dw, dh, this.q)})">${use}</g>`);
    this.stats.images++;
  }

  /** The file: `<defs>` (clips, path glyphs, the images drawn, as PNG data URIs), then the passes in order. */
  async toString(ports: SvgPorts = {}): Promise<string> {
    const images: string[] = [];
    for (let k = 0; k < this.sh.imgs.length; k++) {
      const m = this.sh.imgs[k];
      if (!m.used) continue;
      if (m.layer) images.push(`<g id="i${k}">${m.layer.join("")}</g>`);
      else {
        const png = await encodePng(m.w, m.h, m.px as Uint32Array, ports.deflate);
        images.push(
          `<image id="i${k}" width="${m.w}" height="${m.h}" xlink:href="data:image/png;base64,${base64(png)}"/>`,
        );
      }
    }
    const w = dec(Math.round(this.width * 100));
    const h = dec(Math.round(this.height * 100));
    const scale = this.opts.scale ?? (this.first || 1);
    const back = this.opts.background
      ? `<rect width="${w}" height="${h}"${this.paint(SkinToken.surface, 255, "fill")}/>`
      : "";
    const fonts = (ports.fonts ?? [])
      .map(
        (f) =>
          `@font-face{font-family:"${esc(f.family)}";font-weight:${f.weight};font-style:${f.style};` +
          `src:url(data:${MIME[f.format]};base64,${base64(f.data)}) format("${f.format === "truetype" ? "truetype" : f.format}")}`,
      )
      .join("");
    return (
      `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${dec(Math.round(this.width * scale * 100))}" ` +
      `height="${dec(Math.round(this.height * scale * 100))}" viewBox="0 0 ${w} ${h}">` +
      (this.opts.title ? `<title>${esc(this.opts.title)}</title>` : "") +
      (fonts ? `<style>${fonts}</style>` : "") +
      `<defs>${this.sh.defs.join("")}${images.join("")}</defs>${back}${this.out.join("")}</svg>`
    );
  }

  // ---- numbers, colors and shapes --------------------------------------------------------------------------------

  /** Hundredths as a decimal; a number that is not finite is written 0 and counted. */
  private n(h: number): string {
    if (Number.isFinite(h)) return dec(h);
    this.stats.nonFinite++;
    return "0";
  }
  /** A device px coordinate as CSS px. */
  private p(v: number): string {
    return this.n(Math.round(v * this.q));
  }
  private pt(x: number, y: number): string {
    return `${this.p(x)} ${this.p(y)}`;
  }
  /** A rectangle as a closed subpath: its far edges rounded on their own, so neighbors share edges exactly. */
  private box(x: number, y: number, w: number, h: number): string {
    const q = this.q;
    const x0 = Math.round(x * q);
    const y0 = Math.round(y * q);
    return `M${this.n(x0)} ${this.n(y0)}h${this.n(Math.round((x + w) * q) - x0)}v${this.n(Math.round((y + h) * q) - y0)}h${this.n(x0 - Math.round((x + w) * q))}z`;
  }
  /** ` fill="#rrggbb" fill-opacity="…"` (or stroke): a skin token's CSS color read to hex and its own opacity. */
  private paint(c: number, a: number, attr: "fill" | "stroke"): string {
    let color = hex(c);
    let o = a / 255;
    if (c < NO_COLOR) {
      const [h, k] = cssColor(this.res.token(c));
      color = esc(h);
      o *= k;
    }
    return ` ${attr}="${color}"${o < 1 ? ` ${attr}-opacity="${dec3(o)}"` : ""}`;
  }
  private fontOf(css: string): Font {
    let f = this.fonts.get(css);
    if (!f) {
      const m = CSS_FONT.exec(css);
      f = {
        weight: m?.[1] === "normal" ? "400" : (m?.[1] ?? "400"),
        px: Number(m?.[2] ?? 12),
        family: m?.[3] ?? css,
        box: this.opts.metrics?.(css) ?? null,
      };
      this.fonts.set(css, f);
    }
    return f;
  }
  /**
   * ` dy="…"` placing a non-alphabetic baseline on the alphabetic one, as Canvas2D's em box does (`top` is the box's
   * ascent below the point, `bottom` its descent above it, `middle` half the difference); "" for alphabetic; null when
   * the font has no box (the caller writes `dominant-baseline`). `px` is the font's size in CSS px.
   */
  private dy(font: Font, baseline: number, px: number): string | null {
    const box = font.box;
    if (!box) return null;
    if (baseline === TextBaseline.alphabetic) return "";
    const d =
      baseline === TextBaseline.top
        ? box.ascent
        : baseline === TextBaseline.bottom
          ? -box.descent
          : (box.ascent - box.descent) / 2;
    return ` dy="${this.n(Math.round(d * px * 100))}"`;
  }
  private glyphPath(id: number, rec: GlyphRecord): void {
    if (this.sh.glyphs.has(id)) return;
    this.sh.glyphs.add(id);
    this.sh.defs.push(`<path id="p${id}" d="${esc(rec.path ?? "")}"/>`);
  }

  /** An annular sector: outer arc clockwise, inner arc back; a full ring as two halves each (Canvas2D's `arc`). */
  private sector(cx: number, cy: number, r0: number, r1: number, a0: number, a1: number): string {
    let sweep = a1 - a0;
    const full = sweep >= TAU;
    sweep = full ? TAU : ((sweep % TAU) + TAU) % TAU;
    if (sweep === 0) return "";
    const at = (r: number, a: number) => this.pt(cx + r * Math.cos(a), cy + r * Math.sin(a));
    const R1 = this.p(r1);
    const R0 = this.p(r0);
    if (full) {
      const half = a0 + Math.PI;
      return `M${at(r1, a0)}A${R1} ${R1} 0 1 1 ${at(r1, half)}A${R1} ${R1} 0 1 1 ${at(r1, a0)}L${at(r0, a0)}A${R0} ${R0} 0 1 0 ${at(r0, half)}A${R0} ${R0} 0 1 0 ${at(r0, a0)}Z`;
    }
    const large = sweep > Math.PI ? 1 : 0;
    return `M${at(r1, a0)}A${R1} ${R1} 0 ${large} 1 ${at(r1, a0 + sweep)}L${at(r0, a0 + sweep)}A${R0} ${R0} 0 ${large} 0 ${at(r0, a0)}Z`;
  }
  /** Canvas2D's `arc` as a path verb: a line from the current point to the arc's start, then the arc. */
  private arc(cx: number, cy: number, r: number, a0: number, a1: number, ccw: boolean, has: boolean): string {
    let sweep = ccw ? a0 - a1 : a1 - a0;
    const full = sweep >= TAU;
    sweep = full ? TAU : ((sweep % TAU) + TAU) % TAU;
    const dir = ccw ? -1 : 1;
    const at = (a: number) => this.pt(cx + r * Math.cos(a), cy + r * Math.sin(a));
    const R = this.p(r);
    const to = (s: number, large: number) => `A${R} ${R} 0 ${large} ${ccw ? 0 : 1} ${at(a0 + dir * s)}`;
    const start = `${has ? "L" : "M"}${at(a0)}`;
    if (sweep === 0) return start;
    return full ? `${start}${to(Math.PI, 1)}${to(TAU, 1)}` : start + to(sweep, sweep > Math.PI ? 1 : 0);
  }
  private roundRect(x: number, y: number, w: number, h: number, r: number): string {
    const R = this.p(r);
    const a = (px: number, py: number) => `A${R} ${R} 0 0 1 ${this.pt(px, py)}`;
    return `M${this.pt(x + r, y)}H${this.p(x + w - r)}${a(x + w, y + r)}V${this.p(y + h - r)}${a(x + w - r, y + h)}H${this.p(x + r)}${a(x, y + h - r)}V${this.p(y + r)}${a(x + r, y)}Z`;
  }
  /** A glyph shape as a closed subpath, the geometry Canvas2D's `shape()` builds (seqspec:format#themes, device px). */
  private shape(kind: number, x: number, y: number, r: number): string {
    const P = (px: number, py: number) => this.pt(px, py);
    switch (kind) {
      case 1: // diamond
        return `M${P(x, y - r)}L${P(x + r, y)}L${P(x, y + r)}L${P(x - r, y)}Z`;
      case 2: // triangle
        return `M${P(x, y - r)}L${P(x + r, y + r * 0.8)}L${P(x - r, y + r * 0.8)}Z`;
      case 3: // square
        return this.box(x - r * 0.8, y - r * 0.8, r * 1.6, r * 1.6);
      case 4: {
        // pin: a head over a point
        const hx = (a: number) => P(x + r * 0.65 * Math.cos(a), y - r * 0.35 + r * 0.65 * Math.sin(a));
        const R = this.p(r * 0.65);
        return `M${hx(Math.PI * 0.85)}A${R} ${R} 0 1 1 ${hx(Math.PI * 2.15)}L${P(x, y + r)}Z`;
      }
      case 5: {
        // star
        let d = "";
        for (let k = 0; k < 10; k++) {
          const a = -Math.PI / 2 + (k * Math.PI) / 5;
          const q = k % 2 === 0 ? r : r * 0.45;
          d += `${k === 0 ? "M" : "L"}${P(x + q * Math.cos(a), y + q * Math.sin(a))}`;
        }
        return `${d}Z`;
      }
      case 6: // cross
        return this.box(x - r, y - r * 0.25, 2 * r, r * 0.5) + this.box(x - r * 0.25, y - r, r * 0.5, 2 * r);
      case 7: // tick: one device px wide at any r
        return this.box(x - 0.5, y - r, 1, 2 * r);
      default: {
        // circle
        const c = r * 0.85;
        const R = this.p(c);
        return `M${P(x + c, y)}A${R} ${R} 0 1 1 ${P(x - c, y)}A${R} ${R} 0 1 1 ${P(x + c, y)}Z`;
      }
    }
  }
}

const BASELINE = ["alphabetic", "central", "text-before-edge", "text-after-edge"];
const weight = (f: Font) => (f.weight === "400" ? "" : ` font-weight="${f.weight}"`);
/** Six significant digits, shortest: scales and angles keep relative precision where hundredths would not. */
const sig = (v: number) => String(Number(v.toPrecision(6)));
