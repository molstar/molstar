/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { edt } from '@molstar/core/math/geometry/distance-transform';
import { createTextureImage, type TextureImage } from '@molstar/graphics/gl/renderable/util';
import { TextureChanges } from '@molstar/graphics/gl/renderable/texture-changes';
import { RUNNING_IN_NODEJS } from '@molstar/core/util/nodejs-shims';
import { LRUCache } from '@molstar/core/util/lru-cache';
import { TextMeasure } from '@molstar/core/util/text-measure';

/**
 * Atlases are tens of MB each (43 MB at the default quality), so the cache is keyed by the font (its CSS string) only
 * and keeps the few most recently used. An evicted atlas lives on in the texts built with it; texts built later get a
 * new one. With more fonts than that in use at once, rebuilding all their texts rebuilds atlases too (time, not memory;
 * .v6/designs/labels.md §1).
 */
const FontAtlasCacheCapacity = 4;
const FontAtlasCache = LRUCache.create<FontAtlas>(FontAtlasCacheCapacity);

/** Only the props the atlas depends on, whatever else `props` carries (all of Text's params, say). */
function fontAtlasProps(props: Partial<FontAtlasProps>): FontAtlasProps {
  const d = FontAtlasDefaults;
  return {
    fontFamily: props.fontFamily ?? d.fontFamily,
    fontQuality: props.fontQuality ?? d.fontQuality,
    fontStyle: props.fontStyle ?? d.fontStyle,
    fontVariant: props.fontVariant ?? d.fontVariant,
    fontWeight: props.fontWeight ?? d.fontWeight,
  };
}

function fontAtlasFontSize(fontQuality: number) {
  return 64 * (fontQuality + 1);
}

/** The CSS font the atlas for `props` draws its glyphs with (`ctx.font`), at the atlas's size. */
export function getFontAtlasFont(props: Partial<FontAtlasProps>) {
  const p = fontAtlasProps(props);
  const fonts = FontFamilyFonts[p.fontFamily] ?? p.fontFamily;
  return `${p.fontStyle} ${p.fontVariant} ${p.fontWeight} ${fontAtlasFontSize(p.fontQuality)}px ${fonts}`;
}

export function getFontAtlas(props: Partial<FontAtlasProps>) {
  const p = fontAtlasProps(props);
  const key = getFontAtlasFont(p);
  let atlas = LRUCache.get(FontAtlasCache, key);
  if (!atlas) {
    atlas = new FontAtlas(p);
    LRUCache.set(FontAtlasCache, key, atlas);
    // drawn with fallback glyphs while a web font loads: dropped when it arrives (`loadFont`)
    loadFont(key);
  }
  return atlas;
}

/**
 * The fonts each generic family draws with: a CSS font-family list, tried in order, the generic family last (Prism's
 * theme stacks for sans-serif and monospace; labels.md §2). An entry
 * may name a web font (CSS `@font-face`, a Google Fonts stylesheet, a `FontFace` added to `document.fonts`): it is
 * downloaded only when an atlas needs it (`loadFontAtlasFonts`). Change entries with `setFontFamilyFonts`.
 */
export const FontFamilyFonts: { readonly [F in FontFamily]: string } = {
  'sans-serif': '"Helvetica Neue", "Segoe UI", Helvetica, "Source Sans Pro", Arial, sans-serif',
  monospace: 'ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace',
  serif: 'Georgia, Cambria, "Times New Roman", "Liberation Serif", serif',
  cursive: 'cursive',
};

/**
 * Sets the fonts `family` draws with (a CSS font-family list, e.g. `'Inter, sans-serif'`). Texts built from then on use
 * them; existing texts keep their glyphs until rebuilt, so set the table before the first labels are drawn.
 */
export function setFontFamilyFonts(family: FontFamily, fonts: string) {
  (FontFamilyFonts as { [F in FontFamily]: string })[family] = fonts;
}

/**
 * How long text visuals wait for a web font; past it they draw a fallback, which stays until they rebuild after the font
 * arrived (labels.md §3).
 */
const FontLoadTimeoutMs = 3000;
/** What waiting on a font in flight resolves with: its arrival, or the timeout. */
const FontLoads = new Map<string, Promise<void>>();
/** Fonts whose download failed (or that loaded and still do not match): the fallback stays, nothing to wait for. */
const FailedFonts = new Set<string>();

function isFontLoaded(fontFaces: FontFaceSet, font: string) {
  try {
    return fontFaces.check(font);
  } catch {
    return true; // not a valid CSS font: nothing to load (the atlas draws the generic family)
  }
}

/** A promise while any font of `font` (a CSS font) is a web font not loaded yet; undefined in Node and once loaded or failed (labels.md §3). */
function loadFont(font: string): Promise<void> | undefined {
  const fontFaces = typeof document !== 'undefined' ? document.fonts : undefined;
  if (!fontFaces || FailedFonts.has(font)) return;
  const pending = FontLoads.get(font);
  if (pending) return pending;
  if (isFontLoaded(fontFaces, font)) return;

  // runs as the reaction that settles `arrived`, before anything awaiting it: a visual that waited builds with the
  // font. Any atlas cached for it was drawn before it arrived: while it was in flight, or before the page knew the
  // face at all (`check` is true for a family with no face yet), so it is dropped.
  const settle = () => {
    FontLoads.delete(font);
    clearTimeout(timer);
    if (!isFontLoaded(fontFaces, font)) {
      FailedFonts.add(font);
      return;
    }
    LRUCache.remove(FontAtlasCache, font);
  };
  let timer: ReturnType<typeof setTimeout> | undefined;
  const arrived = fontFaces.load(font).then(settle, settle);
  const waited = Promise.race([
    arrived,
    new Promise<void>((resolve) => (timer = setTimeout(resolve, FontLoadTimeoutMs))),
  ]);
  FontLoads.set(font, waited);
  return waited;
}

/**
 * Downloads the web fonts the atlas for `props` draws with, if any is not loaded yet; undefined when there is nothing
 * to wait for (system fonts, loaded fonts, Node). Await it before building texts so they get the font, not a fallback.
 */
export function loadFontAtlasFonts(props: Partial<FontAtlasProps>): Promise<void> | undefined {
  return loadFont(getFontAtlasFont(props));
}

/**
 * `build()`, after the fonts of the atlas for `props` are loaded (`loadFontAtlasFonts`), or failed, or 3 s passed;
 * synchronous when there is nothing to wait for.
 */
export function withFontAtlasFonts<T>(props: Partial<FontAtlasProps>, build: () => T): T | Promise<T> {
  const loading = loadFontAtlasFonts(props);
  return loading ? loading.then(build) : build();
}

/**
 * Drops every cached atlas and what is known about web fonts (in flight, failed): texts built from then on draw new
 * atlases. For tests, which restore the module state they change, and for an app that changed its fonts.
 */
export function clearFontAtlasCache() {
  while (FontAtlasCache.entries.count > 0) FontAtlasCache.entries.removeFirst();
  FontLoads.clear();
  FailedFonts.clear();
}

export type FontFamily = 'sans-serif' | 'monospace' | 'serif' | 'cursive';
const GenericFamilies = new Set<string>(['sans-serif', 'monospace', 'serif', 'cursive']);
export type FontStyle = 'normal' | 'italic' | 'oblique';
export type FontVariant = 'normal' | 'small-caps';
export type FontWeight = 'normal' | 'bold';

export const FontAtlasParams = {
  fontFamily: PD.Select('sans-serif', [
    ['sans-serif', 'Sans Serif'],
    ['monospace', 'Monospace'],
    ['serif', 'Serif'],
    ['cursive', 'Cursive'],
  ] as [FontFamily, string][]),
  fontQuality: PD.Select(3, [
    [0, 'lower'],
    [1, 'low'],
    [2, 'medium'],
    [3, 'high'],
    [4, 'higher'],
  ]),
  fontStyle: PD.Select('normal', [
    ['normal', 'Normal'],
    ['italic', 'Italic'],
    ['oblique', 'Oblique'],
  ] as [FontStyle, string][]),
  fontVariant: PD.Select('normal', [
    ['normal', 'Normal'],
    ['small-caps', 'Small Caps'],
  ] as [FontVariant, string][]),
  fontWeight: PD.Select('normal', [
    ['normal', 'Normal'],
    ['bold', 'Bold'],
  ] as [FontWeight, string][]),
};
export type FontAtlasParams = typeof FontAtlasParams;
export type FontAtlasProps = PD.Values<FontAtlasParams>;
const FontAtlasDefaults = PD.getDefaultValues(FontAtlasParams);

export type FontAtlasMap = {
  x: number;
  y: number;
  w: number;
  h: number;
  nw: number;
  nh: number; // normalized to lineheight
  /** Where the pen starts in the cell, from its left edge, normalized to lineheight: ink left of the pen fits. */
  left: number;
};

export class FontAtlas {
  readonly props: Readonly<FontAtlasProps>;
  readonly mapped: { [k: string]: FontAtlasMap } = {};
  readonly placeholder: FontAtlasMap;
  readonly texture: TextureImage<Uint8Array>;

  private scratchW = 0;
  private scratchH = 0;
  private scratchLeft = 0;
  private currentX = 0;
  private currentY = 0;
  private readonly scratchData: Uint8Array;

  private readonly cutoff = 0.5;
  readonly buffer: number;
  private readonly radius: number;

  private gridOuter: Float64Array;
  private gridInner: Float64Array;
  private f: Float64Array;
  private d: Float64Array;
  private z: Float64Array;
  private v: Int16Array;

  private scratchContext: CanvasRenderingContext2D;

  readonly lineHeight: number;

  /** The widest glyph cell: an em dash, a wide emoji, with the buffer on both sides. */
  private readonly maxWidth: number;
  private readonly middle: number;
  /** Widths of text in this font, normalized to lineheight (`layout`). */
  private readonly measure: TextMeasure;

  constructor(props: Partial<FontAtlasProps> = {}) {
    const p = fontAtlasProps(props);
    this.props = p;

    // create measurements
    const fontSize = fontAtlasFontSize(p.fontQuality);
    this.buffer = fontSize / 8;
    this.radius = fontSize / 3;
    this.lineHeight = Math.round(fontSize + 2 * this.buffer + this.radius);
    this.maxWidth = Math.round(this.lineHeight * 1.5);

    // create texture (for ~350 characters, a cell 0.75 lineheight wide on average)
    this.texture = createTextureImage(350 * this.lineHeight * Math.round(this.lineHeight * 0.75), 1, Uint8Array);

    // prepare scratch canvas
    this.scratchContext = createCanvasContext(this.maxWidth, this.lineHeight, { willReadFrequently: true })!;

    this.scratchContext.font = getFontAtlasFont(p);
    if (!String(this.scratchContext.font).includes(`${fontSize}px`)) {
      // not a valid CSS font, so the canvas kept its 10px default: draw the generic family instead
      const generic = GenericFamilies.has(p.fontFamily) ? p.fontFamily : 'sans-serif';
      console.warn(`Font atlas: invalid CSS font '${getFontAtlasFont(p)}', drawing ${generic}`);
      this.scratchContext.font = `${p.fontStyle} ${p.fontVariant} ${p.fontWeight} ${fontSize}px ${generic}`;
    }
    this.scratchContext.fillStyle = 'black';
    this.scratchContext.textBaseline = 'middle';
    // kerning on where the canvas has the setting (not node-canvas, which kerns by default)
    if ('fontKerning' in this.scratchContext) this.scratchContext.fontKerning = 'normal';
    const ctx = this.scratchContext;
    const lineHeight = this.lineHeight;
    this.measure = new TextMeasure((text) => ctx.measureText(text).width / lineHeight);

    // SDF scratch values
    this.scratchData = new Uint8Array(this.lineHeight * this.maxWidth);

    // temporary arrays for the distance transform
    this.gridOuter = new Float64Array(this.lineHeight * this.maxWidth);
    this.gridInner = new Float64Array(this.lineHeight * this.maxWidth);
    this.f = new Float64Array(Math.max(this.lineHeight, this.maxWidth));
    this.d = new Float64Array(Math.max(this.lineHeight, this.maxWidth));
    this.z = new Float64Array(Math.max(this.lineHeight, this.maxWidth) + 1);
    this.v = new Int16Array(Math.max(this.lineHeight, this.maxWidth));

    this.middle = Math.ceil(this.lineHeight / 2);

    // replacement Character
    this.placeholder = this.get(String.fromCharCode(0xfffd));
  }

  get(char: string) {
    if (this.mapped[char] === undefined) {
      this.draw(char);

      const { array, width, height } = this.texture;
      const data = this.scratchData;

      if (this.currentX + this.scratchW > width) {
        this.currentX = 0;
        this.currentY += this.scratchH;
      }
      if (this.currentY + this.scratchH > height) {
        console.warn('canvas to small');
        return this.placeholder;
      }

      this.mapped[char] = {
        x: this.currentX,
        y: this.currentY,
        w: this.scratchW,
        h: this.scratchH,
        nw: this.scratchW / this.lineHeight,
        nh: this.scratchH / this.lineHeight,
        left: this.scratchLeft / this.lineHeight,
      };

      for (let y = 0; y < this.scratchH; ++y) {
        for (let x = 0; x < this.scratchW; ++x) {
          array[width * (this.currentY + y) + this.currentX + x] = data[y * this.scratchW + x];
        }
      }

      this.currentX += this.scratchW;
      // written in place: shared textures of the atlas upload these rows again (`.v6/plans/gl/shared-textures.md` §3)
      TextureChanges.writeRows(this.texture, this.currentY, this.currentY + this.scratchH);
    }

    return this.mapped[char];
  }

  /**
   * The graphemes of `text` and where each glyph's pen starts, normalized to lineheight: advances measured in this
   * font with the kerning of each pair (`TextMeasure`). Place a glyph's cell `get(grapheme).left` before its pen.
   */
  layout(text: string) {
    return this.measure.layout(text);
  }

  draw(char: string) {
    const h = this.lineHeight;
    const ctx = this.scratchContext;
    const data = this.scratchData;

    // the cell holds the ink, which may reach left of the pen or right of the advance ('j', 'f', italics), and the
    // buffer for the distance field on both sides
    const m = ctx.measureText(char);
    const inkLeft = Math.max(0, Math.ceil(m.actualBoundingBoxLeft ?? 0));
    const inkRight = Math.max(m.width, m.actualBoundingBoxRight ?? 0);
    const w = Math.min(this.maxWidth, Math.ceil(inkLeft + inkRight + 2 * this.buffer));
    const n = w * h;
    const x = this.buffer + inkLeft;

    ctx.clearRect(0, 0, w, h); // clear scratch area
    ctx.fillText(char, x, this.middle); // draw text
    const imageData = ctx.getImageData(0, 0, w, h);

    for (let i = 0; i < n; i++) {
      const a = imageData.data[i * 4 + 3] / 255; // alpha value
      this.gridOuter[i] = a === 1 ? 0 : a === 0 ? Number.MAX_SAFE_INTEGER : Math.pow(Math.max(0, 0.5 - a), 2);
      this.gridInner[i] = a === 1 ? Number.MAX_SAFE_INTEGER : a === 0 ? 0 : Math.pow(Math.max(0, a - 0.5), 2);
    }

    edt(this.gridOuter, w, h, this.f, this.d, this.v, this.z);
    edt(this.gridInner, w, h, this.f, this.d, this.v, this.z);

    for (let i = 0; i < n; i++) {
      const d = this.gridOuter[i] - this.gridInner[i];
      data[i] = Math.max(0, Math.min(255, Math.round(255 - 255 * (d / this.radius + this.cutoff))));
    }

    this.scratchW = w;
    this.scratchH = h;
    this.scratchLeft = x;
  }
}

/** Type of imported `canvas` module (not using `typeof import('canvas')` to avoid missing types) */
type CanvasModule = any;
let _canvas: CanvasModule | undefined;
function getCanvasModule(): CanvasModule {
  if (!_canvas)
    throw new Error(
      "When running in Node.js and wanting to use Canvas API, call mol-util/data-source's setCanvasModule function first and pass imported `canvas` module to it.",
    );
  return _canvas;
}
/** Set `canvas` module, before using Canvas API functionality in NodeJS. Usage: `setCanvasModule(require('canvas')); // some code `*/
export function setCanvasModule(canvas: CanvasModule) {
  _canvas = canvas;
}
/** Return a newly created canvas context (using a canvas HTML element in browser, canvas module in NodeJS) */
function createCanvasContext(
  width: number,
  height: number,
  options?: CanvasRenderingContext2DSettings,
): CanvasRenderingContext2D | null {
  if (RUNNING_IN_NODEJS) {
    const canvas = getCanvasModule().createCanvas(width, height);
    return canvas.getContext('2d', options) as unknown as CanvasRenderingContext2D;
  } else {
    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;
    return canvas.getContext('2d', options);
  }
}
