/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Canvas2D backend in Chromium (seqspec:viewer#the-render-backend): each batch kind's pixels at DPR 1, 2
// and 3; grouping (style switches at most the distinct colors); a back pass shows only what is presented, stale
// regions under presented ones; copies read the pass's image as it was at its first copy; a recorded frame replayed
// draws the same pixels as the frame drawn directly.

import { describe, expect, it } from 'vitest';
import { Canvas2DBackend, Canvas2DShared } from '../../src/canvas/index.js';
import {
  type DrawResources,
  GlyphBatch,
  GlyphShape,
  NO_COLOR,
  Pass,
  PathBatch,
  PathMode,
  PathVerb,
  RectBatch,
  SectorBatch,
  SHAPE_COUNT,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from '../../src/draw/api.js';
import { RecordingBackend } from '../../src/testing/index.js';
import { drawFrame, Frame, H, PIXELS, W } from '../support/frame.js';

const GLYPHS = [
  { char: 'W', path: null, box: [1, 1] as const, fallback: 0 },
  { char: null, path: 'M0 0 L10 0 L0 10 Z', box: [10, 10] as const, fallback: 0 },
];
const resources: DrawResources = {
  font: () => '20px sans-serif',
  fontPx: () => 20,
  token: (t) => (t === SkinToken.ink ? '#112233' : '#ff00ff'),
  glyph: (id) => GLYPHS[id - SHAPE_COUNT] ?? null,
};

/** A backend over one offscreen canvas per pass (`back`: the passes with a back target). */
function setup(w: number, h: number, back = 0) {
  const canvases = [0, 1].map(() => new OffscreenCanvas(w, h));
  const gs = canvases.map((c) => c.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D);
  const be = new Canvas2DBackend(gs, new Canvas2DShared(resources), back);
  const px = (pass = 0) => gs[pass].getImageData(0, 0, w, h).data;
  /** [r, g, b, a] at (x, y). */
  const at = (x: number, y: number, pass = 0) => [...gs[pass].getImageData(x, y, 1, 1).data];
  return { be, gs, px, at };
}
/** The box of pixels with alpha over 64: [x0, y0, x1, y1). */
function ink(data: Uint8ClampedArray, w: number): number[] {
  const box = [1e9, 1e9, -1, -1];
  for (let i = 3; i < data.length; i += 4)
    if (data[i] > 64) {
      const p = (i - 3) / 4;
      const x = p % w;
      const y = Math.floor(p / w);
      box[0] = Math.min(box[0], x);
      box[1] = Math.min(box[1], y);
      box[2] = Math.max(box[2], x + 1);
      box[3] = Math.max(box[3], y + 1);
    }
  return box;
}
const count = (data: Uint8ClampedArray, rgb: number[]) => {
  let n = 0;
  for (let i = 0; i < data.length; i += 4)
    if (data[i] === rgb[0] && data[i + 1] === rgb[1] && data[i + 2] === rgb[2] && data[i + 3] === 255) n++;
  return n;
};

describe.each([1, 2, 3])('pixels at DPR %i', (d) => {
  it('rects fill and outline whole device px; alpha, skin tokens and NO_COLOR', () => {
    const { be, px, at } = setup(40 * d, 30 * d);
    const b = new RectBatch();
    be.beginPass(Pass.main, 40 * d, 30 * d, d);
    b.begin(0xff0000);
    b.push(2 * d, 3 * d, 4 * d, 5 * d);
    b.begin(0x0000ff, 255, d);
    b.push(10 * d, 10 * d, 8 * d, 6 * d);
    b.begin(0x00ff00, 128);
    b.push(30 * d, 0, d, d);
    b.begin(SkinToken.ink);
    b.push(0, 20 * d, d, d);
    b.begin(NO_COLOR);
    b.push(0, 25 * d, 5 * d, 5 * d);
    be.rects(b);
    be.endPass();
    const data = px();
    expect(count(data, [255, 0, 0])).toBe(20 * d * d);
    expect(at(2 * d, 3 * d)).toEqual([255, 0, 0, 255]);
    expect(at(6 * d, 8 * d - 1)[3]).toBe(0);
    // The outline inside the rect, d px wide.
    expect(count(data, [0, 0, 255])).toBe(48 * d * d - 6 * d * (4 * d));
    expect(at(10 * d, 10 * d)).toEqual([0, 0, 255, 255]);
    expect(at(11 * d, 11 * d)[3]).toBe(0);
    expect(at(30 * d, 0)[3]).toBeGreaterThan(120);
    expect(at(30 * d, 0)[3]).toBeLessThan(136);
    expect(at(0, 20 * d)).toEqual([0x11, 0x22, 0x33, 255]);
    expect(at(d, 26 * d)[3]).toBe(0);
  });

  it('sectors fill the ring between their radii, clockwise from a0 to a1', () => {
    const { be, at } = setup(60 * d, 60 * d);
    const s = new SectorBatch();
    be.beginPass(Pass.main, 60 * d, 60 * d, d);
    s.begin(0xff0000);
    s.push(30 * d, 30 * d, 10 * d, 20 * d, 0, Math.PI / 2);
    be.sectors(s);
    be.endPass();
    const q = (r: number, a: number) =>
      at(Math.round(30 * d + r * d * Math.cos(a)), Math.round(30 * d + r * d * Math.sin(a)));
    expect(q(15, Math.PI / 4)).toEqual([255, 0, 0, 255]);
    expect(q(15, -Math.PI / 4)[3]).toBe(0);
    expect(q(5, Math.PI / 4)[3]).toBe(0);
    expect(q(25, Math.PI / 4)[3]).toBe(0);
  });

  it('paths fill and stroke, with round rects', () => {
    const { be, at } = setup(60 * d, 60 * d);
    const p = new PathBatch();
    be.beginPass(Pass.main, 60 * d, 60 * d, d);
    p.begin(0xff0000);
    p.moveTo(2 * d, 2 * d);
    p.lineTo(20 * d, 2 * d);
    p.lineTo(2 * d, 20 * d);
    p.close();
    p.begin(0x0000ff, 255, PathMode.stroke, 2 * d);
    p.moveTo(0, 40 * d + d / 2);
    p.lineTo(60 * d, 40 * d + d / 2);
    p.begin(0x00ff00);
    const o = p.verb(PathVerb.roundRect);
    p.coords.set([30 * d, 10 * d, 20 * d, 20 * d, 6 * d], o);
    be.paths(p);
    be.endPass();
    expect(at(5 * d, 5 * d)).toEqual([255, 0, 0, 255]);
    expect(at(18 * d, 18 * d)[3]).toBe(0);
    for (let y = 40 * d; y < 41 * d; y++) expect(at(30 * d, y)).toEqual([0, 0, 255, 255]);
    expect(at(30 * d, 40 * d - 1)[3]).toBe(0);
    expect(at(30 * d, 41 * d)[3]).toBe(0);
    expect(at(40 * d, 20 * d)).toEqual([0, 255, 0, 255]);
    expect(at(30 * d, 10 * d)[3]).toBe(0);
  });

  it("glyphs: shapes, a character centered and 2r high in the batch's font, a path in its box", () => {
    const one = (id: number, r: number) => {
      const { be, px } = setup(80 * d, 80 * d);
      const g = new GlyphBatch();
      be.beginPass(Pass.main, 80 * d, 80 * d, d);
      g.begin(0x000000);
      g.push(40 * d, 40 * d, r * d, id);
      be.glyphs(g);
      be.endPass();
      return ink(px(), 80 * d);
    };
    const c = one(GlyphShape.circle, 10);
    expect(c[0]).toBeGreaterThanOrEqual(31 * d - 1);
    expect(c[2]).toBeLessThanOrEqual(49 * d + 1);
    expect(one(GlyphShape.square, 10)).toEqual([32 * d, 32 * d, 48 * d, 48 * d]);
    for (const r of [8, 16]) {
      const w = one(SHAPE_COUNT, r);
      expect(Math.abs((w[0] + w[2]) / 2 - 40 * d)).toBeLessThanOrEqual(1.5 * d);
      expect(Math.abs((w[1] + w[3]) / 2 - 40 * d)).toBeLessThanOrEqual(0.2 * r * d);
      expect(w[3] - w[1]).toBeGreaterThan(r * d * 0.9);
      expect(w[3] - w[1]).toBeLessThanOrEqual(2 * r * d);
    }
    // The triangle path's box, 2r on its larger side, centered: its corners at the box's top-left, top-right and bottom-left.
    expect(one(SHAPE_COUNT + 1, 10)).toEqual([30 * d, 30 * d, 50 * d, 50 * d]);
    const { be, at } = setup(80 * d, 80 * d);
    const g = new GlyphBatch();
    be.beginPass(Pass.main, 80 * d, 80 * d, d);
    g.begin(0xff0000);
    g.push(40 * d, 40 * d, 10 * d, SHAPE_COUNT + 1);
    be.glyphs(g);
    be.endPass();
    expect(at(30 * d + 1, 30 * d + 1)).toEqual([255, 0, 0, 255]);
    expect(at(50 * d - 2, 50 * d - 2)[3]).toBe(0);
  });

  it('text in CSS px scaled by the DPR, placed at its anchor by alignment and baseline, and rotated', () => {
    const run = (align: number, angle: number) => {
      const { be, px } = setup(120 * d, 120 * d);
      const t = new TextBatch();
      be.beginPass(Pass.main, 120 * d, 120 * d, d);
      t.style(0, align, TextBaseline.alphabetic, angle);
      t.begin(0x000000);
      t.push(60 * d, 60 * d, 'MMM');
      be.text(t);
      be.endPass();
      return ink(px(), 120 * d);
    };
    const left = run(TextAlign.left, 0);
    expect(left[0]).toBeGreaterThanOrEqual(60 * d);
    expect(left[0]).toBeLessThan(63 * d);
    expect(left[3]).toBeLessThanOrEqual(60 * d + 1);
    expect(left[3] - left[1]).toBeGreaterThan(10 * d);
    expect(left[3] - left[1]).toBeLessThan(20 * d);
    const center = run(TextAlign.center, 0);
    expect(Math.abs((center[0] + center[2]) / 2 - 60 * d)).toBeLessThanOrEqual(2 * d);
    const right = run(TextAlign.right, 0);
    expect(right[2]).toBeLessThanOrEqual(60 * d + 1);
    // Clockwise a quarter turn: the run reads downward from its anchor, the glyphs' tops facing right.
    const down = run(TextAlign.left, Math.PI / 2);
    expect(down[3] - down[1]).toBeGreaterThan(down[2] - down[0]);
    expect(down[1]).toBeGreaterThanOrEqual(60 * d - 1);
    expect(down[0]).toBeGreaterThanOrEqual(60 * d - 1);
  });

  it("a text item's scale stretches it about its anchor, along and across the run, before the angle (a logo)", () => {
    const run = (sx: number, sy: number, angle = 0) => {
      const { be, px } = setup(160 * d, 160 * d);
      const t = new TextBatch();
      be.beginPass(Pass.main, 160 * d, 160 * d, d);
      t.style(0, TextAlign.left, TextBaseline.alphabetic, angle);
      t.begin(0x000000);
      const i = t.add('M');
      t.x[i] = t.y[i] = 80 * d;
      t.sx[i] = sx;
      t.sy[i] = sy;
      be.text(t);
      be.endPass();
      return ink(px(), 160 * d);
    };
    const [w, h] = [(b: number[]) => b[2] - b[0], (b: number[]) => b[3] - b[1]];
    const plain = run(1, 1);
    const wide = run(3, 0.5);
    expect(w(wide) / w(plain)).toBeCloseTo(3, 0);
    expect(h(wide) / h(plain)).toBeCloseTo(0.5, 0);
    // Standing on its anchor's baseline, from its left edge.
    expect(Math.abs(wide[3] - plain[3])).toBeLessThanOrEqual(1);
    expect(Math.abs(wide[0] - plain[0])).toBeLessThanOrEqual(3 * d);
    // A quarter turn: the stretch follows the run, downward.
    const down = run(3, 0.5, Math.PI / 2);
    expect(h(down) / w(plain)).toBeCloseTo(3, 0);
  });

  it('images from pixels, unsmoothed, and a layer drawn by renderers', () => {
    const { be, at } = setup(40 * d, 40 * d);
    const img = be.images.fromPixels(2, 2, PIXELS);
    const l = be.images.beginLayer(4, 4);
    const r = new RectBatch();
    r.begin(0x00ffff);
    r.push(0, 0, 4, 4);
    l.rects(r);
    const layer = be.images.endLayer();
    be.beginPass(Pass.main, 40 * d, 40 * d, d);
    be.image(img, Float64Array.of(0, 0, 2, 2, 0, 0, 4 * d, 4 * d), false);
    be.image(layer, Float64Array.of(0, 0, 4, 4, 20 * d, 20 * d, 4 * d, 4 * d), false);
    be.endPass();
    expect(at(0, 0)).toEqual([255, 0, 0, 255]);
    expect(at(4 * d - 1, 0)).toEqual([0, 255, 0, 255]);
    expect(at(0, 4 * d - 1)).toEqual([0, 0, 255, 255]);
    expect(at(4 * d - 1, 4 * d - 1)).toEqual([255, 255, 255, 255]);
    expect(at(22 * d, 22 * d)).toEqual([0, 255, 255, 255]);
    be.images.dispose(img);
    be.beginPass(Pass.main, 40 * d, 40 * d, d);
    be.clear(0, 0, 40 * d, 40 * d);
    be.image(img, Float64Array.of(0, 0, 2, 2, 0, 0, 4 * d, 4 * d), false);
    be.endPass();
    expect(at(0, 0)[3]).toBe(0);
  });
});

describe('grouping', () => {
  it("a frame of 1,500 distinct colors (a continuous palette's row) drawn again makes no color string again", () => {
    const { be } = setup(200, 50);
    const b = new RectBatch();
    for (let i = 0; i < 1500; i++) {
      const t = ((i * 7919) % 1500) / 1500;
      b.begin((Math.round(68 + 185 * t) << 16) | (Math.round(1 + 230 * t) << 8) | Math.round(84 - 47 * t));
      b.push(i % 200, 0, 1, 10);
    }
    const frame = () => {
      be.beginPass(Pass.main, 200, 50, 1);
      be.rects(b);
    };
    frame();
    const orig = Number.prototype.toString;
    let made = 0;
    Number.prototype.toString = function (this: number, radix?: number) {
      if (radix === 16) made++;
      return orig.call(this, radix);
    };
    try {
      frame();
    } finally {
      Number.prototype.toString = orig;
    }
    expect(made).toBe(0);
  });

  it('switches style at most once per distinct color, however the items interleave', () => {
    const { be, at } = setup(200, 50);
    const b = new RectBatch();
    be.beginPass(Pass.main, 200, 50, 1);
    for (let i = 0; i < 1000; i++) {
      b.begin(i % 2 === 0 ? 0xff0000 : 0x0000ff);
      b.push(i % 200, Math.floor(i / 200), 1, 1);
    }
    be.styles = 0;
    be.rects(b);
    expect(be.styles).toBe(2);
    expect(at(0, 0)).toEqual([255, 0, 0, 255]);
    expect(at(1, 0)).toEqual([0, 0, 255, 255]);
    b.reset();
    for (let i = 0; i < 1000; i++) {
      b.begin((i * 37) % 50);
      b.push(i % 200, 10 + Math.floor(i / 200), 1, 1);
    }
    be.styles = 0;
    be.rects(b);
    expect(be.styles).toBe(50);
    const t = new TextBatch();
    t.style(0, TextAlign.left, TextBaseline.top);
    for (let i = 0; i < 300; i++) {
      t.begin(i % 3);
      t.push(i % 200, 30, 'x');
    }
    be.styles = 0;
    be.text(t);
    expect(be.styles).toBe(3);
    be.endPass();
  });
});

describe('grouping: the table is a power of two', () => {
  it("a batch of 1,469 items in the logo's twenty colors groups and draws (RN2: it once looped forever)", () => {
    const { be, at } = setup(1500, 4);
    const b = new RectBatch();
    // The colors of a logo over twenty amino acids, as the Taylor scheme gives them.
    const colors = [
      0xffafaf, 0xffc800, 0x6464ff, 0xff0000, 0x00ff00, 0xff00ff, 0xffff00, 0x6464ff, 0xff6600, 0x0066ff, 0xcc00ff,
      0xffcc00, 0x00cc66, 0x990000, 0x009999, 0x666699, 0xcc6600, 0x99ff00, 0x0099ff, 0xff0099,
    ];
    be.beginPass(Pass.main, 1500, 4, 1);
    for (let i = 0; i < 1469; i++) {
      b.begin(colors[(i * 7) % colors.length]);
      b.push(i, 0, 1, 4);
    }
    be.styles = 0;
    be.rects(b);
    expect(be.styles).toBe(colors.length - 1);
    expect(at(0, 0)).toEqual([255, 175, 175, 255]);
    be.endPass();
  });
});

describe('presentation', () => {
  const fill = (be: Canvas2DBackend, color: number, x: number, y: number, w: number, h: number) => {
    const b = new RectBatch();
    b.begin(color);
    b.push(x, y, w, h);
    be.rects(b);
  };

  it('a back pass shows nothing until a region is presented, at endPass, and only that region', () => {
    const { be, at } = setup(40, 20, 1 << Pass.main);
    be.beginPass(Pass.main, 40, 20, 1);
    fill(be, 0xff0000, 0, 0, 40, 20);
    be.endPass();
    expect(at(5, 5)[3]).toBe(0);
    be.beginPass(Pass.main, 40, 20, 1);
    fill(be, 0x0000ff, 0, 0, 40, 20);
    be.present(0, 0, 10, 20);
    expect(at(5, 5)[3]).toBe(0);
    be.endPass();
    expect(at(5, 5)).toEqual([0, 0, 255, 255]);
    expect(at(15, 5)[3]).toBe(0);
  });

  it("stale maps the last presented image to the new view, under that frame's presented regions", () => {
    const { be, at } = setup(40, 20, 1 << Pass.main);
    be.beginPass(Pass.main, 40, 20, 1);
    fill(be, 0xff0000, 0, 0, 10, 20);
    fill(be, 0x00ff00, 10, 0, 10, 20);
    be.present(0, 0, 40, 20);
    be.endPass();
    // A zoom by 2 around x = 0: the old [0, 20) shows over [0, 40) until each entry is drawn whole.
    be.beginPass(Pass.main, 40, 20, 1);
    fill(be, 0x0000ff, 0, 0, 40, 20);
    be.stale(0, 0, 20, 10, 0, 0, 40, 10);
    be.stale(0, 10, 20, 10, 0, 10, 40, 10);
    be.present(0, 10, 40, 10);
    be.endPass();
    expect(at(5, 5)).toEqual([255, 0, 0, 255]);
    expect(at(15, 5)).toEqual([255, 0, 0, 255]);
    expect(at(25, 5)).toEqual([0, 255, 0, 255]);
    expect(at(25, 15)).toEqual([0, 0, 255, 255]);
  });

  it('a pass without a back target draws at once and ignores present and stale', () => {
    const { be, at } = setup(40, 20, 1 << Pass.main);
    be.beginPass(Pass.overlay, 40, 20, 1);
    fill(be, 0xff0000, 0, 0, 10, 10);
    expect(at(5, 5, 1)).toEqual([255, 0, 0, 255]);
    be.stale(0, 0, 10, 10, 20, 0, 10, 10);
    be.present(0, 0, 40, 20);
    be.endPass();
    expect(at(25, 5, 1)[3]).toBe(0);
  });

  it("copies read the pass's image as it was at its first copy", () => {
    const { be, at } = setup(20, 30);
    be.beginPass(Pass.main, 20, 30, 1);
    fill(be, 0xff0000, 0, 0, 20, 5);
    fill(be, 0x00ff00, 0, 5, 20, 5);
    fill(be, 0x0000ff, 0, 10, 20, 20);
    be.endPass();
    be.beginPass(Pass.main, 20, 30, 1);
    be.copy(0, 5, 0, 0, 20, 10);
    be.copy(0, 5, 0, 10, 20, 20);
    be.endPass();
    expect(at(5, 2)[3]).toBe(0);
    expect(at(5, 7)).toEqual([255, 0, 0, 255]);
    expect(at(5, 12)).toEqual([0, 255, 0, 255]);
    expect(at(5, 17)).toEqual([0, 0, 255, 255]);
  });

  it('a clip holds the draws inside it, and the style after it is set again', () => {
    const { be, at } = setup(20, 20);
    be.beginPass(Pass.main, 20, 20, 1);
    fill(be, 0xff0000, 0, 0, 2, 2);
    be.pushClip(5, 5, 5, 5);
    fill(be, 0x00ff00, 0, 0, 20, 20);
    be.popClip();
    fill(be, 0x00ff00, 15, 15, 2, 2);
    be.endPass();
    expect(at(7, 7)).toEqual([0, 255, 0, 255]);
    expect(at(12, 12)[3]).toBe(0);
    expect(at(15, 15)).toEqual([0, 255, 0, 255]);
  });
});

describe('recording', () => {
  it('a recorded frame replayed into Canvas2D draws the pixels the frame drew directly', () => {
    const direct = setup(W, H, 1 << Pass.main);
    const replayed = setup(W, H, 1 << Pass.main);
    const f = new Frame(2000, 400);
    for (const d of [1, 2]) {
      drawFrame(direct.be, f, 0, direct.be.images.fromPixels(2, 2, PIXELS), d);
      const rec = new RecordingBackend();
      drawFrame(rec, f, 0, rec.images.fromPixels(2, 2, PIXELS), d);
      rec.replay(replayed.be);
      const a = direct.px();
      const b = replayed.px();
      expect(ink(a, W)).toEqual([0, 0, W, H - 1]);
      let diff = 0;
      for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) diff++;
      expect(diff).toBe(0);
    }
  });
});
