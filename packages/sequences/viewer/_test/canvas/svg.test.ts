/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SVG backend against Canvas2D in Chromium (seqspec:viewer#the-render-backend): the same batches drawn by both, the
// SVG rasterized as an image, compared by pixels: text anchoring and baselines first (alignment, baseline, angle and
// per-item scale of a run), then every batch kind's coverage. Rectangles are exact at DPR 1; curves and text differ by
// antialiasing and hinting.

import { describe, expect, it } from 'vitest';
import { probeFonts } from '../../scripts/font-probe.js';
import { Canvas2DBackend, Canvas2DShared } from '../../src/canvas/index.js';
import {
  type DrawResources,
  GlyphBatch,
  Pass,
  PathBatch,
  PathMode,
  PathVerb,
  RectBatch,
  type RenderBackend,
  SectorBatch,
  SHAPE_COUNT,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from '../../src/draw/api.js';
import { type FontMetrics, metricsMeasurer, SvgBackend } from '../../src/svg/index.js';

const GLYPHS = [
  { char: 'W', path: null, box: [1, 1] as const, fallback: 0 },
  { char: null, path: 'M0 0 L10 0 L0 10 Z', box: [10, 10] as const, fallback: 0 },
];
const resources: DrawResources = {
  font: (id) => (id === 1 ? '600 16px serif' : '400 20px sans-serif'),
  fontPx: (id) => (id === 1 ? 16 : 20),
  token: (t) => (t === SkinToken.ink ? '#112233' : 'rgb(43 127 255 / 0.5)'),
  glyph: (id) => GLYPHS[id - SHAPE_COUNT] ?? null,
};

/** The em boxes of this browser's fonts, read as the committed tables' generator reads them. */
const live = metricsMeasurer(
  probeFonts({ families: { 'sans-serif': 'sans-serif', serif: 'serif' }, weights: [400, 600] }).map(
    (p): FontMetrics => ({ ...p, unitsPerEm: 1000, advances: Uint16Array.from(p.advances) }),
  ),
);
const metrics = (font: string) => live.metrics(font);

/** The batches drawn on Canvas2D and as SVG, both as [w × h × 4] pixels at `dpr`. */
async function both(draw: (be: RenderBackend) => void, w: number, h: number, dpr = 1, baselines = true) {
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  const canvas = new Canvas2DBackend([g], new Canvas2DShared(resources), 0);
  canvas.beginPass(Pass.main, w, h, dpr);
  draw(canvas);
  canvas.endPass();
  const svg = new SvgBackend(resources, { metrics: baselines ? metrics : undefined });
  svg.beginPass(Pass.main, w, h, dpr);
  draw(svg);
  svg.endPass();
  const url = URL.createObjectURL(new Blob([await svg.toString()], { type: 'image/svg+xml' }));
  const img = new Image();
  img.src = url;
  await img.decode();
  const s = new OffscreenCanvas(w, h);
  const sg = s.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  sg.drawImage(img, 0, 0, w, h);
  URL.revokeObjectURL(url);
  return { canvas: g.getImageData(0, 0, w, h).data, svg: sg.getImageData(0, 0, w, h).data, text: await svg.toString() };
}
/** The box of pixels with alpha over 64: [x0, y0, x1, y1). */
function ink(data: Uint8ClampedArray, w: number): number[] {
  const box = [1e9, 1e9, -1, -1];
  for (let i = 3; i < data.length; i += 4)
    if (data[i] > 64) {
      const p = (i - 3) / 4;
      box[0] = Math.min(box[0], p % w);
      box[1] = Math.min(box[1], Math.floor(p / w));
      box[2] = Math.max(box[2], (p % w) + 1);
      box[3] = Math.max(box[3], Math.floor(p / w) + 1);
    }
  return box;
}
const coverage = (d: Uint8ClampedArray) => {
  let a = 0;
  for (let i = 3; i < d.length; i += 4) a += d[i];
  return a / 255;
};
const differing = (a: Uint8ClampedArray, b: Uint8ClampedArray) => {
  let n = 0;
  for (let i = 0; i < a.length; i += 4)
    if (Math.abs(a[i + 3] - b[i + 3]) > 8 || (a[i + 3] > 8 && Math.abs(a[i] - b[i]) > 8)) n++;
  return n;
};
const near = (a: number[], b: number[], tol: number) => a.every((v, i) => Math.abs(v - b[i]) <= tol);

function run(
  be: RenderBackend,
  align: number,
  baseline: number,
  o: { angle?: number; sx?: number; sy?: number; font?: number } = {},
) {
  const t = new TextBatch();
  t.style(o.font ?? 0, align, baseline, o.angle ?? 0);
  t.begin(0xcc0000);
  const k = t.add('Hxgé');
  t.x[k] = 80;
  t.y[k] = 60;
  t.sx[k] = o.sx ?? 1;
  t.sy[k] = o.sy ?? 1;
  be.text(t);
}

describe('text placed as Canvas2D places it', () => {
  const alignments = [TextAlign.left, TextAlign.center, TextAlign.right];
  const baselines = [TextBaseline.alphabetic, TextBaseline.middle, TextBaseline.top, TextBaseline.bottom];
  for (const dpr of [1, 2])
    for (const a of alignments)
      for (const b of baselines)
        it(`alignment ${a}, baseline ${b}, DPR ${dpr}: the same ink box`, async () => {
          const w = 160 * dpr;
          const { canvas, svg } = await both((be) => run(be, a, b), w, 120 * dpr, dpr);
          const [c, s] = [ink(canvas, w), ink(svg, w)];
          // Anchoring: every edge within a device px. The SVG is written with `dy` from this browser's own em box (the
          // tables' generator reads it the same way), so `middle`, `top` and `bottom` are Canvas2D's and not the
          // `dominant-baseline` ones, 0.05 to 0.1 em off (the next test).
          expect(near(c, s, 1), `canvas ${c} svg ${s}`).toBe(true);
        });

  it("without the em box the baselines are the reader's own (dominant-baseline): within 0.12 em, and dy is what closes it", async () => {
    // 20 px sans-serif: `central`, `text-before-edge` and `text-after-edge` sit 0.05, 0.1 and 0.05 em off (SV1a).
    for (const b of [TextBaseline.middle, TextBaseline.top, TextBaseline.bottom]) {
      const [off, on] = [
        await both((be) => run(be, TextAlign.left, b), 160, 120, 1, false),
        await both((be) => run(be, TextAlign.left, b), 160, 120),
      ];
      const [c, s, d] = [ink(off.canvas, 160), ink(off.svg, 160), ink(on.svg, 160)];
      expect(near(c, s, 0.12 * 20 + 1), `canvas ${c} dominant-baseline ${s}`).toBe(true);
      expect(near(c, d, 1), `canvas ${c} dy ${d}`).toBe(true);
      expect(off.text).toContain('dominant-baseline');
      expect(on.text).not.toContain('dominant-baseline');
    }
  });

  it('an angle and a per-item scale: the same box', async () => {
    for (const o of [{ angle: Math.PI / 2 }, { sx: 2, sy: 0.5 }, { angle: -Math.PI / 6, sx: 1.5, sy: 1.2 }]) {
      const { canvas, svg } = await both((be) => run(be, TextAlign.center, TextBaseline.middle, o), 160, 120);
      const [c, s] = [ink(canvas, 160), ink(svg, 160)];
      expect(near(c, s, 1), `${JSON.stringify(o)}: canvas ${c} svg ${s}`).toBe(true);
    }
  });

  it("a second font (weight, size, family) and its glyph in the batch's font", async () => {
    const { canvas, svg } = await both(
      (be) => {
        run(be, TextAlign.left, TextBaseline.alphabetic, { font: 1 });
        const g = new GlyphBatch();
        g.font = 1;
        g.begin(0x0000cc);
        g.push(40, 100, 9, SHAPE_COUNT);
        be.glyphs(g);
      },
      160,
      120,
    );
    expect(near(ink(canvas, 160), ink(svg, 160), 2)).toBe(true);
    expect(Math.abs(coverage(canvas) - coverage(svg)) / coverage(canvas)).toBeLessThan(0.15);
  });
});

describe('every batch kind covers what Canvas2D covers', () => {
  it('rects: exact, fills, alpha, an outline and a skin token', async () => {
    const { canvas, svg } = await both(
      (be) => {
        const b = new RectBatch();
        b.begin(0xff0000);
        b.push(2, 3, 20, 10);
        b.begin(0x00aa00, 128);
        b.push(10, 8, 20, 10);
        b.begin(0x0000ff, 255, 2);
        b.push(40, 10, 30, 20);
        b.begin(SkinToken.selection, 255);
        b.push(5, 30, 50, 8);
        be.rects(b);
      },
      80,
      50,
    );
    expect(differing(canvas, svg)).toBe(0);
  });

  it('sectors, paths and glyph shapes: the same coverage', async () => {
    const draw = (be: RenderBackend) => {
      const s = new SectorBatch();
      s.begin(0x884400, 200);
      s.push(40, 40, 10, 30, 0, 2);
      s.push(110, 40, 12, 30, 1, 1 + 2 * Math.PI);
      s.push(180, 40, 0, 30, 4, 1);
      be.sectors(s);
      const p = new PathBatch();
      p.begin(0x00aa44, 255, PathMode.fill);
      p.moveTo(10, 90);
      p.lineTo(40, 90);
      p.quadTo(50, 120, 20, 110);
      p.close();
      p.begin(0x00aa44, 255, PathMode.fill);
      p.moveTo(60, 100);
      p.lineTo(80, 100);
      p.lineTo(70, 120);
      p.close();
      p.begin(0x2244cc, 255, PathMode.stroke, 4);
      p.moveTo(100, 90);
      let o = p.verb(PathVerb.arc);
      p.coords.set([120, 100, 12, 0.5, 4, 1], o);
      o = p.verb(PathVerb.cubic);
      p.coords.set([150, 80, 160, 120, 190, 100], o);
      p.begin(0xcc2244, 255, PathMode.fill);
      o = p.verb(PathVerb.roundRect);
      p.coords.set([200, 90, 30, 20, 6], o);
      be.paths(p);
      const g = new GlyphBatch();
      g.font = 0;
      g.begin(0x222222);
      for (let k = 0; k < SHAPE_COUNT + 2; k++) g.push(16 + k * 26, 150, 9, k);
      be.glyphs(g);
    };
    const { canvas, svg } = await both(draw, 280, 170);
    // Antialiasing of curves differs (measured: 546 of 47,600 pixels, coverage within 0.6%); a missing shape differs by far more.
    expect(coverage(canvas)).toBeGreaterThan(6000);
    expect(differing(canvas, svg)).toBeLessThan(0.015 * 280 * 170);
    expect(Math.abs(coverage(canvas) - coverage(svg)) / coverage(canvas)).toBeLessThan(0.02);
  });

  it('images: smoothed or not, whole or a source rectangle', async () => {
    const px = new Uint32Array(4 * 4);
    for (let i = 0; i < 16; i++)
      px[i] = new Uint32Array(new Uint8Array([(i * 16) & 255, 255 - i * 16, 100, 255]).buffer)[0];
    const c = new OffscreenCanvas(60, 40);
    const g = c.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
    const shared = new Canvas2DShared(resources);
    const canvas = new Canvas2DBackend([g], shared, 0);
    const svg = new SvgBackend(resources);
    const r = new Float64Array([0, 0, 4, 4, 2, 2, 24, 24]);
    const r2 = new Float64Array([1, 1, 2, 2, 34, 4, 20, 20]);
    for (const be of [canvas, svg]) {
      be.beginPass(Pass.main, 60, 40, 1);
      const id = be.images.fromPixels(4, 4, px);
      be.image(id, r, false);
      be.image(id, r2, false);
      be.endPass();
    }
    const url = URL.createObjectURL(new Blob([await svg.toString()], { type: 'image/svg+xml' }));
    const img = new Image();
    img.src = url;
    await img.decode();
    const s = new OffscreenCanvas(60, 40);
    const sg = s.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
    sg.drawImage(img, 0, 0);
    expect(differing(g.getImageData(0, 0, 60, 40).data, sg.getImageData(0, 0, 60, 40).data)).toBeLessThan(40);
  });
});
