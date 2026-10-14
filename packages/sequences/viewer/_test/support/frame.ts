/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A frame through a backend, as the composer and renderers will draw one (seqspec:viewer#the-render-backend): a back
// pass, a pan's copy after the first frame, a clip, every batch kind (thousands of rects and text runs), an image, the
// presented region. Batches are owned here and reused by every frame, as a renderer module owns them. Doubles go
// through `add()`'s slot.

import {
  GlyphBatch,
  GlyphShape,
  Pass,
  PathBatch,
  PathMode,
  PathVerb,
  RectBatch,
  type RenderBackend,
  SectorBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from '../../src/draw/api.js';

export const W = 400;
export const H = 120;
const COLORS = [0xe41a1c, 0x377eb8, 0x4daf4a, 0x984ea3, 0xff7f00, 0xffff33, 0xa65628, SkinToken.ink];
/** One-character strings, made once: a text run holds a reference. */
export const CHARS = Array.from({ length: 128 }, (_, i) => String.fromCharCode(i));

export class Frame {
  readonly rects = new RectBatch();
  readonly sectors = new SectorBatch();
  readonly paths = new PathBatch();
  readonly glyphs = new GlyphBatch();
  readonly text = new TextBatch();
  readonly at = new Float64Array(8);
  readonly nRects: number;
  readonly nText: number;
  constructor(nRects: number, nText: number) {
    this.nRects = nRects;
    this.nText = nText;
  }
}

/** Frame `step` of a pan by one device px per frame; `img`: an image id of the backend's. */
export function drawFrame(be: RenderBackend, f: Frame, step: number, img: number, dpr = 1): void {
  be.beginPass(Pass.main, W, H, dpr);
  if (step > 0) be.copy(-1, 0, 0, 0, W, H);
  be.pushClip(0, 0, W, H);
  const r = f.rects;
  r.reset();
  for (let i = 0; i < f.nRects; i++) {
    if ((i & 63) === 0) r.begin(COLORS[(i >> 6) & 7], i & 64 ? 255 : 160);
    r.push((i * 3 + step) % W, (i >> 7) % 100, 2, 10);
  }
  r.begin(SkinToken.focusRing, 255, 1);
  r.push(1, 1, W - 2, H - 2);
  be.rects(r);
  const t = f.text;
  t.reset();
  t.style(0, TextAlign.center, TextBaseline.middle);
  for (let i = 0; i < f.nText; i++) {
    if ((i & 31) === 0) t.begin(COLORS[(i >> 5) & 7]);
    const k = t.add(CHARS[65 + (i % 26)]);
    t.x[k] = ((i * 7 + step) % W) + 0.5;
    t.y[k] = 20 + (i % 9) * 10;
  }
  be.text(t);
  const g = f.glyphs;
  g.reset();
  g.begin(0x222222);
  for (let i = 0; i < 64; i++) {
    const k = g.add(i % 10 < 8 ? i % 8 : GlyphShape.star);
    g.x[k] = 6 + i * 6;
    g.y[k] = 110;
    g.radius[k] = 2.5;
  }
  be.glyphs(g);
  const s = f.sectors;
  s.reset();
  s.begin(0x377eb8, 200);
  const k = s.add();
  s.cx[k] = 60;
  s.cy[k] = 60;
  s.r0[k] = 20;
  s.r1[k] = 30;
  s.a0[k] = 0.25 + step * 0.01;
  s.a1[k] = 2;
  be.sectors(s);
  const p = f.paths;
  p.reset();
  p.begin(0x4daf4a, 255, PathMode.stroke, 3);
  p.moveTo(0, 50.5);
  p.quadTo(W / 2, 10.5, W, 50.5);
  p.begin(0xa65628, 128);
  const o = p.verb(PathVerb.roundRect);
  p.coords[o] = 300;
  p.coords[o + 1] = 70;
  p.coords[o + 2] = 40;
  p.coords[o + 3] = 20;
  p.coords[o + 4] = 6;
  be.paths(p);
  const at = f.at;
  at[0] = at[1] = 0;
  at[2] = at[3] = 2;
  at[4] = 200;
  at[5] = 80;
  at[6] = at[7] = 20;
  be.image(img, at, false);
  be.popClip();
  be.present(0, 0, W, H);
  be.endPass();
}

/** A 2 × 2 image: red, green, blue, white. */
export const PIXELS = new Uint32Array(
  new Uint8Array([255, 0, 0, 255, 0, 255, 0, 255, 0, 0, 255, 255, 255, 255, 255, 255]).buffer,
);
