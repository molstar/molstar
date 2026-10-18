/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An animated zoom's text work, counted (seqspec:viewer#text-measuring, seqspec:viewer#the-render-backend): the
// navigator's double click on case 14 (a protein's residues over its coding DNA) and case 15 (four designs' residues
// placed on the wild type through their mappings), each letter size a zoom passes through. A zoom that reached a new
// letter size made some 80 glyph atlas cells in one frame, each drawn at once from a page rasterized on the CPU, so
// each draw after a make uploaded the page to the GPU target: a 50 to 110 ms frame (David's jitter, 2026-10-13). The
// atlas uploads each page at most once a pass now (new cells show from the next), and after one zoom in and out
// nothing is measured again: no `measureText`, no new font. Positive controls: cells drawn from at once (the old
// behavior) upload a page per new cell; the text caches dropped every frame measure again.

import { describe, expect, it } from 'vitest';
import { case14, case15 } from '../../../core/_test/format/cases.js';
import type { Canvas2DBackend } from '../../src/canvas/index.js';
import type { Viewport } from '../../src/viewport.js';
import { close, type Opened, opened } from './open.js';

/** Pages the atlas holds (backend.ts `PAGES`): the most uploads a pass can need. */
const PAGES = 4;
/**
 * The sections at DPR 2, the letters' cap above the header labels' 12 px (`maxLetterPx`), so a zoom passes through
 * more letter sizes, each new to the atlas, as it did before the cap.
 */
const OPEN = (height: number) => ({ width: 1344, height, dpr: 2, render: { maxLetterPx: 16 } });
/** A 250 ms animation at 60 Hz, and a few frames over. */
const FRAMES = 20;

function atlasOf(t: Opened) {
  const a = (t.surface as unknown as { backend: Canvas2DBackend }).backend.images.atlas;
  if (!a) throw new Error('no atlas');
  return a;
}
const fontsOf = (t: Opened) => (t.ctx.text as unknown as { fonts: unknown[] }).fonts.length;

/** `measureText` calls over both context classes (the measurer's is an OffscreenCanvas's). */
function countMeasure() {
  const protos = [CanvasRenderingContext2D.prototype, OffscreenCanvasRenderingContext2D.prototype] as const;
  const originals = protos.map((p) => p.measureText);
  const counter = {
    calls: 0,
    restore() {
      for (let i = 0; i < protos.length; i++) protos[i].measureText = originals[i];
    },
  };
  protos.forEach((p, i) => {
    p.measureText = function (this: unknown, text: string) {
      counter.calls++;
      return (originals[i] as (t: string) => TextMetrics).call(this, text);
    } as typeof p.measureText;
  });
  return counter;
}

/**
 * The navigator's double click (gestures.ts `toggleZoom`): about 60 positions across the data area centered at `x`,
 * then the full view, each animated; `each` runs after every frame. The most uploads a frame made, and the frames
 * whose cell differed from the last (the animation's steps).
 */
function zoomInOut(t: Opened, x: number, each?: () => void) {
  const vp = t.ctx.viewports.get(t.section) as unknown as Viewport;
  const atlas = atlasOf(t);
  let most = 0;
  let steps = 0;
  let cell = 0;
  const run = () => {
    for (let f = 0; f < FRAMES; f++) {
      const before = atlas.uploads;
      t.advance(16);
      t.surface.renderNow();
      each?.();
      most = Math.max(most, atlas.uploads - before);
      if (vp.live.current.cell !== cell) steps++;
      cell = vp.live.current.cell;
    }
  };
  const s = vp.live.current;
  vp.zoomCentered(Math.max(vp.props.minLetterCell, s.width / 60), x, { animate: true, origin: 'pointer' });
  run();
  vp.fit({ animate: true, origin: 'pointer' });
  run();
  return { most, steps };
}

const CASES = [
  ['case 14', case14, 300],
  ['case 15', case15, 400],
] as const;

describe("an animated zoom's text work", () => {
  for (const [name, doc, height] of CASES) {
    it(`${name}: each atlas page uploaded at most once a pass, from cold`, async () => {
      const t = await opened(doc(), OPEN(height));
      const r = zoomInOut(t, 600);
      const a = atlasOf(t);
      console.log(`${name} cold: ${JSON.stringify({ ...r, made: a.made, uploads: a.uploads })}`);
      expect(a.made, 'cells made by the zoom').toBeGreaterThan(20);
      expect(r.steps, "the animation's steps").toBeGreaterThan(20);
      expect(r.most).toBeLessThanOrEqual(PAGES);
      close(t);
    });
    it(`${name}: positive control, cells drawn from at once upload a page per new cell`, async () => {
      const t = await opened(doc(), OPEN(height));
      atlasOf(t).prepass = false;
      const r = zoomInOut(t, 600);
      console.log(`${name} undeferred: ${JSON.stringify(r)}`);
      expect(r.most).toBeGreaterThan(PAGES);
      close(t);
    });
    it(`${name}: after one zoom in and out, nothing is measured again and no font is made`, async () => {
      const t = await opened(doc(), OPEN(height));
      zoomInOut(t, 600);
      const fonts = fontsOf(t);
      const made = atlasOf(t).made;
      const c = countMeasure();
      try {
        zoomInOut(t, 600);
        zoomInOut(t, 600);
      } finally {
        c.restore();
      }
      expect(c.calls, 'measureText calls').toBe(0);
      expect(fontsOf(t), 'fonts').toBe(fonts);
      expect(atlasOf(t).made - made, 'atlas cells made').toBeLessThanOrEqual(8);
      close(t);
    });
    it(`${name}: positive control, the text caches dropped each frame measure again`, async () => {
      const t = await opened(doc(), OPEN(height));
      zoomInOut(t, 600);
      const c = countMeasure();
      try {
        zoomInOut(t, 600, () => t.ctx.text.invalidate());
      } finally {
        c.restore();
      }
      expect(c.calls).toBeGreaterThan(FRAMES * 50);
      close(t);
    });
  }
});
