/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Canvas2D backend on the live path allocates nothing per frame (seqspec:viewer#the-render-backend): pans, each a
// back pass with a copy, 4,000 rects in 8 colors, 100 text runs, glyphs, a sector, paths, an image, a present; and a
// frame of 1,500 distinct colors (a continuous palette's row: each color's string made once, kept while it is drawn),
// drawn again and again. The pans over 1,000 frames and 4,000, so a few hundred bytes a frame show; the colors over 100
// and 400 (a string re-made per color a frame is about 30 KB a frame). The positive control: the pans allocating
// garbage a frame, sized to the run (`leaking`), measured alike, are caught.
// Chrome collects on its own while it draws, though the JS heap does not grow: about every 2,500,000 `drawImage` calls
// (2 s here) whatever their source (an OffscreenCanvas, an ImageBitmap, a canvas), and within the colors' 10,000 frames
// of 1,500 `fillStyle` changes. A run that meets a collection is refused, so with 1,000 text runs a frame from the
// glyph atlas the pans' runs of 4,000 frames, and the colors' at 1,000 frames, failed as "a GC ran" (2026-10-14, the
// closing pass): the pans draw 100 text runs (the text path still measured), the colors run shorter. The bindings box a
// fractional argument to `lineTo`, `arc`, `rect`, `drawImage` and `lineWidth` (12 bytes each; not to `moveTo`,
// `fillRect` or `setTransform`): the frame's sector and paths took up to 120 bytes a frame there, below what the runs
// catch (todo.md, Costs).

import { expect, it } from "vitest";
import { Canvas2DBackend, Canvas2DShared } from "../../src/canvas/index.js";
import { Pass, RectBatch, SkinToken } from "../../src/draw/api.js";
import { drawFrame, Frame, H, PIXELS, W } from "../support/frame.js";
import { controlCaught, heapScaling, isFlat, leaking } from "./heap-browser.js";

const N = 1_000;
/** The colors' run: shorter, as Chrome collects within longer ones (above). */
const COLORS_N = 100;

/** The backend over a W x H canvas, and a pan's step: one frame of the shared fixture. */
function panning() {
  const canvas = document.createElement("canvas");
  canvas.width = W;
  canvas.height = H;
  const shared = new Canvas2DShared({
    font: () => "12px sans-serif",
    fontPx: () => 12,
    token: (t) => (t === SkinToken.ink ? "#222" : "#06c"),
    glyph: () => null,
  });
  const be = new Canvas2DBackend([canvas.getContext("2d")], shared, 1 << Pass.main);
  const f = new Frame(4000, 100);
  const img = shared.fromPixels(2, 2, PIXELS);
  return { be, step: (i: number) => drawFrame(be, f, i, img) };
}

it("pans through the Canvas2D backend are flat", () => {
  const r = heapScaling(N, panning().step);
  expect(isFlat(r), JSON.stringify(r)).toBe(true);
}, 120_000);

it("a frame of 1,500 distinct colors drawn again and again is flat (each color's string made once)", () => {
  const { be } = panning();
  const b = new RectBatch();
  for (let i = 0; i < 1500; i++) {
    const t = ((i * 7919) % 1500) / 1500;
    b.begin((Math.round(68 + 185 * t) << 16) | (Math.round(1 + 230 * t) << 8) | Math.round(84 - 47 * t));
    b.push(i % W, 0, 1, 10);
  }
  const r = heapScaling(COLORS_N, () => {
    be.beginPass(Pass.main, W, H, 1);
    be.rects(b);
    be.endPass();
  });
  expect(isFlat(r), JSON.stringify(r)).toBe(true);
}, 120_000);

it("the same pans allocating garbage sized to the run a frame are caught (positive control)", () => {
  expect(controlCaught(() => heapScaling(N, leaking(panning().step, N)))).toBe(true);
}, 120_000);
