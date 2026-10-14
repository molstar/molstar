/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The glyph atlas of the Canvas2D backend (seqspec:viewer#the-render-backend; CS2 letters) in Chromium: letters drawn
// from its cells equal `fillText`'s, within the tolerance pinned here, at DPR 1 and 2, in the light and dark skins, for
// every character the letters renderers draw (printable ASCII: the protein and nucleic alphabets, gaps, lower case,
// atom symbols' first letters) in their fonts (mono 6 to 14 px, regular and bold; ui 7 to 12 px bold), at whole, half
// and odd anchors and in each alignment and baseline. What the atlas leaves to `fillText` (more than one character, a
// scale, an angle) draws as before; a font load and a skin's new token colors drop its cells; a full atlas evicts.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { Canvas2DBackend, Canvas2DShared } from "../../src/canvas/index.ts";
import {
  type DrawResources,
  NO_COLOR,
  Pass,
  RectBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "../../src/draw/api.ts";
import { Dirty } from "../../src/draw/index.ts";
import { DefaultSkin } from "../../src/index.ts";
import { canvasesOf, close, diff as diffPixels, type Opened, opened, pixelsOf } from "./open.ts";

const MONO = DefaultSkin.fontMono;
const UI = DefaultSkin.fontUi;
/** The tokens the letters draw in: the light skin's, and skin.css's dark block. */
const SKINS = {
  light: { surface: DefaultSkin.surface, ink: DefaultSkin.ink, inkQuiet: DefaultSkin.inkQuiet },
  dark: { surface: "rgb(24 26 30)", ink: "rgb(230 232 235)", inkQuiet: "rgb(156 163 175)" },
};
type SkinName = keyof typeof SKINS;
/** The fonts, by id: [CSS font, px]. */
const FONTS: [string, number][] = [];
for (let px = 6; px <= 14; px++) for (const w of [400, 600]) FONTS.push([`${w} ${px}px ${MONO}`, px]);
for (let px = 7; px <= 12; px++) FONTS.push([`600 ${px}px ${UI}`, px]);
const CHARS = Array.from({ length: 94 }, (_, i) => String.fromCharCode(33 + i));
/** A letter's inks: the skin's ink and quiet ink, the inks `cellInk` picks on a cell, a theme's color; translucent. */
const OPAQUE: [number, number][] = [
  [SkinToken.ink, 255],
  [SkinToken.inkQuiet, 255],
  [0xffffff, 255],
  [0x1b1f24, 255],
  [0x2f7fd0, 255],
];
const TRANSLUCENT: [number, number][] = [
  [SkinToken.ink, 128],
  [0x2f7fd0, 90],
];
/** Anchors' fractions (device px): letters sit on whole and half pixels; 0.37 and 0.62 go to fillText. */
const FX = [0, 0.5, 0.25, 0.75, 0.37, 0.5, 0];
const FY = [0.5, 0, 0.25, 0.75, 0.5, 0.62, 0.5];
const STYLES: [number, number][] = [
  [TextAlign.center, TextBaseline.alphabetic],
  [TextAlign.center, TextBaseline.middle],
  [TextAlign.left, TextBaseline.alphabetic],
  [TextAlign.right, TextBaseline.top],
  [TextAlign.center, TextBaseline.bottom],
];

function resources(skin: SkinName, version = { v: 0 }): DrawResources & { version: { v: number } } {
  const t = SKINS[skin];
  return {
    version,
    font: (id) => FONTS[id][0],
    fontPx: (id) => FONTS[id][1],
    token: (k) => (k === SkinToken.ink ? t.ink : k === SkinToken.inkQuiet ? t.inkQuiet : t.surface),
    glyph: () => null,
  };
}

/** A backend drawing into a visible canvas through a back target, as the surface's main pass. */
function target(w: number, h: number, shared: Canvas2DShared) {
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  const be = new Canvas2DBackend([g], shared, 1 << Pass.main);
  return { be, px: () => g.getImageData(0, 0, w, h).data };
}

const COLS = 24;
/** Every character of font f in each ink and style on the skin's surface, the cells `pitch` device px apart. */
function frame(be: Canvas2DBackend, dpr: number, f: number, inks: [number, number][], w: number, pitch: number) {
  const h = Math.ceil((STYLES.length * inks.length * CHARS.length) / COLS) * pitch;
  const bg = new RectBatch(1);
  const t = new TextBatch(1024);
  be.beginPass(Pass.main, w, h, dpr);
  bg.begin(SkinToken.surface);
  bg.push(0, 0, w, h);
  be.rects(bg);
  let k = 0;
  for (const [align, base] of STYLES) {
    t.reset();
    t.style(f, align, base);
    for (const [ink, alpha] of inks)
      for (const ch of CHARS) {
        t.begin(ink, alpha);
        const i = t.add(ch);
        t.x[i] = (k % COLS) * pitch + pitch / 2 + FX[k % FX.length];
        t.y[i] = Math.floor(k / COLS) * pitch + pitch / 2 + FY[Math.floor(k / FX.length) % FY.length];
        k++;
      }
    be.text(t);
  }
  be.present(0, 0, w, h);
  be.endPass();
  return k;
}

interface Diff {
  pixels: number;
  differ: number;
  max: number;
}
/** Pixels whose channels differ, and the largest difference. */
function diff(a: Uint8ClampedArray, b: Uint8ClampedArray, into: Diff) {
  for (let i = 0; i < a.length; i += 4) {
    let m = 0;
    for (let c = 0; c < 4; c++) m = Math.max(m, Math.abs(a[i + c] - b[i + c]));
    if (m > 0) into.differ++;
    into.max = Math.max(into.max, m);
  }
  into.pixels += a.length / 4;
}

/**
 * The tolerance, measured (2026-10-11, the browser suite's Chrome, every case below): opaque letters equal fillText's,
 * pixel for pixel; translucent ones differ by at most 1 in a channel (a cell's 8-bit premultiplied ink blended at the
 * item's alpha rounds apart from the glyph's coverage blended at once). Anchors at other fractions (here 0.37, 0.62)
 * were up to 76 apart through the atlas (the baseline's offset rounds them to another glyph position): fillText's.
 */
const TOLERANCE = 1;

describe("letters from the atlas equal fillText's", () => {
  for (const dpr of [1, 2])
    for (const skin of ["light", "dark"] as const)
      it(`at DPR ${dpr}, ${skin}`, () => {
        const opaque: Diff = { pixels: 0, differ: 0, max: 0 };
        const translucent: Diff = { pixels: 0, differ: 0, max: 0 };
        const atlas = new Canvas2DShared(resources(skin));
        const plain = new Canvas2DShared(resources(skin), { atlas: false });
        for (let f = 0; f < FONTS.length; f++) {
          const pitch = Math.ceil(2 * FONTS[f][1] * dpr) + 4;
          for (const [inks, into] of [
            [OPAQUE, opaque],
            [TRANSLUCENT, translucent],
          ] as const) {
            const w = COLS * pitch;
            const h = Math.ceil((STYLES.length * inks.length * CHARS.length) / COLS) * pitch;
            const a = target(w, h, atlas);
            const b = target(w, h, plain);
            frame(a.be, dpr, f, inks, w, pitch);
            frame(b.be, dpr, f, inks, w, pitch);
            diff(a.px(), b.px(), into);
          }
        }
        console.log(`atlas DPR ${dpr} ${skin}: ${JSON.stringify({ opaque, translucent })} cells ${atlas.atlas?.made}`);
        expect(atlas.atlas?.made).toBeGreaterThan(5000);
        expect(opaque.max, "opaque letters").toBeLessThanOrEqual(TOLERANCE);
        expect(translucent.max, "translucent letters").toBeLessThanOrEqual(TOLERANCE);
      }, 120_000);
});

describe("the atlas's bounds and invalidation", () => {
  it("leaves text of more than one character, scaled or rotated items to fillText", () => {
    const shared = new Canvas2DShared(resources("light"));
    const { be } = target(200, 100, shared);
    const t = new TextBatch(4);
    be.beginPass(Pass.main, 200, 100, 1);
    t.style(8, TextAlign.center, TextBaseline.middle);
    t.begin(SkinToken.ink);
    t.push(20, 20, "AB");
    const i = t.add("C");
    t.x[i] = 60;
    t.y[i] = 20;
    t.sx[i] = 2;
    be.text(t);
    t.reset();
    t.style(8, TextAlign.center, TextBaseline.middle, Math.PI / 2);
    t.push(100, 50, "D");
    be.text(t);
    be.endPass();
    expect(shared.atlas?.made).toBe(0);
    t.reset();
    t.style(8, TextAlign.center, TextBaseline.middle);
    t.push(100, 50, "D");
    t.push(120, 50, "D");
    be.beginPass(Pass.main, 200, 100, 1);
    be.text(t);
    be.endPass();
    expect(shared.atlas?.made).toBe(1);
    expect(t.color[0]).not.toBe(NO_COLOR);
  });
  it("drops its cells on a font load and on a token's new color; a DPR is part of the key", () => {
    const res = resources("light");
    const shared = new Canvas2DShared(res, { fonts: () => res.version.v });
    const { be } = target(200, 100, shared);
    const t = new TextBatch(4);
    const draw = (dpr: number) => {
      be.beginPass(Pass.main, 200, 100, dpr);
      be.text(t);
      be.endPass();
    };
    t.style(8, TextAlign.center, TextBaseline.middle);
    t.begin(SkinToken.ink);
    t.push(20, 20.5, "A");
    draw(1);
    draw(1);
    expect(shared.atlas?.made).toBe(1);
    draw(2);
    expect(shared.atlas?.made).toBe(2);
    res.version.v++;
    draw(1);
    expect(shared.atlas?.made).toBe(3);
    draw(1);
    expect(shared.atlas?.made).toBe(3);
    // A skin whose ink is another color: the token's cells go.
    const dark = resources("dark");
    (res as { token: DrawResources["token"] }).token = dark.token;
    draw(1);
    expect(shared.atlas?.made).toBe(4);
  });
  it("a full atlas evicts the page drawn least recently, never one the batch drew; its letters stay right", () => {
    const shared = new Canvas2DShared(resources("light"));
    const plain = new Canvas2DShared(resources("light"), { atlas: false });
    const [w, h] = [1280, 500];
    const a = target(w, h, shared);
    const b = target(w, h, plain);
    const bg = new RectBatch(1);
    const t = new TextBatch(256);
    // 14 px bold mono at DPR 3: 96 px cells, 100 to a page, 400 in all; 250 inks a batch overflow them each round.
    for (let round = 0; round < 24; round++) {
      t.reset();
      t.style(17, TextAlign.center, TextBaseline.middle);
      for (let k = 0; k < 250; k++) {
        t.begin(((round * 250 + k) * 2654435761) & 0xffffff);
        const i = t.add(CHARS[k % CHARS.length]);
        t.x[i] = 25 + (k % 25) * 50;
        t.y[i] = 25.5 + Math.floor(k / 25) * 48;
      }
      for (const s of [a, b]) {
        s.be.beginPass(Pass.main, w, h, 3);
        bg.reset();
        bg.begin(SkinToken.surface);
        bg.push(0, 0, w, h);
        s.be.rects(bg);
        s.be.text(t);
        s.be.present(0, 0, w, h);
        s.be.endPass();
      }
    }
    expect(shared.atlas?.evicted).toBeGreaterThan(0);
    const d: Diff = { pixels: 0, differ: 0, max: 0 };
    diff(a.px(), b.px(), d);
    expect(d.max).toBeLessThanOrEqual(TOLERANCE);
  });
});

describe("the surface's `atlas` option (SurfaceOptions.atlas)", () => {
  /** A DNA sequence of 24 letters at 20 px a cell, with a ruler. */
  const doc = () => {
    const b = createSeqViewSpec({ title: "gene" });
    const s = b.data.sequence({ ref: "s", alphabet: "dna", value: "GAATTCGCCACCATGGCTTCCAAG" });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: s });
    v.track({ label: "" }).ruler();
    v.track({ label: "DNA" }).letters();
    return b.getDocument();
  };
  /** `fillText` calls over the two context classes (the page's canvases and the atlas's and back targets' offscreen ones). */
  function countFillText() {
    const protos = [CanvasRenderingContext2D.prototype, OffscreenCanvasRenderingContext2D.prototype] as const;
    const originals = protos.map((p) => p.fillText);
    const counter = {
      calls: 0,
      restore() {
        for (let i = 0; i < protos.length; i++) protos[i].fillText = originals[i];
      },
    };
    protos.forEach((p, i) => {
      p.fillText = function (this: unknown, ...a: [string, number, number, number?]) {
        counter.calls++;
        return (originals[i] as (...x: unknown[]) => void).apply(this, a);
      } as typeof p.fillText;
    });
    return counter;
  }
  /** Whole frames drawn after the first, `fillText` calls counted over them. */
  function framesFillText(t: Opened, frames: number) {
    const c = countFillText();
    try {
      for (let k = 0; k < frames; k++) {
        t.surface.composer.redraw();
        t.surface.invalidate(Dirty.all);
        t.surface.renderNow();
      }
    } finally {
      c.restore();
    }
    return c.calls / frames;
  }

  it("true (the default) draws a frame's letters from cells, false through fillText; the pixels are alike", async () => {
    const on = await opened(doc(), { width: 480, height: 120 });
    const off = await opened(doc(), { width: 480, height: 120, atlas: false });
    const withAtlas = framesFillText(on, 3);
    const without = framesFillText(off, 3);
    // 24 letters a frame, plus the ruler's labels and the row label, which both draw by fillText.
    expect(without - withAtlas).toBeGreaterThanOrEqual(24);
    expect(withAtlas).toBeLessThan(without);
    const [a] = canvasesOf(on);
    const [b] = canvasesOf(off);
    expect(diffPixels(pixelsOf(a), pixelsOf(b), a.width, 1).count).toBe(0);
    close(on);
    close(off);
  });
});
