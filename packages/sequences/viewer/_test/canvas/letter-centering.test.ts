/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Letters centered on their cells, by pixels (seqspec:viewer#the-render-backend, text placement): a letter on the
// alphabetic baseline `letterLines` gives, from the font's ink heights the canvas measurer reads at the device size
// (`TextManager.ink`), inks a capital centered on its cell, and a lower-case letter's x-height, to half a device px: in
// Canvas2D (from the glyph atlas and by `fillText`) and in the SVG rasterized, at the letters' sizes, cells of the
// letters', the rows' and an odd height, DPR 1 and 2. A `middle` baseline drew capitals 1 to 2.4 device px high
// (ui-monospace, 10 to 14 px). Which renderers place their letters so is the recordings' (_test/letter-lines.test.ts).

import { describe, expect, it } from "vitest";
import { Canvas2DBackend, Canvas2DShared, canvasMeasurer } from "../../src/canvas/index.ts";
import {
  type DrawResources,
  LETTER_LINE,
  letterLines,
  Pass,
  type RenderBackend,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "../../src/draw/api.ts";
import { DefaultSkin } from "../../src/index.ts";
import { SvgBackend } from "../../src/svg/index.ts";
import { TextManager } from "../../src/text.ts";

/**
 * Capitals, and lower-case letters of the x-height alone (no ascender or descender: their ink is the x-height's). Not
 * "W": its pointed feet and diagonals at 8 px, DPR 1, measure 0.6 device px off in an odd cell.
 */
const CHARS = ["M", "S", "D", "E", "H", "x", "a", "e", "s"];
/** A column per letter (CSS px). */
const COLUMN = 32;
/** Cells: [top, height] in CSS px, as the letters (22, inset 2), the rows (22, inset 1) and a short track (15, inset 2). */
const CELLS: [number, number][] = [
  [4, 18],
  [3, 20],
  [6.5, 11],
];

function textManager() {
  const tm = new TextManager(canvasMeasurer(), () => 64);
  tm.setFamilies(DefaultSkin.fontMono, DefaultSkin.fontUi);
  return tm;
}
const resources = (tm: TextManager): DrawResources => ({
  font: (id) => tm.fontString(id),
  fontPx: (id) => Number(/([\d.]+)px/.exec(tm.fontString(id))?.[1] ?? 0),
  token: () => "#000",
  glyph: () => null,
});

/** Each letter on the line `letterLines` gives for the band [y0, y1) device px, a column each. */
function letters(be: RenderBackend, tm: TextManager, font: number, dpr: number, y0: number, y1: number) {
  const lines = new Float64Array(2);
  letterLines({ text: tm, dpr }, font, y0, y1, lines);
  const t = new TextBatch();
  t.style(font, TextAlign.center, TextBaseline.alphabetic);
  t.begin(SkinToken.ink);
  for (const [i, ch] of CHARS.entries()) {
    const k = t.add(ch);
    t.x[k] = Math.round((i + 0.5) * COLUMN * dpr);
    t.y[k] = lines[LETTER_LINE[ch.charCodeAt(0)]];
  }
  be.text(t);
}

/** The pixels of the letters drawn by Canvas2D (`atlas` or `fillText`) and of the SVG rasterized: w × h × 4. */
async function drawn(tm: TextManager, font: number, dpr: number, y0: number, y1: number) {
  const w = CHARS.length * COLUMN * dpr;
  const h = 40 * dpr;
  const canvas = (atlas: boolean) => {
    const c = new OffscreenCanvas(w, h);
    const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
    const be = new Canvas2DBackend([g], new Canvas2DShared(resources(tm), { atlas }), 0);
    be.beginPass(Pass.main, w, h, dpr);
    letters(be, tm, font, dpr, y0, y1);
    be.endPass();
    return g.getImageData(0, 0, w, h).data;
  };
  const svg = new SvgBackend(resources(tm));
  svg.beginPass(Pass.main, w, h, dpr);
  letters(svg, tm, font, dpr, y0, y1);
  svg.endPass();
  const url = URL.createObjectURL(new Blob([await svg.toString()], { type: "image/svg+xml" }));
  const img = new Image();
  img.src = url;
  await img.decode();
  const s = new OffscreenCanvas(w, h);
  const sg = s.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  sg.drawImage(img, 0, 0, w, h);
  URL.revokeObjectURL(url);
  return { w, h, atlas: canvas(true), fillText: canvas(false), svg: sg.getImageData(0, 0, w, h).data };
}

/**
 * The ink's top and bottom in [x0, x1) (device px): the first and last rows with ink, each inking the share of a px its
 * coverage is of its column's most (a stem a px wide covers less than a px across, whatever its height).
 */
function inkY(d: Uint8ClampedArray, w: number, h: number, x0: number, x1: number): [number, number] {
  const most = new Float64Array(x1 - x0);
  for (let y = 0; y < h; y++)
    for (let x = x0; x < x1; x++) most[x - x0] = Math.max(most[x - x0], d[(y * w + x) * 4 + 3]);
  const share = (y: number) => {
    let s = 0;
    for (let x = x0; x < x1; x++) if (most[x - x0] > 0) s = Math.max(s, d[(y * w + x) * 4 + 3] / most[x - x0]);
    return s;
  };
  let t = 0;
  while (t < h && share(t) === 0) t++;
  let b = h - 1;
  while (b > t && share(b) === 0) b--;
  return [t + 1 - share(t), b + share(b)];
}

describe("letters centered on their cells", () => {
  const fonts: [string, number, 400 | 600][] = [
    ["mono", 8, 400],
    ["mono", 10, 400],
    ["mono", 12, 400],
    ["mono", 14, 400],
    ["mono", 14, 600],
    ["ui", 10, 600],
  ];
  for (const dpr of [1, 2])
    for (const [role, px, weight] of fonts)
      it(`${role} ${px} px ${weight}, DPR ${dpr}: a capital's ink and a lower-case letter's x-height on the cell's center`, async () => {
        const tm = textManager();
        const font = tm.font(role as "mono" | "ui", px, weight);
        const off: string[] = [];
        for (const [top, height] of CELLS) {
          const y0 = Math.round(top * dpr);
          const y1 = Math.round((top + height) * dpr);
          const center = (y0 + y1) / 2;
          const r = await drawn(tm, font, dpr, y0, y1);
          for (const [i, ch] of CHARS.entries()) {
            const x0 = i * COLUMN * dpr;
            const x1 = (i + 1) * COLUMN * dpr;
            for (const [by, d] of [
              ["atlas", r.atlas],
              ["fillText", r.fillText],
              ["svg", r.svg],
            ] as const) {
              const [a, b] = inkY(d, r.w, r.h, x0, x1);
              const mid = (a + b) / 2;
              if (Math.abs(mid - center) > 0.5)
                off.push(`${ch} ${by}, cell ${top} + ${height}: ${(mid - center).toFixed(2)} device px`);
            }
          }
          // The atlas draws what fillText draws.
          expect(r.atlas).toEqual(r.fillText);
        }
        expect(off).toEqual([]);
      });
});
