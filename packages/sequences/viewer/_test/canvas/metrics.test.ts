/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `metricsMeasurer`'s committed tables against this browser's fonts (SV1b; seqspec:viewer#the-render-backend,
// seqspec:viewer#text-measuring): the tables are the fonts of the machine that generated them
// (scripts/font-metrics.ts), so where this browser's differ the comparison says so and skips; where they match, the
// tables equal a fresh probe, and the measurer agrees with canvas `measureText`.

import { describe, expect, it } from "vitest";
import { FONT_FAMILIES, FONT_WEIGHTS, probeFonts } from "../../scripts/font-probe.ts";
import { DefaultFontMetrics, metricsMeasurer } from "../../src/svg/index.ts";

const live = probeFonts({ families: FONT_FAMILIES, weights: FONT_WEIGHTS });
const same = (p: (typeof live)[number]) => {
  const t = DefaultFontMetrics.find((m) => m.family === p.family && m.weight === p.weight);
  return t?.advances.every((a, i) => a === p.advances[i]) ?? false;
};
const m = metricsMeasurer(DefaultFontMetrics);
const g = new OffscreenCanvas(10, 10).getContext("2d") as OffscreenCanvasRenderingContext2D;

describe("the committed tables against this browser", () => {
  for (const p of live)
    it(`${p.family} ${p.weight}: the same advances, em box, cap and x-height as a fresh probe`, (ctx) => {
      if (!same(p)) return ctx.skip(`this browser's ${p.family} ${p.weight} is not the font the tables were made from`);
      const t = DefaultFontMetrics.find((x) => x.family === p.family && x.weight === p.weight);
      expect([t?.fallback, t?.ascent, t?.descent, t?.cap, t?.xHeight]).toEqual([
        p.fallback,
        p.ascent,
        p.descent,
        p.cap,
        p.xHeight,
      ]);
    });

  it("measure() is canvas measureText to within kerning and an optical-size font's looseness, at the skin's sizes", (ctx) => {
    if (!live.every(same)) return ctx.skip("this browser's fonts are not the ones the tables were made from");
    const strings = [
      "Hamburgefonstiv",
      "chr12:1,234,567-1,299,999",
      "Zinc ligands (UniProt)",
      "ABCDEFGHIKLMNPQRSTVWY",
      "Tetramerization 324-356",
    ];
    for (const stack of ["ui-sans-serif, system-ui, sans-serif", "ui-monospace, SFMono-Regular, Menlo, monospace"])
      for (const [weight, px] of [
        [400, 10],
        [400, 12],
        [600, 13],
      ]) {
        const font = `${weight} ${px}px ${stack}`;
        g.font = font;
        for (const s of strings) {
          const want = g.measureText(s).width;
          expect(Math.abs(m.measure(font, s) - want) / want, `${font}: ${s}`).toBeLessThan(0.03);
        }
      }
  });

  it("the em box is where this browser puts top, middle and bottom: the ink of an H sits on the baseline dy gives", (ctx) => {
    if (!live.every(same)) return ctx.skip("this browser's fonts are not the ones the tables were made from");
    const px = 40;
    for (const stack of ["ui-sans-serif, system-ui, sans-serif", "serif"]) {
      const font = `400 ${px}px ${stack}`;
      const box = m.metrics(font);
      const c = new OffscreenCanvas(100, 200);
      const x = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
      const foot = (baseline: CanvasTextBaseline) => {
        x.clearRect(0, 0, 100, 200);
        x.font = font;
        x.textBaseline = baseline;
        x.fillText("H", 10, 100);
        const d = x.getImageData(0, 0, 100, 200).data;
        for (let y = 199; y >= 0; y--) for (let i = 0; i < 100; i++) if (d[(y * 100 + i) * 4 + 3] > 127) return y + 1;
        return 0;
      };
      const base = foot("alphabetic");
      expect(foot("top") - base, stack).toBeCloseTo((box?.ascent ?? 0) * px, 0);
      expect(base - foot("bottom"), stack).toBeCloseTo((box?.descent ?? 0) * px, 0);
      expect(foot("middle") - base, stack).toBeCloseTo((((box?.ascent ?? 0) - (box?.descent ?? 0)) / 2) * px, 0);
    }
  });
});
