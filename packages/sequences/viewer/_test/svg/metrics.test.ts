/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `metricsMeasurer` and its committed tables (SV1b; seqspec:viewer#the-render-backend, seqspec:viewer#text-measuring)
// in Node: widths are the advances' sums, the font string picks the table the skin's stacks resolve to, the em box
// comes back in em, and `dy` follows from it. Whether the tables match a browser's fonts is checked in Chromium
// (_test/canvas/metrics.test.ts).

import { describe, expect, it } from "vitest";
import { GlyphBatch, Pass, SHAPE_COUNT, TextAlign, TextBaseline, TextBatch } from "../../src/draw/api.ts";
import { DefaultSkin } from "../../src/index.ts";
import { DefaultFontMetrics, type FontMetrics, hasMetrics, metricsMeasurer, SvgBackend } from "../../src/svg/index.ts";
import { monospaceMeasurer } from "../../src/text.ts";
import { all, parseXml } from "./xml.ts";

const m = metricsMeasurer(DefaultFontMetrics);
const table = (family: string, weight = 400) =>
  DefaultFontMetrics.find((t) => t.family === family && t.weight === weight) as FontMetrics;

describe("the committed tables", () => {
  it("each family and weight: 95 advances, a fallback, an em box of one em", () => {
    expect(DefaultFontMetrics.map((t) => `${t.family} ${t.weight}`).sort()).toEqual(
      ["monospace", "sans-serif", "serif", "ui-monospace", "ui-sans-serif"]
        .flatMap((f) => [`${f} 400`, `${f} 600`])
        .sort(),
    );
    for (const t of DefaultFontMetrics) {
      expect(t.unitsPerEm).toBe(1000);
      expect(t.advances.length).toBe(95);
      expect(Math.min(...t.advances), `${t.family} ${t.weight}`).toBeGreaterThan(0);
      expect(t.fallback).toBeGreaterThan(300);
      // Canvas2D's top and bottom baselines are the em box's edges: they are one em apart (to a pixel at 1,000 px).
      expect(Math.abs(t.ascent + t.descent - 1000)).toBeLessThanOrEqual(2);
      expect(t.ascent).toBeGreaterThan(t.descent);
    }
  });
  it("monospace families are monospace; a bold face is no narrower", () => {
    for (const f of ["monospace", "ui-monospace"]) expect(new Set(table(f).advances).size).toBe(1);
    for (const f of ["ui-sans-serif", "sans-serif", "serif"]) {
      const [a, b] = [table(f, 400), table(f, 600)];
      expect(b.advances.reduce((s, v) => s + v, 0)).toBeGreaterThanOrEqual(a.advances.reduce((s, v) => s + v, 0));
    }
  });
});

describe("widths", () => {
  it("a character is its advance, a string their sum, in the font's size", () => {
    const t = table("sans-serif");
    expect(m.measure("400 1000px sans-serif", "H")).toBeCloseTo(t.advances["H".charCodeAt(0) - 32], 6);
    expect(m.measure("400 12px sans-serif", "")).toBe(0);
    let sum = 0;
    for (const c of "Hello, world") sum += t.advances[c.charCodeAt(0) - 32];
    expect(m.measure("400 12px sans-serif", "Hello, world")).toBeCloseTo((sum / 1000) * 12, 9);
    expect(m.measure("400 24px sans-serif", "Hello, world")).toBeCloseTo(
      2 * m.measure("400 12px sans-serif", "Hello, world"),
      9,
    );
  });
  it("any other code point is the fallback (a surrogate pair once); controls too", () => {
    const f = table("sans-serif").fallback;
    expect(m.measure("400 1000px sans-serif", "é")).toBe(f);
    expect(m.measure("400 1000px sans-serif", "日本")).toBe(2 * f);
    expect(m.measure("400 1000px sans-serif", "\u{1F600}")).toBe(f);
    expect(m.measure("400 1000px sans-serif", "\t")).toBe(f);
  });
  it("a string that is not a font measures 0, and has no box", () => {
    expect(m.measure("bold sans-serif", "abc")).toBe(0);
    expect(m.metrics("bold sans-serif")).toBeNull();
  });
});

describe("the font string picks the table", () => {
  const widthOf = (font: string) => m.measure(font, "Hamburgefonstiv");
  const want = (family: string, weight = 400) =>
    (metricsMeasurer([table(family, weight)]).measure("400 1000px x", "Hamburgefonstiv") / 1000) * 12;
  it("the skin's own stacks: its first family", () => {
    expect(widthOf(`400 12px ${DefaultSkin.fontUi}`)).toBeCloseTo(want("ui-sans-serif"), 9);
    expect(widthOf(`400 12px ${DefaultSkin.fontMono}`)).toBeCloseTo(want("ui-monospace"), 9);
  });
  it("system-ui and the system monospace fonts name those tables; generic names the generic ones", () => {
    expect(widthOf("400 12px system-ui, sans-serif")).toBeCloseTo(want("ui-sans-serif"), 9);
    expect(widthOf("400 12px Menlo, monospace")).toBeCloseTo(want("ui-monospace"), 9);
    expect(widthOf("400 12px Consolas, 'Courier New', monospace")).toBeCloseTo(want("monospace"), 9);
    expect(widthOf("400 12px Georgia, serif")).toBeCloseTo(want("serif"), 9);
    expect(widthOf('400 12px "Helvetica Neue", Arial')).toBeCloseTo(want("sans-serif"), 9);
  });
  it("an unknown family goes on to the next, and a stack of unknown ones is the first table's", () => {
    expect(widthOf('400 12px "Inter", monospace')).toBeCloseTo(want("monospace"), 9);
    expect(widthOf('400 12px "Inter", "Nobody"')).toBeCloseTo(want("ui-sans-serif"), 9);
  });
  it("the nearest weight; bold and a leading style word read", () => {
    expect(widthOf("400 12px serif")).toBeCloseTo(want("serif", 400), 9);
    expect(widthOf("700 12px serif")).toBeCloseTo(want("serif", 600), 9);
    expect(widthOf("bold 12px serif")).toBeCloseTo(want("serif", 600), 9);
    expect(widthOf("italic 300 12px serif")).toBeCloseTo(want("serif", 400), 9);
    expect(widthOf("normal 12px serif")).toBeCloseTo(want("serif", 400), 9);
  });
  it("host tables: only theirs are used", () => {
    const own: FontMetrics = {
      family: "Mine",
      weight: 400,
      unitsPerEm: 100,
      advances: new Uint16Array(95).fill(50),
      fallback: 80,
      ascent: 90,
      descent: 10,
    };
    const mine = metricsMeasurer([own]);
    expect(mine.measure("400 10px Mine", "abc")).toBe(15);
    expect(mine.measure("400 10px sans-serif", "ab é")).toBe(5 + 5 + 5 + 8);
    expect(mine.metrics("400 10px Mine")).toEqual({ ascent: 0.9, descent: 0.1 });
    expect(() => metricsMeasurer([])).toThrow(/no tables/);
  });
});

describe("the em box", () => {
  it("the table's ascent, descent, cap and x-height in em, and a measurer without one is not a MetricsMeasurer", () => {
    expect(m.metrics("400 12px ui-sans-serif")).toEqual({
      ascent: table("ui-sans-serif").ascent / 1000,
      descent: table("ui-sans-serif").descent / 1000,
      cap: (table("ui-sans-serif").cap ?? 0) / 1000,
      xHeight: (table("ui-sans-serif").xHeight ?? 0) / 1000,
    });
    expect(hasMetrics(m)).toBe(true);
    expect(hasMetrics(monospaceMeasurer)).toBe(false);
  });
});

describe("baselines as dy", () => {
  const FONT = '400 20px "x", sans-serif';
  const resources = { font: () => FONT, fontPx: () => 20, token: () => "red", glyph: () => null };
  async function written(baseline: number, metrics?: (f: string) => { ascent: number; descent: number } | null) {
    const be = new SvgBackend(resources, { metrics });
    be.beginPass(Pass.main, 200, 100, 1);
    const t = new TextBatch();
    t.style(0, TextAlign.left, baseline, 0);
    t.begin(0);
    t.x[t.add("Hi")] = 10;
    be.text(t);
    be.endPass();
    return await be.toString();
  }
  it("middle, top and bottom are explicit offsets from the alphabetic baseline, in CSS px, and no dominant-baseline", async () => {
    const box = { ascent: 0.8, descent: 0.2 };
    const dy = async (b: number) => all(parseXml(await written(b, () => box)), "text")[0].attrs.dy;
    expect(await dy(TextBaseline.alphabetic)).toBeUndefined();
    expect(await dy(TextBaseline.top)).toBe("16"); // 0.8 em of 20 px below the point
    expect(await dy(TextBaseline.bottom)).toBe("-4");
    expect(await dy(TextBaseline.middle)).toBe("6"); // (0.8 - 0.2) / 2 em
    for (const b of [TextBaseline.middle, TextBaseline.top, TextBaseline.bottom])
      expect(await written(b, () => box)).not.toContain("dominant-baseline");
  });
  it("without a box the baseline is dominant-baseline, as in SV1a", async () => {
    expect(await written(TextBaseline.middle)).toContain('dominant-baseline="central"');
    expect(await written(TextBaseline.top, () => null)).toContain('dominant-baseline="text-before-edge"');
    expect(await written(TextBaseline.bottom, () => null)).toContain('dominant-baseline="text-after-edge"');
  });
  it("a character glyph is centered by dy on its 2r-high font, a rotated item's dy is in its own axes", async () => {
    const box = { ascent: 0.8, descent: 0.2 };
    const res = { ...resources, glyph: () => ({ char: "W", path: null, box: [1, 1] as const, fallback: 0 }) };
    const be = new SvgBackend(res, { metrics: () => box });
    be.beginPass(Pass.main, 200, 100, 2);
    const g = new GlyphBatch();
    g.font = 0;
    g.begin(0);
    g.push(40, 40, 10, SHAPE_COUNT); // 2r = 20 device px = 10 CSS px at DPR 2
    be.glyphs(g);
    const t = new TextBatch();
    t.style(0, TextAlign.center, TextBaseline.top, Math.PI / 2);
    t.begin(0);
    const k = t.add("Hi");
    t.x[k] = 60;
    t.y[k] = 60;
    t.sx[k] = 2;
    be.text(t);
    be.endPass();
    const svg = parseXml(await be.toString());
    const [glyph, rotated] = all(svg, "text");
    expect(glyph.attrs.dy).toBe("3"); // (0.8 - 0.2) / 2 em of 10 CSS px
    expect(glyph.attrs["font-size"]).toBe("10");
    expect(all(svg, "g").some((e) => e.attrs["dominant-baseline"])).toBe(false);
    // The offset is along the text's own (rotated, scaled) axes: 0.8 em of 20 px font, set in the item's local frame.
    expect(rotated.attrs.transform).toMatch(/rotate\(90\) scale\(2 1\)$/);
    expect(rotated.attrs.dy).toBe("16");
  });
});
