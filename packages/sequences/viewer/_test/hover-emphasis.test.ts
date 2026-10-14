/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover emphasis (HV1; seqspec:interaction#hover) on the recording backend: the hover ink's contrast over a grid of
// colors in both skins; per renderer, the outline drawn for the hovered element alone, in the hover ink of its drawn
// color, in both skins, every stroke inside the item's own drawn box (David, 2026-10-12: never over its edge), and a
// glyph's ring only inside its cell; the keyboard cursor's item; the opt-out (`hoverEmphasis`, the style key
// `hover_emphasis`); and no emphasis in an export or in the main pass.

import { describe, expect, it } from "vitest";
import { contrastRatio, HoverContrast } from "../src/draw/hover.ts";
import { Dirty, HOVER_WIDTH, hoverInk, hoverWidth, Pass } from "../src/draw/index.ts";
import { DarkSkin, DefaultSkin, type Skin } from "../src/index.ts";
import { renderSvg } from "../src/svg/index.ts";
import { recording } from "./composer-fixtures.ts";
import { entryOf } from "./draw-fixtures.ts";
import { type Box, BoxKind, hoverDoc, passBoxes, within } from "./hover-fixtures.ts";
import { logoDoc } from "./logo-fixtures.ts";
import { lit } from "./overlay-fixtures.ts";

const lum = (c: number) => {
  const l = (x: number) => (x <= 10.31475 ? x / 3294.6 : ((x / 255 + 0.055) / 1.055) ** 2.4);
  return 0.2126 * l((c >> 16) & 255) + 0.7152 * l((c >> 8) & 255) + 0.0722 * l(c & 255);
};
const SURFACE = { light: 0xffffff, dark: 0x181a1e };

describe("the hover ink", () => {
  /** A 6 × 6 × 6 grid of opaque colors, and every tenth of them at a third, two thirds alpha. */
  const grid: number[] = [];
  for (let r = 0; r < 6; r++)
    for (let g = 0; g < 6; g++) for (let b = 0; b < 6; b++) grid.push(((r * 51) << 16) | ((g * 51) << 8) | (b * 51));
  const keys = [...grid, ...grid.filter((_, i) => i % 10 === 0).flatMap((c) => [c | (170 << 24), c | (85 << 24)])];
  const seen = (key: number, bg: number) => {
    const a = (255 - (key >>> 24)) / 255;
    const mix = (s: number) => Math.round(((key >> s) & 255) * a + ((bg >> s) & 255) * (1 - a)) << s;
    return mix(16) | mix(8) | mix(0);
  };

  it.each([
    ["light", DefaultSkin],
    ["dark", DarkSkin],
  ] as const)(
    "%s: 3:1 against the surface always, 1.5:1 against the fill wherever black or white allows it",
    (name, skin) => {
      const bg = SURFACE[name];
      const far = name === "light" ? 0 : 0xffffff;
      const onBg: number[] = [];
      const onFill: number[] = [];
      let reversed = 0;
      for (const key of keys) {
        const ink = hoverInk(skin, key);
        const fill = seen(key, bg);
        onBg.push(contrastRatio(lum(ink), lum(bg)));
        const cf = contrastRatio(lum(ink), lum(fill));
        onFill.push(cf);
        expect(onBg.at(-1), `#${key.toString(16)}`).toBeGreaterThanOrEqual(HoverContrast.background);
        // 1.5:1 against the fill always: away from the surface where the far end allows it (darker on a light skin,
        // lighter on a dark one), else the other way (a near-black bar on light gets a lighter outline inside it).
        expect(cf, `#${key.toString(16)}`).toBeGreaterThanOrEqual(HoverContrast.fill);
        const reach = contrastRatio(lum(far), lum(fill));
        const darker = lum(ink) <= lum(fill);
        expect(darker, `#${key.toString(16)}`).toBe((name === "light") === reach >= HoverContrast.fill);
        if (reach < HoverContrast.fill) reversed++;
      }
      const sorted = (a: number[]) => [...a].sort((x, y) => x - y);
      const q = (a: number[], p: number) => sorted(a)[Math.floor(p * (a.length - 1))].toFixed(2);
      console.log(
        `HV1 hover ink, ${name}, ${keys.length} colors: against the surface min ${q(onBg, 0)}, median ${q(onBg, 0.5)}; ` +
          `against the fill min ${q(onFill, 0)}, 10th percentile ${q(onFill, 0.1)}, median ${q(onFill, 0.5)}; ` +
          `${reversed} toward the surface`,
      );
    },
  );

  it("is deterministic, cached per color and skin, and keeps a hue: a blue's ink is a darker blue, a lighter one on dark", () => {
    expect(hoverInk(DefaultSkin, 0x4878d0)).toBe(hoverInk({ ...DefaultSkin } as Skin, 0x4878d0));
    const ink = hoverInk(DefaultSkin, 0x4878d0);
    expect(ink & 255).toBeGreaterThan((ink >> 16) & 255);
    expect(lum(ink)).toBeLessThan(lum(0x4878d0));
    const dark = hoverInk(DarkSkin, 0x4878d0);
    expect(dark & 255).toBeGreaterThan((dark >> 16) & 255);
    expect(lum(dark)).toBeGreaterThan(lum(0x4878d0));
    // A skin token resolves through the skin's own color: the ink (near black) gets a lighter ink on light.
    expect(lum(hoverInk(DefaultSkin, -4))).toBeGreaterThan(lum(0x1b1f24));
    expect(hoverWidth(2)).toBe(3);
    expect(hoverWidth(1)).toBe(2);
    expect(HOVER_WIDTH).toBe(1.5);
  });
});

type T = Awaited<ReturnType<typeof lit>>;
/** A context showing the hover document at DPR 2, in `skin`, zoomed to `cell` CSS px a position from position `from`. */
async function hoverAt(
  skin: "light" | "dark",
  cell?: number,
  style?: Record<string, unknown>,
  from = 0,
  n = 100,
  lineColor?: unknown,
) {
  const t = await lit(hoverDoc(n, { style, lineColor }), { width: 1000, height: 600, dpr: 2 });
  t.ctx.behaviors.get("hoverLabel")?.remove();
  if (skin === "dark") t.ctx.setSkin(DarkSkin);
  if (cell) t.vp?.setView(cell, from * cell, "api");
  await t.ctx.graph.update();
  return t;
}
/** What the overlay draws for the hover at a track's position (y: the entry's own, or `dy` CSS px off its top). */
function hovered(t: T, track: string, index: number, dy?: number) {
  const l = t.layout(t.section);
  const x = t.xOf(t.section, 0, index);
  const y = dy === undefined ? t.yOf(t.section, track) : t.yOf(t.section, track) - entryHeight(t, track) / 2 + dy;
  t.hover(x, y);
  const on = passBoxes(t.overlay(), Pass.overlay);
  t.ctx.setRenderProps({ hoverEmphasis: false });
  const off = passBoxes(t.overlay(), Pass.overlay);
  t.ctx.setRenderProps({ hoverEmphasis: true });
  const key = (b: Box) => b.join(",");
  const rest = new Map<string, number>();
  for (const b of off) rest.set(key(b), (rest.get(key(b)) ?? 0) + 1);
  // What the emphasis added: the boxes of the frame with it less those of the frame without it.
  const added = on.filter((b) => {
    const n = rest.get(key(b)) ?? 0;
    if (n > 0) rest.set(key(b), n - 1);
    return n === 0;
  });
  const main = passBoxes(t.frame(recording(), Dirty.tracks), Pass.main);
  return { added, main, l, x: x * 2, y: y * 2 };
}
function entryHeight(t: T, track: string): number {
  const l = t.layout(t.section);
  const k = t.trackIndex(track);
  for (let e = 0; e < l.entryCount; e++) if (l.entryTrack[e] === k) return l.entryHeight[e];
  throw new Error(track);
}
const skinOf = (s: string) => (s === "dark" ? DarkSkin : DefaultSkin);
/** Every added box inside a main-pass box of its source color (the ink of that color): the item's own pixels. */
function inside(added: Box[], main: Box[], skin: Skin) {
  for (const a of added) {
    const own = main.filter((m) => m[5] !== BoxKind.outline && hoverInk(skin, m[4]) === a[4] && within(a, m, 0.01));
    expect(own.length, `box ${a.join(",")} lies within no drawn box of its color`).toBeGreaterThan(0);
  }
}

describe.each(["light", "dark"] as const)("per renderer, %s", (skin) => {
  it("ranges: each fragment of the hovered item outlined inside its bar, not the item over the same positions", async () => {
    const t = await hoverAt(skin);
    // g1's second fragment [30, 45) alone at 42, on lane 0 (g2 [14, 40) is on lane 1).
    const { added, main } = hovered(t, "t.ranges", 42, 6);
    expect(added.length).toBe(2);
    expect(added.every((b) => b[5] === BoxKind.outline)).toBe(true);
    inside(added, main, skinOf(skin));
    // The two fragments: [10, 20) and [30, 45).
    const xs = added.map((b) => b[0]).sort((a, b) => a - b);
    expect(xs[0]).toBe(Math.round((t.ctx.layout.axis(t.section)?.x(0, 10 as never) ?? 0) * 2));
  });

  it("points: the hovered glyph outlined inset at the fit (its ring would leave its cell), ringed zoomed in", async () => {
    const t = await hoverAt(skin);
    const { added, main } = hovered(t, "t.points", 50);
    expect(added.map((b) => b[5]).sort()).toEqual([BoxKind.glyph, BoxKind.glyph]);
    inside(added.slice(0, 1), main, skinOf(skin));
    const z = await hoverAt(skin, 24, undefined, 40);
    const r = hovered(z, "t.points", 50);
    expect(r.added.map((b) => b[5])).toEqual([BoxKind.ring]);
    // Inside the position's cell (24 CSS px) and the entry's band.
    const cell = 24 * 2;
    expect(r.added[0][2] - r.added[0][0]).toBeLessThanOrEqual(cell);
    expect(Math.abs((r.added[0][0] + r.added[0][2]) / 2 - r.x)).toBeLessThan(cell / 2);
  });

  it("variants: the hovered marker alone, of the two stacked at its position", async () => {
    const t = await hoverAt(skin, 24);
    // A marker of 6 px in a pitch of 13: no room for a ring, so the marker outlined inset (its glyph in the ink, its
    // fill inside); one marker, though two are stacked at 12.
    const { added, main } = hovered(t, "t.variants", 12, 2 + 6);
    expect(added.map((b) => b[5])).toEqual([BoxKind.glyph, BoxKind.glyph]);
    inside(added.slice(0, 1), main, skinOf(skin));
    expect(within(added[1], added[0])).toBe(true);
  });

  it("secondary: a helix's box outlined inside, a strand's arrow inset in its own shape", async () => {
    const t = await hoverAt(skin);
    const helix = hovered(t, "t.secondary", 8);
    expect(helix.added.map((b) => b[5])).toEqual([BoxKind.outline]);
    inside(helix.added, helix.main, skinOf(skin));
    const strand = hovered(t, "t.secondary", 28);
    expect(strand.added.map((b) => b[5])).toEqual([BoxKind.shape, BoxKind.shape]);
    // Both inside the strand's own arrow as the main pass drew it.
    const arrow = strand.main.filter((m) => m[5] === BoxKind.shape && within(strand.added[0], m, 0.01));
    expect(arrow.length).toBe(1);
    expect(within(strand.added[1], arrow[0])).toBe(true);
  });

  it("values: the hovered position's bar outlined inside it; a line's point ringed only inside its column", async () => {
    const t = await hoverAt(skin);
    const bar = hovered(t, "t.values", 33);
    expect(bar.added.map((b) => b[5])).toEqual([BoxKind.outline]);
    inside(bar.added, bar.main, skinOf(skin));
    // A ring inside the position's column, at the fit (about 8 CSS px) and zoomed; none where the column is narrower.
    const ax = t.ctx.layout.axis(t.section);
    const column = [Math.round((ax?.x(0, 33 as never) ?? 0) * 2), Math.round((ax?.x(0, 34 as never) ?? 0) * 2)];
    const line = hovered(t, "t.line", 33);
    expect(line.added.map((b) => b[5])).toEqual([BoxKind.ring]);
    expect(line.added[0][0]).toBeGreaterThanOrEqual(column[0]);
    expect(line.added[0][2]).toBeLessThanOrEqual(column[1]);
    const thin = await hoverAt(skin, undefined, undefined, 0, 400);
    expect(hovered(thin, "t.line", 133).added).toEqual([]);
  });

  it("values: a banded line's ring in the hovered position's band (seqspec:format#layer-types)", async () => {
    const colors = [0xff7d45, 0xffdb13, 0x65cbf3, 0x0053d6];
    const hex = (c: number) => `#${c.toString(16).padStart(6, "0")}`;
    const palette = {
      kind: "discrete",
      mode: "absolute",
      colors: colors.map((c, k) => [hex(c), [0, 0.5, 0.7, 0.9][k]]),
    };
    const t = await hoverAt(skin, undefined, undefined, 0, 100, {
      theme: "field_palette",
      params: { field: "value", palette },
    });
    // Position 33's value is 0.83 (the third band), 34's 0.34 (the first): `hoverDoc`'s 0.2 + ((i * 31) % 80) / 100.
    for (const [i, band] of [
      [33, 2],
      [34, 0],
    ]) {
      const { added } = hovered(t, "t.line", i);
      expect(added.map((b) => b[5])).toEqual([BoxKind.ring]);
      expect(added[0][4], `position ${i}`).toBe(hoverInk(skinOf(skin), colors[band]));
    }
  });

  it("matrix: the hovered cell, its position's columns by its band", async () => {
    const t = await hoverAt(skin, 12);
    const h = entryHeight(t, "t.matrix");
    const { added } = hovered(t, "t.matrix", 7, h / 6);
    expect(added.length).toBe(1);
    const ax = t.ctx.layout.axis(t.section);
    const x0 = Math.round((ax?.x(0, 7 as never) ?? 0) * 2);
    const x1 = Math.round((ax?.x(0, 8 as never) ?? 0) * 2);
    expect([added[0][0], added[0][2]]).toEqual([x0, x1]);
    // Band 0 of 3: the entry's top third.
    const top = (t.yOf(t.section, "t.matrix") - h / 2) * 2;
    expect(added[0][1]).toBeGreaterThanOrEqual(Math.round(top));
    expect(added[0][3]).toBeLessThanOrEqual(Math.round(top + (2 * h) / 3) + 1);
  });

  it("logo: nothing over its stacked letters (David, 2026-10-12), while the profile's matrix below outlines its cell", async () => {
    const t = await lit(logoDoc(), { width: 1000, height: 600, dpr: 2 });
    if (skin === "dark") t.ctx.setSkin(DarkSkin);
    await t.ctx.graph.update();
    const e = entryOf(t.ctx, "Logo");
    t.hover(t.xOf(t.section, 0, 5), e.top + e.height / 2);
    const emphasis = (b: Box) => b[5] === BoxKind.outline || b[5] === BoxKind.ring || b[5] === BoxKind.glyph;
    expect(passBoxes(t.overlay(), Pass.overlay).filter(emphasis)).toEqual([]);
    const m = entryOf(t.ctx, "Matrix");
    t.hover(t.xOf(t.section, 0, 5), m.top + 2);
    const cell = passBoxes(t.overlay(), Pass.overlay).filter((b) => b[5] === BoxKind.outline);
    expect(cell.length).toBe(1);
    expect(cell[0][1]).toBeGreaterThanOrEqual(Math.round(m.top * 2));
  });

  it("letters: nothing (the hover band is its cell's column; an outline would sit over the letters)", async () => {
    const t = await hoverAt(skin, 12);
    expect(hovered(t, "t.letters", 7).added).toEqual([]);
  });
});

describe("what the emphasis follows", () => {
  it("the keyboard cursor's item: the hover it sets outlines the item under the cursor", async () => {
    const t = await hoverAt("light");
    const x = t.xOf(t.section, 0, 42);
    const y = t.yOf(t.section, "t.ranges") - entryHeight(t, "t.ranges") / 2 + 6;
    const hit = (await import("../src/index.ts")).Hit.create();
    expect(t.ctx.layout.hitTest(t.section, x, y, hit)).toBe(true);
    t.ctx.hover.setFromHit(hit, "keyboard");
    const on = passBoxes(t.overlay(), Pass.overlay).filter((b) => b[5] === BoxKind.outline);
    expect(on.length).toBeGreaterThanOrEqual(2);
  });

  it("an external hover (3D) on the same positions outlines nothing", async () => {
    const t = await hoverAt("light");
    const before = passBoxes(t.overlay(), Pass.overlay).length;
    t.ctx.hover.setLoci({ kind: "seq", seq: { a: [[30, 45]] } } as never, "3d");
    const boxes = passBoxes(t.overlay(), Pass.overlay);
    expect(boxes.filter((b) => b[5] === BoxKind.outline || b[5] === BoxKind.ring).length).toBe(0);
    expect(boxes.length).toBeGreaterThanOrEqual(before);
  });

  it("off: the render prop and the style key `hover_emphasis: false`", async () => {
    const t = await hoverAt("light", undefined, { hover_emphasis: false });
    expect(t.ctx.layout.props(t.section).hoverEmphasis).toBe(false);
    const { added } = hovered(t, "t.ranges", 42, 6);
    expect(added).toEqual([]);
    const x = t.xOf(t.section, 0, 42);
    t.hover(x, t.yOf(t.section, "t.ranges") - entryHeight(t, "t.ranges") / 2 + 6);
    expect(passBoxes(t.overlay(), Pass.overlay).filter((b) => b[5] === BoxKind.outline)).toEqual([]);
  });

  it("never in the main pass, nor in an export: `renderSvg` without `hover` draws none, with the overlay or without", async () => {
    const t = await hoverAt("light");
    const x = t.xOf(t.section, 0, 42);
    t.hover(x, t.yOf(t.section, "t.ranges") - entryHeight(t, "t.ranges") / 2 + 6);
    const ink = hoverInk(DefaultSkin, passBoxes(t.overlay(), Pass.overlay).find((b) => b[5] === 1)?.[4] ?? 0);
    expect(passBoxes(t.frame(recording(), Dirty.tracks), Pass.main).some((b) => b[5] === BoxKind.outline)).toBe(false);
    const hex = (c: number) => `#${c.toString(16).padStart(6, "0")}`;
    const outline = hex(passBoxes(t.overlay(), Pass.overlay).find((b) => b[5] === BoxKind.outline)?.[4] ?? -1);
    for (const overlay of [false, true]) {
      const svg = await renderSvg(t.ctx, { width: 1000, overlay });
      expect(svg.includes(outline), `overlay ${overlay}`).toBe(false);
    }
    // With `hover`, the host asks for the hover: its element's emphasis is part of it.
    expect((await renderSvg(t.ctx, { width: 1000, overlay: true, hover: true })).includes(outline)).toBe(true);
    expect(ink).toBeGreaterThanOrEqual(0);
  });
});
