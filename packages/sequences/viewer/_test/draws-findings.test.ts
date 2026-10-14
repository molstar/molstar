/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN1's findings on the recording backend (K1a, K2, K4, K5, L2): K1a's column oracle for cells (per device column the
// drawn color against a brute force at three zooms and two DPRs; a `mean` theme's palette(mean) over dense values and
// over a row field; a categorical theme's plurality class; more than 255 colors through the raw path; either block
// order), letters with their theme's alpha, the ink on a cell, every location filled, element symbols past barium and
// the letters LOD per block; the ruler reading the numbering live, its ticks and breaks on a filtered axis against a
// brute force, `visibleRange` in closed form against the linear walk, `numberText` bounded; emphasis clamped to each
// wrapped line, past 64 items in view as density, a ligand's label kept and atoms' symbols.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { Numbering } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { case3 } from "../../core/_test/format/cases.js";
import { SkinToken } from "../src/draw/backend.js";
import { visibleRange } from "../src/draw/cells.js";
import { numberText } from "../src/draw/numbers.js";
import type { DrawView } from "../src/draw/view.js";
import { type Theme, ThemeLocation } from "../src/index.js";
import { type LettersState, type PerPosition, rowValues } from "../src/prepare.js";
import { PYRAMID_BASE, valuePyramid } from "../src/pyramids.js";
import { CellInk } from "../src/themes.js";
import { set } from "./axis-fixtures.js";
import { columns, drawn, entryOf, layerOf, meanDoc, record, rects, texts } from "./draw-fixtures.js";

type Opened = Awaited<ReturnType<typeof drawn>>;
async function zoom(t: Opened, section: string, cell: number, x = 0) {
  t.ctx.viewports.get(section)?.setView(cell, x, "api");
  await t.ctx.graph.update();
}
const prepared = (t: Opened, k: number) =>
  t.ctx.graph.assembled.value?.view?.prepared[k] as { state: unknown; themes: readonly (Theme | null)[] };

/** Colors within one step per channel (0xRRGGBB); -2: none, equal only to none. */
const near = (a: number, b: number) =>
  a === b || (a >= 0 && b >= 0 && [16, 8, 0].every((s) => Math.abs(((a >> s) & 255) - ((b >> s) & 255)) <= 1));

/** Per device column of block b's view, the shown range it holds: [x, k0, k1]. */
function columnsOf(v: DrawView, b: number): [number, number, number][] {
  const out: [number, number, number][] = [];
  const range = new Int32Array(2);
  visibleRange(v, b, range);
  for (let k = range[0]; k < range[1]; ) {
    const x = v.xd(b, k);
    let k1 = k + 1;
    while (k1 < range[1] && v.xd(b, k1) === x) k1++;
    out.push([x, k, k1]);
    k = k1;
  }
  return out;
}

/**
 * The oracle for a `mean` theme (seqspec:sources#levels-and-the-zoom-gate): per device column the bins of the level its
 * positions hold whose
 * start lies in it (raw positions below the base), their values' mean by brute force, through the theme.
 */
/**
 * The block indices a device column of shown [k0, k1) reduces: the positions of the bins at its level whose start lies
 * in it (a level from the positions per column and the pyramid's `levels`; raw below the base: its own).
 */
function reduced(v: DrawView, b: number, k0: number, k1: number, levels: number): [number, number] {
  const ax = v.axis;
  if (!ax) throw new Error("no axis");
  const per = 1 / (v.cell * v.span[b] * v.dpr);
  const L =
    per >= 2 ** PYRAMID_BASE && levels > 0 ? Math.min(Math.floor(Math.log2(per)), PYRAMID_BASE + levels - 1) : -1;
  const base = ax.base[b];
  const [i0, i1] = [ax.indexOf(b, k0 as never), ax.indexOf(b, (k1 - 1) as never) + 1];
  if (L < 0) return [i0, i1];
  const size = 2 ** L;
  return [
    base + Math.ceil((i0 - base) / size) * size,
    Math.min(base + ax.length[b], base + Math.ceil((i1 - base) / size) * size),
  ];
}
/** The oracle for a `mean` theme: per device column its reduced positions' mean by brute force, through the theme. */
function meanOracle(v: DrawView, t: Theme, b: number, values: Float64Array, levels: number) {
  const ax = v.axis;
  if (!ax) throw new Error("no axis");
  const base = ax.base[b];
  const loc = ThemeLocation.create();
  const finite = values.filter((x) => !Number.isNaN(x));
  const [lo, hi] = [Math.min(...finite.slice(0, 1e5)), Math.max(...finite.slice(0, 1e5))];
  return columnsOf(v, b).map(([x, k0, k1]) => {
    const [i0, i1] = reduced(v, b, k0, k1, levels);
    let [sum, n] = [0, 0];
    for (let i = i0; i < i1; i++) if (!Number.isNaN(values[i - base])) [sum, n] = [sum + values[i - base], n + 1];
    if (n === 0) return [x, -2];
    ThemeLocation.at(loc, ax, v.layer, b, ax.indexOf(b, k0 as never));
    loc.annotation = loc.row = -1;
    // The mean of values lies within their range (a sum of equal values may pass it by an ulp).
    loc.value = Math.min(hi, Math.max(lo, sum / n));
    return [x, t.value(loc)];
  });
}

describe("K1a's oracle: a `mean` theme draws palette(mean) per device column", () => {
  for (const dpr of [1, 2])
    for (const per of [6, 40, 300])
      it(`dense values and a row field under Viridis, ${per} positions a column at DPR ${dpr}, either block order`, async () => {
        for (const reversed of [false, true]) {
          const doc = meanDoc(40_000, 2);
          if (reversed) {
            const view = doc.root.children?.find((c) => c.kind === "view");
            const section = view?.children?.[0];
            const blocks = section?.children?.filter((c) => c.kind === "block") ?? [];
            blocks.reverse();
            if (section?.children)
              section.children = [...blocks, ...section.children.filter((c) => c.kind !== "block")];
          }
          const t = await drawn(doc, { width: 1000, height: 420, dpr });
          const { section } = layerOf(t.ctx, "Shade");
          await zoom(t, section, 1 / (per * dpr));
          for (const label of ["Shade", "Scores"]) {
            const { k } = layerOf(t.ctx, label);
            const p = prepared(t, k);
            const s = p.state as PerPosition;
            const theme = p.themes[0] as Theme;
            expect(theme.reduction).toBe("mean");
            const { rec, frame } = record(t.ctx, label);
            const v = frame.v;
            const e = entryOf(t.ctx, label);
            const y = Math.round((e.top + e.height / 2) * dpr);
            let compared = 0;
            for (let b = 0; b < 2; b++) {
              const values =
                label === "Shade"
                  ? (s.value[b] as Float64Array)
                  : rowValues(s, b, theme.rowValue as (a: number, r: number) => number);
              const levels = valuePyramid(values, false).sum.length;
              const expected = meanOracle(v, theme, b, values, levels);
              const got = columns(rec, 0, 4000, y);
              for (const [x, color] of expected) {
                if (x < Math.round(e.left * dpr) || x >= Math.round(e.right * dpr)) continue;
                // The pyramid sums in another order than the brute force: a channel may round the other way.
                const where = `${label} block ${b} column ${x}${reversed ? " (reversed)" : ""}`;
                expect(near(got[x], color), `${where}: ${got[x].toString(16)} for ${color.toString(16)}`).toBe(true);
                compared++;
              }
            }
            expect(compared).toBeGreaterThan(200);
          }
          t.ctx.dispose();
        }
      });
});

/** A `color` layer over `categories` ranges of a 40,000-residue sequence, one category each (a categorical palette). */
function categoriesDoc(categories: number, n = 40_000): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  // Runs of one letter, so a letter holds most of a column too.
  const runs = Array.from({ length: n }, (_, i) => "ACDEFGHIKLMNPQRSTVWY"[Math.floor(i / 350) % 20]).join("");
  const a = b.data.sequence({ ref: "a", value: runs });
  const rows = Array.from({ length: categories }, (_, r) => r);
  const value = {
    start: rows.map((r) => Math.floor((r * n) / categories)),
    stop: rows.map((r) => Math.floor(((r + 1) * n) / categories)),
    category: rows.map((r) => `c${r}`),
  };
  const d = b.data.annotation({ ref: "d", geometry: "ranges", target: a, value } as never);
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: a });
  v.track({ label: "Letters" }).letters();
  const palette = { kind: "categorical", colors: "Category10", repeat_color_list: false };
  const many = categories > 255;
  const colors = many
    ? { kind: "categorical", colors: rows.map((r) => `#${(r * 40503).toString(16).padStart(6, "0").slice(-6)}`) }
    : palette;
  v.track({ label: "Codes" }).color({
    data: d,
    color: { theme: "field_palette", params: { field: "category", palette: colors } },
  } as never);
  return b.getDocument();
}

describe("K1a's oracle: a `classes` theme draws the plurality class per device column", () => {
  /** Per column its reduced positions' colors by brute force: the strict majority, or -1 (none holds most). */
  function majority(v: DrawView, color: (b: number, i: number) => number, b: number, levels: number) {
    return columnsOf(v, b).map(([x, k0, k1]) => {
      const [i0, i1] = reduced(v, b, k0, k1, levels);
      const n = new Map<number, number>();
      for (let i = i0; i < i1; i++) {
        const c = color(b, i);
        if (c >= 0) n.set(c, (n.get(c) ?? 0) + 1);
      }
      const [best, count] = [...n.entries()].sort((p, q) => q[1] - p[1])[0] ?? [-1, 0];
      return [x, count * 2 > i1 - i0 ? best : -1] as const;
    });
  }
  for (const dpr of [1, 2])
    for (const [categories, per] of [
      [40, 6],
      [40, 40],
      [300, 6],
      [300, 40],
    ] as const)
      it(`${categories} categories (${categories > 255 ? "raw" : "pyramid"}), letters too, ${per} positions a column at DPR ${dpr}`, async () => {
        const t = await drawn(categoriesDoc(categories), { width: 1000, height: 300, dpr });
        const { section } = layerOf(t.ctx, "Codes");
        await zoom(t, section, 1 / (per * dpr));
        for (const label of ["Codes", "Letters"]) {
          const { k } = layerOf(t.ctx, label);
          const p = prepared(t, k);
          const theme = p.themes[0] as Theme;
          const { rec, frame } = record(t.ctx, label);
          const v = frame.v;
          const ax = v.axis;
          if (!ax) throw new Error("no axis");
          const loc = ThemeLocation.create();
          const color = (b: number, i: number) => {
            ThemeLocation.at(loc, ax, k, b, i);
            if (label === "Letters") {
              loc.letter = ((p.state as LettersState).codes[b] as Uint8Array)[i - ax.base[b]];
              return theme.value(loc);
            }
            const s = p.state as PerPosition;
            loc.annotation = s.annotation[b]?.[i] ?? s.only[b];
            loc.row = s.row[b]?.[i] ?? -1;
            return loc.row < 0 ? -1 : theme.value(loc);
          };
          const e = entryOf(t.ctx, label);
          const got = columns(rec, 0, 4000, Math.round((e.top + e.height / 2) * dpr));
          let [held, equal] = [0, 0];
          // Past 255 colors a block reads raw positions; else the class pyramid's levels.
          const raw = label === "Codes" && categories > 255;
          const levels = raw ? 0 : valuePyramid(new Float64Array(ax.length[0]), false).sum.length;
          for (const [x, best] of majority(v, color, 0, levels)) {
            if (best < 0 || x < Math.round(e.left * dpr) || x >= Math.round(e.right * dpr)) continue;
            held++;
            if (got[x] === best) equal++;
          }
          expect(held, label).toBeGreaterThan(100);
          // Every column a class holds most of draws it (the tournament is exact there).
          expect(equal, `${label}: ${equal} of ${held}`).toBe(held);
        }
        t.ctx.dispose();
      });
});

describe("letters", () => {
  it("cells carry the theme's alpha byte (case 3's conservation shading on a row's letters)", async () => {
    const doc = case3();
    const view = doc.root.children?.find((c) => c.kind === "view");
    const section = view?.children?.[0];
    const color = { theme: "conservation_shading", params: { data: "cons", threshold: 0.4 } };
    const track = {
      kind: "track",
      ref: "family.shaded",
      params: { label: "Shaded" },
      children: [{ kind: "layer", ref: "family.shaded-l", params: { type: "letters", row: 0, color } }],
    };
    section?.children?.splice(1, 0, track as never);
    const t = await drawn(doc, { width: 1000, height: 600 });
    await zoom(t, layerOf(t.ctx, "Shaded").section, 10);
    const alphas = rects(record(t.ctx, "Shaded").rec).map((r) => r[5]);
    expect(alphas.some((a) => a > 0 && a < 255)).toBe(true);
    expect(alphas.some((a) => a === 255)).toBe(true);
  });

  it("dark ink on light cells, light on dark", async () => {
    for (const [fill, ink] of [
      ["#ffff00", CellInk.dark],
      ["#000080", CellInk.light],
    ] as const) {
      const b = createSeqViewSpec();
      b.data.sequence({ ref: "a", value: "MKVLATG" });
      b.view({ ref: "v", label: "V" })
        .track({ label: "Letters" })
        .letters({ color: { theme: "uniform", params: { color: fill } } } as never);
      const t = await drawn(b.getDocument(), { width: 600, height: 200 });
      await zoom(t, layerOf(t.ctx, "Letters").section, 16);
      const ts = texts(record(t.ctx, "Letters").rec);
      expect(ts.length).toBe(7);
      expect(ts.every((d) => d.color === ink)).toBe(true);
    }
  });

  it("atoms by element symbol past barium, drawn only where the block's widest symbol fits", async () => {
    const b = createSeqViewSpec();
    const atoms = ["C", "N", "FE", "W", "PT", "AU"];
    const lig = b.data.ligand({
      ref: "lig",
      value: Array.from({ length: 60 }, (_, i) => atoms[i % 6]),
      chem_comp_ids: ["X"],
    });
    b.view({ ref: "v", label: "V", ligand_form: "atoms" }).block({ ligand: lig }).track({ label: "Atoms" }).letters();
    const t = await drawn(b.getDocument(), { width: 600, height: 200 });
    const { section } = layerOf(t.ctx, "Atoms");
    await zoom(t, section, 30);
    expect(
      texts(record(t.ctx, "Atoms").rec)
        .map((d) => d.text)
        .slice(0, 6),
    ).toEqual(["C", "N", "Fe", "W", "Pt", "Au"]);
    // At 8 px a cell a two-letter symbol does not fit: cells only (the old LOD drew them overlapping).
    await zoom(t, section, 8);
    expect(texts(record(t.ctx, "Atoms").rec)).toEqual([]);
  });

  it("a ligand's atoms draw as cells in their element's color at every zoom, symbols where they fit (case 12's heme)", async () => {
    const b = createSeqViewSpec();
    // Heme: Fe, four N, carbons, four O.
    const heme = ["FE", "N", "N", "N", "N", ...Array.from({ length: 30 }, () => "C"), "O", "O", "O", "O"];
    const lig = b.data.ligand({ ref: "lig", value: heme, chem_comp_ids: ["HEM"] });
    b.view({ ref: "v", label: "V", ligand_form: "atoms" }).block({ ligand: lig }).track({ label: "Atoms" }).letters();
    const t = await drawn(b.getDocument(), { width: 600, height: 200 });
    const { section } = layerOf(t.ctx, "Atoms");
    // A pixel, a few px and a symbol's width a cell: always the cells (an element's color, the Fe's apart from the
    // C's), the symbols only at the last.
    for (const [cell, symbols] of [
      [0.5, false],
      [3, false],
      [8, false],
      [30, true],
    ] as const) {
      await zoom(t, section, cell);
      const { rec } = record(t.ctx, "Atoms");
      const colors = new Set(rects(rec).map((x) => x[4]));
      expect(colors.size, `cell ${cell}: element colors`).toBeGreaterThan(2);
      expect(rects(rec).every((x) => x[4] >= 0)).toBe(true);
      expect(texts(rec).length > 0, `cell ${cell}: symbols`).toBe(symbols);
    }
  });

  it("a `color` layer by unit colors every covered cell (the location's unit filled)", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(4) });
    const d = b.data.annotation({ ref: "d", geometry: "ranges", target: a, value: { start: [5], stop: [15] } });
    b.view({ ref: "v", label: "V" })
      .track({ label: "Units" })
      .color({ data: d, color: { theme: "by_unit" } } as never);
    const t = await drawn(b.getDocument(), { width: 600, height: 200 });
    await zoom(t, layerOf(t.ctx, "Units").section, 10);
    const r = rects(record(t.ctx, "Units").rec);
    expect(r.length).toBe(10);
    expect(r.every((x) => x[4] >= 0)).toBe(true);
  });
});

describe("the ruler", () => {
  it("a numbering switch relabels the ticks with no prepare rerun", async () => {
    const b = createSeqViewSpec();
    b.data
      .sequence({ ref: "a", value: "MKVLATGHEK".repeat(10) })
      .numbering({ name: "auth", method: "offset", first: 101 });
    const v = b.view({ ref: "v", label: "V" });
    v.track({ label: "" }).ruler();
    const t = await drawn(b.getDocument(), { width: 800, height: 200 });
    const { k, section } = layerOf(t.ctx, "");
    await zoom(t, section, 6);
    const state = prepared(t, k).state;
    const before = texts(record(t.ctx, "").rec).map((d) => d.text);
    expect((await t.ctx.applyPatch(set(section, { numbering: "auth" }))).ok).toBe(true);
    await t.ctx.graph.update();
    const after = texts(record(t.ctx, "").rec).map((d) => d.text);
    expect(prepared(t, k).state).toBe(state);
    expect(before.slice(0, 3)).toEqual(["10", "20", "30"]);
    expect(after.slice(0, 3)).toEqual(["110", "120", "130"]);
  });

  /** 2,000 residues in a `labels` numbering (with gaps and insertion codes), filtered to three runs. */
  function filtered(): SeqViewSpecDoc {
    const b = createSeqViewSpec();
    const n = 2000;
    const labels = Array.from({ length: n }, (_, i) => (i % 97 === 50 ? null : i % 211 === 7 ? `${i + 1}A` : i + 1));
    b.data
      .sequence({ ref: "a", value: "MKVLATGHEK".repeat(n / 10) })
      .numbering({ name: "auth", method: "labels", labels });
    const keep = b.data.annotation({
      ref: "keep",
      geometry: "ranges",
      target: "a",
      value: { start: [100, 500, 900], stop: [300, 520, 1400] },
    });
    const v = b.view({ ref: "v", label: "V", numbering: "auth", positions: [keep] } as never);
    v.track({ label: "" }).ruler();
    return b.getDocument();
  }
  it.each([0.5, 2, 9])(
    "on a filtered axis at %f px, ticks and breaks equal a brute force over every position",
    async (cell) => {
      const t = await drawn(filtered(), { width: 1000, height: 200 });
      const { section } = layerOf(t.ctx, "");
      await zoom(t, section, cell);
      const { rec, frame } = record(t.ctx, "");
      const v = frame.v;
      const ax = v.axis;
      if (!ax) throw new Error("no axis");
      const dpr = v.dpr;
      const r = rects(rec);
      const ticks = r.filter((x) => x[4] === SkinToken.inkQuiet).map((x) => x[0]);
      const breaks = r.filter((x) => x[4] === SkinToken.inkFaint).map((x) => x[0]);
      // Brute force: the labels of every shown position in view, the step as the draw's rule finds it.
      const n = v.numbering[0];
      const font = t.ctx.text.font("ui", 11);
      const strings = t.ctx.graph.assembled.value?.view?.set.strings;
      const labelOf = (k: number) => {
        const text = n && strings ? Numbering.format(n, ax.indexOf(0, k as never) as never, strings) : "";
        return /^-?\d+$/.test(text) ? Number(text) : Number.NaN;
      };
      const last = labelOf(ax.shown[0] - 1);
      const room = t.ctx.text.numberWidth(font, last) + 14;
      let step = 1e12;
      for (let p = 1; p < 1e12; p *= 10) {
        const c = [p, 2 * p, ...(p >= 10 ? [2.5 * p] : []), 5 * p].find((s) => s * v.cell >= room);
        if (c !== undefined) {
          step = c;
          break;
        }
      }
      const wantTicks: number[] = [];
      const wantBreaks: number[] = [];
      for (let k = v.visLo[0]; k < v.visHi[0]; k++) {
        const brk = k > 0 && ax.breakAt(0, k as never);
        if (brk) wantBreaks.push(v.xd(0, k));
        const m = labelOf(k);
        if (!Number.isNaN(m) && (m % step === 0 || brk)) wantTicks.push(Math.round(v.xd(0, k) + (v.cell * dpr) / 2));
      }
      expect(breaks).toEqual(wantBreaks);
      expect([...new Set(ticks)].sort((p, q) => p - q)).toEqual([...new Set(wantTicks)].sort((p, q) => p - q));
      expect(wantTicks.length).toBeGreaterThan(2);
      expect(dpr).toBe(1);
    },
  );

  it("number strings come from one bounded cache, a collision replacing its slot", () => {
    const a = numberText(12345);
    expect(numberText(12345)).toBe(a);
    for (let n = 0; n < 20_000; n++) numberText(n);
    expect(numberText(12345)).toBe("12345");
    expect(numberText(12, "GFP")).toBe("GFP 12");
    expect(numberText(12)).toBe("12");
  });
});

describe("visibleRange in closed form", () => {
  it("equals the old linear walk to the device column's edges, at random views, zooms and DPRs", async () => {
    const t = await drawn(meanDoc(40_000, 2), { width: 1000, height: 300 });
    const { section } = layerOf(t.ctx, "Shade");
    let seed = 5;
    const rnd = () => {
      seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
      return seed / 2 ** 32;
    };
    const out = new Int32Array(2);
    let checked = 0;
    for (const dpr of [1, 1.5, 2, 3]) {
      for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 300, dpr);
      for (let s = 0; s < 40; s++) {
        const cell = 10 ** (-4 + 4 * rnd());
        await zoom(t, section, cell, rnd() * 80_000 * cell);
        const v = record(t.ctx, "Shade").frame.v;
        for (let b = 0; b < 2; b++) {
          visibleRange(v, b, out);
          let lo = v.visLo[b];
          let hi = v.visHi[b];
          if (lo < hi && v.cell * v.span[b] * v.dpr < 1) {
            const c0 = v.xd(b, lo);
            while (lo > 0 && v.xd(b, lo - 1) === c0) lo--;
            const c1 = v.xd(b, hi - 1);
            while (hi < (v.axis?.shown[b] ?? 0) && v.xd(b, hi) === c1) hi++;
          }
          expect([out[0], out[1]], `cell ${cell} dpr ${dpr} block ${b}`).toEqual(lo < hi ? [lo, hi] : [0, 0]);
          checked++;
        }
      }
    }
    expect(checked).toBe(320);
  });
});

describe("emphasis", () => {
  it("wrapped, each line's emphasized letters are its own, clamped to its positions", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(8) });
    const d = b.data.annotation({ ref: "d", geometry: "ranges", target: a, value: { start: [25], stop: [40] } });
    const v = b.view({ ref: "v", label: "V", layout: "text" });
    const track = v.track({ label: "Letters" });
    track.letters();
    track.emphasis({ data: d });
    // About 30 residues a line: the range wraps.
    const t = await drawn(b.getDocument(), { width: 300, height: 600, render: { wrappedGroupLabels: "none" } });
    const { rec, frame } = record(t.ctx, "Letters", undefined, false);
    expect(rec.textItems.n).toBeGreaterThan(0);
    const l = t.ctx.layout.get(layerOf(t.ctx, "Letters").section)?.current;
    expect((l?.lines.count ?? 0) > 1).toBe(true);
    // Emphasis is the track's second layer: its own frame per line.
    const k = layerOf(t.ctx, "Letters").k + 1;
    const entry = t.ctx.registry.get("renderers", "emphasis");
    const { LayerFrame, RecordingBackend } = await import("../src/testing/index.js");
    const f = new LayerFrame(t.ctx, entry as never, layerOf(t.ctx, "Letters").section, k);
    const r = new RecordingBackend();
    f.draw(r);
    const letters = texts(r);
    expect(letters.length).toBe(15);
    // Every underline lies on one line's positions: within the data area, never past a line's last letter.
    const unders = rects(r).filter((x) => x[4] === SkinToken.accent);
    expect(unders.length).toBe(15);
    expect(frame.v.geometry).toBe(1);
  });

  it("past 64 items in view below the letters LOD, density: at most a bar per device column and color", async () => {
    const t = await drawn(meanDoc(40_000), { width: 1000, height: 300 });
    const { section } = layerOf(t.ctx, "Marks");
    for (const cell of [0.05, 0.4]) {
      await zoom(t, section, cell);
      const r = rects(record(t.ctx, "Marks").rec);
      const xs = new Set(r.map((x) => x[0]));
      expect(r.length, `cell ${cell}`).toBeLessThanOrEqual(xs.size * 8);
      expect(xs.size, `cell ${cell}`).toBeGreaterThan(100);
      expect(
        r.every((x) => x[2] === 1),
        `cell ${cell}: one device px each`,
      ).toBe(true);
    }
    // 64 items or fewer in view: each its cell.
    await zoom(t, section, 3);
    const items = rects(record(t.ctx, "Marks").rec);
    expect(items.length).toBeLessThanOrEqual(64);
    expect(items.every((x) => x[2] >= 3)).toBe(true);
    // At letters: the 7th positions' letters, bold, over their cleared cells.
    await zoom(t, section, 16);
    const { rec } = record(t.ctx, "Marks");
    expect(texts(rec).length).toBeGreaterThan(5);
    expect(rects(rec).some((x) => x[4] === SkinToken.surface)).toBe(true);
  });

  it("over a ligand's one position its label stays (an underline only); over atoms, their symbols", async () => {
    for (const form of ["position", "atoms"] as const) {
      const b = createSeqViewSpec();
      const lig = b.data.ligand({ ref: "lig", value: ["C", "N", "FE", "O"], chem_comp_ids: ["HEM"] });
      const d = b.data.annotation({ ref: "d", geometry: "ranges", target: lig, value: { start: [0], stop: [4] } });
      const track = b.view({ ref: "v", label: "V", ligand_form: form }).block({ ligand: lig }).track({ label: "Lig" });
      track.letters();
      track.emphasis({ data: d });
      const t = await drawn(b.getDocument(), { width: 600, height: 200 });
      await zoom(t, layerOf(t.ctx, "Lig").section, 30);
      const { LayerFrame, RecordingBackend } = await import("../src/testing/index.js");
      const k = layerOf(t.ctx, "Lig").k + 1;
      const f = new LayerFrame(
        t.ctx,
        t.ctx.registry.get("renderers", "emphasis") as never,
        layerOf(t.ctx, "Lig").section,
        k,
      );
      const r = new RecordingBackend();
      f.draw(r);
      if (form === "position") {
        expect(texts(r)).toEqual([]);
        expect(rects(r).filter((x) => x[4] === SkinToken.surface)).toEqual([]);
        expect(rects(r).length).toBeGreaterThan(0);
      } else expect(texts(r).map((x) => x.text)).toEqual(["C", "N", "Fe", "O"]);
    }
  });
});
