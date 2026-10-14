/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Text measuring (LY1; seqspec:viewer#text-measuring) with the fixture measurer: label widths and fits, graphemes,
// numbers by arithmetic, wrapped runs, a resize without re-preparing, a font change, the wrapped groups' labels in
// each placement, none on a ligand and none on a group's continuation line. The heap guards
// are layout.heap.ts's.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { case1 } from "../../core/_test/format/cases.js";
import type { ShownIndex } from "../src/axis.js";
import { EntryKind, Hit, type RenderProps } from "../src/index.js";
import { GROUP_LABEL_PAD, groupLabelBand, groupLabelPx, layoutRuns, lineTable, type PreparedRun } from "../src/runs.js";
import { fixtureMeasurer, fixtureWidth } from "../src/testing/index.js";

import { fitCount, fitted, fontPx } from "../src/text.js";
import { set } from "./axis-fixtures.js";
import { opened } from "./layout-fixtures.js";

const ui = (text: string, px = 12) =>
  [...text].reduce((w, c) => w + fixtureWidth(c.codePointAt(0) as number, false), 0) * px;

describe("labels", async () => {
  const t = await opened(case1());
  const text = t.ctx.text;
  const f = text.font("ui", 12);

  it("measures once per font; cuts to a column with an ellipsis, the cut string built once", () => {
    const s = "Homo sapiens (human) hemoglobin subunit alpha";
    const p = text.label(s, f);
    expect(p.width).toBeCloseTo(ui(s));
    expect(text.label(s, f)).toBe(p);
    const ellipsis = text.label("…", f).width;
    const n = fitCount(p, 120, ellipsis);
    const cut = fitted(p, n);
    expect(cut.endsWith("…")).toBe(true);
    expect(ui(cut)).toBeLessThanOrEqual(120);
    expect(ui(fitted(p, n + 1))).toBeGreaterThan(120);
    expect(fitted(p, n)).toBe(cut);
    expect(fitted(p, p.ends.length)).toBe(s);
    expect(fitCount(p, 1000, ellipsis)).toBe(s.length);
    expect(fitCount(p, 1, ellipsis)).toBe(0);
  });

  it("cuts on grapheme boundaries: combining marks and an emoji sequence stay whole", () => {
    const s = "Café 👩‍🔬 lab";
    const p = text.label(s, f);
    const starts = [...new Intl.Segmenter(undefined, { granularity: "grapheme" }).segment(s)].map((g) => g.index);
    expect(p.ends.length).toBe(starts.length);
    expect(Array.from(p.offsets ?? [])).toEqual([...starts, s.length]);
    for (let k = 0; k < starts.length; k++) {
      expect(fitted(p, k)).toBe(`${s.slice(0, starts[k])}…`);
      expect(p.ends[k]).toBeCloseTo(ui(s.slice(0, starts[k + 1] ?? s.length)));
    }
    // "é" is one grapheme of two code units; the emoji one of five.
    expect(fitted(p, 4)).toBe("Café…");
  });

  it("numbers by arithmetic over digit advances; letter cells per alphabet", () => {
    const a = text.advances(f);
    const w = (s: string) => [...s].reduce((x, c) => x + a[c.charCodeAt(0)], 0);
    expect(text.numberWidth(f, -1234)).toBe(Math.ceil(w("-1234")));
    expect(text.numberWidth(f, 0)).toBe(Math.ceil(w("0")));
    expect(text.numberWidth(f, 27, "A".charCodeAt(0))).toBe(Math.ceil(w("27A")));
    const mono = text.font("mono", 15);
    expect(text.cell(mono, "protein")).toBe(9);
    expect(text.cell(f, "dna")).toBeCloseTo(0.65 * 12);
    expect(fontPx(text.fontString(mono))).toBe(15);
  });
});

describe("wrapped runs", () => {
  const run = (widths: number[], ends: number[], starts: number[], block = 0): PreparedRun => ({
    block,
    count: ends.length,
    segments: widths.length,
    widths: Float32Array.from(widths),
    positionEnds: null,
    startShown: Int32Array.from(starts),
    advance: ends[0],
    inline: false,
    labels: widths.map(() => ""),
    labelWidths: new Float32Array(widths.length),
  });

  it("fills lines greedily, each block from a new line; a group wider than a line is cut at a position", () => {
    const ten = Array.from({ length: 10 }, (_, k) => (k + 1) * 10);
    const a = run([50, 50], [...ten], [0, 5]);
    const b = run([100], ten, [0], 1);
    const out = lineTable(1);
    // 50 + 8 + 50 = 108 fits 110; b's 100 does not fit 60: cut at 6 positions per line.
    expect(layoutRuns([a, b], 110, 8, out)).toBe(2);
    expect(layoutRuns([a, b], 60, 8, out)).toBe(4);
    expect([0, 1, 2, 3].map((k) => [out.block[k], out.start[k], out.stop[k]])).toEqual([
      [0, 0, 5],
      [0, 5, 10],
      [1, 0, 6],
      [1, 6, 10],
    ]);
  });

  it("a CASE_1 text view resized from 900 to 600 px reruns only layoutRuns; a font change rebuilds once", async () => {
    const t = await opened(case1(), { size: { "text.main": [900, 400] } });
    const stats = t.layouts.stats;
    const lines = () => t.layout("text.main").lines.count;
    const wide = lines();
    const [prepared, laid] = [stats.prepareRuns, stats.layoutRuns];
    expect(prepared).toBe(1);
    t.vp("text.main").setSize(600, 400, 1);
    expect(lines()).toBeGreaterThan(wide);
    expect([stats.prepareRuns, stats.layoutRuns]).toEqual([prepared, laid + 1]);
    // Every line is the one block's, from its start.
    const l = t.layout("text.main").lines;
    expect([l.block[0], l.start[0], l.stop[l.count - 1]]).toEqual([0, 0, 393]);

    const label = t.ctx.text.label("p53", t.ctx.text.font("ui", 12));
    const builds = stats.builds;
    const version = t.ctx.text.version.value;
    const m = fixtureMeasurer();
    t.ctx.setTextMeasurer(m);
    expect(t.ctx.text.version.value).toBe(version + 1);
    expect(t.ctx.text.label("p53", t.ctx.text.font("ui", 12))).not.toBe(label);
    expect(m.calls).toBeGreaterThan(0);
    t.layout("text.main");
    t.layout("text.main");
    t.scheduler.advance(16);
    const sections = t.ctx.viewports.sections.value.length;
    expect(stats.builds - builds).toBe(sections);
    expect(stats.prepareRuns).toBe(prepared + 1);
    // A patch that keeps the axis (a layer's color) prepares nothing again.
    const ref = t.entries("text.main")[1][1];
    await t.ctx.applyPatch(set(ref, { hidden: true }));
    t.layout("text.main");
    expect(stats.prepareRuns).toBe(prepared + 1);
  });
});

describe("group labels: above each group's end, inline before its start, or none", () => {
  // 30 residues in groups of 10 under a ruler, three numberings: an offset, labels with insertion codes and unlabeled
  // positions.
  const doc = (numbering?: string) => {
    const b = createSeqViewSpec();
    const kabat = Array.from({ length: 30 }, (_, i) =>
      i < 10 ? i + 1 : i === 10 ? "10A" : i === 19 ? "18A" : i < 20 ? i : null,
    );
    b.data
      .sequence({ ref: "a", value: "MKVLATGHEK".repeat(3) })
      .numbering({ name: "auth", method: "offset", first: 101 })
      .numbering({ name: "kabat", method: "labels", labels: kabat });
    const v = b.view({ ref: "v", label: "V", layout: "text", ...(numbering ? { numbering } : {}) });
    v.track({ ref: "v.ruler" }).ruler();
    v.track({ ref: "v.letters" }).letters();
    return b.getDocument();
  };
  type Mode = RenderProps["wrappedGroupLabels"];
  // The data's left: the gutter (its lines' numbers) takes its 56 px only with no group labels.
  // The text's inset (8 px), then the gutter where it draws the lines' numbers.
  const INSET = 8;
  const left = (m: Mode) => INSET + (m === "none" ? 56 : 0);
  const L = left("aboveEnd");
  const mode = (m: Mode, d?: number) => ({
    // At the 15 px these widths are for (the default is the letters' cap, `maxLetterPx`), lines touching.
    render: { wrappedGroupLabels: m, wrappedFontPx: 15, wrappedGap: 0 },
    ...(d === undefined ? {} : { size: { "v.main": [left(m) + d + 8, 400] as const } }),
  });
  // Mono 15 px: 9 px cells, 90 px a group; the labels in mono 9 px (0.6 em a digit).
  const sup = (text: string) => text.length * 0.6 * groupLabelPx(15) + GROUP_LABEL_PAD;
  const T = EntryKind.track;
  const runOf = (t: Awaited<ReturnType<typeof opened>>) => {
    t.layout("v.main");
    return t.ctx.layout.axis("v.main")?.runs[0] as PreparedRun;
  };
  const labels = (t: Awaited<ReturnType<typeof opened>>) => {
    const run = runOf(t);
    return [[...run.labels], Array.from(run.labelWidths, (w) => w > 0)];
  };
  const hitAt = (t: Awaited<ReturnType<typeof opened>>, entry: number) => {
    const l = t.layout("v.main");
    const y = l.scrolledTop + l.entryTop[entry] + 4;
    return (x: number) => {
      const h = Hit.create();
      t.ctx.layout.hitTest("v.main", x, y, h);
      return h.shown;
    };
  };

  describe("above-end (the default)", () => {
    it("each group's last label, following the numbering; insertion codes as they are, none where unlabeled", async () => {
      const t = await opened(doc());
      expect(t.ctx.layout.props("v.main").wrappedGroupLabels).toBe("aboveEnd");
      expect(labels(t)).toEqual([
        ["10", "20", "30"],
        [true, true, true],
      ]);
      await t.ctx.applyPatch(set("v.main", { numbering: "auth" }));
      expect(labels(t)).toEqual([
        ["110", "120", "130"],
        [true, true, true],
      ]);
      await t.ctx.applyPatch(set("v.main", { numbering: "kabat" }));
      expect(labels(t)).toEqual([
        ["10", "18A", ""],
        [true, true, false],
      ]);
    });

    it("take no room: wrapping and hits are as without labels", async () => {
      const d = 2 * 90 + 8;
      const t = await opened(doc(), mode("aboveEnd", d));
      const run = runOf(t);
      // No per-position array without inline labels; a position's end is (k + 1) · advance.
      expect([run.inline, Array.from(run.widths), run.positionEnds, run.advance]).toEqual([
        false,
        [90, 90, 90],
        null,
        9,
      ]);
      const l = t.layout("v.main").lines;
      expect([l.count, l.start[1]]).toEqual([2, 20]);
      const view = t.ctx.layout.axis("v.main");
      expect(view && [view.x(0, 0 as ShownIndex), view.x(0, 10 as ShownIndex)]).toEqual([L, L + 98]);
      const at = hitAt(t, 1);
      expect([at(L + 0.5), at(L + 90 + 4), at(L + 98 + 0.5)]).toEqual([0, -1, 10]);
    });

    it("a band above each line holds them, in place of the ruler's (which draws nothing wrapped)", async () => {
      const t = await opened(doc(), mode("aboveEnd", 140));
      const band = groupLabelBand(15);
      expect(band).toBe(14);
      expect(t.entries("v.main").slice(0, 4)).toEqual([
        [EntryKind.labels, "block 0", 0, INSET, band],
        [T, "v.letters", 0, INSET + band, 24],
        [EntryKind.labels, "block 0", 1, INSET + band + 24, band],
        [T, "v.letters", 1, INSET + 2 * band + 24, 24],
      ]);
      // The band hits nothing.
      expect(hitAt(t, 0)(L + 85)).toBe(-1);
    });
  });

  describe("inlineStart", () => {
    it("each group's first label, following the numbering", async () => {
      const t = await opened(doc(), mode("inlineStart"));
      expect(labels(t)).toEqual([
        ["1", "11", "21"],
        [true, true, true],
      ]);
      await t.ctx.applyPatch(set("v.main", { numbering: "auth" }));
      expect(labels(t)).toEqual([
        ["101", "111", "121"],
        [true, true, true],
      ]);
      await t.ctx.applyPatch(set("v.main", { numbering: "kabat" }));
      expect(labels(t)).toEqual([
        ["1", "10A", ""],
        [true, true, false],
      ]);
    });

    it("line breaking counts their widths: two groups that fit without them wrap with them; no band", async () => {
      const d = 2 * 90 + 8;
      const on = await opened(doc(), mode("inlineStart", d));
      const run = runOf(on);
      expect(Array.from(run.widths)).toEqual(["1", "11", "21"].map((s) => Math.fround(sup(s) + 90)));
      const l = on.layout("v.main").lines;
      expect([l.count, l.start[1], l.start[2]]).toEqual([3, 10, 20]);
      expect(on.entries("v.main").slice(0, 2)).toEqual([
        [T, "v.letters", 0, INSET, 24],
        [T, "v.letters", 1, INSET + 24, 24],
      ]);
      // Room for both labels: two groups on a line again.
      const fits = await opened(doc(), mode("inlineStart", Math.ceil(d + sup("1") + sup("11"))));
      expect(fits.layout("v.main").lines.count).toBe(2);
    });

    it("hit testing maps x past a label to its letter; the label itself hits nothing", async () => {
      const t = await opened(doc(), mode("inlineStart", 1000));
      const view = t.ctx.layout.axis("v.main");
      const at = hitAt(t, 0);
      const s1 = sup("1");
      const s11 = sup("11");
      // The letters start after "1"; group 1 after the gap and "11".
      expect(view?.x(0, 0 as ShownIndex)).toBeCloseTo(L + s1, 4);
      expect(view?.x(0, 10 as ShownIndex)).toBeCloseTo(L + s1 + 90 + 8 + s11, 4);
      expect([at(L + s1 / 2), at(L + s1 + 0.5), at(L + s1 + 9 * 9 + 4)]).toEqual([-1, 0, 9]);
      const g1 = L + s1 + 90 + 8;
      expect([at(g1 + s11 / 2), at(g1 + s11 + 0.5), at(g1 + s11 + 9 + 4)]).toEqual([-1, 10, 11]);
    });
  });

  it("inline-start: a group cut onto a new line has no label on its continuation; x and hits start there", async () => {
    // 60 px: group 0 ("1" and ten letters, 96 px) is cut; the line after the cut starts mid-group, with no label.
    const t = await opened(doc(), mode("inlineStart", 60 - 8));
    const l = t.layout("v.main").lines;
    const s = l.start[1];
    expect([l.block[1], s > 0 && s < 10]).toEqual([0, true]);
    const view = t.ctx.layout.axis("v.main");
    if (!view) throw new Error("no view");
    view.lineStart = s;
    expect(view.x(0, s as ShownIndex)).toBe(L);
    expect(hitAt(t, 1)(L + 0.5)).toBe(s);
    // The next group's label is on its own first position, wherever its line starts.
    const next = Array.from({ length: l.count }, (_, n) => l.start[n]).indexOf(10);
    view.lineStart = 10;
    expect(next > 0 && view.x(0, 10 as ShownIndex) > L).toBe(true);
  });

  it("none on a ligand, in either placement", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEKMK" });
    b.data.ligand({ ref: "lig", value: ["C", "C", "N", "O"], chem_comp_ids: ["XYZ"] });
    const v = b.view({ ref: "v", label: "V", layout: "text" });
    v.block({ sequence: "a" }).block({ ligand: "lig" } as never);
    v.track({ ref: "v.letters" }).letters();
    for (const m of ["aboveEnd", "inlineStart"] as const) {
      const t = await opened(b.getDocument(), mode(m));
      t.layout("v.main");
      const runs = t.ctx.layout.axis("v.main")?.runs ?? [];
      expect([runs.length, [...runs[0].labels], [...runs[1].labels], runs[1].labelWidths[0]]).toEqual([
        2,
        m === "aboveEnd" ? ["10", "12"] : ["1", "11"],
        [""],
        0,
      ]);
    }
  });

  it("none: no labels and no band; switching the mode prepares the runs again", async () => {
    const t = await opened(doc(), mode("none", 1000));
    const run = runOf(t);
    expect([run.inline, [...run.labels], Array.from(run.widths)]).toEqual([false, ["", "", ""], [90, 90, 90]]);
    expect(t.entries("v.main")).toEqual([[T, "v.letters", 0, INSET, 24]]);
    const prepared = t.layouts.stats.prepareRuns;
    t.ctx.setRenderProps({ wrappedGroupLabels: "aboveEnd" });
    expect(labels(t)[0]).toEqual(["10", "20", "30"]);
    expect(t.layouts.stats.prepareRuns).toBe(prepared + 1);
  });
});
