/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AX1 (seqspec:model#pieces-and-layout-maps; seqspec:format#blocks): pieces in the format through the composer. A
// reversed piece's whole frame is the forward one mirrored within each block, pass by pass (the tracks, the overlay's
// marks for a selection,
// a highlight and the hover, the navigator and its marks); the ruler labels descend; a boundary sits on its position's
// right; hits read the piece; the SVG export of a reversed piece.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import type { ShownIndex } from "../src/axis.js";
import { Dirty, Pass } from "../src/draw/index.js";
import { Hit, HoverKind } from "../src/index.js";
import { renderSvg } from "../src/svg/index.js";
import { type MirrorItem, mirrorItems, type RecordingBackend, unmatched } from "../src/testing/index.js";
import { recording } from "./composer-fixtures.js";
import { directed } from "./conformance-fixtures.js";
import { entryOf } from "./draw-fixtures.js";
import { family, lit, rowLoci, seqLoci } from "./overlay-fixtures.js";
import { allTypes } from "./prepare-fixtures.js";
import { all, parseXml } from "./svg/xml.js";

type Lit = Awaited<ReturnType<typeof lit>>;

/** Per block its extent in a pass (device px): the tracks' and overlay's on the axis, the navigator's on its map. */
function extents(t: Lit, pass: number): number[][] {
  const v = t.ctx.layout.axis(t.section);
  const ax = v?.axis;
  if (!v || !ax) return [];
  const nav = pass === Pass.nav || pass === Pass.navOverlay;
  return Array.from({ length: ax.count }, (_, b) =>
    nav
      ? [Math.round(v.navX0s[b] * v.dpr), Math.round(v.navX1s[b] * v.dpr)]
      : [Math.round(v.x(b, 0 as ShownIndex) * v.dpr), Math.round(v.x(b, ax.shown[b] as ShownIndex) * v.dpr)],
  );
}
/** The center x of block b's index i as shown (CSS px). */
function xOf(t: Lit, b: number, i: number): number {
  const v = t.ctx.layout.axis(t.section);
  const k = v?.shownOf(b, i) ?? (-1 as ShownIndex);
  return v ? (v.x(b, k) + v.x(b, (k + 1) as ShownIndex)) / 2 : 0;
}
/**
 * The passes compared: the tracks and the overlay; zoomed, the navigator's overlay too (its window, marks and grips; at
 * the fit the window spans both blocks, mirrored as one). The navigator's map is an image of the tracks in miniature,
 * drawn by the renderers the kit checks.
 */
const PASSES = { main: Pass.main, overlay: Pass.overlay, navOverlay: Pass.navOverlay };

/** One protein of 120 residues with a ruler, letters, features, values and boundaries: one block, to zoom in on. */
function oneBlock() {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "a", label: "A", value: "MKVLATGHEKGSWQERTYIP".repeat(6) });
  const f = b.data.annotation({
    ref: "f",
    geometry: "ranges",
    target: a,
    value: { start: [3, 30, 70], stop: [12, 31, 95], label: ["one", "two", "three"] },
  } as never);
  const cut = b.data.annotation({ ref: "cut", geometry: "ranges", target: a, value: { boundary: [8, 20] } } as never);
  const vals = b.data.annotation({
    ref: "vals",
    geometry: "values",
    target: a,
    value: { value: Array.from({ length: 120 }, (_, i) => ((i * 31) % 100) / 100) },
  } as never);
  const v = b.view({ ref: "v", label: "V", height: "box" });
  v.block({ sequence: a });
  v.track({ label: "" }).ruler();
  v.track({ label: "Letters" }).letters();
  v.track({ label: "Features" }).ranges({ data: f, labels: true } as never);
  v.track({ label: "Values" }).values({ data: vals } as never);
  v.track({ label: "Boundaries" }).ranges({ data: cut } as never);
  return b.getDocument();
}

interface Case {
  name: string;
  doc: () => SeqViewSpecDoc;
  marks: (t: Lit) => void;
  /** Block, index and track label hovered. */
  hover: readonly [number, number, string];
  /** Zoomed to a cell: forward at the start, reversed at the end; mirrored across the data area and the map. */
  cell?: number;
}
const CASES: Case[] = [
  {
    name: "every layer type on two chains at the fit: a selection, a highlight and the hover",
    doc: () => allTypes({ length: 50 }),
    marks: (t) => {
      t.ctx.selection.set(seqLoci(t.set, "a", 5, 12), "api");
      t.ctx.highlight.set("host", seqLoci(t.set, "c", 20, 26));
    },
    hover: [0, 30, "Letters"],
  },
  {
    name: "an alignment's rows at the fit, a row's residues selected",
    doc: () => family(),
    marks: (t) => t.ctx.selection.set(rowLoci(t.set, 1, 1, 4), "api"),
    hover: [0, 4, "Rows"],
  },
  {
    name: "one block zoomed in: the forward start against the reversed end, and the navigator",
    doc: oneBlock,
    marks: (t) => t.ctx.selection.set(seqLoci(t.set, "a", 5, 12), "api"),
    hover: [0, 8, "Letters"],
    cell: 16,
  },
];

describe("a reversed piece's whole frame is the forward one mirrored, pass by pass", () => {
  for (const c of CASES)
    it(c.name, async () => {
      const frames: { t: Lit; rec: RecordingBackend }[] = [];
      for (const direction of [1, -1] as const) {
        const t = await lit(directed(c.doc(), direction), { width: 1000, height: 420 });
        if (c.cell) {
          t.vp?.setView(c.cell, 0, "api");
          if (direction < 0) t.vp?.setView(c.cell, 1e9, "api");
          await t.ctx.graph.update();
        }
        c.marks(t);
        const [b, i, label] = c.hover;
        const e = entryOf(t.ctx, label);
        t.hover(xOf(t, b, i), e.top + e.height / 2);
        expect(t.ctx.hover.live.current.kind).not.toBe(HoverKind.none);
        frames.push({ t, rec: t.frame(recording(), Dirty.all) });
      }
      const [fwd, rev] = frames;
      for (const [name, pass] of Object.entries(PASSES)) {
        if (pass === Pass.navOverlay && !c.cell) continue;
        const adv = (font: number) => fwd.t.ctx.text.advances(font);
        const whole = c.cell !== undefined;
        const area = (t: Lit) => {
          const v = t.ctx.layout.axis(t.section);
          return v ? [[Math.round(v.left), Math.round(v.right)]] : [];
        };
        const ends = (t: Lit) => (whole && pass !== Pass.navOverlay ? area(t) : extents(t, pass));
        const a = mirrorItems(fwd.rec, ends(fwd.t), true, adv, 1, pass);
        const b = mirrorItems(rev.rec, ends(rev.t), false, adv, 1, pass);
        const show = (x: MirrorItem[]) => x.slice(0, 3).map((q) => `${q.kind} ${q.x}+${q.w} @${q.y} ${q.key}`);
        expect(a.length, name).toBeGreaterThan(0);
        expect({ pass: name, unmatched: show(unmatched(a, b)), extra: show(unmatched(b, a)) }).toEqual({
          pass: name,
          unmatched: [],
          extra: [],
        });
      }
    });
});

/** 30 residues (`SEQ`) with a ruler and letters; `range` and `direction` on its one block. */
const SEQ = "MKVLATGHEKGSWQERTYIPACDNLLFRSV";
function short(o: { range?: [number, number]; direction?: 1 | -1 } = {}) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "a", label: "A", value: SEQ });
  const v = b.view({ ref: "v", label: "V", height: "box" });
  v.block({ sequence: a, ...o });
  v.track({ label: "" }).ruler();
  v.track({ label: "Letters" }).letters();
  return b.getDocument();
}
/** Text drawn in the main pass, left to right: the ruler's numbers and the letters. */
function textsByX(rec: RecordingBackend) {
  const t = rec.textItems;
  const out: [number, string][] = [];
  for (let i = 0; i < t.n; i++) out.push([t.x[i], t.text[i]]);
  return out.sort((a, b) => a[0] - b[0]).map((x) => x[1]);
}

describe("pieces in the format, through the viewer", () => {
  it("a ranged reversed piece: its letters from the range's last down, the ruler's numbers descending", async () => {
    const t = await lit(short({ range: [4, 24], direction: -1 }), { width: 900, height: 200 });
    const texts = textsByX(t.frame(recording(), Dirty.all));
    expect(texts.filter((s) => /^[A-Z]$/.test(s)).join("")).toBe([...SEQ.slice(4, 24)].reverse().join(""));
    const numbers = texts.filter((s) => /^\d+$/.test(s)).map(Number);
    expect(numbers.length).toBeGreaterThan(2);
    expect(numbers).toEqual([...numbers].sort((a, b) => b - a));
    expect(Math.max(...numbers)).toBeLessThanOrEqual(24);
    expect(Math.min(...numbers)).toBeGreaterThanOrEqual(5);
  });

  it("hits read the piece: each position's center gives its index; the hover label names the piece hovered", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", label: "A", value: SEQ });
    const v = b.view({ ref: "v", label: "V", height: "box" });
    v.block({ sequence: a, range: [2, 10], direction: -1, label: "first" });
    v.block({ sequence: a, range: [20, 28], label: "second" });
    v.track({ label: "Letters" }).letters();
    const t = await lit(b.getDocument(), { width: 900, height: 200 });
    const ax = t.ctx.layout.axis(t.section);
    if (!ax?.axis) throw new Error("no layout");
    const e = entryOf(t.ctx, "Letters");
    const y = e.top + e.height / 2;
    const hit = Hit.create();
    const read = (blk: number) =>
      Array.from({ length: ax.axis?.shown[blk] ?? 0 }, (_, k) => {
        const x = (ax.x(blk, k as ShownIndex) + ax.x(blk, (k + 1) as ShownIndex)) / 2;
        t.ctx.layout.hitTest(t.section, x, y, hit);
        return hit.kind === HoverKind.seq && hit.block === blk ? hit.index : -1;
      });
    expect(read(0)).toEqual([9, 8, 7, 6, 5, 4, 3, 2]);
    expect(read(1)).toEqual([20, 21, 22, 23, 24, 25, 26, 27]);
    t.hover(xOf(t, 1, 22), y);
    expect(t.ctx.labels.compute(t.ctx.hover.live.current)?.title).toBe(`second ${SEQ[22]}23`);
    t.hover(xOf(t, 0, 5), y);
    expect(t.ctx.labels.compute(t.ctx.hover.live.current)?.title).toBe(`first ${SEQ[5]}6`);
  });

  it("the hover label says which copy of a source residue a position is, past the first, either way (Ub 45 (copy 2))", async () => {
    const UB = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    for (const direction of [1, -1] as const) {
      const b = createSeqViewSpec();
      b.data.sequence({ ref: "ub", value: UB }).numbering({ name: "uniprot", method: "offset" });
      b.data.sequence({ ref: "diub", value: UB + UB });
      const segments = [
        [0, 0, 76],
        [76, 0, 76],
      ];
      b.data.mapping({ ref: "diub-ub", from: "diub", to: "ub", method: "segments", segments } as never);
      const v = b.view({ ref: "v", label: "V", numbering: "uniprot", direction });
      v.block({ sequence: "diub", label: "Ub2" }).track({ label: "Letters" }).letters();
      const t = await lit(b.getDocument(), { width: 1600, height: 200 });
      t.vp?.setView(10, direction > 0 ? 0 : 1e9, "api");
      await t.ctx.graph.update();
      const e = entryOf(t.ctx, "Letters");
      const y = e.top + e.height / 2;
      const title = (i: number) => {
        t.hover(xOf(t, 0, i), y);
        return t.ctx.labels.compute(t.ctx.hover.live.current)?.title;
      };
      expect([title(44), title(120)]).toEqual([`Ub2 ${UB[44]}45`, `Ub2 ${UB[44]}45 (copy 2)`]);
    }
  });

  it("the SVG export of a reversed piece: its letters and numbers as drawn", async () => {
    const xml = await renderSvg(short({ direction: -1 }), { width: 900 });
    const texts = all(parseXml(xml), "text")
      .map((e) => [Number(e.attrs.x), e.text] as const)
      .sort((p, q) => p[0] - q[0])
      .map((p) => p[1]);
    expect(texts.filter((s) => /^[A-Z]$/.test(s)).join("")).toBe([...SEQ].reverse().join(""));
    const numbers = texts.filter((s) => /^\d+$/.test(s)).map(Number);
    expect(numbers.length).toBeGreaterThan(2);
    expect(numbers).toEqual([...numbers].sort((p, q) => q - p));
  });
});
