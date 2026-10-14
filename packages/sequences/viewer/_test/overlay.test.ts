/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's overlay on the recording backend (CS1b; seqspec:interaction#the-skin; K3): what the channels light,
// each block's marks within its own lane when stacked; row marks across header, data columns and data, a selected row's
// band over its gaps; the heads band at `headerWidth` 0 with the hovered row across the data columns, the data at half
// strength; one `paintMark`, so a selection's border is as wide wrapped as linear; the keyboard cursor's ring after an
// arrow key; renderers' `overlay` members and the navigator's marks; highlights of style `band` in the main pass; a
// tile arrival redrawing only its columns, in its own section only; and the ruler's copies ("48 (copy 2)").

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type Loci, OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Dirty, Pass, SectionComposer, SkinToken } from "../src/draw/index.ts";
import type { DrawView } from "../src/draw/view.ts";
import { HoverKind, KernelRenderers, type LayerRenderer, MarkBit, SeqViewerConfig } from "../src/index.ts";
import type { OverlayMarks } from "../src/marks.ts";
import type { Streams, StreamView } from "../src/streams.ts";
import { Op } from "../src/testing/index.ts";
import { COPY_BACK, composed, recording, textsOf } from "./composer-fixtures.ts";
import { driven } from "./input-fixtures.ts";
import { lazyLoaders } from "./layout-fixtures.ts";
import { family, lit, passRects, rowLoci, seqLoci, twoBlocks } from "./overlay-fixtures.ts";
import { allTypes } from "./prepare-fixtures.ts";

const S = SkinToken;
/** Rects of one color: [x, y, w, h, color, alpha, stroke]. */
const of = (rects: number[][], color: number, alpha = -1) =>
  rects.filter((r) => r[4] === color && (alpha < 0 || r[5] === alpha));

describe("position marks", () => {
  it("stacked: a block's marks stay within its own lane, never across the other block's", async () => {
    const t = await lit(twoBlocks({ stacked: true }), { width: 800, height: 400 });
    t.ctx.selection.set(seqLoci(t.set, "a", 10, 13), "api");
    t.ctx.highlight.set("host", seqLoci(t.set, "c", 4, 6));
    const r = passRects(t.overlay(), Pass.overlay);
    const l = t.layout(t.section);
    const lane = (b: number) => {
      const [e0, e1] = [l.blockEntries[2 * b], l.blockEntries[2 * b + 1] - 1];
      return [l.scrolledTop + l.entryTop[e0], l.scrolledTop + l.entryTop[e1] + l.entryHeight[e1]].map(Math.round);
    };
    const [a, c] = [lane(0), lane(1)];
    expect(a[1]).toBeLessThanOrEqual(c[0]);
    const sel = of(r, S.selection);
    const hl = of(r, S.highlight);
    expect([sel.length, hl.length]).toEqual([1, 1]);
    expect([sel[0][1], sel[0][1] + sel[0][3]]).toEqual(a);
    expect([hl[0][1], hl[0][1] + hl[0][3]]).toEqual(c);
    // Over the cells' pixels: residues 10 to 12 of block a.
    const ax = t.ctx.layout.axis(t.section);
    expect(sel[0][0]).toBe(Math.round(ax?.x(0, 10 as never) ?? 0));
    // Its sides outlined, one device px each.
    expect(of(r, S.selectionBorder).map((x) => [x[2], x[3]])).toEqual([
      [1, a[1] - a[0]],
      [1, a[1] - a[0]],
    ]);
  });

  it("side by side: a column band from the navigator's bottom to the last entry drawn, the hover over it", async () => {
    const t = await lit(twoBlocks(), { width: 900, height: 400 });
    t.ctx.selection.set(seqLoci(t.set, "c", 3, 4), "api");
    const l = t.layout(t.section);
    t.hover(t.xOf(t.section, 1, 3), l.scrolledTop + l.entryTop[l.pinnedCount] + 4);
    const r = passRects(t.overlay(), Pass.overlay);
    const [sel] = of(r, S.selection);
    let bottom = 0;
    for (let e = 0; e < l.entryCount; e++)
      bottom = Math.max(bottom, (e < l.pinnedCount ? 0 : l.scrolledTop) + l.entryTop[e] + l.entryHeight[e]);
    expect([sel[1], sel[1] + sel[3]]).toEqual([Math.round(l.navHeight), Math.round(Math.min(bottom, l.height))]);
    expect(of(r, S.hover).length).toBeGreaterThan(0);
  });

  it("one `paintMark`: a selection's border as wide wrapped as linear, at DPR 2", async () => {
    const widths = async (wrapped: boolean) => {
      const t = await lit(twoBlocks({ wrapped }), { width: 600, height: 600, dpr: 2 });
      t.vp?.setView(12, 0, "api");
      t.ctx.selection.set(seqLoci(t.set, "a", 5, 9), "api");
      return of(passRects(t.overlay(), Pass.overlay), S.selectionBorder).map((x) => x[2]);
    };
    const [linear, wrapped] = [await widths(false), await widths(true)];
    expect(linear.length).toBeGreaterThan(0);
    expect(wrapped.length).toBeGreaterThan(0);
    expect(new Set([...linear, ...wrapped])).toEqual(new Set([2]));
  });
});

describe("row marks", () => {
  /** Row d's band (device px, DPR 1) in a family's rows. */
  const rowY = (t: Awaited<ReturnType<typeof lit>>, d: number) => {
    const l = t.layout(t.section);
    const top = l.scrolledTop + l.entryTop[l.rows.entry] - (t.vp?.live.current.scrollY ?? 0) + d * l.rows.height;
    return [Math.round(top), Math.round(top + l.rows.height)];
  };

  it("a selected row: its label and data columns lit, its band over its gap, no column lit", async () => {
    const t = await lit(family(), { width: 800, height: 300 });
    t.vp?.setView(30, 0, "api");
    t.ctx.selection.set(rowLoci(t.set, 0, 0, 5), "api");
    const r = passRects(t.overlay(), Pass.overlay);
    const l = t.layout(t.section);
    const [y0, y1] = rowY(t, 0);
    const sel = of(r, S.selection);
    // Header and data columns: from 0 to the data area; the data: columns 0 to 6, the gap at 2 included.
    const ax = t.ctx.layout.axis(t.section);
    const x = (k: number) => Math.round(ax?.x(0, k as never) ?? 0);
    expect(sel.map((s) => [s[0], s[1], s[2], s[3]])).toEqual([
      [0, y0, Math.round(l.dataLeft), y1 - y0],
      [x(0), y0, x(6) - x(0) - 1, y1 - y0],
    ]);
    // The cells outlined as a box.
    expect(of(r, S.selectionBorder).map((s) => s[6])).toEqual([1]);
  });

  it("aligned cells: rectangles on their rows, a gap's too; the hover lights the cell's column and its row", async () => {
    const t = await lit(family(), { width: 800, height: 300 });
    t.vp?.setView(30, 0, "api");
    const cells = { kind: "aligned-loci", set: t.set, alignment: 0, elements: [{ row: 2, columns: 2 }] };
    t.ctx.highlight.set("host", {
      ...cells,
      elements: [{ row: 2, columns: OrderedSet.ofRange(1, 4) }],
    } as unknown as Loci);
    const [y0, y1] = rowY(t, 2);
    const r = passRects(t.overlay(), Pass.overlay);
    expect(of(r, S.highlight).map((s) => [s[1], s[3]])).toEqual([[y0, y1 - y0]]);
    // The hover on row 1, column 4: its row across (the data at half strength), its column over every row.
    const [h0, h1] = rowY(t, 1);
    t.hover(t.xOf(t.section, 0, 4), (h0 + h1) / 2);
    const o = passRects(t.overlay(), Pass.overlay);
    const l = t.layout(t.section);
    expect(of(o, S.hover, 255).filter((s) => s[1] === h0 && s[0] === 0)[0]?.[2]).toBe(Math.round(l.dataLeft));
    expect(of(o, S.hover, 128).map((s) => [s[0], s[1], s[3]])).toEqual([[Math.round(l.dataLeft), h0, h1 - h0]]);
    expect(of(o, S.hover, 255).some((s) => s[3] > 3 * (h1 - h0))).toBe(true);
  });

  it("heads at `headerWidth` 0: the heads band drawn, the hovered row across the data columns", async () => {
    const t = await lit(family({ bare: true }), { width: 800, height: 300, render: { headerWidth: 0 } });
    const l = t.layout(t.section);
    expect(l.columnCount).toBe(1);
    expect(textsOf(t.frame())).toContain("Score");
    const [h0, h1] = rowY(t, 1);
    t.hover(t.xOf(t.section, 0, 3), (h0 + h1) / 2);
    const o = passRects(t.overlay(), Pass.overlay);
    const band = of(o, S.hover, 255).find((s) => s[0] === 0 && s[1] === h0);
    expect(band?.[2]).toBe(Math.round(l.dataLeft));
    expect(band?.[2]).toBeGreaterThanOrEqual(Math.round(l.columnX[0] + l.columnWidth[0]));
  });

  it("heads on a row of their own (`ownRow`): row marks and the hovered row on their rows, none in the heads' row", async () => {
    const t = await lit(family(), { width: 800, height: 300, render: { dataColumnHeads: "ownRow" } });
    t.vp?.setView(30, 0, "api");
    const l = t.layout(t.section);
    expect([l.pinnedCount, l.headsTop, l.scrolledTop - l.headsTop]).toEqual([
      1,
      l.entryTop[0] + l.entryHeight[0],
      t.ctx.layout.props(t.section).rowHeight,
    ]);
    t.ctx.selection.set(rowLoci(t.set, 0, 0, 5), "api");
    const [y0, y1] = rowY(t, 0);
    expect(y0).toBe(Math.round(l.scrolledTop));
    const sel = of(passRects(t.overlay(), Pass.overlay), S.selection);
    expect(sel.map((s) => [s[0], s[1], s[3]])).toEqual([
      [0, y0, y1 - y0],
      [Math.round(t.ctx.layout.axis(t.section)?.x(0, 0 as never) ?? 0), y0, y1 - y0],
    ]);
    // The hovered row 1 across header and data columns; a point in the heads' row hovers no row.
    const [h0, h1] = rowY(t, 1);
    t.hover(t.xOf(t.section, 0, 4), (h0 + h1) / 2);
    const o = passRects(t.overlay(), Pass.overlay);
    expect(of(o, S.hover, 255).filter((s) => s[1] === h0 && s[0] === 0)[0]?.[2]).toBe(Math.round(l.dataLeft));
    t.hover(l.columnX[0] + 4, l.headsTop + l.headsHeight / 2);
    const none = passRects(t.overlay(), Pass.overlay);
    expect(of(none, S.hover).filter((s) => s[0] === 0)).toEqual([]);
  });
});

describe("the ring, renderers' overlays and the navigator's marks", () => {
  it("after an arrow key, the ring at the cursor's cell (`cursorRect`, `--seq-focus-ring`)", async () => {
    const t = await driven(twoBlocks(), { renderers: KernelRenderers as LayerRenderer[], size: { v: [800, 300] } });
    const section = t.ctx.viewports.sections.value[0];
    const c = new SectionComposer(t.ctx, section);
    const rec = recording();
    c.frame(rec, Dirty.overlay);
    expect(of(passRects(rec, Pass.overlay), S.focusRing)).toEqual([]);
    t.focus(section);
    t.key(section, "ArrowRight");
    const cur = t.ctx.ui.keyCursor.value;
    expect(cur?.index).toBe(1);
    rec.reset();
    c.frame(rec, Dirty.overlay);
    const [ring] = of(passRects(rec, Pass.overlay), S.focusRing);
    const ax = t.ctx.layout.axis(section);
    expect(ring[6]).toBe(2);
    expect(ring[0]).toBe(Math.round(ax?.x(0, 1 as never) ?? 0) - 2);
  });

  it("a renderer's `overlay` member per entry, reading the section's marks; the navigator draws them", async () => {
    const seen: number[] = [];
    const letters = (KernelRenderers as LayerRenderer[]).find((r) => r.name === "letters") as LayerRenderer;
    const spy = { ...letters, overlay: (v: DrawView) => void seen.push(v.marks.bitsAt(0, 10)) } as LayerRenderer;
    const t = await lit(twoBlocks(), { width: 800, height: 300, replace: { letters: spy } });
    t.ctx.selection.set(seqLoci(t.set, "a", 10, 11), "api");
    const r = t.frame(recording(), Dirty.all);
    expect(seen).toEqual([MarkBit.selection]);
    expect(of(passRects(r, Pass.navOverlay), S.selectionBorder).length).toBeGreaterThan(0);
  });
});

describe("the hovered item", () => {
  it("is the marks' for its layer only (`OverlayMarks.item`, `itemAnnotation`): the pairs overlay lights its own pair", async () => {
    const t = await lit(allTypes(), { width: 900, height: 600 });
    const l = t.layout(t.section);
    const m = () => (t.composer as unknown as { channels: { marks(s: number): OverlayMarks } }).channels.marks(0);
    t.overlay();
    const h = t.ctx.hover.live.current;
    for (let e = l.pinnedCount; e < l.entryCount && h.kind !== HoverKind.item; e++)
      for (let x = l.dataLeft + 2; x < l.dataRight && h.kind !== HoverKind.item; x += 3)
        t.hover(x, l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2);
    expect(h.kind).toBe(HoverKind.item);
    expect([m().item(h.layer), m().item(h.layer + 1)]).toEqual([h.item.item, -1]);
    // With its annotation: a layer drawn from two annotations has an item of that id in each.
    expect(h.item.annotation).toBeGreaterThanOrEqual(0);
    expect([m().itemAnnotation(h.layer), m().itemAnnotation(h.layer + 1)]).toEqual([h.item.annotation, -1]);
    t.ctx.hover.clear("pointer");
    expect([m().item(h.layer), m().itemAnnotation(h.layer)]).toEqual([-1, -1]);
  });
});

describe("the main pass", () => {
  it("a highlight of style `band`: a full-height band behind the entries, the tracks redrawn on its change", async () => {
    const t = await lit(twoBlocks(), { width: 800, height: 300 });
    const be = recording(COPY_BACK);
    t.frame(be);
    expect(t.frame(be).opCount).toBe(0);
    t.ctx.highlight.set("site", seqLoci(t.set, "a", 20, 30), { style: "band" });
    const r = t.frame(be);
    const l = t.layout(t.section);
    const [band] = of(passRects(r, Pass.main), S.highlight);
    expect([band[1], band[3]]).toEqual([Math.round(l.navHeight), Math.round(l.height) - Math.round(l.navHeight)]);
    // Drawn before the entries: the first rects of the frame.
    expect(r.rectItems.color[0]).toBe(S.highlight);
    t.ctx.highlight.clear("site");
    expect(of(passRects(t.frame(be), Pass.main), S.highlight)).toEqual([]);
  });

  it("a stream's arrival redraws only the columns of the regions that changed", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(200) });
    const vals = b.data.annotation({ ref: "vals", geometry: "values", target: a, value: { value: [0.5, 1, 0.2] } });
    const g = b.view({ ref: "g", label: "G" });
    const one = g.section({ ref: "g.one" }).block({ sequence: a });
    one.track({ label: "V" }).values({ ref: "g.vals", data: vals });
    g.section({ ref: "g.two" }).block({ sequence: a }).track({ label: "L" }).letters();
    const t = await composed(b.getDocument(), { section: "g.one", width: 1000 });
    t.vp?.setView(4, 0, "api");
    const other = new SectionComposer(t.ctx, "g.two");
    // A layer held by region, whose changes are columns [40, 60): `changed` as the stream manager's.
    const fake = {
      version: 0,
      changed: (since: number, out: Int32Array) => {
        out[0] = 40;
        out[1] = 60;
        return since === fake.version - 1 ? 1 : -1;
      },
    } as unknown as StreamView & { version: number };
    const get = t.ctx.streams.get.bind(t.ctx.streams);
    (t.ctx.streams as { get: typeof get }).get = (ref: string) => (ref === "g.vals" ? fake : get(ref));
    const [be, quiet] = [recording(COPY_BACK), recording(COPY_BACK)];
    t.frame(be);
    other.frame(quiet, Dirty.tracks);
    fake.version++;
    const r = t.frame(be);
    quiet.reset();
    other.frame(quiet, Dirty.tracks);
    const ax = t.ctx.layout.axis(t.section);
    const x = (k: number) => Math.round(ax?.x(0, k as never) ?? 0);
    const k = r.ops.findIndex((_, i) => i % 4 === 0 && r.ops[i] === Op.clear) / 4;
    const at = r.ops[4 * k + 3];
    expect([r.count(Op.clear), quiet.opCount]).toEqual([1, 0]);
    expect([r.args[at], r.args[at + 2]]).toEqual([x(40), x(60) - x(40)]);
  });

  it("a real tile arrival redraws the rows within its columns, and nothing in another section", async () => {
    const uri = "https://example.org/lazy";
    const b = createSeqViewSpec();
    const source = { uri, protocol: "seqstar-tiles/1" };
    const aln = b.data.alignment({ ref: "lazy", row_count: 400, column_count: 2000, source } as never);
    const g = b.view({ ref: "g", label: "G", style: "bare", height: "box" });
    g.section({ ref: "g.rows" }).block({ alignment: aln }).track({ label: "Rows" }).rows();
    const s = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(50) });
    g.section({ ref: "g.seq" }).block({ sequence: s }).track({ label: "L" }).letters();
    const t = await composed(b.getDocument(), { section: "g.rows", loaders: lazyLoaders(uri, 400, 2000) });
    // No margin: with one the tiles past the view arrive before the pan, which is what a margin is for (ST2).
    t.ctx.config.set(SeqViewerConfig.Stream.Margin, 0);
    t.vp?.setView(4, 0, "api");
    const other = new SectionComposer(t.ctx, "g.seq");
    const [be, quiet] = [recording(COPY_BACK), recording(COPY_BACK)];
    t.frame(be);
    other.frame(quiet, Dirty.tracks);
    expect(t.frame(be).opCount).toBe(0);
    const arrive = async () => {
      (t.ctx.streams as Streams).pass();
      await new Promise((ok) => setTimeout(ok, 0));
      await new Promise((ok) => setTimeout(ok, 0));
    };
    await arrive();
    // The first arrival holds the row labels' chunk: the band is redrawn whole, its labels with it.
    const labelled = t.frame(be);
    const l = t.layout(t.section);
    const clearsFrom = (r: typeof labelled) =>
      Array.from({ length: r.opCount }, (_, k) => k)
        .filter((k) => r.ops[4 * k] === Op.clear)
        .map((k) => r.args[r.ops[4 * k + 3]]);
    expect(Math.min(...clearsFrom(labelled))).toBe(0);
    // Past the first tile's columns, only tiles arrive: the rows within their columns.
    t.vp?.setView(4, 600, "api");
    t.frame(be);
    other.frame(quiet, Dirty.tracks);
    await arrive();
    const r = t.frame(be);
    quiet.reset();
    other.frame(quiet, Dirty.tracks);
    expect([r.count(Op.clear) > 0, r.count(Op.copy), quiet.opCount]).toEqual([true, 0, 0]);
    // Each clear within the data area and the arrived tiles' columns (256 a tile at the finest level).
    const ax = t.ctx.layout.axis(t.section);
    for (let k = 0; k < r.opCount; k++) {
      if (r.ops[4 * k] !== Op.clear) continue;
      const at = r.ops[4 * k + 3];
      expect(r.args[at]).toBeGreaterThanOrEqual(
        Math.max(Math.round(l.dataLeft), Math.round(ax?.x(0, 256 as never) ?? 0)),
      );
      expect(r.args[at] + r.args[at + 2]).toBeLessThanOrEqual(Math.round(l.dataRight));
      expect(r.args[at] + r.args[at + 2]).toBeLessThanOrEqual(Math.round(ax?.x(0, 512 as never) ?? 0));
    }
  });

  it("the ruler labels a copy past the first: di-ubiquitin's second Ub 45 is “45 (copy 2)”", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "ub", value: "MQIFVKTLTG".repeat(8).slice(0, 76) });
    b.data.sequence({ ref: "diub", value: "MQIFVKTLTG".repeat(16).slice(0, 152) });
    const segments = [
      [0, 0, 76],
      [76, 0, 76],
    ];
    b.data.mapping({ ref: "diub-ub", from: "diub", to: "ub", method: "segments", segments } as never);
    const v = b.view({ ref: "v", label: "V" });
    v.block({ sequence: "diub" });
    v.track({ label: "" }).ruler({ of: "ub" } as never);
    const t = await composed(b.getDocument(), { width: 1600 });
    t.vp?.setView(10, 0, "api");
    const texts = textsOf(t.frame());
    expect(texts).toContain("40");
    // The second copy's labels are wider: thinned to every second tick of 5.
    expect(texts).toContain("45 (copy 2)");
    expect(texts.filter((s) => s.includes("copy 1"))).toEqual([]);
  });
});
