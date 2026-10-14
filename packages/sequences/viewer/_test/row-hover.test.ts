/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `hoverRow` behavior on the recording backend (seqspec:interaction#hover): a hover on a track's data draws its
// row band on the overlay at half a header hover's strength, a half over the header and a quarter over the data (the track's box, a split row's, a matrix's
// hovered band alone, none on the ruler), under the hover's own marks; `band: "off"`, the strengths, the behavior
// removed; on a stacked axis past a block's end, its lane's; none for a collapsed clade's wedge row (no tree hover);
// none in an SVG export that does not draw the hover. And it is a drawing only: the hover's kind, header description, loci, size and version are what they are
// without it, so no reader (a 3D link, a host) takes a data hover for a header's.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loci } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Pass, SkinToken } from "../src/draw/index.ts";
import { HoverKind, Inputs, lociSize, PointerType, type RowHoverParams } from "../src/index.ts";
import { renderSvg } from "../src/svg/index.ts";
import { FAMILY_SECTION, familyWithTree } from "./families.ts";
import { lit, passRects, twoBlocks } from "./overlay-fixtures.ts";

/** Letters, a ruler, ranges, values, a split track and a three-state matrix on one chain. */
function doc() {
  const b = createSeqViewSpec({ title: "Row hover" });
  const n = 40;
  const a = b.data.sequence({ ref: "a", value: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQGNPEEAL" });
  const ann = (ref: string, geometry: string, value: unknown, more: object = {}) =>
    b.data.annotation({ ref, geometry, target: a, value, ...more } as never);
  const rg = ann("rg", "ranges", { start: [3, 12], stop: [8, 20], label: ["p", "q"], kind: ["x", "y"] });
  const index = Array.from({ length: n }, (_, i) => i);
  const vv = ann("vv", "values", { index, value: index.map((i) => 0.2 + ((i * 31) % 80) / 100) });
  const mx = ann(
    "mx",
    "matrix",
    { index, value: index.map((i) => [(i % 7) / 7, (i % 11) / 11, (i % 3) / 3]) },
    { states: ["x", "y", "z"] },
  );
  const v = b.view({ ref: "v", label: "V", height: "box" } as never);
  v.block({ sequence: a });
  v.track({ ref: "t.ruler", label: "" }).ruler();
  v.track({ ref: "t.letters", label: "Letters" }).letters();
  v.track({ ref: "t.ranges", label: "Ranges" }).ranges({ data: rg });
  v.track({ ref: "t.values", label: "Values" }).values({ data: vv });
  v.track({ ref: "t.split", label: "Split", expanded: true } as never).ranges({ data: rg, split: "kind" } as never);
  v.track({ ref: "t.matrix", label: "Matrix" }).matrix({ data: mx });
  return b.getDocument();
}
const open = () => lit(doc(), { width: 800, height: 600 });
type T = Awaited<ReturnType<typeof open>>;

/** The band's rects on the overlay: [y, h, alpha] of its header part (x 0) and of its data part (x dataLeft). */
function band(t: T) {
  const o = passRects(t.overlay(), Pass.overlay).filter((r) => r[4] === SkinToken.hover);
  const left = Math.round(t.layout(t.section).dataLeft);
  const head = o.filter((r) => r[0] === 0 && r[2] === left).map((r) => [r[1], r[3], r[5]]);
  const data = o.filter((r) => r[0] === left && r[3] === head[0]?.[1]).map((r) => [r[1], r[3], r[5]]);
  return { head, data };
}
/** Entry e's visible box of track `ref` (its `row`th split row) as [y, h], device px at DPR 1. */
function boxOf(t: T, ref: string, row = -1) {
  const l = t.layout(t.section);
  const k = t.trackIndex(ref);
  for (let e = 0; e < l.entryCount; e++)
    if (l.entryTrack[e] === k && (row < 0 || l.entryRow[e] === row)) {
      const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e];
      return [Math.round(top), Math.round(top + l.entryHeight[e]) - Math.round(top)];
    }
  throw new Error(`no entry for ${ref}`);
}

describe("a hover on a track's data draws its row band (the hoverRow behavior)", () => {
  it("the track's box: half strength over the header, a quarter over the data, on an item and off one", async () => {
    const t = await open();
    for (const i of [5, 25]) {
      t.hover(t.xOf(t.section, 0, i), t.yOf(t.section, "t.ranges"));
      const h = t.ctx.hover.live.current;
      expect(h.kind).toBe(i === 5 ? HoverKind.item : HoverKind.seq);
      const [y, hh] = boxOf(t, "t.ranges");
      expect(band(t)).toEqual({ head: [[y, hh, 128]], data: [[y, hh, 64]] });
    }
    // Values too; a split row's band is its row's.
    t.hover(t.xOf(t.section, 0, 10), t.yOf(t.section, "t.values"));
    expect(band(t).head).toEqual([[...boxOf(t, "t.values"), 128]]);
    const l = t.layout(t.section);
    const k = t.trackIndex("t.split");
    let row = -1;
    for (let e = 0; e < l.entryCount; e++) if (l.entryTrack[e] === k && l.entryRow[e] >= 0) row = l.entryRow[e];
    expect(row).toBeGreaterThanOrEqual(0);
    t.hover(t.xOf(t.section, 0, 30), t.yOf(t.section, "t.split", row));
    expect(t.ctx.hover.live.current.split).toBe(row);
    expect(band(t).head).toEqual([[...boxOf(t, "t.split", row), 128]]);
  });

  it("a matrix lights its hovered band alone (its renderer's rowHover); the ruler none", async () => {
    const t = await open();
    const [y, hh] = boxOf(t, "t.matrix");
    const states: number[] = [];
    for (let dy = 2; dy < hh; dy += 6) {
      t.hover(t.xOf(t.section, 0, 6), y + dy);
      const s = t.ctx.hover.live.current.state;
      if (s < 0 || states.includes(s)) continue;
      states.push(s);
      const { head, data } = band(t);
      expect(head.length).toBe(1);
      expect(data).toEqual([[head[0][0], head[0][1], 64]]);
      // The band's rows hold the point, within the matrix's box and short of it.
      expect(head[0][0]).toBeLessThanOrEqual(y + dy);
      expect(head[0][0] + head[0][1]).toBeGreaterThan(y + dy);
      expect(head[0][0]).toBeGreaterThanOrEqual(y);
      expect(head[0][1]).toBeLessThan(hh);
    }
    expect(states.length).toBe(3);
    t.hover(t.xOf(t.section, 0, 6), t.yOf(t.section, "t.ruler"));
    expect(t.ctx.hover.live.current.track).toBeGreaterThanOrEqual(0);
    expect(band(t).head).toEqual([]);
  });

  it("`band: off`, the strengths, and the behavior removed; a header's hover keeps its own band", async () => {
    const t = await open();
    const at = () => t.hover(t.xOf(t.section, 0, 25), t.yOf(t.section, "t.values"));
    const p = t.ctx.behaviors.get<RowHoverParams>("hoverRow");
    p?.setParams({ band: "off" });
    at();
    expect(band(t).head).toEqual([]);
    p?.setParams({ band: "entry", header: 1, data: 0.75 });
    const [y, hh] = boxOf(t, "t.values");
    expect(band(t)).toEqual({ head: [[y, hh, 255]], data: [[y, hh, 191]] });
    p?.remove();
    expect(band(t).head).toEqual([]);
    // The header's hover is the headerHover behavior's: its band at full and half strength, as before.
    t.hover(10, t.yOf(t.section, "t.values"));
    expect(t.ctx.hover.live.current.kind).toBe(HoverKind.header);
    expect(band(t)).toEqual({ head: [[y, hh, 255]], data: [[y, hh, 128]] });
  });
});

describe("the band is a drawing only: a data hover is never a header's", () => {
  /** What a reader sees of the hover: kind, header description, loci (without its set), size, version bumps. */
  const seen = async (band: boolean) => {
    const t = await open();
    if (!band) t.ctx.behaviors.get("hoverRow")?.remove();
    const loci = (l: Loci | null) =>
      JSON.stringify(l, (k, v) => (k === "set" ? undefined : ArrayBuffer.isView(v) ? Array.from(v as never) : v));
    const out: unknown[] = [];
    for (const [track, i] of [
      ["t.ranges", 5],
      ["t.ranges", 25],
      ["t.values", 10],
      ["t.letters", 3],
      ["t.matrix", 6],
      ["t.ruler", 8],
    ] as const) {
      const v0 = t.ctx.hover.live.version.value;
      t.hover(t.xOf(t.section, 0, i), t.yOf(t.section, track));
      t.overlay();
      const h = t.ctx.hover.live.current;
      const l = t.ctx.hover.toLoci() as Loci | null;
      expect(h.kind, track).not.toBe(HoverKind.header);
      expect(h.header.on, track).toBe(false);
      // The element alone: an item, never the track's items (two on the ranges track, as its header gives).
      expect(lociSize(l)?.items ?? 0, track).toBeLessThanOrEqual(1);
      out.push([track, i, h.kind, h.region, h.track, h.header.on, t.ctx.hover.live.version.value - v0]);
      out.push([loci(l), lociSize(l), t.ctx.hover.size()]);
    }
    return out;
  };

  it("kind, loci, size and version bumps are the same with the band as without it", async () => {
    expect(await seen(true)).toEqual(await seen(false));
    // The control: the ranges track's header gives both its items.
    const t = await open();
    t.hover(10, t.yOf(t.section, "t.ranges"));
    expect([t.ctx.hover.live.current.kind, lociSize(t.ctx.hover.toLoci() as Loci | null)?.items]).toEqual([
      HoverKind.header,
      2,
    ]);
  });
});

describe("where the band is not, and where it still is", () => {
  it("a collapsed clade's wedge row: its data hover draws no band and no tree hover; its header does", async () => {
    const t = await lit(familyWithTree(40, 10), { width: 1000, height: 400, section: FAMILY_SECTION });
    t.ctx.ui.toggleClade("family", 2);
    await t.ctx.graph.update();
    const l = t.layout(t.section);
    const top = l.scrolledTop + l.entryTop[l.rows.entry];
    let wedge = -1;
    for (let d = 0; d < l.rows.count && wedge < 0; d++) if (l.rows.node?.[d] === 2) wedge = d;
    expect(wedge).toBeGreaterThanOrEqual(0);
    const y = top + (wedge + 0.5) * l.rows.height;
    t.hover(l.dataLeft + 40, y);
    const h = t.ctx.hover.live.current;
    expect([h.kind, h.lane, h.treeNode]).toEqual([HoverKind.none, -1, 2]);
    const o = passRects(t.overlay(), Pass.overlay).filter((r) => r[4] === SkinToken.hover && r[0] === 0);
    expect(o).toEqual([]);
    // The control: the wedge's header hovers its node, which lights the clade's row.
    t.hover(l.headerWidth - 4, y);
    expect(t.ctx.hover.live.current.kind).toBe(HoverKind.header);
    expect(passRects(t.overlay(), Pass.overlay).some((r) => r[4] === SkinToken.hover && r[0] === 0)).toBe(true);
  });

  it("stacked, past a block's end (no block hovered): the band of the lane under the pointer", async () => {
    const t = await lit(twoBlocks({ stacked: true }), { width: 800, height: 400 });
    const l = t.layout(t.section);
    // Block c (90 residues) is the shorter lane: a point past its end, on its letters.
    let e = 0;
    while (l.entryBlock[e] !== 1 || l.entryTrack[e] < 0) e++;
    const p = Inputs.pointer();
    p.section = t.ctx.layout.sectionIndex(t.section);
    p.type = PointerType.move;
    p.x = t.xOf(t.section, 1, 89) + 40;
    p.y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
    t.ctx.input.pointer(p);
    const h = t.ctx.hover.live.current;
    expect([h.block, h.track >= 0]).toEqual([-1, true]);
    const top = Math.round(l.scrolledTop + l.entryTop[e]);
    const lit0 = passRects(t.overlay(), Pass.overlay).filter((r) => r[4] === SkinToken.hover && r[0] === 0);
    expect(lit0.map((r) => r[1])).toEqual([top]);
  });

  it("an SVG export with the overlay and not the hover draws no band; with the hover, it does", async () => {
    const t = await open();
    const opts = { width: 800, overlay: true };
    const before = await renderSvg(t.ctx, opts);
    t.hover(t.xOf(t.section, 0, 25), t.yOf(t.section, "t.values"));
    expect(band(t).head.length).toBe(1);
    expect(await renderSvg(t.ctx, opts)).toBe(before);
    expect(await renderSvg(t.ctx, { ...opts, hover: true })).not.toBe(before);
  });
});
