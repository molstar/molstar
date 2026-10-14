/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot (seqspec:format#layer-types `dot_plot`, seqspec:interaction#hover, seqspec:interaction#selection): its
// grid drawn as tiles over the plot and its pairs as segments, y showing the section's window down the track; a hover
// on a cell lights both positions and labels the partner and the value; a drag selects a rectangle as two ranges; the
// overlay draws lit runs as rows; the SVG export draws what the screen does.

import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { SkinToken } from "../src/draw/backend.ts";
import { Dirty, Pass } from "../src/draw/index.ts";
import { Hit, KernelRenderers, type LayerRenderer, MarkBit, type SeqLociAny } from "../src/index.ts";
import { renderSvg } from "../src/svg/index.ts";
import { Op } from "../src/testing/index.ts";
import { recording } from "./composer-fixtures.ts";
import { directed } from "./conformance-fixtures.ts";
import { dotPlotDoc, SEQ } from "./dot-plot-fixtures.ts";
import { driven } from "./input-fixtures.ts";
import { lit, passRects } from "./overlay-fixtures.ts";
import { all, parseXml } from "./svg/xml.ts";

const SECTION = "v.s";
const N = SEQ.length;

/** The plot entry's box (surface CSS px) and the data area's edges. */
function plotBox(ctx: Awaited<ReturnType<typeof driven>>["ctx"]) {
  const l = ctx.layout.get(SECTION)?.current;
  const v = ctx.layout.axis(SECTION);
  const a = ctx.graph.assembled.value;
  if (!l || !v || !a) throw new Error("no layout");
  const track = a.view?.resolved.tracks.ref.findIndex((r) => a.set.strings.get(r) === "t-plot") ?? -1;
  const e = Array.from({ length: l.entryCount }, (_, k) => k).find((k) => l.entryTrack[k] === track) ?? -1;
  const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - v.scrollY;
  return { top, height: l.entryHeight[e], left: v.left, right: v.right, v };
}
/** The surface point of cell (i, j): i's center across, j's down. */
function cellAt(t: Awaited<ReturnType<typeof driven>>, i: number, j: number): [number, number] {
  const b = plotBox(t.ctx);
  const x = t.xOf(SECTION, 0, i);
  const y = b.top + ((t.xOf(SECTION, 0, j) - b.left) / (b.right - b.left)) * b.height;
  return [x, y];
}
function indices(loci: SeqLociAny | null): number[] {
  const out: number[] = [];
  const els = (loci as unknown as { elements: { indices: OrderedSet }[] } | null)?.elements ?? [];
  for (const e of els) OrderedSet.forEach(e.indices, (x) => out.push(x), null);
  return out.sort((a, b) => a - b);
}

describe("the dot plot's draw", () => {
  it("draws its grid as tiles over the plot, y the section's window down the track, and its pairs as segments", async () => {
    const t = await lit(dotPlotDoc(), { width: 900, height: 600 });
    const r = t.frame(recording(), Dirty.all);
    const images: number[][] = [];
    for (let k = 0; k < r.opCount; k++)
      if (r.ops[4 * k] === Op.image) images.push(Array.from(r.args.subarray(r.ops[4 * k + 3], r.ops[4 * k + 3] + 8)));
    // One tile at level 0: every shown position a sample, the whole window over the plot's box.
    expect(images.length).toBe(1);
    const [sx, sy, sw, sh, dx, dy, dw, dh] = images[0];
    expect([sx, sy, sw, sh]).toEqual([0, 0, N, N]);
    const l = t.ctx.layout.get(t.section)?.current;
    const v = t.ctx.layout.axis(t.section);
    expect(Math.abs(dx - (v?.left ?? 0))).toBeLessThanOrEqual(1);
    expect(dw).toBeGreaterThan(0);
    expect(dh).toBeGreaterThan(0);
    expect(dy).toBeGreaterThanOrEqual(l?.scrolledTop ?? 0);
    // The repeat's matches at offset 30, both ways: segments (paths), none on the diagonal.
    expect(r.pathItems.n).toBeGreaterThanOrEqual(2);
  });
});

describe("the dot plot's interaction", () => {
  const behaviors = { renderers: KernelRenderers as LayerRenderer[] };

  it("a hover on a cell names its partner under y, lights both and labels the partner and the value", async () => {
    const t = await driven(dotPlotDoc(), behaviors);
    const [x, y] = cellAt(t, 5, 35);
    const hit = Hit.create();
    expect(t.ctx.layout.hitTest(SECTION, x, y, hit)).toBe(true);
    expect([hit.index, hit.pairIndex]).toEqual([5, 35]);
    t.pointer("move", SECTION, x, y);
    const h = t.ctx.hover.live.current;
    expect([h.index, h.pairIndex]).toEqual([5, 35]);
    expect(indices(t.ctx.hover.toLoci())).toEqual([5, 35]);
    const label = t.ctx.labels.current.value;
    expect(label?.lines.map((l) => [l.key, l.value])).toEqual(
      expect.arrayContaining([
        ["with", expect.stringContaining(`${SEQ[35]}36`)],
        ["Distance", "30"],
      ]),
    );
    // Another row of the same column is another element.
    const [, y2] = cellAt(t, 5, 36);
    t.pointer("move", SECTION, x, y2);
    expect(t.ctx.hover.live.current.pairIndex).toBe(36);
  });

  it("a drag selects a rectangle as two ranges, held as the pending selection while it lasts", async () => {
    const t = await driven(dotPlotDoc(), behaviors);
    const [x0, y0] = cellAt(t, 2, 30);
    const [x1, y1] = cellAt(t, 6, 34);
    t.pointer("down", SECTION, x0, y0);
    t.pointer("move", SECTION, (x0 + x1) / 2, (y0 + y1) / 2);
    t.pointer("move", SECTION, x1, y1);
    const pd = t.ctx.selection.pending.current;
    expect([pd.active, pd.startShown, pd.endShown, pd.yStartShown, pd.yEndShown]).toEqual([true, 2, 6, 30, 34]);
    t.pointer("up", SECTION, x1, y1);
    expect(t.ctx.selection.pending.current.active).toBe(false);
    const sel = t.ctx.selection.current.value as unknown as { seq: SeqLociAny };
    expect(indices(sel.seq)).toEqual([2, 3, 4, 5, 6, 30, 31, 32, 33, 34]);
    // Shift adds a rectangle; the view never pans.
    const before = t.vp(SECTION).live.current.scrollX;
    const [x2, y2] = cellAt(t, 10, 12);
    const [x3, y3] = cellAt(t, 11, 13);
    t.drag(SECTION, x2, y2, x3, y3, { shift: true });
    expect(indices((t.ctx.selection.current.value as unknown as { seq: SeqLociAny }).seq)).toEqual([
      2, 3, 4, 5, 6, 10, 11, 12, 13, 30, 31, 32, 33, 34,
    ]);
    expect(t.vp(SECTION).live.current.scrollX).toBe(before);
  });

  it("the cursor over a plot's data is a crosshair, kept for the drag; a track's data elsewhere grabs", async () => {
    const t = await driven(dotPlotDoc(), behaviors);
    const [x, y] = cellAt(t, 5, 35);
    t.pointer("move", SECTION, x, y);
    expect(t.ctx.ui.cursor(SECTION).value).toBe("crosshair");
    t.pointer("down", SECTION, x, y);
    t.pointer("move", SECTION, x + 40, y + 40);
    expect(t.ctx.ui.cursor(SECTION).value).toBe("crosshair");
    t.pointer("up", SECTION, x + 40, y + 40);
    t.pointer("move", SECTION, x, plotBox(t.ctx).top - 2);
    expect(t.ctx.ui.cursor(SECTION).value).not.toBe("crosshair");
  });

  it("a click on a cell selects both positions", async () => {
    const t = await driven(dotPlotDoc(), behaviors);
    const [x, y] = cellAt(t, 8, 20);
    t.click(SECTION, x, y);
    t.advance(1000);
    expect(indices((t.ctx.selection.current.value as unknown as { seq: SeqLociAny }).seq)).toEqual([8, 20]);
  });
});

describe("the dot plot's overlay", () => {
  it("draws every lit run as a row across the plot, the hover's crosshair and the hovered cell", async () => {
    const t = await lit(dotPlotDoc(), { width: 900, height: 600 });
    const v = t.ctx.layout.axis(t.section);
    const l = t.ctx.layout.get(t.section)?.current;
    if (!v || !l) throw new Error("no layout");
    t.frame(recording(), Dirty.all);
    t.ctx.selection.set(
      { kind: "seq-loci", set: t.set, elements: [{ unit: 0, indices: OrderedSet.ofRange(30, 34) }] } as never,
      "api",
    );
    const rects = passRects(t.overlay(), Pass.overlay);
    const left = Math.floor(v.left);
    const width = Math.ceil(v.right) - left;
    // The selection's column (the composer's) and its row across the plot (the dot plot's), with its edges.
    const rows = rects.filter((r) => r[0] === left && r[2] === width && r[4] === SkinToken.selection);
    expect(rows.length).toBe(1);
    expect(rects.filter((r) => r[4] === SkinToken.selectionBorder).length).toBeGreaterThanOrEqual(4);
    expect(MarkBit.selection).toBe(1);
  });
});

describe("the dot plot's SVG export", () => {
  it("draws the grid as an image and the pairs as paths, as the screen does", async () => {
    const xml = await renderSvg(dotPlotDoc(), { width: 900 });
    const doc = parseXml(xml);
    expect(all(doc, "image").length).toBe(1);
    expect(all(doc, "path").length).toBeGreaterThanOrEqual(2);
  });
});

describe("the dot plot on a reversed piece", () => {
  it("mirrors both ways: a cell is hit where its positions are shown, across and down", async () => {
    const t = await driven(directed(dotPlotDoc(), -1), { renderers: KernelRenderers as LayerRenderer[] });
    const [x, y] = cellAt(t, 5, 35);
    const hit = Hit.create();
    expect(t.ctx.layout.hitTest(SECTION, x, y, hit)).toBe(true);
    expect([hit.index, hit.pairIndex]).toEqual([5, 35]);
    const b = plotBox(t.ctx);
    // Position 35 lies left of 5 when reversed, so its row is above 5's.
    expect(y).toBeLessThan(b.top + ((x - b.left) / (b.right - b.left)) * b.height);
  });
});
