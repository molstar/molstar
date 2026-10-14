/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's edges (seqspec:format#layer-types `dot_plot`): grids sized to what their rows place, the cell cap and
// its report, the pyramid's means against the cells' at levels past one (sizes no power of two, missing values), two
// grids on one block pair, a plot after another layer of its track, stacked lanes, a drag across blocks and a cancelled
// one, tiles rebuilt under a new theme and built a few a frame on a progressive backend, and a pair of single
// positions drawn as a dot at every zoom.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { OrderedSet } from "@molstar/seq-core/model";
import { afterEach, describe, expect, it } from "vitest";
import { SkinToken } from "../src/draw/backend.ts";
import { DotPlotLimits, type DotPlotState, DotPlotVisits, dotPlotSample, dotPlotValue } from "../src/draw/dot-plot.ts";
import { BackendCaps, Dirty, Pass } from "../src/draw/index.ts";
import { Hit, KernelRenderers, type LayerRenderer } from "../src/index.ts";
import { Op, type RecordingBackend } from "../src/testing/index.ts";
import { recording } from "./composer-fixtures.ts";
import { dotPlotDoc } from "./dot-plot-fixtures.ts";
import { driven } from "./input-fixtures.ts";
import { lit, passRects } from "./overlay-fixtures.ts";

const renderers = { renderers: KernelRenderers as LayerRenderer[] };
const cells = DotPlotLimits.cells;
afterEach(() => {
  DotPlotLimits.cells = cells;
});

/** Two sequences side by side (a, b) and grids over them, a plot track whose first layer is `first` (dot_plot, ranges). */
function twoBlocks(o: {
  na?: number;
  nb?: number;
  grids?: Record<string, unknown>[];
  first?: "ranges";
  stacked?: boolean;
  range?: [number, number];
}) {
  const [na, nb] = [o.na ?? 37, o.nb ?? 53];
  const b = createSeqViewSpec({ title: "two" });
  const a = b.data.sequence({ ref: "a", value: "M".repeat(na) });
  const q = b.data.sequence({ ref: "b", value: "K".repeat(nb) });
  const r = b.data.annotation({ ref: "r", geometry: "ranges", target: a, value: { start: [1], stop: [3] } });
  const refs = (o.grids ?? []).map((value, k) =>
    b.data.annotation({ ref: `g${k}`, label: `G${k}`, geometry: "grid", value }),
  );
  const sec = b.view({ ref: "v", label: "V" }).section({ ref: "v.s", arrangement: o.stacked ? "stacked" : "side" });
  sec.block({ sequence: a, ...(o.range ? { range: o.range } : {}) });
  sec.block({ sequence: q });
  const t = sec.track({ ref: "t-plot", label: "Plot" });
  if (o.first === "ranges") t.ranges({ data: r });
  t.dot_plot({ ref: "l-plot", data: refs.length ? refs : "r" });
  return b.getDocument();
}
function prepared(t: { ctx: { graph: { assembled: { value: unknown } } } }, layer = "l-plot"): DotPlotState {
  const m = (
    t.ctx.graph.assembled.value as {
      view: {
        resolved: { layers: { ref: Int32Array } };
        set: { strings: { get(i: number): string } };
        prepared: { state: unknown }[];
      };
    }
  ).view;
  const k = Array.from(m.resolved.layers.ref).findIndex((x) => m.set.strings.get(x) === layer);
  return m.prepared[k].state as DotPlotState;
}
function images(r: RecordingBackend) {
  let n = 0;
  for (let k = 0; k < r.opCount; k++) if (r.ops[4 * k] === Op.image) n++;
  return n;
}

describe("grids", () => {
  it("are sized to what their rows place, within the blocks' pieces, float32; a row that is the whole grid is kept as given", async () => {
    // a[5, 15) × b[10, 30) on 37 × 53 residues; a's block shows [8, 30) only.
    const value = Array.from({ length: 10 * 20 }, (_, x) => x);
    const t = await lit(
      twoBlocks({
        grids: [
          { target: ["a"], start: [5], stop: [15], b_target: ["b"], b_start: [10], b_stop: [30], value: [value] },
        ],
        range: [8, 30],
      }),
    );
    const s = prepared(t);
    expect(s.grids.map((g) => [g.bx, g.by, g.x0, g.nx, g.y0, g.ny])).toEqual([[0, 1, 8, 7, 10, 20]]);
    expect(s.grids[0].values).toBeInstanceOf(Float32Array);
    // a 9, b 12: the row's (9 − 5) × 20 + (12 − 10).
    expect(dotPlotValue(s, 0, 9, 1, 12)).toBe(4 * 20 + 2);
  });

  it("past the cap a block pair is left out and reported, with its cells, what it allocates and the limit", async () => {
    DotPlotLimits.cells = 1000;
    const t = await lit(twoBlocks({ grids: [{ target: ["a"], b_target: ["b"], value: [Array(37 * 53).fill(1)] }] }));
    expect(prepared(t).grids).toEqual([]);
    const issues = t.ctx.issues.current.value.filter((i) => i.code === "svs.render.omitted");
    // A copy (1,961) and its pyramid (19 × 27 + 10 × 14 + 5 × 7 + 3 × 4 + 2 × 2 + 1 = 705).
    expect(issues.map((i) => i.message)).toEqual([
      expect.stringContaining(
        "is 37 × 53 = 1961 cells, and a copy of its cells and its mean pyramid would take 2666 floats, past the 1000",
      ),
    ]);
  });

  it("a Float32Array row that is the whole grid is referenced: only its pyramid counts against the cap", async () => {
    DotPlotLimits.cells = 1000;
    const row = new Float32Array(37 * 53).fill(2);
    const t = await lit(twoBlocks({ grids: [{ target: ["a"], b_target: ["b"], value: [row] }] }));
    const s = prepared(t);
    expect(s.grids.length).toBe(1);
    expect(s.grids[0].values).toBe(row);
    expect(t.ctx.issues.current.value.filter((i) => i.code === "svs.render.omitted")).toEqual([]);
  });

  it("two grids on one block pair are one grid, the later row over the earlier", async () => {
    const all = { target: ["a"], b_target: ["b"], value: [Array(37 * 53).fill(1)] };
    const part = {
      target: ["a"],
      start: [0],
      stop: [2],
      b_target: ["b"],
      b_start: [0],
      b_stop: [2],
      value: [[7, 7, 7, 7]],
    };
    const s = prepared(await lit(twoBlocks({ grids: [all, part] })));
    expect(s.grids.length).toBe(1);
    expect([dotPlotValue(s, 0, 1, 1, 1), dotPlotValue(s, 0, 5, 1, 5)]).toEqual([7, 1]);
  });

  it("the pyramid's sample is the mean of its shown cells, missing ones left out, at levels 1 to 5 on sizes no power of two", async () => {
    const rnd = (k: number) => ((k * 2654435761) % 1000) / 10;
    const value = Array.from({ length: 37 * 53 }, (_, k) => (k % 7 === 3 ? null : rnd(k)));
    const t = await lit(twoBlocks({ grids: [{ target: ["a"], b_target: ["b"], value: [value] }] }));
    const s = prepared(t);
    const ax = t.ctx.layout.axis(t.section)?.axis;
    if (!ax) throw new Error("no axis");
    const g = s.grids[0];
    expect(g.levels?.length).toBeGreaterThan(4);
    for (let L = 1; L <= 6; L++) {
      const step = 2 ** L;
      for (let sx = 0; sx * step < 37; sx++)
        for (let sy = 0; sy * step < 53; sy++) {
          let [sum, n] = [0, 0];
          for (let i = sx * step; i < Math.min(37, (sx + 1) * step); i++)
            for (let j = sy * step; j < Math.min(53, (sy + 1) * step); j++) {
              const v = value[i * 53 + j];
              if (v !== null) [sum, n] = [sum + Math.fround(v), n + 1];
            }
          const want = n > 0 ? sum / n : Number.NaN;
          const got = dotPlotSample(ax, g, L, sx, sy);
          if (Number.isNaN(want)) expect(got).toBeNaN();
          else expect(Math.abs(got - want)).toBeLessThan(1e-3);
        }
    }
  });
});

describe("the overlay and the hit", () => {
  it("a plot after another layer of its track still draws its rows (the track's first dot plot)", async () => {
    const value = Array(37 * 53).fill(1);
    const t = await lit(twoBlocks({ first: "ranges", grids: [{ target: ["a"], b_target: ["b"], value: [value] }] }), {
      width: 900,
      height: 600,
    });
    expect(prepared(t).first).toBe(true);
    t.frame(recording(), Dirty.all);
    t.ctx.selection.set(
      { kind: "seq-loci", set: t.set, elements: [{ unit: 1, indices: OrderedSet.ofRange(10, 14) }] } as never,
      "api",
    );
    const v = t.ctx.layout.axis(t.section);
    const left = Math.floor(v?.left ?? 0);
    const rows = passRects(t.overlay(), Pass.overlay).filter((r) => r[0] === left && r[4] === SkinToken.selection);
    expect(rows.length).toBe(1);
  });

  it("on a stacked section each lane plots its own block, y within it; the hit names the partner's block", async () => {
    const value = Array(37 * 37).fill(2);
    const t = await driven(
      twoBlocks({ stacked: true, na: 37, nb: 37, grids: [{ target: ["b"], b_target: ["b"], value: [value] }] }),
      { ...renderers, size: { "v.s": [1000, 1600] } },
    );
    const l = t.layout("v.s");
    const v = t.ctx.layout.axis("v.s");
    if (!v) throw new Error("no axis");
    // The second lane's plot entry.
    const e =
      Array.from({ length: l.entryCount }, (_, k) => k)
        .filter((k) => l.entryBlock[k] === 1)
        .pop() ?? -1;
    const top = l.scrolledTop + l.entryTop[e] - v.scrollY;
    const x = t.xOf("v.s", 1, 5);
    const y = top + ((t.xOf("v.s", 1, 20) - v.left) / (v.right - v.left)) * l.entryHeight[e];
    const hit = Hit.create();
    t.ctx.layout.hitTest("v.s", x, y, hit);
    expect([hit.block, hit.index, hit.pairBlock, hit.pairIndex]).toEqual([1, 5, 1, 20]);
  });

  it("a drag across into the next block stops at its own block's end; a cancelled drag selects nothing", async () => {
    const value = Array(37 * 53).fill(1);
    const t = await driven(twoBlocks({ grids: [{ target: ["a"], b_target: ["b"], value: [value] }] }), renderers);
    const l = t.layout("v.s");
    const v = t.ctx.layout.axis("v.s");
    if (!v) throw new Error("no axis");
    const e = l.entryCount - 1;
    const top = l.scrolledTop + l.entryTop[e] - v.scrollY;
    const yOf = (k: number) => top + ((t.xOf("v.s", 1, k) - v.left) / (v.right - v.left)) * l.entryHeight[e];
    t.pointer("down", "v.s", t.xOf("v.s", 0, 30), yOf(10));
    t.pointer("move", "v.s", t.xOf("v.s", 1, 4), yOf(12));
    const pd = t.ctx.selection.pending.current;
    expect([pd.block, pd.startShown, pd.endShown, pd.yBlock, pd.yEndShown]).toEqual([0, 30, 36, 1, 12]);
    t.pointer("cancel", "v.s", t.xOf("v.s", 1, 4), yOf(12));
    expect(t.ctx.selection.pending.current.active).toBe(false);
    expect((t.ctx.selection.current.value as unknown as { seq: unknown }).seq ?? null).toBe(null);
  });
});

describe("tiles", () => {
  it("are built again under a new theme", async () => {
    const t = await lit(dotPlotDoc(), { width: 900, height: 600 });
    t.frame(recording(), Dirty.all);
    const before = DotPlotVisits.tiles;
    t.frame(recording(), Dirty.all);
    expect(DotPlotVisits.tiles).toBe(before);
    await t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "l-dist", params: { color: { theme: "uniform", params: { color: "#336699" } } } }],
    });
    t.frame(recording(), Dirty.all);
    expect(DotPlotVisits.tiles).toBeGreaterThan(before);
  });

  it("a progressive backend builds a few a frame, the coarser level meanwhile, and asks for another frame", async () => {
    const seq = "MKTAYIAKQR".repeat(240);
    const t = await lit(dotPlotDoc(seq, true), { width: 900, height: 600, dpr: 2 });
    // The fit: 2,400 residues over about 800 px at DPR 2, level 1, some 25 tiles.
    const be = recording(BackendCaps.Progressive);
    t.frame(be, Dirty.all);
    expect(t.composer.again).toBe(true);
    expect(images(be)).toBeGreaterThan(0);
    // Frames until every tile in view is built: then no more asked.
    let frames = 1;
    while (t.composer.again && frames < 50) {
      t.composer.redraw();
      t.frame(be, Dirty.all);
      frames++;
    }
    expect([t.composer.again, frames > 1]).toEqual([false, true]);
  });
});

describe("pairs", () => {
  it("a pair of single positions is a dot at every zoom: a circle in its cell (about 20 × 8 px), never a stroke", async () => {
    const b = createSeqViewSpec({ title: "dots" });
    const s = b.data.sequence({ ref: "s", value: "MKTAYIAKQR".repeat(4) });
    const p = b.data.annotation({ ref: "p", geometry: "pairs", target: s, value: { start: [5], b_start: [20] } });
    b.view({ ref: "v", label: "V" })
      .section({ ref: "v.s" })
      .block({ sequence: s })
      .track({ label: "P" })
      .dot_plot({ data: p, aspect_ratio: null, color: { theme: "uniform", params: { color: "#d1495b" } } });
    const t = await lit(b.getDocument() as SeqViewSpecDoc, { width: 900, height: 600 });
    // 40 positions over the data area at 20 px a cell: the plot's 320 px a window of about 8 px rows.
    t.vp?.setView(20, 0, "api");
    const r = t.frame(recording(), Dirty.all);
    const red = 0xd1495b;
    expect(Array.from(r.pathItems.color.subarray(0, r.pathItems.n)).filter((c) => c === red)).toEqual([]);
    expect(Array.from(r.glyphItems.color.subarray(0, r.glyphItems.n)).filter((c) => c === red).length).toBe(1);
  });
});
