/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's hot paths allocate nothing per frame after warm-up (RN1's pattern): the shared density (points stacked and
// majority, raw and from the count pyramid; ranges' coverage; a block counted through its items beside one through its
// pyramid) through the test entries `ranges`, `points` and `variants` will become; `pairs`' draw and `tiny` through the
// kit at each of its zooms, and its hit over a grid of points; the matrix's tiles at the kit's zooms. The positive
// control goes through the kit's `heap` hook on raw points: at each frame the same warm-up and measurement of the frame
// allocating an object each (`leaking`) must catch it. Each case loads the viewer afresh (`vi.resetModules`;
// draws.heap.ts says why: Maglev and the module runner).

import { describe, expect, it, vi } from "vitest";
import type { LayerRenderer } from "../src/index.ts";
import { controlCaught, describeScaling, FRAME_N, FRAME_WARM, heapScaling, isFlat, leaking } from "./heap.ts";

const N = FRAME_N;
const WARM = FRAME_WARM;
/**
 * The guards' measurement: warmed by count, then `N` frames and `4N`, measured twice when the first is not flat: the first run after the kit switches contexts can grow once by about
 * 3 MB (measured: the short run only, the long run and every later run flat; none in one context), which no frame
 * repeats. A frame that allocates fails both.
 */
const measure = (frame: (i: number) => void) => {
  for (let i = 0; i < WARM; i++) frame(i);
  const r = heapScaling(N, frame);
  return isFlat(r) ? r : heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};
/** The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. */
const control = (frame: () => void) => {
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error("the leaking frame was not caught");
};

type Case = { entry: "points" | "ranges" | "pairs" | "matrix"; doc: string; mixed?: boolean };
const CASES: Case[] = [
  { entry: "points", doc: "points, raw" },
  { entry: "points", doc: "points, pyramid, majority" },
  { entry: "points", doc: "points, pyramid, majority", mixed: true },
  { entry: "ranges", doc: "ranges, raw" },
  { entry: "ranges", doc: "ranges, pyramid" },
  { entry: "pairs", doc: "pairs" },
  { entry: "matrix", doc: "matrix" },
];

async function run(c: Case, heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import("../src/index.ts");
  const { rendererConformance } = await import("../src/testing/index.ts");
  const { opened } = await import("./layout-fixtures.ts");
  const fx = await import("./rn2-fixtures.ts");
  const { asText } = await import("./conformance-fixtures.ts");
  const { OrderedSet } = await import("@molstar/seq-core/model");
  const make = {
    "points, raw": () => fx.densityDoc(),
    "points, pyramid, majority": () => fx.densityDoc({ n: 200_000, blocks: 2, items: 20_000, density: "majority" }),
    "ranges, raw": () => fx.densityDoc({ type: "ranges", items: 6000, len: 300 }),
    "ranges, pyramid": () => fx.densityDoc({ n: 400_000, type: "ranges", items: 20_000, len: 3000 }),
    pairs: () => fx.pairsDoc({ count: 2000 }),
    matrix: () => fx.wideMatrix(),
  }[c.doc] as () => ReturnType<typeof fx.densityDoc>;
  const small =
    c.entry === "pairs"
      ? () => fx.pairsDoc({ n: 2000, count: 300 })
      : () =>
          fx.densityDoc({
            n: 2000,
            items: 300,
            blocks: c.mixed ? 2 : 1,
            type: c.entry === "ranges" ? "ranges" : "points",
          });
  const all =
    c.entry === "points" || c.entry === "ranges" ? fx.withDensity(c.entry) : ([...KernelRenderers] as LayerRenderer[]);
  const entry = all.find((r) => r.name === c.entry) as LayerRenderer;
  fx.forced.positions = c.mixed ? [null, OrderedSet.ofRange(0, 200_000)] : null;
  try {
    return await rendererConformance(entry, {
      open: async ({ geometry, entry: has }) =>
        // Wrapped reads items (no pyramid there): a small document, so the text layout's lines stay few.
        (
          await opened(geometry === "wrapped" ? asText(c.entry === "matrix" ? make() : small()) : make(), {
            renderers: has ? all : all.filter((r) => r !== entry),
          })
        ).ctx,
      heap,
      items: false,
    });
  } finally {
    fx.forced.positions = null;
  }
}

describe("RN2's hot paths allocate nothing after warm-up", () => {
  for (const c of CASES)
    it(`${c.entry} on ${c.doc}${c.mixed ? ", a block through its items" : ""}`, async () => {
      const heap = (await run(c)).filter((x) => x.check === "heap");
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
    }, 240_000);

  it("points, raw, an object allocated a frame (the positive control)", async () => {
    const heap = (await run({ entry: "points", doc: "points, raw" }, control)).filter((x) => x.check === "heap");
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
  }, 240_000);

  it("pairs: the hit over a grid of fractional points beside one long pair", async () => {
    vi.resetModules();
    const { KernelRenderers, Hit } = await import("../src/index.ts");
    const { opened } = await import("./layout-fixtures.ts");
    const { layerOf, entryOf } = await import("./draw-fixtures.ts");
    const { pairsDoc } = await import("./rn2-fixtures.ts");
    const t = await opened(pairsDoc(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { section } = layerOf(t.ctx, "Contacts");
    t.ctx.viewports.get(section)?.setView(4, 4 * 15_000, "api");
    await t.ctx.graph.update();
    const e = entryOf(t.ctx, "Contacts");
    const h = Hit.create();
    const xs = Float64Array.from({ length: 64 }, (_, i) => e.left + 3.3 + i * ((e.right - e.left - 8) / 64));
    const ys = Float64Array.from({ length: 8 }, (_, i) => e.top + 0.7 + i * (e.height / 8));
    guard(() => {
      for (let i = 0; i < xs.length; i++)
        for (let j = 0; j < ys.length; j++) t.ctx.layout.hitTest(section, xs[i], ys[j], h);
    });
  }, 240_000);

  it("pairs: the overlay's lit pairs through the channels' runs and the end index", async () => {
    vi.resetModules();
    const { KernelRenderers } = await import("../src/index.ts");
    const { litPairs, MarkBit } = await import("../src/draw/api.ts");
    const { MarkSet } = await import("../src/marks.ts");
    const { opened } = await import("./layout-fixtures.ts");
    const { layerOf } = await import("./draw-fixtures.ts");
    const { pairsDoc } = await import("./rn2-fixtures.ts");
    const t = await opened(pairsDoc(), { renderers: [...KernelRenderers] as LayerRenderer[] });
    const { k, section } = layerOf(t.ctx, "Contacts");
    const st = t.ctx.graph.assembled.value?.view?.prepared[k].state as never;
    const v = t.ctx.layout.axis(section);
    if (!v) throw new Error("no axis");
    const set = new MarkSet();
    set.reset(1);
    set.begin(0);
    for (let i = 0; i < 40; i++) set.add(i * 400, i * 400 + 30);
    set.close(0, 20_000);
    const marks = {
      bitsAt: () => 0,
      any: () => false,
      runs: () => set,
      primary: () => -1,
      item: () => -1,
      itemAnnotation: () => -1,
      element: () => false,
    };
    const out = new Int32Array(5000);
    guard(() => litPairs(v, st, marks, MarkBit.selection | MarkBit.hover, out));
  }, 240_000);
});
