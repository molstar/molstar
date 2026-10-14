/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's feature, pair-overlay and logo draws through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): a frame of `ranges`, `points`, `variants`, `secondary` and `logo` allocates
// nothing after warm-up, at each of the kit's zooms and below the row's minimum, linear and wrapped, on documents that
// reach their paths (every feature layer; boundaries pinned over a track's letters; thousands of items through their density, raw and from the count pyramid;
// transcripts with hairlines;
// 8,000 alignment columns for the logo's thin path, raw and from its pyramid); the `pairs` overlay in both of its modes over many lit ends.
// The positive control goes through the kit's `heap` hook on ranges over the features document: at each frame the same
// warm-up and measurement of the frame allocating an object each (`leaking`) must catch it.
// Each case loads the viewer afresh (`vi.resetModules`; draws.heap.ts says why: Maglev and the module runner).

import { describe, expect, it, vi } from "vitest";
import type { LayerRenderer } from "../src/index.js";
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from "./heap.js";

const N = FRAME_N;
/**
 * The guards' measurement: warmed until the frame stops allocating (`warmFrame`: a tiny draw climbs to TurboFan over
 * tens of thousands), then `N` frames and `4N` measured, again when not flat.
 */
const measure = (step: (i: number) => void) => {
  let i = 0;
  const frame = () => step(i++);
  warmFrame(frame);
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

type Case = { entry: "ranges" | "points" | "variants" | "secondary" | "logo"; doc: string };
const CASES: Case[] = [
  { entry: "ranges", doc: "features" },
  { entry: "ranges", doc: "boundaries over letters" },
  { entry: "ranges", doc: "transcripts" },
  { entry: "ranges", doc: "density, raw" },
  { entry: "ranges", doc: "density, pyramid" },
  { entry: "points", doc: "features" },
  { entry: "points", doc: "density, raw" },
  { entry: "points", doc: "density, pyramid" },
  { entry: "variants", doc: "features" },
  { entry: "variants", doc: "scale" },
  { entry: "secondary", doc: "features" },
  { entry: "logo", doc: "logo" },
  { entry: "logo", doc: "wide logo" },
];

async function run(c: Case, heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import("../src/index.js");
  const { rendererConformance } = await import("../src/testing/index.js");
  const { asText } = await import("./conformance-fixtures.js");
  const { opened } = await import("./layout-fixtures.js");
  const { scaleDoc, transcripts } = await import("./prepare-fixtures.js");
  const fx = await import("./rn2-fixtures.js");
  const mine = await import("./rn2b-fixtures.js");
  const { pinsDoc } = await import("./boundary-pins-fixtures.js");
  const make = {
    features: () => mine.featuresDoc({ n: 300 }),
    "boundaries over letters": () => pinsDoc(),
    transcripts: () => transcripts(),
    "density, raw": () => fx.densityDoc({ type: c.entry === "ranges" ? "ranges" : "points", items: 6000, len: 300 }),
    "density, pyramid": () =>
      fx.densityDoc({ n: 400_000, type: c.entry === "ranges" ? "ranges" : "points", items: 20_000, len: 3000 }),
    scale: () => scaleDoc({ length: 20_000 }),
    logo: () => mine.logoDoc(),
    "wide logo": () => mine.logoDoc(12, 8000),
  }[c.doc] as () => ReturnType<typeof mine.featuresDoc>;
  const all = [...KernelRenderers] as LayerRenderer[];
  const entry = all.find((r) => r.name === c.entry) as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) =>
      (
        await opened(geometry === "wrapped" ? asText(make()) : make(), {
          renderers: has ? all : all.filter((r) => r !== entry),
        })
      ).ctx,
    heap,
    items: false,
  });
}

describe("RN2's feature and logo frames allocate nothing after warm-up", () => {
  for (const c of CASES)
    it(`${c.entry} on ${c.doc}`, async () => {
      const heap = (await run(c)).filter((x) => x.check === "heap");
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
    }, 240_000);

  it("ranges on features, an object allocated a frame (the positive control)", async () => {
    const heap = (await run({ entry: "ranges", doc: "features" }, control)).filter((x) => x.check === "heap");
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
  }, 240_000);
});

describe("the pairs overlay allocates nothing after warm-up", () => {
  it.each(["always", "hover", "hover, an item hovered"])(
    "arcs %s, many lit ends",
    async (mode) => {
      const arcs = mode.startsWith("hover") ? "hover" : "always";
      vi.resetModules();
      const { KernelRenderers } = await import("../src/index.js");
      const { opened } = await import("./layout-fixtures.js");
      const { layerOf } = await import("./draw-fixtures.js");
      const { RecordingBackend } = await import("../src/testing/index.js");
      const { pairsDoc } = await import("./rn2-fixtures.js");
      const { lit, MarkBit, overlayFrame } = await import("./rn2b-fixtures.js");
      const t = await opened(pairsDoc({ arcs }), { renderers: [...KernelRenderers] as LayerRenderer[] });
      const { section } = layerOf(t.ctx, "Contacts");
      t.ctx.viewports.get(section)?.setView(0.05, 0, "api");
      await t.ctx.graph.update();
      const runs: [number, number][] = [];
      for (let i = 0; i < 40; i++) runs.push([i * 400, i * 400 + 30]);
      const marks = lit(MarkBit.hover, runs, 1, 20_000) as ReturnType<typeof lit> & {
        hoverLayer: number;
        hoverAnnotation: number;
        hoverItem: number;
      };
      if (mode !== arcs) {
        // The hovered pair alone for the hover (`marks.item` of its annotation), through its index by item.
        const { k } = layerOf(t.ctx, "Contacts");
        const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as { annotation: Int32Array };
        marks.hoverLayer = k;
        marks.hoverAnnotation = s.annotation[0];
        marks.hoverItem = 7;
      }
      const frame = overlayFrame(t.ctx, "Contacts", marks);
      const rec = new RecordingBackend();
      frame.draw(rec);
      expect(rec.pathItems.n).toBeGreaterThan(0);
      guard(() => frame.draw(rec));
    },
    240_000,
  );
});
