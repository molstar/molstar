/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's frames and hits allocate nothing after warm-up (seqspec:viewer#extension-points): the conformance
// kit's guard (`heapScaling` per frame, at the kit's three zooms, `tiny` too) on case 6, on 20 poses × 8 kinds
// (stripes) and at a short track (compact), on mixed data and the absent property (the tiles beside the rows), and on
// 20,000 positions zoomed out (merged cells), case 6 and the wide one reversed (AX1); and a million hits over every row
// and stripe, compact rows and a block without the property, each row and each kind hit somewhere first. Each case
// loads the viewer afresh (`vi.resetModules`): a frame built in this module would box from the second case on (047,
// RN3). The positive control, last: the kit's first frame on the 20 poses allocating an object each (`leaking`), warmed
// and measured alike, is caught.

import { describe, expect, it, vi } from "vitest";
import { controlCaught, leaking, warmFrame } from "../heap.ts";
import { describeScaling, heapScaling, isFlat } from "./heap.ts";

const N = 1_000;
/** Warmed until the frame stops allocating (`warmFrame`: the wide scene's tiny draw climbed past 2,000 frames). */
const measure = (frame: () => void) => {
  warmFrame(frame);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};

const SCENES = ["case 6", "many", "compact", "mixed", "absent", "wide", "case 6 reversed", "wide reversed"] as const;
type Name = (typeof SCENES)[number];

async function load() {
  vi.resetModules();
  const fx = await import("./fixtures.ts");
  const kit = await import("../../src/testing/index.ts");
  const viewer = await import("../../src/index.ts");
  const draws = await import("../conformance-fixtures.ts");
  const layout = await import("../layout-fixtures.ts");
  const cases = await import("../../../core/_test/format/cases.ts");
  const src = await import("../../src/interactions/viewer/index.ts");
  return { fx, kit, viewer, draws, layout, cases, src };
}

describe("the fingerprint's frames allocate nothing after warm-up, at every zoom the kit draws", () => {
  for (const name of SCENES)
    it(name, async () => {
      const { fx, kit, viewer, draws, layout, cases, src } = await load();
      const scene = {
        "case 6": () => ({ doc: cases.case6(), section: "poses-view.main", view: "poses-view" }),
        many: () => ({ ...fx.manyScene(), view: "v" }),
        compact: () => ({ ...fx.manyScene(20, 8, 1), view: "v" }),
        mixed: () => ({ ...fx.mixedScene(), view: "v" }),
        absent: () => ({ ...fx.absentScene(), view: "v" }),
        wide: () => ({ ...fx.wideScene(), view: "v" }),
        // AX1: every piece reversed.
        "case 6 reversed": () => ({
          doc: draws.directed(cases.case6(), -1),
          section: "poses-view.main",
          view: "poses-view",
        }),
        "wide reversed": () => ({ ...fx.wideScene(), doc: draws.directed(fx.wideScene().doc, -1), view: "v" }),
      }[name as Name]();
      const entry = src.FingerprintRenderer as never;
      const checks = await kit.rendererConformance(entry, {
        open: async ({ geometry, entry: has }) =>
          (
            await layout.opened(geometry === "wrapped" ? draws.asText(scene.doc) : scene.doc, {
              renderers: [...viewer.KernelRenderers],
              view: scene.view,
              registry: has
                ? [src.InteractionsViewerRegistry, fx.ShortFingerprint]
                : [{ svs: src.InteractionsViewerRegistry.svs }],
            })
          ).ctx,
        heap: guard,
        items: false,
      });
      const heap = checks.filter((c) => c.check === "heap");
      expect(heap.length).toBeGreaterThanOrEqual(6);
      expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
    }, 240_000);
});

describe("the hit allocates nothing: a million over every row and stripe", () => {
  it.each(["case 6", "compact", "mixed"] as const)(
    "%s",
    async (name) => {
      const { fx, layout, viewer, cases, src } = await load();
      const scene = {
        "case 6": () => ({
          doc: cases.case6(),
          section: "poses-view.main",
          track: "poses-view.fingerprint",
          view: "poses-view",
        }),
        compact: () => ({ ...fx.manyScene(20, 8, 1), view: "v" }),
        mixed: () => ({ ...fx.mixedScene(), view: "v" }),
      }[name]();
      const t = await layout.opened(scene.doc, {
        renderers: [...viewer.KernelRenderers],
        view: scene.view,
        registry: [src.InteractionsViewerRegistry, fx.ShortFingerprint],
      });
      for (const ref of t.ctx.viewports.sections.value) t.ctx.viewports.get(ref)?.setSize(1000, 600, 2);
      await t.ctx.graph.update();
      const l = t.layout(scene.section);
      const e = Array.from(l.entryTrack.subarray(0, l.entryCount)).indexOf(t.trackIndex(scene.track));
      const top = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e];
      const points: number[] = [];
      for (let x = l.dataLeft; x < l.dataRight; x += 3)
        for (let y = top - 2; y < top + l.entryHeight[e] + 2; y += 0.5) points.push(x, y);
      const hit = viewer.Hit.create();
      // Every row (state) and every kind (part) is hit somewhere.
      const states = new Set<number>();
      const parts = new Set<number>();
      const blocks = new Set<number>();
      for (let k = 0; k < points.length; k += 2) {
        t.ctx.layout.hitTest(scene.section, points[k], points[k + 1], hit);
        states.add(hit.state);
        parts.add(hit.part);
        if (hit.state >= 0) blocks.add(hit.block);
      }
      // Mixed: the block without the property (the third) is hit as bands too.
      if (name === "mixed") expect([...blocks].sort()).toEqual([0, 1, 2]);
      const s = t.ctx.graph.assembled.value?.view?.prepared.find((p) => p.renderer === "interactions:fingerprint")
        ?.state as {
        matrix: { states: number };
        kinds: Uint8Array;
      };
      for (let q = 0; q < s.matrix.states; q++) expect(states.has(q)).toBe(true);
      for (const k of s.kinds) if (name !== "mixed") expect(parts.has(k)).toBe(true);
      expect(states.has(-1)).toBe(true);
      const n = points.length / 2;
      const step = (i: number) => {
        const k = (i * 7) % n;
        t.ctx.layout.hitTest(scene.section, points[2 * k], points[2 * k + 1], hit);
      };
      // Warmed past the hit's tier-up (compact: about 1.4 bytes a hit for 1.5 million hits, then none), a thousand a frame.
      let w = 0;
      warmFrame(() => {
        for (let q = 0; q < 1_000; q++) step(w++);
      });
      const r = heapScaling(250_000, step);
      expect(isFlat(r), describeScaling(r, 250_000)).toBe(true);
    },
    240_000,
  );
});

describe("the positive control", () => {
  it("the kit's first frame on the 20 poses, allocating an object each, is caught", async () => {
    const { fx, kit, viewer, draws, layout, src } = await load();
    const scene = fx.manyScene();
    const caught: boolean[] = [];
    await kit.rendererConformance(src.FingerprintRenderer as never, {
      open: async ({ geometry, entry: has }) =>
        (
          await layout.opened(geometry === "wrapped" ? draws.asText(scene.doc) : scene.doc, {
            renderers: [...viewer.KernelRenderers],
            view: "v",
            registry: has
              ? [src.InteractionsViewerRegistry, fx.ShortFingerprint]
              : [{ svs: src.InteractionsViewerRegistry.svs }],
          })
        ).ctx,
      heap: (frame) => {
        if (caught.length === 0) caught.push(controlCaught(() => measure(leaking(frame) as () => void)));
      },
      items: false,
    });
    expect(caught).toEqual([true]);
  }, 240_000);
});
