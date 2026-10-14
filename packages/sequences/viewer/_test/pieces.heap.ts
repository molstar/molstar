/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AX1's reversed draws through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): every built-in renderer's frame on reversed and ranged pieces (every section
// `direction: -1`, every block of a sequence a `range`) allocates nothing after warm-up, at each of the kit's zooms and
// below the row's minimum, linear and wrapped, on documents that reach their paths: every type on two chains, 40,000
// positions zoomed out (classes, pyramids), dense means, thousands of features through their density, an alignment's
// rows and logo, a matrix and a miniature. The run conversions are the axis's (`runToShown`, `shownToRuns`),
// allocation-free in either direction. Then the overlay's hover frames (overlay.heap.ts's) on reversed pieces. The
// positive control goes through the kit's `heap` hook on reversed letters over every type: at each frame the same
// warm-up and measurement of the frame allocating an object each (`leaking`) must catch it.
// Each case loads the viewer afresh (`vi.resetModules`; draws.heap.ts says why: Maglev and the module runner).

import { describe, expect, it, vi } from "vitest";
import type { LayerRenderer } from "../src/index.js";
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  HEAP_DEEP,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from "./heap.js";

const N = FRAME_N;
const WARM = FRAME_WARM;
/** The kit's guards' measurement: warmed by count, then `N` frames and `4N` measured, again when not flat. */
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

// A reversed piece goes through the same code for every renderer of a kind (the axis's run conversions), so one
// renderer is kept for each of the three: letters, values and color (values), ranges (ranges, points, variants,
// secondary); the renderers with a code path of their own stay. The forward cases of every renderer are in
// draws.heap.ts and rn2-draws.heap.ts. Dropped on 2026-10-12 (the tier's run time): letters on 40,000 positions,
// color, ruler, emphasis, variants, secondary, ranges and points on their density documents; `HEAP_DEEP` runs them.
const CASES: [string, string][] = [
  ["letters", "all types"],
  ["values", "40,000 positions"],
  ["pairs", "all types"],
  ["ranges", "all types"],
  ["rows", "family"],
  ["logo", "family"],
  ["matrix", "matrix"],
  ["seqv:map", "miniature"],
  ...((HEAP_DEEP
    ? [
        ["letters", "40,000 positions"],
        ["color", "means"],
        ["ruler", "all types"],
        ["emphasis", "all types"],
        ["ranges", "density"],
        ["points", "density"],
        ["variants", "all types"],
        ["secondary", "all types"],
      ]
    : []) as [string, string][]),
];

async function run(name: string, doc: string, heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import("../src/index.js");
  const { rendererConformance } = await import("../src/testing/index.js");
  const { asText, directed, matrixDoc, miniatureDoc } = await import("./conformance-fixtures.js");
  const { meanDoc } = await import("./draw-fixtures.js");
  const { opened } = await import("./layout-fixtures.js");
  const { allTypes, scaleDoc } = await import("./prepare-fixtures.js");
  const { densityDoc } = await import("./rn2-fixtures.js");
  const { case3 } = await import("../../core/_test/format/cases.js");
  const make = {
    "all types": () => allTypes(),
    "40,000 positions": () => scaleDoc({ length: 20_000 }),
    means: () => meanDoc(10_000),
    density: () => densityDoc({ type: name === "ranges" ? "ranges" : "points", items: 6000, len: 300 }),
    family: case3,
    matrix: matrixDoc,
    miniature: () => miniatureDoc(),
  }[doc] as () => ReturnType<typeof allTypes>;
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === name) as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      const d = directed(make(), -1, true);
      return (
        await opened(geometry === "wrapped" ? asText(d) : d, {
          renderers: has ? [...all] : all.filter((r) => r !== entry),
        })
      ).ctx;
    },
    heap,
    items: false,
  });
}

describe("AX1: reversed and ranged pieces' frames allocate nothing after warm-up", () => {
  for (const [name, doc] of CASES)
    it(`${name} on ${doc}, reversed`, async () => {
      const heap = (await run(name, doc)).filter((c) => c.check === "heap");
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
    }, 240_000);

  it("letters on all types, reversed, an object allocated a frame (the positive control)", async () => {
    const heap = (await run("letters", "all types", control)).filter((c) => c.check === "heap");
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 240_000);
});

describe("AX1: the overlay's hover frames on reversed pieces allocate nothing after warm-up", () => {
  it.each(["scale", "rows"] as const)(
    "hover frames over %s, reversed, with a selection and a highlight",
    async (kind) => {
      vi.resetModules();
      const { hoverFrames } = await import("./overlay-fixtures.js");
      const { frame, hovers } = await hoverFrames(kind, -1);
      const t0 = performance.now();
      for (let i = 0; i < WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
      // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
      let settled = 0;
      settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
      let i = 0;
      const r = heapScaling(N, () => frame(i++));
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(hovers()).toBeGreaterThan(WARM);
    },
    120_000,
  );
});
