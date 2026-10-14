/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LM1's letter maps through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): the map lookup is per drawn letter and per cell, so `letters` (its own, `of`
// another polymer's, and 20,000 positions into classes) and `rows` (their tiles and letters) with a complement map, and
// `emphasis` redrawing its track's mapped letters, forward and as a reverse complement, linear and wrapped, allocate
// nothing after warm-up. The positive control goes through the kit's `heap` hook on letters over every type: at each
// frame the same warm-up and measurement of the frame allocating an object each (`leaking`) must catch it. Each case
// loads the viewer afresh (`vi.resetModules`; draws.heap.ts says why).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it, vi } from "vitest";
import type { LayerRenderer } from "../src/index.js";
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_WARM,
  HEAP_DEEP,
  heapScaling,
  isFlat,
  leaking,
} from "./heap.js";

const N = FRAME_N;
const WARM = FRAME_WARM;
/** The guards' measurement: warmed by count, then `N` frames and `4N` measured, again when not flat. */
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

const CASES: [string, string][] = [
  ["letters", "all types"],
  ["letters", "20,000 positions"],
  ["letters", "gene"],
  ["rows", "family"],
  ["emphasis", "emphasized letters"],
];

/** 2,000 nucleotides, every third emphasized, over the letters of their own track. */
function emphasized() {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", label: "S", alphabet: "dna", value: "ACGTTGCAGGCT".repeat(167) });
  const start = Array.from({ length: 600 }, (_, k) => 3 * k);
  const f = b.data.annotation({
    ref: "f",
    geometry: "ranges",
    target: s,
    value: { start, stop: start.map((x) => x + 2) },
  });
  b.view({ ref: "v", label: "V" }).track({ label: "S" }).letters().emphasis({ data: f });
  return b.getDocument();
}

async function run(name: string, doc: string, reversed: boolean, heap = guard) {
  vi.resetModules();
  const { KernelRenderers } = await import("../src/index.js");
  const { rendererConformance } = await import("../src/testing/index.js");
  const { asText, directed, mapped } = await import("./conformance-fixtures.js");
  const { opened } = await import("./layout-fixtures.js");
  const { allTypes, gene, scaleDoc } = await import("./prepare-fixtures.js");
  const { case3 } = await import("../../core/_test/format/cases.js");
  const make = {
    "all types": () => allTypes(),
    "20,000 positions": () => scaleDoc({ length: 20_000 }),
    gene: () => gene(),
    family: case3,
    "emphasized letters": emphasized,
  }[doc] as () => ReturnType<typeof allTypes>;
  const all = KernelRenderers as LayerRenderer[];
  const entry = all.find((r) => r.name === name) as LayerRenderer;
  return rendererConformance(entry, {
    open: async ({ geometry, entry: has }) => {
      const d = reversed ? directed(mapped(make()), -1, true) : mapped(make());
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

// The reverse complement goes through the same map lookups for every letters document and for emphasis, so it is kept
// on one letters document and on rows (dropped on 2026-10-12: letters on 20,000 positions and on gene, emphasis;
// `HEAP_DEEP` runs them).
const REVERSED = new Set(["letters on all types", "rows on family"]);

describe("LM1: frames drawn through a letter map allocate nothing after warm-up", () => {
  for (const reversed of [false, true])
    for (const [name, doc] of CASES)
      if (!reversed || HEAP_DEEP || REVERSED.has(`${name} on ${doc}`))
        it(`${name} on ${doc}, ${reversed ? "as a reverse complement" : "complemented"}`, async () => {
          const heap = (await run(name, doc, reversed)).filter((c) => c.check === "heap");
          expect(heap.length).toBeGreaterThan(0);
          expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
        }, 240_000);

  it("letters on all types, an object allocated a frame (the positive control)", async () => {
    const heap = (await run("letters", "all types", false, control)).filter((c) => c.check === "heap");
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 240_000);
});
