/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// XR1b's heap guard through the conformance kit (seqspec:viewer#extension-points): a frame of `org.example:stars` (a
// star per item through `mark`, colored by a field provider's value) allocates nothing after warm-up, at each of the
// kit's zooms and on rings (CG1b: drawn through the ring map). The positive control goes through the kit's `heap` hook:
// at each of those frames the same warm-up and measurement of the frame allocating an object each (`leaking`) must
// catch it.

import { expect, it } from "vitest";
import type { LayerRenderer } from "../src/draw/api.js";
import { KernelRenderers } from "../src/index.js";
import { DefaultPresets } from "../src/presets.js";
import { MapRenderer } from "../src/renderers/map.js";
import { RowsRenderer } from "../src/renderers/rows.js";
import { rendererConformance, StarsFormat, StarsRenderer } from "../src/testing/index.js";
import { miniatureDoc, ringed, starsDoc } from "./conformance-fixtures.js";
import { familyWithTree } from "./families.js";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, warmFrame } from "./heap.js";
import { opened } from "./layout-fixtures.js";

const N = 2_000;
/**
 * The guard's measurement: warmed until the frame stops allocating (`warmFrame`: `seqv:map`'s layer at 2.08 px climbed
 * past 2,000 frames), then `N` frames and `4N` measured.
 */
const measure = (step: (i: number) => void) => {
  let i = 0;
  const frame = () => step(i++);
  warmFrame(frame);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};
/** The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. */
const control = (frame: () => void) => {
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error("the leaking frame was not caught");
};
const run = (entry: LayerRenderer, heap = guard) =>
  rendererConformance(entry, {
    open: async ({ geometry, entry: has }) =>
      (
        await opened(geometry === "circular" ? ringed(starsDoc()) : starsDoc(), {
          renderers: [...KernelRenderers, ...(has ? [entry] : [])],
          registry: [{ svs: StarsFormat }, { presets: DefaultPresets.filter((p) => p.name === "circular") }],
        })
      ).ctx,
    heap,
  });

it("a stars frame allocates nothing after warm-up, at every zoom and on rings", async () => {
  const heap = (await run(StarsRenderer as LayerRenderer)).filter((c) => c.check === "heap");
  expect(heap.filter((c) => c.geometry === "circular").length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 120_000);

it("a stars frame allocating an object each is caught, at every zoom and on rings (positive control)", async () => {
  const heap = (await run(StarsRenderer as LayerRenderer, control)).filter((c) => c.check === "heap");
  expect(heap.filter((c) => c.geometry === "circular").length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 120_000);

// RN3's entries through the kit with its heap guard: the rows (single rows at every zoom the kit draws, then the dense
// band) and `seqv:map` (its image built once). Hits and legends are conformance.test.ts's.
for (const [entry, doc] of [
  [RowsRenderer, () => familyWithTree(120, 12, { columns: true, sites: true })],
  [MapRenderer, miniatureDoc],
] as const)
  it(`${entry.name} allocates nothing per frame after warm-up, at every zoom the kit draws`, async () => {
    const checks = await rendererConformance(entry as LayerRenderer, {
      open: async ({ entry: has }) =>
        (await opened(doc(), { renderers: has ? [...KernelRenderers] : KernelRenderers.filter((r) => r !== entry) }))
          .ctx,
      heap: guard,
    });
    const heap = checks.filter((c) => c.check === "heap");
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 120_000);
