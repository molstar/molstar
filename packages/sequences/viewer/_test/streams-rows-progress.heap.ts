/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The 2-D skeleton allocates nothing per frame (seqspec:sources#levels-and-the-zoom-gate; seqspec:viewer#performance):
// `progressRows`, the overlay over a lazy alignment's rows (3,000 × 400 through the tiles protocol), live, DPR 2, with
// every tile pending (skeletons shimmering, the clock advancing a frame each step) and with every tile failed (error
// bands). Each case on fresh modules (`vi.resetModules`); warmed by count and by time, then settled. The positive
// control: the pending frames allocating an object each (`leaking`), warmed and measured alike, are caught.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loaders } from "@molstar/seq-core/format";
import { describe, expect, it, vi } from "vitest";
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from "./heap.js";

const N = 1_000;
const URI = "https://example.org/lazy";

/** The guard's measurement: warmed by count and by time, settled (`settleFrame`), then `N` frames and `4N`. */
function measure(frame: (i: number) => void) {
  const t0 = performance.now();
  let i = 0;
  for (; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  settleFrame(() => frame(i++), FRAME_SETTLE_MS);
  return heapScaling(N, () => frame(i++));
}

/** The lazy rows composed live with every tile pending (its answers never come) or failed, and an overlay frame. */
async function rowsFrame(c: "pending" | "failed") {
  vi.resetModules();
  const { tilesStream } = await import("@molstar/seq-core/stream/tiles");
  const { tileServer } = await import("@molstar/seq-core/stream/tiles/server");
  const { composed, recording } = await import("./composer-fixtures.js");
  const { Dirty } = await import("../src/draw/composer.js");
  const { SkinToken } = await import("../src/draw/backend.js");
  const b = createSeqViewSpec({ title: "lazy" });
  const aln = b.data.alignment({
    ref: "lazy",
    row_count: 3000,
    column_count: 400,
    source: { uri: URI, protocol: "seqstar-tiles/1" },
  });
  b.view({ ref: "v", label: "V", style: "msa", height: "box" })
    .block({ alignment: aln })
    .track({ label: "Rows" })
    .rows();
  // The reference server, but for its tiles: never answered (pending) or refused (failed).
  const server = tileServer({ rowCount: 3000, columnCount: 400 });
  const loaders: Loaders = {
    load: () => Promise.reject(new Error("none")),
    stream: tilesStream((u) =>
      u === URI
        ? (path, signal) =>
            !path.startsWith("tiles/")
              ? server.transport(path, signal)
              : c === "pending"
                ? new Promise<Uint8Array>(() => {})
                : Promise.reject(new Error("503 from the fixture"))
        : null,
    ),
  };
  const t = await composed(b.getDocument(), { width: 1000, height: 300, dpr: 2, loaders });
  t.composer.live = true;
  // The scheduler's clock, read through an array (a double set on a closure would be boxed).
  const clock = new Float64Array(1);
  clock[0] = 1000;
  (t.scheduler as { now: () => number }).now = () => clock[0];
  t.vp?.setView(1.5, 0, "api");
  await t.ctx.graph.update();
  const be = recording();
  t.frame(be, Dirty.tracks);
  // The view's tiles asked (and refused, when failing), then drawn again.
  (t.ctx.streams as unknown as { pass(): void }).pass();
  for (let k = 0; k < 4; k++) await new Promise((ok) => setTimeout(ok, 0));
  t.frame(be, Dirty.tracks);
  const frame = () => {
    // The shimmer's clock moves a frame; failed bands do not animate.
    if (c === "pending") clock[0] += 16;
    t.frame(be, Dirty.overlay);
  };
  const shown = (token: number) => {
    let n = 0;
    for (let k = 0; k < be.rectItems.n; k++) if (be.rectItems.color[k] === token) n++;
    return n;
  };
  const token = c === "pending" ? SkinToken.pending : SkinToken.conflict;
  return { t, frame, drawn: () => shown(token) };
}

describe("the 2-D skeleton and error bands allocate nothing per frame", () => {
  for (const c of ["pending", "failed"] as const)
    it(`every tile ${c}`, async () => {
      const { t, frame, drawn } = await rowsFrame(c);
      frame();
      // Not vacuous: the shades drawn; a pending skeleton shimmers.
      expect(drawn()).toBeGreaterThan(0);
      expect(t.composer.animating).toBe(c === "pending");
      const r = measure(frame);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      t.ctx.dispose();
    }, 120_000);

  it("the same pending frames allocating an object each are caught (positive control)", async () => {
    const { t, frame } = await rowsFrame("pending");
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
    t.ctx.dispose();
  }, 120_000);
});
