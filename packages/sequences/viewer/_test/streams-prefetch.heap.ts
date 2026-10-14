/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST6's heap guards (seqspec:sources#the-viewers-side; the memory rule): a pass with prefetch on allocates nothing once
// what it asks
// is held: the lead past the margin as a pan alternates sides, the next level as a zoom alternates in and out, and the
// snapshots' focus regions (two hints, held). A pass that allocates an object a step fails the same guard (the positive
// control). Fresh modules a case (`vi.resetModules`), in a file of its own.

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
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

async function prefetchFrame(motion: "pan" | "zoom") {
  const { createSeqViewSpec } = await import("@molstar/seq-core/builder");
  const { GRCH38, hg38Streams } = await import("../../core/_test/format/stream-provider.js");
  const { composed } = await import("./composer-fixtures.js");
  const { GENOME } = await import("./stream-fixtures.js");
  const AT = 100_000_000;
  const b = createSeqViewSpec({ title: "hg38" });
  const genome = b.data.source({ ref: "genome", stream: "genome", data_kind: "letters" });
  const hg38 = b.data.complex({ ref: "hg38", source_ref: genome, label: "GRCh38" });
  b.data.source({ ref: "phylop_src", stream: "phylop", data_kind: "values", coordinate_system: "hg38" });
  const phylop = b.data.annotation({ ref: "phylop", geometry: "values", source_ref: "phylop_src" });
  const v = b.view({ ref: "browser", label: "hg38" });
  v.block({ complex: hg38 });
  v.focus({ ref: "browser.focus", target: "hg38/chr2", start: AT, stop: AT + 10_000 });
  v.track({ label: "Ruler" }).ruler();
  v.track({ label: "phyloP" }).values({ data: phylop });
  for (const [ref, at] of [
    ["a", AT],
    ["b", 150_000_000],
    ["c", 20_000_000],
  ] as const)
    b.snapshot({
      ref,
      duration_ms: 1000,
      transition_ms: 0,
      ops: [{ op: "set", ref: "browser.focus", params: { start: at, stop: at + 10_000 } }],
    });
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3) });
  const t = await composed(b.getDocument(), {
    width: 1000,
    height: 200,
    dpr: 2,
    loaders: { load: () => Promise.reject(new Error("none")), stream: s.stream },
  });
  const m = t.ctx.streams as import("../src/streams.js").Streams;
  const vp = t.vp;
  if (!vp) throw new Error("no viewport");
  await t.ctx.snapshots.apply("a");
  vp.zoomTo(1, AT, AT + 10_000);
  const flush = () => new Promise((ok) => setTimeout(ok, 0));
  // The cell at the first view and at a view a third as wide (a zoom in), set through the viewport's own verb.
  const wide = vp.live.current.cell;
  const narrow = wide * 3;
  const x = vp.live.current.scrollX;
  for (let k = 0; k < 8; k++) {
    for (const [dx, cell] of [
      [9, wide],
      [-9, narrow],
    ] as const) {
      if (motion === "pan") vp.panBy(dx, 0, "pointer");
      else vp.setView(cell, x, "api");
      t.layout(GENOME);
      m.pass();
      await flush();
    }
  }
  let i = 0;
  return {
    manager: m,
    asked: () => s.phylop.requests.length,
    prefetched: () => s.phylop.requests.filter((r) => r.priority === "prefetch").length,
    pass: () => {
      // A step of the motion (small ints: a double passed to a call is boxed), and its pass.
      if (motion === "pan") vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, "pointer");
      else vp.setView((i++ & 1) === 0 ? wide : narrow, x, "api");
      m.pass();
    },
    dispose: () => t.ctx.dispose(),
  };
}

describe("a pass with prefetch on allocates nothing once what it asks is held", () => {
  for (const motion of ["pan", "zoom"] as const)
    it(`a pass as the view ${motion === "pan" ? "pans to and fro" : "zooms in and out"}, hints held`, async () => {
      vi.resetModules();
      const g = await prefetchFrame(motion);
      // Not vacuous: prefetches were asked (the lead or the zoom's level, and the snapshots' regions).
      expect(g.prefetched()).toBeGreaterThan(motion === "pan" ? 3 : 1);
      const asked = g.asked();
      expect(guard(g.pass)).toBe("");
      expect(g.asked()).toBe(asked);
      g.dispose();
    }, 120_000);

  it("the same pass allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const g = await prefetchFrame("pan");
    expect(controlCaught(() => measure(leaking(g.pass) as () => void))).toBe(true);
    g.dispose();
  }, 120_000);
});
