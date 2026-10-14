/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The folded genes' heap guard (the lane fold, ST4b follow-up): a pan's pass and composer frame over the demo's GENCODE
// slices at TP53 (about 30 transcripts folded into two lanes) allocate nothing once the regions are held, the fold's
// `fold(fit)` and the folded items' read included. Not vacuous: the folded items are what is drawn (two lanes, TP53
// labeled once), and nothing is asked while measured. A frame allocating an object a step fails the same guard (the
// positive control).

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
const WARM = 6 * FRAME_WARM;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

async function foldedFrame() {
  const { composed, recording } = await import("./composer-fixtures.js");
  const { Dirty } = await import("../src/draw/composer.js");
  const { demoSources } = await import("./jbrowse-sources.js");
  const { JBROWSE_SECTION, jbrowseDoc } = await import("./jbrowse-story-fixtures.js");
  const src = demoSources("nav");
  const t = await composed(jbrowseDoc("chr17", 7_661_779, 7_687_538), {
    width: 1000,
    height: 300,
    dpr: 2,
    loaders: src.loaders,
  });
  const m = t.ctx.streams as import("../src/streams.js").Streams;
  const vp = t.vp;
  if (!vp) throw new Error("no viewport");
  t.layout(JBROWSE_SECTION);
  let done = false;
  void m.settled().then(() => {
    done = true;
  });
  for (let k = 0; k < 2000 && !done; k++) {
    t.scheduler.advance(16);
    await new Promise((ok) => setTimeout(ok, 0));
  }
  const be = recording();
  let i = 0;
  const frame = () => {
    vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, "pointer");
    m.pass();
    t.frame(be, Dirty.all);
  };
  // Both pan positions' regions asked before measuring: a region asked at one position and not at the other stays in
  // flight within `Stream.CancelGraceMs` (VT-02), never asked and cancelled every other frame, so nothing is asked
  // while measured.
  for (let k = 0; k < 8; k++) {
    frame();
    await new Promise((ok) => setTimeout(ok, 0));
  }
  return { frame, backend: be, dispose: () => t.ctx.dispose(), requests: () => m.stats.current.requests };
}

describe("the folded genes' pans allocate nothing per frame once their regions are held", () => {
  it("a pass and a frame, TP53 folded", async () => {
    vi.resetModules();
    const g = await foldedFrame();
    g.frame();
    const names = g.backend.textItems.text.slice(0, g.backend.textItems.n);
    expect(names.filter((x) => x === "TP53").length).toBe(1);
    const asked = g.requests();
    expect(guard(g.frame)).toBe("");
    expect(g.requests()).toBe(asked);
    g.dispose();
  }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const g = await foldedFrame();
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
    g.dispose();
  }, 120_000);
});
