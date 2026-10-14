/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's frames allocate nothing after warm-up (CS1a; seqspec:viewer#the-render-backend): pans and scrolls
// through the copy and the strips' draws at three zooms (every type), range labels sticking at the view's edges (two
// zooms), a wrapped
// section's scroll frames (letters, emphasis, each group label mode, the gutter's numbers with `none`), and whole frames
// and pans over 400,000 positions zoomed out (cells, values' columns, the ruler), on a recording backend that reports
// `Copy` and `Back` (the ops logged, not applied).
// The overlay's guards (hover frames with a selection, a highlight and row loci) are CS1b's. The positive control runs
// the first case's pans allocating an object a frame (`leaking`) through the same warm-up, settle and measurement,
// which must catch it.
//
// Each case loads the viewer afresh (`vi.resetModules`), as RN1's draw guards do: in one process V8 can leave a draw in
// its middle tier after other documents deoptimized it, where vitest's module runner boxes doubles (RN1).

import { describe, expect, it, vi } from "vitest";
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
} from "./heap.ts";

const N = FRAME_N;
/** A copied frame is cheap, and the draws compile concurrently while frames run: warmed by count and by time. */
const WARM = FRAME_WARM;
const WARM_MS = FRAME_WARM_MS;
/** The guards' measurement: warmed by count and by time, settled (`settleFrame`), then `N` frames and `4N` measured. */
function measure(frame: (i: number) => void) {
  const t0 = performance.now();
  let i = 0;
  for (; i < WARM || performance.now() - t0 < WARM_MS; i++) frame(i);
  settleFrame(() => frame(i++), FRAME_SETTLE_MS);
  return heapScaling(N, () => frame(i++));
}
async function guard(frame: (i: number) => void) {
  const r = measure(frame);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
}

async function load() {
  vi.resetModules();
  const fx = await import("./composer-fixtures.ts");
  const docs = await import("./prepare-fixtures.ts");
  const { createSeqViewSpec } = await import("@molstar/seq-core/builder");
  return { ...fx, ...docs, createSeqViewSpec };
}

describe("composer frames allocate nothing after warm-up", () => {
  it.each([13, 3, 0.4])(
    "pans and scrolls through the copy at cell %f (every type)",
    async (cell) => {
      const { composed, allTypes, recording, COPY_BACK } = await load();
      const t = await composed(allTypes({ length: 2000 }), { width: 1000, height: 300, dpr: 2 });
      const be = recording(COPY_BACK);
      t.vp?.setView(cell, 200, "api");
      t.frame(be);
      const copies = t.composer.stats.copies;
      await guard((i) => {
        const k = i & 3;
        t.vp?.panBy(k === 0 ? 7 : k === 1 ? -7 : 0, k === 2 ? 11 : k === 3 ? -11 : 0, "pointer");
        t.frame(be);
      });
      expect(t.composer.stats.copies - copies).toBeGreaterThan(2 * WARM);
    },
    120_000,
  );

  it("the same pans allocating an object a frame are caught (positive control)", async () => {
    const { composed, allTypes, recording, COPY_BACK } = await load();
    const t = await composed(allTypes({ length: 2000 }), { width: 1000, height: 300, dpr: 2 });
    const be = recording(COPY_BACK);
    t.vp?.setView(13, 200, "api");
    t.frame(be);
    const frame = (i: number) => {
      const k = i & 3;
      t.vp?.panBy(k === 0 ? 7 : k === 1 ? -7 : 0, k === 2 ? 11 : k === 3 ? -11 : 0, "pointer");
      t.frame(be);
    };
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 120_000);

  it.each(["band", "rulerRail"])(
    "every pass through a pan, the %s navigator's window lit",
    async (nav) => {
      const { composed, allTypes, recording, COPY_BACK } = await load();
      const preset = { name: `lab:${nav}`, geometry: "linear", navigator: nav } as const;
      const t = await composed(allTypes({ length: 2000 }), {
        width: 1000,
        height: 300,
        dpr: 2,
        registry: [{ presets: [preset] }],
      });
      t.ctx.ui.setStyle(t.section, preset.name);
      await t.ctx.graph.update();
      const be = recording(COPY_BACK);
      t.vp?.setView(3, 200, "api");
      t.frame(be, 7);
      await guard((i) => {
        t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, "pointer");
        t.frame(be, 7);
      });
    },
    120_000,
  );

  it.each(["aboveEnd", "inlineStart", "none"] as const)(
    "a wrapped section's scroll frames with %s group labels",
    async (mode) => {
      const { composed, recording, COPY_BACK, createSeqViewSpec } = await load();
      const b = createSeqViewSpec();
      const a = b.data
        .sequence({ ref: "a", value: "MKVLATGHEKWRSPQIDNCY".repeat(150) })
        .numbering({ name: "auth", method: "offset", first: 1001 });
      const starts = Array.from({ length: 60 }, (_, k) => k * 50);
      const f = b.data.annotation({
        ref: "f",
        geometry: "ranges",
        target: a,
        value: { start: starts, stop: starts.map((s) => s + 30) },
      });
      const v = b.view({ ref: "v", label: "V", layout: "text", numbering: "auth", height: "box" });
      v.block({ sequence: a });
      const track = v.track({ label: "L" });
      track.letters();
      track.emphasis({ data: f });
      const t = await composed(b.getDocument(), {
        width: 900,
        height: 400,
        dpr: 2,
        render: { wrappedGroupLabels: mode },
      });
      const be = recording(COPY_BACK);
      t.vp?.scrollTo(0, 2000, "api");
      t.frame(be);
      await guard((i) => {
        t.vp?.panBy(0, (i & 1) === 0 ? 37 : -37, "pointer");
        t.frame(be);
      });
    },
    120_000,
  );

  it.each([13, 4])(
    "pans with range labels centered, sticking at the view's edges and cut, at cell %f (seqspec:sources#levels-and-the-zoom-gate)",
    async (cell) => {
      const { composed, recording, COPY_BACK, createSeqViewSpec } = await load();
      const b = createSeqViewSpec();
      const a = b.data.sequence({ ref: "a", value: "MKVLATGHEKWRSPQIDNCY".repeat(30) });
      const r = b.data.annotation({
        ref: "r",
        geometry: "ranges",
        target: a,
        value: {
          start: [20, 250, 300, 330, 420],
          stop: [260, 290, 380, 336, 590],
          label: ["Kinase domain", "Linker", "SH2 domain", "Site", "C-terminal tail"],
        },
      });
      const v = b.view({ ref: "v", label: "V" });
      v.block({ sequence: a });
      v.track({ label: "Features" }).ranges({ data: r, labels: true });
      const t = await composed(b.getDocument(), { width: 1000, height: 200, dpr: 2 });
      const be = recording(COPY_BACK);
      t.vp?.setView(cell, cell * 200, "api");
      t.frame(be);
      const copies = t.composer.stats.copies;
      // Out and back over 160 px: labels move with their ranges, stick, and are cut as their ranges leave the view.
      await guard((i) => {
        t.vp?.panBy((i & 32) === 0 ? 5 : -5, 0, "pointer");
        t.frame(be);
      });
      expect(t.composer.stats.copies - copies).toBeGreaterThan(2 * WARM);
    },
    120_000,
  );

  it.each(HEAP_DEEP ? [0.0019, 0.03, 0.4] : [0.0019, 0.4])(
    "whole frames and pans over 400,000 positions at cell %f",
    async (cell) => {
      const { composed, scaleDoc, recording, COPY_BACK } = await load();
      const t = await composed(scaleDoc(), { width: 900, height: 420, dpr: 2 });
      const be = recording(COPY_BACK);
      t.vp?.setView(cell, 0, "api");
      t.frame(be);
      await guard((i) => {
        if (i % 3 === 0) t.composer.redraw();
        else t.vp?.panBy(i % 2 === 0 ? 37 : -37, 0, "pointer");
        t.frame(be);
      });
    },
    120_000,
  );
});
