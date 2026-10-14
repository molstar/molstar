/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A data hover's row band (the `hoverRow` behavior; seqspec:interaction#hover, seqspec:viewer#performance) allocates
// nothing per hover frame after warm-up: the pointer sweeping down the hover document's tracks (letters, ranges, points,
// values, a matrix's bands through its renderer's `rowHover`) at DPR 2, the element changing on most moves, an overlay
// frame per move on the recording backend. Counted over the measured frames alone, not vacuous: the element changed on
// most moves, and most frames drew the band (its header part at x 0 in the hover tint). The positive control: the same
// frames allocating an object each (`leaking`), through the same warm-up, settle and measurement, are caught.
//
// Its own file, the viewer loaded afresh (`vi.resetModules`): cases in one file share V8's feedback.

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
} from "./heap.ts";

const N = 1_000;

/**
 * The guard's measurement: warmed by count and by time, settled (`settleFrame`), `warmed` called, then `N` frames and
 * `4N` measured.
 */
function measure(frame: (i: number) => void, warmed = () => {}) {
  const t0 = performance.now();
  let i = 0;
  for (; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  settleFrame(() => frame(i++), FRAME_SETTLE_MS);
  warmed();
  return heapScaling(N, () => frame(i++));
}

/** The hover document lit at DPR 2, and a frame: the pointer on a track's row, an overlay frame, the band counted. */
async function setup() {
  vi.resetModules();
  const { Dirty, SkinToken } = await import("../src/draw/index.ts");
  const { Inputs, PointerType } = await import("../src/index.ts");
  const { COPY_BACK, recording } = await import("./composer-fixtures.ts");
  const { hoverDoc } = await import("./hover-fixtures.ts");
  const { lit } = await import("./overlay-fixtures.ts");
  const t = await lit(hoverDoc(100), { width: 1000, height: 600, dpr: 2 });
  // The overlay alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get("hoverLabel")?.remove();
  t.vp?.setView(12, 0, "api");
  await t.ctx.graph.update();
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const tracks = ["t.letters", "t.ranges", "t.points", "t.values", "t.matrix"];
  const ys: number[] = [];
  for (const ref of tracks) ys.push(t.yOf(t.section, ref));
  // A matrix's three bands, each its own row band.
  const m = ys[ys.length - 1];
  ys.push(m - 12, m + 12);
  // Eight positions against seven rows (coprime): the sweep meets every pair.
  const xs = [3, 4, 5, 6, 7, 8, 9, 10].map((i) => t.xOf(t.section, 0, i));
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const HOVER = SkinToken.hover;
  let drawn = 0;
  const frame = (i: number) => {
    p.x = xs[i % xs.length];
    p.y = ys[i % ys.length];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    const r = be.rectItems;
    let band = 0;
    for (let j = 0; j < r.n && band === 0; j++) if (r.x[j] === 0 && r.color[j] === HOVER) band++;
    drawn += band;
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover, drawn: () => drawn };
}

describe("a data hover's row band allocates nothing after warm-up", () => {
  it("a hover sweep down the tracks, the band drawn", async () => {
    const { frame, hovers, drawn } = await setup();
    let hovers0 = 0;
    let drawn0 = 0;
    let measured = 0;
    const r = measure(
      (i) => {
        frame(i);
        measured++;
      },
      () => {
        hovers0 = hovers();
        drawn0 = drawn();
        measured = 0;
      },
    );
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // Not vacuous, over the measured frames: the element changed on most moves, and most drew the band.
    expect(measured).toBeGreaterThanOrEqual(5 * N);
    expect(hovers() - hovers0).toBeGreaterThan(0.8 * measured);
    expect(drawn() - drawn0).toBeGreaterThan(0.8 * measured);
  }, 120_000);

  it("the same frames allocating an object each are caught (positive control)", async () => {
    const { frame } = await setup();
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 120_000);
});
