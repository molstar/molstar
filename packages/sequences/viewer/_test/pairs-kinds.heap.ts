/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames over cells that hold several kinds of contact allocate nothing after warm-up: the split cells' glyphs and "+n"
// markers through pans, and the overlay lighting the hovered glyph's kind (a hit and a hover per frame). The positive
// control runs the same frames allocating an object each (`leaking`) through the same warm-up and measurement, which must
// catch it. The case loads the viewer afresh (`vi.resetModules`), as the composer's guards do.

import { describe, expect, it, vi } from "vitest";
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from "./heap.ts";

async function setup() {
  vi.resetModules();
  const { lit } = await import("./overlay-fixtures.ts");
  const { recording, COPY_BACK } = await import("./composer-fixtures.ts");
  const { layerOf } = await import("./draw-fixtures.ts");
  const { createSeqViewSpec } = await import("@molstar/seq-core/builder");
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: "a", value: "MKVLATGHEKWRSPQIDNCY".repeat(100) });
  // Hubs every 12 residues: one with 3 kinds, the next with 11 (more than fit: a marker), partners on the next cells.
  const start: number[] = [];
  const bStart: number[] = [];
  const category: string[] = [];
  for (let hub = 0; hub < 150; hub++) {
    const kinds = hub % 2 === 0 ? 3 : 11;
    for (let k = 0; k < kinds; k++) {
      start.push(hub * 12);
      bStart.push(hub * 12 + 1 + (k % 8));
      category.push(`k${k}`);
    }
  }
  const pairs = b.data.annotation({
    ref: "pp",
    geometry: "pairs",
    target: a,
    value: { start, b_target: start.map(() => "a"), b_start: bStart, category },
  } as never);
  const v = b.view({ ref: "v", label: "V", height: "box" }).block({ sequence: a });
  v.track({ label: "Letters" }).letters();
  v.track({ label: "Contacts" }).pairs({
    data: pairs,
    arcs: "hover",
    glyph: { theme: "field_map", params: { map: { k0: "diamond", k1: "circle", k2: "square" } } },
  } as never);
  const t = await lit(b.getDocument(), { width: 1000, height: 300, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(4, 0, "api");
  t.frame(be, 7);
  // The hover label's strings are one allocation per hovered element by design (IN2's guard): the overlay alone.
  t.ctx.behaviors.get("hoverLabel")?.remove();
  const { k } = layerOf(t.ctx, "Contacts");
  const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as { frame: Float64Array; cellKinds: Int32Array };
  return { t, be, s };
}

/** The guard's measurement: warmed by count and by time, then `heapScaling` over `n` frames and `4n`. */
function run(n: number, frame: (i: number) => void) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  let i = 0;
  return heapScaling(n, () => frame(i++));
}

describe("frames over cells of several kinds", () => {
  it("pans, a hover on each glyph and its overlay allocate nothing; an object allocated per frame is caught", async () => {
    const { t, be, s } = await setup();
    expect(s.cellKinds.some((n) => n > 1)).toBe(true);
    const into = (await import("./composer-fixtures.ts")).recording();
    // The glyphs of the first hubs, from a draw of the layer.
    const { record, layerOf } = await import("./draw-fixtures.ts");
    const glyphs = Array.from(record(t.ctx, "Contacts").rec.glyphItems.x.subarray(0, 24), (x) => x / 2);
    expect(layerOf(t.ctx, "Contacts").section).toBe(t.section);
    const y = s.frame[0];
    const frame = (i: number) => {
      t.vp?.panBy((i & 1) === 0 ? 5 : -5, 0, "pointer");
      t.frame(be, 7);
      t.hover(glyphs[i % glyphs.length], y);
      t.overlay(into);
    };
    const r = run(FRAME_N, frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);
});
