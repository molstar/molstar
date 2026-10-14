/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG2's ring draws allocate nothing per frame after warm-up (seqspec:viewer#the-render-backend): whole frames of a
// circular section with items' labels outside the ring (their leaders and texts from the kept placement), arrowheads by
// a direction field
// and pairs as chords, the pointer moving round the rings (polar hits and the hover). The positive control: the same
// steps allocating an object each (`leaking`), through the same warm-up and measurement, are caught. Each case loads the
// viewer afresh (`vi.resetModules`).

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

/**
 * The guards' measurement: the shared warm-up (2,000 frames and 3 s, restored after the trim: at 1 s, 240 positions
 * failed most runs, about 45 to 150 bytes a step, flat with HEAP_DEEP, the pointer still meeting new items, positions
 * and rings), then `FRAME_N` steps and 4× measured.
 */
function measure(step: (i: number) => void) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) step(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => step(settled++), FRAME_SETTLE_MS);
  let i = 0;
  return heapScaling(FRAME_N, () => step(i++));
}

/** The ring composed over `n` positions, and a step: a pointer move round the rings and a whole frame. */
async function setup(n: number) {
  vi.resetModules();
  const { composed, recording } = await import("./composer-fixtures.ts");
  const { createSeqViewSpec } = await import("@molstar/seq-core/builder");
  const { DefaultPresets } = await import("../src/presets.ts");
  const { Inputs, PointerType } = await import("../src/index.ts");
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "ACGTACGTTA".repeat(n / 10), alphabet: "dna", cyclic: true });
  const at = (x: number) => Math.floor(x * n);
  const value = {
    start: [at(0.8), at(1 / 6), at(0.55)],
    stop: [at(1 / 8), at(1 / 2), at(0.7)],
    label: ["across", "forward", "reverse"],
    strand: ["+", "+", "-"],
  };
  const f = b.data.annotation({ ref: "f", geometry: "ranges", target: s, value } as never);
  const cuts = { boundary: [at(0.3), at(0.31), at(0.32), at(0.9)], label: ["A", "B", "C", "D"] };
  const c = b.data.annotation({ ref: "c", geometry: "ranges", target: s, value: cuts } as never);
  const pv = { start: [at(0.1), at(0.4)], b_target: ["s", "s"], b_start: [at(0.6), at(0.45)] };
  const pp = b.data.annotation({ ref: "pp", geometry: "pairs", target: s, value: pv } as never);
  const v = b.view({ ref: "v", label: "V", style: "circular" }).block({ sequence: s });
  v.track({ label: "" }).ruler();
  v.track({ label: "F", height: 2 }).ranges({ data: f, labels: true, direction_field: "strand" } as never);
  v.track({ label: "C" }).ranges({ data: c, labels: true } as never);
  v.track({ label: "P" }).pairs({ data: pp, display: "both" } as never);
  const registry = [{ presets: DefaultPresets.filter((p) => p.name === "circular") }];
  const t = await composed(b.getDocument(), { width: 600, height: 600, dpr: 2, registry });
  const l = t.ctx.layout.get(t.section)?.current;
  if (!l) throw new Error("no layout");
  const be = recording();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get("hoverLabel")?.remove();
  // Marks across the origin, and the keyboard cursor on the letters (set once: a cursor move publishes an object).
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const step = (i: number) => {
    // Round the rings, in and out across the entries: a polar hit and the hover each move.
    const a = (i * 0.37) % (2 * Math.PI);
    const r = l.ringR * (0.5 + 0.5 * ((i * 0.13) % 1));
    p.x = l.ringCx + r * Math.cos(a);
    p.y = l.ringCy + r * Math.sin(a);
    t.ctx.input.pointer(p);
    t.composer.redraw();
    t.frame(be, 7);
  };
  return { step, be };
}

describe("CG2's ring frames allocate nothing after warm-up", () => {
  it.each([240, 120_000])(
    "labels with leaders, arrowheads and chords, polar hovers round the rings, %i positions",
    async (n) => {
      const { step, be } = await setup(n);
      const r = measure(step);
      expect(be.textItems.n).toBeGreaterThan(0);
      expect(be.pathItems.n).toBeGreaterThan(0);
      expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    },
    120_000,
  );

  it("the same steps allocating an object each are caught (positive control)", async () => {
    const { step } = await setup(240);
    expect(controlCaught(() => measure(leaking(step)))).toBe(true);
  }, 120_000);
});
