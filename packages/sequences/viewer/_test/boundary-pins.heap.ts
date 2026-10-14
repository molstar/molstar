/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Boundaries pinned over their track's letters allocate nothing per hover frame after warm-up
// (seqspec:viewer#performance): the pointer moving across a design row's pinned boundaries and the letters beside them,
// each move a hit (the pin's reach, the residue beside it) and an overlay frame outlining the hovered pin. Positive
// controls: most frames hovered a boundary, and the same frames allocating an object each (`leaking`), through the same
// warm-up and measurement, are caught.

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
} from "./heap.js";

async function setup() {
  vi.resetModules();
  const { lit } = await import("./overlay-fixtures.js");
  const { recording, COPY_BACK } = await import("./composer-fixtures.js");
  const { pinsDoc, S } = await import("./boundary-pins-fixtures.js");
  const { Dirty } = await import("../src/draw/index.js");
  const { HoverKind, Inputs, PointerType } = await import("../src/index.js");
  const ITEM = HoverKind.item;
  const t = await lit(pinsDoc(), { width: 900, height: 400, dpr: 2 });
  t.ctx.behaviors.get("hoverLabel")?.remove();
  t.vp?.setView(16, 0, "api");
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  const a = t.ctx.layout.axis(S);
  const at = a ? a.x(0, a.shownOf(0, 20)) : 0;
  // On the middle boundary and a few px off it (its pin), then the residues beside it; on both rows.
  const xs = [at, at - 2, at + 3, at - 8, at + 8, at + 1];
  const ys = [t.yOf(S, "t-of"), t.yOf(S, "t-own")];
  let boundary = 0;
  let measuring = false;
  const frame = (i: number) => {
    p.x = xs[i % xs.length];
    p.y = ys[(i / xs.length) & 1];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
    if (measuring && t.ctx.hover.live.current.kind === ITEM) boundary++;
  };
  return {
    frame,
    measure: () => {
      measuring = true;
    },
    boundary: () => boundary,
  };
}

/** The guard's measurement: warmed by count and by time, then `heapScaling` over `n` frames and `4n`. */
function run(n: number, frame: (i: number) => void, measure = () => {}) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  measure();
  let i = 0;
  return heapScaling(n, () => frame(i++));
}

describe("pinned boundaries over a hover sweep", () => {
  it("allocate nothing; the pins hover; an object allocated per frame is caught", async () => {
    const { frame, measure, boundary } = await setup();
    const r = run(FRAME_N, frame, measure);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(boundary()).toBeGreaterThan(FRAME_N);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);
});
