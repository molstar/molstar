/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the `follow` behavior (LK1; seqspec:viewer#performance), on `heapScaling`: a leader with two
// followers (a protein through the codon step and the gene reversed), driven at input rate by a drag and by a Ctrl
// wheel zoom; each step moves all three. Fractional coordinates, as an engine's; every default behavior subscribed. The
// positive control: the wheel's steps allocating an object each (`leaking`), primed and measured alike, are caught.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { Inputs, Modifier, PointerType } from "../src/index.js";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.js";
import { driven } from "./input-fixtures.js";

const N = 50_000;
let sink = 0;

/** Case 14's gene, a second copy of its two exons' worth of codons, and three sections: gene, protein, gene reversed. */
function doc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const coding = "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCCATTGGTACTTCAACAGCCAGGTGGAG".repeat(20);
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: `GAATTCGCCACC${coding}TAAGCGGCC` });
  const protein = translate(b.data, { ref: "protein", from: gene, cds: [[12, 12 + coding.length]] });
  const v = b.view({ ref: "v", label: "V" });
  v.section({ ref: "v.dna" }).block({ sequence: gene }).track().letters();
  v.section({ ref: "v.protein", follow: "v.dna" }).block({ sequence: protein }).track().letters();
  v.section({ ref: "v.rev", follow: "v.dna", direction: -1 }).block({ sequence: gene }).track().letters();
  return b.getDocument();
}

function prime(step: (i: number) => void) {
  const gc = (globalThis as { gc?: () => void }).gc;
  let k = 0;
  for (let chunk = 0, quiet = 0; chunk < 100 && (quiet < 3 || chunk < 15); chunk++) {
    gc?.();
    const before = process.memoryUsage().heapUsed;
    for (let j = 0; j < 20_000; j++) step(k++);
    quiet = process.memoryUsage().heapUsed - before < 2048 ? quiet + 1 : 0;
  }
}

/** The guard's measurement: primed, then `N` and `4N` steps. */
function measure(step: (i: number) => void) {
  prime(step);
  return heapScaling(N, step);
}

/** A Ctrl wheel zoom on the leader at x 400.5, in and out by turns: its step and the three viewports. */
async function wheeling() {
  const t = await driven(doc(), { labels: false });
  const [dna, protein, rev] = ["v.dna", "v.protein", "v.rev"].map((r) => t.vp(r));
  dna.zoomTo(0, 300, 400);
  const w = Inputs.wheel();
  w.section = t.ctx.layout.sectionIndex("v.dna");
  w.y = t.layout("v.dna").scrolledTop + 8.5;
  w.modifiers = Modifier.ctrl;
  w.deltaX = 0;
  const step = (i: number) => {
    w.x = 400.5 + (i & 7);
    w.deltaY = i & 1 ? -3.25 : 3.25;
    sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
  };
  return { step, dna, protein, rev };
}

describe("follow at input rate", () => {
  it("a drag on the leader moves both followers each step: flat", async () => {
    const t = await driven(doc(), { labels: false });
    const [dna, protein, rev] = ["v.dna", "v.protein", "v.rev"].map((r) => t.vp(r));
    dna.zoomTo(0, 300, 360);
    const y = t.layout("v.dna").scrolledTop + 8.5;
    const e = Inputs.pointer();
    e.section = t.ctx.layout.sectionIndex("v.dna");
    t.pointer("down", "v.dna", 500.5, y);
    t.pointer("move", "v.dna", 490.5, y);
    expect(t.ctx.ui.drag.value?.owner).toBe("dragPan");
    e.type = PointerType.move;
    e.buttons = 1;
    e.button = -1;
    e.y = y;
    const [p0, r0] = [protein.live.version.value, rev.live.version.value];
    const step = (i: number) => {
      // Whole px: the leader's scroll snaps to device px, so a fractional step would drift it to an end over the run.
      e.x = 490.5 + (i & 15) * 3;
      t.ctx.input.pointer(e);
    };
    const r = measure(step);
    t.pointer("up", "v.dna", 500.5, y);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(protein.live.version.value - p0).toBeGreaterThan(N);
    expect(rev.live.version.value - r0).toBeGreaterThan(N);
    expect(protein.live.current.origin).toBe("behavior");
  }, 120_000);

  it("a Ctrl wheel zoom on the leader, in and out, zooms both followers each step: flat", async () => {
    const { step, dna, protein, rev } = await wheeling();
    const [p0, r0] = [protein.live.version.value, rev.live.version.value];
    const r = measure(step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(protein.live.version.value - p0).toBeGreaterThan(N);
    expect(rev.live.version.value - r0).toBeGreaterThan(N);
    expect(rev.live.current.cell).toBe(dna.live.current.cell);
  }, 120_000);

  it("the same wheel steps allocating an object each are caught (positive control)", async () => {
    const { step } = await wheeling();
    expect(controlCaught(() => measure(leaking(step)))).toBe(true);
  }, 120_000);
});
