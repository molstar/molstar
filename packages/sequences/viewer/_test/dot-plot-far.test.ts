/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-12: a dot plot's segments far off the view never wrap into it. At 28 px a position and DPR 32, a position
// 4,793,490 positions away is 2^32 device px away, which an int32 device x wrapped back on screen.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { pendingSpan } from "../src/draw/dot-plot.js";
import { recording } from "./composer-fixtures.js";
import { lit } from "./overlay-fixtures.js";

const N = 4_800_000;

describe("VL-12: a dot plot at letter zoom on a long sequence", () => {
  it("draws nothing of pairs 2^32 device px off the view, and its own pair in view", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", alphabet: "dna", value: "ACGT".repeat(N / 4) });
    // The view from `here − 2`: position 0 is 2^32 − 2,944 device px left of it, so an int32 x put it 2,944 px in.
    const here = 4_793_489;
    // A dot at [0, 1) on both axes (2^32 px left of the view, and up), one in view, and a diagonal from 100,000
    // positions back through the view (cut to it).
    const pairs = b.data.annotation({
      ref: "pp",
      geometry: "pairs",
      value: {
        target: ["s", "s", "s"],
        start: [0, here, here - 100_000],
        stop: [1, here + 1, here + 3],
        b_target: ["s", "s", "s"],
        b_start: [0, here, here - 100_000],
        b_stop: [1, here + 1, here + 3],
      },
    });
    b.view({ ref: "v", label: "V" })
      .section({ ref: "v.s" })
      .block({ sequence: s })
      .track({ label: "Plot" })
      .dot_plot({ data: pairs });
    const t = await lit(b.getDocument(), { width: 300, height: 300, dpr: 32 });
    t.vp?.setView(28, (here - 2) * 28, "api");
    expect(t.vp?.live.current.cell).toBe(28);
    const r = t.frame(recording());
    const g = r.glyphItems;
    // The one dot in view, its circle; none for the far one.
    expect(g.n).toBe(1);
    // The diagonal, cut to the view: every coordinate within a few device px of the data area.
    const p = r.pathItems;
    const coords = Array.from(p.coords.subarray(0, p.nc));
    expect(p.n).toBe(1);
    const l = t.layout("v.s");
    for (let i = 0; i < coords.length; i += 2) {
      expect(coords[i]).toBeGreaterThanOrEqual(l.dataLeft * 32 - 70);
      expect(coords[i]).toBeLessThanOrEqual(l.dataRight * 32 + 70);
    }
  }, 60_000);
});

describe("LC-10: a drag's rectangle cut to a position past the view each way", () => {
  it("both ends clamped; wholly off one side, none", () => {
    const out = new Int32Array(2);
    // View [100, 200): a drag inside, across either edge, and wholly off either side.
    expect([pendingSpan(100, 200, 150, 120, out), ...out]).toEqual([true, 120, 151]);
    expect([pendingSpan(100, 200, 10, 500, out), ...out]).toEqual([true, 99, 201]);
    expect(pendingSpan(100, 200, 10, 40, out)).toBe(false);
    expect(pendingSpan(100, 200, 300, 260, out)).toBe(false);
  });
});
