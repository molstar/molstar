/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-07 (seqspec:viewer#extension-points, `stable`): a stable renderer's draw of a pan's uncovered strip equals a full
// frame's in that strip. Emphasis decides items or density over the whole view, not the strip's few items.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { SkinToken } from "../src/draw/backend.ts";
import { Dirty, SectionComposer } from "../src/draw/index.ts";
import { COPY_BACK, composed, opNames, recording } from "./composer-fixtures.ts";

function doc() {
  const b = createSeqViewSpec();
  const n = 4000;
  const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(n / 10) });
  const start = Array.from({ length: n / 7 }, (_, i) => i * 7);
  const d = b.data.annotation({
    ref: "d",
    geometry: "ranges",
    target: a,
    value: { start, stop: start.map((s) => s + 1) },
  });
  const v = b.view({ ref: "v", label: "V", height: "box" });
  v.block({ sequence: a });
  v.track({ label: "E" }).emphasis({ data: d });
  return b.getDocument();
}

/** The accent rects at x in [x0, x1) (device px) as [x, y, w, h], sorted. */
const accent = (r: ReturnType<typeof recording>, x0: number, x1: number) => {
  const b = r.rectItems;
  const out: number[][] = [];
  for (let i = 0; i < b.n; i++)
    if (b.color[i] === SkinToken.accent && b.x[i] >= x0 && b.x[i] < x1) out.push([b.x[i], b.y[i], b.w[i], b.h[i]]);
  return out.sort((p, q) => p[0] - q[0] || p[1] - q[1]);
};

describe("VL-07: a stable draw's strip equals the full frame there", () => {
  it("emphasis at 2 px a cell: the pan's strip draws density, as the full frame does", async () => {
    const t = await composed(doc(), { width: 1400, height: 300, dpr: 1 });
    t.vp?.setView(2, 200, "api");
    const be = recording(COPY_BACK);
    t.frame(be);
    t.vp?.panBy(30, 0, "api");
    t.frame(be);
    expect(opNames(be)).toContain("copy");
    // The strip uncovered at the data area's right.
    const right = Math.round(t.layout("v.main").dataRight);
    const strip = accent(be, right - 30, right);
    // The same view drawn whole, by a fresh composer.
    const fresh = new SectionComposer(t.ctx, t.section);
    const whole = recording(COPY_BACK);
    fresh.frame(whole, Dirty.all);
    fresh.dispose();
    const full = accent(whole, right - 30, right);
    expect(strip.length).toBeGreaterThan(0);
    expect(strip).toEqual(full);
  });
});
