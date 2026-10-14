/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CM-10's viewer side (seqspec:seqql#labels-and-indices): a SeqQL range across a cyclic system's origin
// (`(range ring 38-3)`) selected in the viewer lights both of its pieces at the block's two ends, and on a block shown
// from origin 30, where those positions run on, one stretch (its two runs side by side).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { Q } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";
import { Pass, SkinToken } from "../src/draw/index.ts";
import { lit, passRects } from "./overlay-fixtures.ts";

async function ring(origin: number) {
  const b = createSeqViewSpec();
  const r = b.data.sequence({ ref: "ring", value: "ACGTACGTTA".repeat(4), alphabet: "dna", cyclic: true });
  b.data.computed({ ref: "wrap", provider: "query", label: "Wrap", query: Q.read("(range ring 38-3)") });
  b.view({ ref: "v", label: "V" })
    .block({ sequence: r, ...(origin ? { origin } : {}) })
    .track({ label: "L" })
    .letters();
  const t = await lit(b.getDocument(), { width: 900, height: 200 });
  const loci = t.ctx.refs.toLoci("wrap");
  expect(loci).not.toBeNull();
  t.ctx.selection.set(loci, "api");
  const ax = t.ctx.layout.axis(t.section);
  if (!ax) throw new Error("no axis");
  const x = (k: number) => Math.round(ax.x(0, k as never) * ax.dpr);
  const rects = passRects(t.overlay(), Pass.overlay)
    .filter((q) => q[4] === SkinToken.selection)
    .map((q) => [q[0], q[0] + q[2]]);
  return { x, rects, shownOf: (i: number) => ax.shownOf(0, i) };
}

describe("CM-10, the viewer: a range across a cyclic origin", () => {
  it("from origin 0: lit at the block's two ends, positions 37 to 39 and 0 to 2", async () => {
    const { x, rects } = await ring(0);
    // Each piece its cells, its right edge a device px in (the selection's band).
    expect(rects.sort((a, b) => a[0] - b[0])).toEqual([
      [x(0), x(3) - 1],
      [x(37), x(40) - 1],
    ]);
  });

  it("from origin 30: one stretch, shown 7 to 12 (its two runs side by side)", async () => {
    const { x, rects, shownOf } = await ring(30);
    expect([shownOf(37), shownOf(2)]).toEqual([7, 12]);
    // A run each, drawn as its own rect: a device px seam at the origin inside the stretch, accepted (LC-12; the
    // selection's runs are the model's, two across a cyclic origin).
    expect(rects.sort((a, b) => a[0] - b[0])).toEqual([
      [x(7), x(10) - 1],
      [x(10), x(13) - 1],
    ]);
  });
});
