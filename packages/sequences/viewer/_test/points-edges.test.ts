/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-05: a variant over several positions (a deletion) on edge geometry: partly scrolled out, on a reversed piece and
// across a rotated cyclic block's seam. Its line is a piece per shown run in view, its marker the middle of its longest
// shown run, and the hit finds it there (draw and hit by one center).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { Hit } from "../src/draw/api.ts";
import { recording } from "./composer-fixtures.ts";
import { lit } from "./overlay-fixtures.ts";

const S = "v.main";

async function variants(block: Record<string, unknown>, start: number[], stop: number[], scroll: number) {
  const b = createSeqViewSpec();
  const p = b.data.sequence({ ref: "p", value: "ACGTACGTTA".repeat(20), alphabet: "dna", cyclic: true });
  const v = b.data.annotation({ ref: "vars", geometry: "ranges", target: p, value: { start, stop } });
  const sec = b.view({ ref: "v", label: "V" }).section({ ref: S });
  sec.block({ sequence: p, ...block } as never);
  sec.track({ ref: "t-variants", label: "Variants", height: 2 }).variants({ data: v });
  const t = await lit(b.getDocument(), { width: 900, height: 400 });
  t.vp?.setView(20, scroll, "api");
  const be = t.frame(recording());
  const ax = t.ctx.layout.axis(S);
  if (!ax) throw new Error("no axis");
  const r = be.rectItems;
  const lines = Array.from({ length: r.n }, (_, i) => [r.x[i], r.y[i], r.w[i], r.h[i]]).filter((q) => q[3] <= 2);
  const g = be.glyphItems;
  const glyphs = Array.from({ length: g.n }, (_, i) => [g.x[i], g.y[i]]);
  /** A cell's left edge (CSS px). */
  const left = (i: number) => ax.x(0, ax.shownOf(0, i));
  const right = (i: number) => ax.x(0, (ax.shownOf(0, i) + 1) as never);
  const hit = Hit.create();
  /** The item hit at a glyph's center (device px), -1: none. */
  const at = (x: number, y: number) => (t.ctx.layout.hitTest(S, x, y, hit) ? hit.item : -1);
  return { t, ax, lines, glyphs, left, right, at };
}

describe("VL-05: a deletion's line and marker on edge geometry; the hit at the marker", () => {
  it("partly scrolled out past the right edge: its line from its first cell to the edge, its marker on what is shown", async () => {
    const { ax, lines, glyphs, left, at } = await variants({}, [2, 10], [5, 150], 0);
    const edge = left(ax.visHi[0]);
    // Item 0 over its three cells; item 1 from cell 10 to the data area's right edge (before VL-05: a stray 1 px tick
    // and no line).
    expect(lines.map((q) => [q[0], q[2]])).toEqual([
      [left(2), left(5) - left(2)],
      [left(10), edge - left(10)],
    ]);
    expect(glyphs.map((g) => g[0])).toEqual([(left(2) + left(5)) / 2, (left(10) + edge) / 2]);
    expect(glyphs.map(([x, y]) => at(x, y))).toEqual([0, 1]);
  });

  it("a reversed piece: one line over its cells, the marker at their middle", async () => {
    const { lines, glyphs, left, right, at } = await variants({ direction: -1 }, [10], [16], 3600);
    // Reversed: position 15 on the left, 10 on the right (before VL-05: a 1 px tick).
    expect(lines.map((q) => [q[0], q[2]])).toEqual([[left(15), right(10) - left(15)]]);
    expect(glyphs.map((g) => g[0])).toEqual([(left(15) + right(10)) / 2]);
    expect(at(glyphs[0][0], glyphs[0][1])).toBe(0);
  });

  it("across a rotated cyclic block's seam: a piece at each end, the marker on the piece in view", async () => {
    // Origin 20: [20, 22) at the block's start, [17, 20) at its end.
    const start = await variants({ origin: 20 }, [17], [22], 0);
    expect(start.lines.map((q) => [q[0], q[2]])).toEqual([[start.left(20), start.right(21) - start.left(20)]]);
    expect(start.glyphs.map((g) => g[0])).toEqual([(start.left(20) + start.right(21)) / 2]);
    expect(start.at(start.glyphs[0][0], start.glyphs[0][1])).toBe(0);
    const end = await variants({ origin: 20 }, [17], [22], 1e9);
    expect(end.lines.map((q) => [q[0], q[2]])).toEqual([[end.left(17), end.right(19) - end.left(17)]]);
    expect(end.glyphs.map((g) => g[0])).toEqual([(end.left(17) + end.right(19)) / 2]);
    expect(end.at(end.glyphs[0][0], end.glyphs[0][1])).toBe(0);
  });
});
