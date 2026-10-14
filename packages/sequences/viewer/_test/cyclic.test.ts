/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Cyclic coordinate systems in the viewer (CY1; seqspec:model#cyclic-systems): a feature across the origin is drawn as
// its two runs with no hairline round the outside of its arc, and a focus across the origin centers its shown runs.
// CG1b: range args with start > stop (a link's range across the origin) name both runs.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Pass, SkinToken } from "../src/draw/index.js";
import { HoverKind } from "../src/index.js";
import type { Features } from "../src/prepare.js";
import { open, vp } from "./axis-fixtures.js";
import { layer, opened, preparedRenderers } from "./layout-fixtures.js";
import { lit, passRects } from "./overlay-fixtures.js";

/**
 * 40 nucleotides, cyclic or not: `a` is one item in three fragments, [36, 4) written across the origin and [20, 24);
 * `b` is one row across the origin, [30, 2).
 */
function ring(cyclic: boolean, view?: (v: ReturnType<ReturnType<typeof createSeqViewSpec>["view"]>) => void) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "ACGTACGTTA".repeat(4), alphabet: "dna", cyclic });
  const value = cyclic
    ? { id: ["a", "a", "b"], start: [36, 20, 30], stop: [4, 24, 2] }
    : { id: ["a", "a", "a", "b", "b"], start: [36, 0, 20, 30, 0], stop: [40, 4, 24, 40, 2] };
  const f = b.data.annotation({ ref: "f", geometry: "ranges", target: s, value });
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: s, ...(cyclic ? { origin: 30 } : {}) })
    .track({ label: "F" })
    .ranges({ data: f });
  view?.(v);
  return b.getDocument();
}

describe("cyclic systems in the viewer (CY1)", () => {
  it("joins an item's fragments within its arc, never round the outside of it; a linear one's every gap", async () => {
    const joins = async (cyclic: boolean) => {
      const t = await opened(ring(cyclic), { renderers: preparedRenderers() });
      const f = layer(t.ctx, "F").state as Features;
      const out: number[][] = [];
      for (let k = 0; k < f.count; k++) if (f.joinFrom[k] >= 0) out.push([f.stop[f.joinFrom[k]], f.start[k]]);
      return out.sort((x, y) => x[0] - y[0]);
    };
    // a's arc is [20, 4) round the origin: the gap [24, 36) is in it, [4, 20) is outside; b's [2, 30) is outside.
    expect(await joins(true)).toEqual([[24, 36]]);
    expect(await joins(false)).toEqual([
      [2, 30],
      [4, 20],
      [24, 36],
    ]);
  });

  it("a focus across the origin centers its shown runs (the piece from origin 30: [36, 4) is shown [6, 14))", async () => {
    const { ctx } = await open(ring(true, (v) => v.focus({ target: "s", start: 36, stop: 4 })));
    const s = vp(ctx, ctx.viewports.sections.value[0]).live.current;
    // 8 positions in 160 px: a cell of 20, centered.
    expect([s.cell, s.scrollX, s.visLo[0], s.visHi[0]]).toEqual([20, 120, 6, 14]);
  });

  it("the hovered item lights both its shown runs in a section with an origin (from 2: [36, 4) is [34, 40) and [0, 2))", async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "ACGTACGTTA".repeat(4), alphabet: "dna", cyclic: true });
    const f = b.data.annotation({ ref: "f", geometry: "ranges", target: s, value: { start: [36], stop: [4] } });
    b.view({ ref: "v", label: "V" }).block({ sequence: s, origin: 2 }).track({ label: "F" }).ranges({ data: f });
    const t = await lit(b.getDocument(), { width: 600, height: 200 });
    // The item's runs alone: not the row band its hover draws (the hoverRow behavior's, the same tint).
    t.ctx.behaviors.get("hoverRow")?.remove();
    const l = t.layout(t.section);
    const ax = t.ctx.layout.axis(t.section);
    if (!ax) throw new Error("no axis");
    t.overlay();
    const h = t.ctx.hover.live.current;
    const e = l.entryCount - 1;
    // Shown 36 is index 38, in the item's first run.
    t.hover(ax.x(0, 36 as never) + 2, l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2);
    expect([h.kind, h.index]).toEqual([HoverKind.item, 38]);
    const spans = passRects(t.overlay(), Pass.overlay)
      .filter((r) => r[4] === SkinToken.hover)
      .map((r) => [r[0], r[0] + r[2]])
      .sort((a, c) => a[0] - c[0])
      // Rects that meet (across a cell gap) are one lit run: [36, 40) and [0, 2) are shown side by side.
      .reduce<number[][]>((acc, r) => {
        const last = acc[acc.length - 1];
        if (last && r[0] <= last[1] + 1) last[1] = Math.max(last[1], r[1]);
        else acc.push([...r]);
        return acc;
      }, []);
    const at = (k: number) => Math.round(ax.x(0, k as never));
    // Both runs, each from its first cell's left to its last cell's right (less the cell gap), and nothing else.
    expect(spans.map(([x0, x1]) => [x0, x1 - at(0) > 0 && Math.abs(x1 - at(2)) <= 1 ? at(2) : x1])).toEqual([
      [at(0), at(2)],
      [at(34), spans[1]?.[1]],
    ]);
    expect(Math.abs((spans[1]?.[1] ?? 0) - at(40))).toBeLessThanOrEqual(1);
  });

  it("range args across the origin (a link's start > stop, or beg after end) name both runs; on a linear one nothing", async () => {
    const runs = (l: unknown) =>
      ((l as { elements?: { indices: OrderedSet }[] } | null)?.elements ?? []).map((e) =>
        Array.from({ length: OrderedSet.size(e.indices) }, (_, i) => OrderedSet.getAt(e.indices, i)),
      );
    const t = await opened(ring(true), { renderers: preparedRenderers() });
    expect(runs(t.ctx.refs.toLoci("s", { start: 36, stop: 4 }))).toEqual([[0, 1, 2, 3, 36, 37, 38, 39]]);
    // Labels in the ordinal numbering (1-based): 37 … 4 is the same range.
    expect(runs(t.ctx.refs.toLoci("s", { beg: "37", end: "4", numbering: "ordinal" }))).toEqual(
      runs(t.ctx.refs.toLoci("s", { start: 36, stop: 4 })),
    );
    // A link's select: the selection holds both runs of the unit.
    await t.ctx.commands.run(
      { command: "seqv:select", params: { target: { refs: ["s"], args: { start: 36, stop: 4 } }, mode: "set" } },
      "host:links",
    );
    const sel = t.ctx.selection.current.value.seq;
    expect(sel ? runs(sel).flat().length : 0).toBe(8);
    const linear = await opened(ring(false), { renderers: preparedRenderers() });
    expect(linear.ctx.refs.toLoci("s", { start: 36, stop: 4 })).toBeNull();
  });
});
