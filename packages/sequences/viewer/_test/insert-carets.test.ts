/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Carets (seqspec:model#projection, seqspec:interaction#hover): an insert lit on one side of a partial map lights the
// boundary it sits at on the other, a caret, not a position; a stretch a design deletes lights the junction it leaves;
// a boundary item lights a caret where it is; beyond the map's pairs, nothing. The composer draws a caret as a bar on
// the edge between two cells.

import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Dirty } from "../src/draw/index.js";
import { ChannelMarks, MarkBit } from "../src/index.js";
import type { MarkSet } from "../src/marks.js";
import { recording } from "./composer-fixtures.js";
import { caretsDoc as doc } from "./insert-carets-fixtures.js";
import { lit } from "./overlay-fixtures.js";

/** A channel's runs and carets per section: [runs as [s, t], carets]. */
function read(ch: ChannelMarks, bit: number, sections: number) {
  const out: [number[][], number[]][] = [];
  for (let s = 0; s < sections; s++) {
    const m = ch.marks(s)?.runs(bit) as MarkSet;
    const runs: number[][] = [];
    for (let r = m.from(0); r < m.to(0); r++) runs.push([m.start[r], m.stop[r]]);
    out.push([runs, Array.from(m.caretAt.subarray(0, m.carets))]);
  }
  return out;
}

async function open() {
  const t = await lit(doc(), { width: 900, height: 400 });
  const set = t.set;
  const unit = (key: string) => set.units.key.findIndex((k) => set.strings.get(k) === key);
  return { ...t, ch: new ChannelMarks(t.ctx), unit };
}

describe("carets", () => {
  it("an insert lights a caret on the other side, a deleted stretch the junction it leaves; beyond the pairs nothing", async () => {
    const t = await open();
    t.ctx.selection.set(
      { kind: "seq-loci", set: t.set, elements: [{ unit: t.unit("d/unit"), indices: range(20, 24) }] } as never,
      "api",
    );
    expect(read(t.ch, MarkBit.selection, 2)).toEqual([
      [[], [20]],
      [[[20, 24]], []],
    ]);
    // WT 30–32, deleted in the design: the junction between design 32 and 33 (WT 29 and 33).
    t.ctx.selection.set(
      { kind: "seq-loci", set: t.set, elements: [{ unit: t.unit("wt/unit"), indices: range(29, 32) }] } as never,
      "api",
    );
    expect(read(t.ch, MarkBit.selection, 2)).toEqual([
      [[[29, 32]], []],
      [[], [33]],
    ]);
    // The tag, beyond the pairs: nothing on the wild type.
    t.ctx.selection.set(
      { kind: "seq-loci", set: t.set, elements: [{ unit: t.unit("d/unit"), indices: range(41, 44) }] } as never,
      "api",
    );
    expect(read(t.ch, MarkBit.selection, 2)[0]).toEqual([[], []]);
  });

  it("a boundary item lights a caret where it is, and through the map on the other side", async () => {
    const t = await open();
    const a = t.set.annotations.findIndex((x) => x && t.set.strings.get(x.key) === "marks");
    t.ctx.selection.set(
      { kind: "item-loci", set: t.set, elements: [{ annotation: a, items: range(0, 1) }] } as never,
      "api",
    );
    expect(read(t.ch, MarkBit.selection, 2)).toEqual([
      [[], [10]],
      [[], [10]],
    ]);
  });

  it("the composer draws a caret as a bar on the edge between two cells", async () => {
    const t = await open();
    t.ctx.selection.set(
      { kind: "seq-loci", set: t.set, elements: [{ unit: t.unit("d/unit"), indices: range(20, 24) }] } as never,
      "api",
    );
    const be = t.frame(recording(), Dirty.all);
    const ax = t.ctx.layout.axis(t.section);
    const edge = Math.round(ax?.x(0, ax.shownOf(0, 20)) ?? 0);
    const rects = Array.from({ length: be.rectItems.n }, (_, i) => [be.rectItems.x[i], be.rectItems.w[i]]);
    expect(rects.some(([x, w]) => w === 3 && Math.abs(x + 1 - edge) <= 1)).toBe(true);
  });
});

const range = (s: number, t: number) => OrderedSet.ofRange(s, t);
