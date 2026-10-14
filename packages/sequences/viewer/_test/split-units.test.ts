/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A split track whose rows sit on units of their own (HR1): each row knows its unit,
// a highlight or hover on one unit marks its row alone (the section's other tracks as before), a run on the shared
// sequence every row, and a position hovered on a row is that unit's residue.

import { CoordinateSystemKind as K, OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Pass, SkinToken } from "../src/draw/index.ts";
import { ChannelMarks, Hit, HoverKind, SeqViewerConfig } from "../src/index.ts";
import type { SplitFeatures } from "../src/layout.ts";
import type { MarkSet } from "../src/marks.ts";
import { opened, preparedRenderers } from "./layout-fixtures.ts";
import { lit, passRects } from "./overlay-fixtures.ts";
import { CHAIN_AT, splitUnitsDoc } from "./split-units-fixtures.ts";

const S = "v.main";

async function open() {
  const t = await opened(splitUnitsDoc(), { renderers: preparedRenderers(), size: { [S]: [1000, 400] } });
  const set = () => t.ctx.graph.assembled.value?.set as NonNullable<typeof t.ctx.graph.assembled.value>["set"];
  const unit = (name: string) => {
    const s = { kind: 0, id: -1, row: -1 };
    set().byName(name, s);
    return s.id;
  };
  const ch = new ChannelMarks(t.ctx);
  const runs = (m: MarkSet | undefined) => {
    const out: number[][] = [];
    for (let r = m?.from(0) ?? 0; m && r < m.to(0); r++) out.push([m.start[r], m.stop[r]]);
    return out;
  };
  /** A highlight on unit `key`'s [s, t), then what the section, and the rows on A, B and C, light. */
  const lit = (key: string, s: number, e: number) => {
    const elements = [{ unit: unit(key), indices: OrderedSet.ofRange(s, e) }];
    t.ctx.highlight.set("test", { kind: "seq-loci", set: set(), elements } as never);
    const m = ch.marks(0);
    return {
      all: runs(m?.highlight),
      rows: ["A", "B", "C"].map((u) => runs(m?.ownOf(unit(u))?.highlight)),
    };
  };
  return { ...t, set, unit, ch, lit, runs };
}

describe("split rows on units of their own", () => {
  it("each row's unit, from its data's target", async () => {
    const t = await open();
    const k = t.ctx.graph.assembled.value?.view?.resolved.tracks.split[t.trackIndex("t-chains")] ?? -1;
    const state = t.ctx.graph.assembled.value?.view?.prepared[k]?.state as SplitFeatures;
    expect(state.split?.labels).toEqual(["A", "B", "C"]);
    expect(Array.from(state.split?.unit ?? [])).toEqual(["A", "B", "C"].map(t.unit));
    t.ch.sync();
    expect(Array.from(t.ch.marks(0)?.ownUnits ?? [])).toEqual(["A", "B", "C"].map(t.unit));
  });

  it("a run on unit B marks B's row alone (and the reference's tracks); one on the reference every row", async () => {
    const t = await open();
    const B = [[CHAIN_AT + 4, CHAIN_AT + 8]];
    expect(t.lit("B", 4, 8)).toEqual({ all: B, rows: [[], B, []] });
    expect(t.lit("ref/unit", 30, 34)).toEqual({ all: [[30, 34]], rows: [[[30, 34]], [[30, 34]], [[30, 34]]] });
  });

  it("the hover: a position on C's row is C's residue, lighting C's row alone; an item its own unit's", async () => {
    const t = await open();
    const hit = Hit.create();
    // Reference 47 is chain residue 37: C observed [5, 35) has no item there, so the hit is a position.
    t.ctx.layout.hitTest(S, t.xOf(S, 0, CHAIN_AT + 37), t.yOf(S, "t-chains", 2), hit);
    expect([hit.kind, hit.split, hit.seq.unit, hit.seq.index, hit.system.kind, hit.system.id]).toEqual([
      HoverKind.seq,
      2,
      t.unit("C"),
      37,
      K.unit,
      t.unit("C"),
    ]);
    t.ctx.hover.setFromHit(hit, "pointer");
    const m = t.ch.marks(0);
    expect(t.runs(m?.hover)).toEqual([[47, 48]]);
    expect(["A", "B", "C"].map((u) => t.runs(m?.ownOf(t.unit(u))?.hover))).toEqual([[], [], [[47, 48]]]);
    // On B's row over its item, `UI.OwnRowHover` "item": the item, on B.
    t.ctx.config.set(SeqViewerConfig.UI.OwnRowHover, "item");
    t.ctx.layout.hitTest(S, t.xOf(S, 0, CHAIN_AT + 5), t.yOf(S, "t-chains", 1), hit);
    expect([hit.kind, hit.split]).toEqual([HoverKind.item, 1]);
    t.ctx.hover.setFromHit(hit, "pointer");
    expect(["A", "B", "C"].map((u) => t.runs(t.ch.marks(0)?.ownOf(t.unit(u))?.hover).length)).toEqual([0, 2, 0]);
  });

  it("a position on the reference's own track stays the reference's", async () => {
    const t = await open();
    const hit = Hit.create();
    t.ctx.layout.hitTest(S, t.xOf(S, 0, 20), t.yOf(S, "t-letters"), hit);
    expect([hit.kind, hit.seq.unit, hit.seq.index]).toEqual([HoverKind.seq, t.unit("ref/unit"), 20]);
  });
});

describe("drawn: a split row on its own unit takes that unit's marks over its band", () => {
  it("a highlight on B: the column over every entry but A's and C's rows; one on the reference over all", async () => {
    const t = await lit(splitUnitsDoc(), { width: 1000, height: 400 });
    const l = t.layout(t.section);
    const band = (row: number) => {
      const y = t.yOf(t.section, "t-chains", row);
      const e = Array.from({ length: l.entryCount }, (_, k) => k).find(
        (k) => l.entryRow[k] === row && l.entryTrack[k] === t.trackIndex("t-chains"),
      ) as number;
      return [Math.round(y - l.entryHeight[e] / 2), Math.round(y + l.entryHeight[e] / 2)];
    };
    const [a, b, c] = [0, 1, 2].map(band);
    const ys = (unit: string, s: number, e: number) => {
      const u = t.set.units.key.findIndex((k) => t.set.strings.get(k) === unit);
      const elements = [{ unit: u, indices: OrderedSet.ofRange(s, e) }];
      t.ctx.highlight.set("host", { kind: "seq-loci", set: t.set, elements } as never);
      const rects = passRects(t.overlay(), Pass.overlay).filter((r) => r[4] === SkinToken.highlight);
      return rects.map((r) => [r[1], r[1] + r[3]]);
    };
    const onB = ys("B", 4, 8);
    expect(onB.some(([s, e]) => s <= b[0] && e >= b[1])).toBe(true);
    for (const row of [a, c]) expect(onB.filter(([s, e]) => s < row[1] && e > row[0])).toEqual([]);
    const onRef = ys("ref/unit", 30, 34);
    // Every row lit: the segments tile the column with no gap.
    const sorted = onRef.sort((x, y) => x[0] - y[0]);
    for (let k = 1; k < sorted.length; k++) expect(sorted[k][0]).toBe(sorted[k - 1][1]);
    expect(sorted[0][0]).toBeLessThanOrEqual(a[0]);
    expect(sorted[sorted.length - 1][1]).toBeGreaterThanOrEqual(c[1]);
  });
});
