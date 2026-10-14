/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the channels' per-move path (IN3; seqspec:interaction#channels: after a route's first use nothing
// allocates; CP2's `spans.heap`), on RB0b's stepped `heapScaling`: a hover element change per step, set through the
// hover manager, then the channels' rebuild (closed, projected into every section's hover set, its primary) and the
// overlay's reads; and a million element changes over 10,000 distinct alignment rows and their positions (a row's
// residue on its own columns read without the row's cached map, which allocates on a miss). The positive control: each
// step allocating an object too (`leaking`), measured alike, is caught.

import { CoordinateSystemKind as K, type SeqSet } from "@molstar/seq-core/model";
import { beforeAll, describe, expect, it } from "vitest";
import { gene, rows } from "../../core/_test/format/project-fixtures.ts";
import { ChannelMarks, Hit, HitRegion, HoverKind, MarkBit } from "../src/index.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";
import { shown } from "./projection-fixtures.ts";

let sink = 0;
let g: Awaited<ReturnType<typeof shown>>;
let r: Awaited<ReturnType<typeof shown>>;
const ROWS = 10_000;

beforeAll(async () => {
  [g, r] = await Promise.all([shown(gene({ chromosome: true })), shown(rows(ROWS))]);
  // The channels alone: the hover label's churn per element change is IN2's guard (input.heap.ts), not this one.
  for (const t of [g, r]) t.ctx.behaviors.get("hoverLabel")?.remove();
}, 60_000);

/** A hit on a residue of `unit`, reused. */
function seqHit(set: SeqSet, unit: number) {
  const hit = Hit.create();
  Object.assign(hit, { section: 0, region: HitRegion.data, kind: HoverKind.seq, block: 0 });
  Object.assign(hit.seq, { set, unit });
  Object.assign(hit.system, { kind: K.unit, id: unit, row: -1 });
  return hit;
}

/** A nucleotide or a residue hovered by turns, the channels synced, a section's marks read. */
function geneStep(ch: ChannelMarks) {
  const [gu, pu] = [g.unit("gene/unit"), g.unit("protein/unit")];
  const [nt, aa] = [seqHit(g.set, gu), seqHit(g.set, pu)];
  return (i: number) => {
    const h = i & 1 ? aa : nt;
    const x = i & 1 ? i % 20 : i % 120;
    h.index = h.shown = x;
    h.seq.index = x as never;
    g.ctx.hover.setFromHit(h, "pointer");
    ch.sync();
    const m = ch.marks(i % 3);
    sink += (m?.bitsAt(0, x) ?? 0) + (m?.primary(0) ?? 0) + (m?.any(0, 40, 80, MarkBit.hover) ? 1 : 0);
  };
}

describe("the channels' hover allocates nothing per element change", () => {
  it("a nucleotide or a residue hovered, closed, projected into every section with its primary, and read", () => {
    const ch = new ChannelMarks(g.ctx);
    const result = heapScaling(50_000, geneStep(ch));
    expect(ch.rebuilds.hover).toBeGreaterThan(200_000);
    expect(sink).not.toBe(0);
    expect(isFlat(result), describeScaling(result, 50_000)).toBe(true);
  });

  it(`a million element changes over ${ROWS} distinct rows and positions of an alignment`, () => {
    const ch = new ChannelMarks(r.ctx);
    const hit = Hit.create();
    Object.assign(hit, { section: 0, region: HitRegion.data, kind: HoverKind.aligned, block: 0 });
    Object.assign(hit.aligned, { set: r.set, alignment: 0 });
    Object.assign(hit.system, { kind: K.row, id: 0, row: 0 });
    const result = heapScaling(250_000, (i) => {
      const row = (i * 7919) % ROWS; // every row, in no order
      const column = i % 5;
      hit.dataRow = hit.system.row = row;
      hit.index = hit.shown = column;
      hit.aligned.row = row as never;
      hit.aligned.column = column as never;
      r.ctx.hover.setFromHit(hit, "pointer");
      ch.sync();
      sink += ch.marks(0)?.hover.count ?? 0;
    });
    expect(ch.rebuilds.hover).toBeGreaterThan(1_000_000);
    expect(sink).not.toBe(0);
    expect(isFlat(result), describeScaling(result, 250_000)).toBe(true);
  });

  it("the same steps allocating an object each are caught (positive control)", () => {
    const ch = new ChannelMarks(g.ctx);
    expect(controlCaught(() => heapScaling(50_000, leaking(geneStep(ch))))).toBe(true);
  });
});
