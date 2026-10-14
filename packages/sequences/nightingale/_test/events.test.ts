/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What `mount` decides, without a DOM (events.ts and the axis edges): the position and ranges an event names, the
// arbiter, the pointer's hover deduped before it is located, the per-manager highlight dedupe (`lit`) and the merge of
// fetched values; the converter onto the axis (`AxisRuns`, `AxisSink`) and back (`extentOf`, `locate`) for blocks shown
// whole, in a `positions` subset, a ligand, reversed and across a cyclic origin.

import { layoutOf } from "@molstar/seq-core/format";
import { OrderedSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { AxisRuns, AxisSink, BlockRuns, extentOf, type Located, locate } from "../src/axis.js";
import {
  arbitrate,
  HoverResult,
  type HoverTarget,
  hoverAt,
  mergeSeries,
  newLight,
  nextHighlight,
  noteLit,
  positionOf,
  rangesOf,
} from "../src/events.js";
import type { AxisBlock, AxisMeta } from "../src/index.js";
import { axisView, type LabelBox, type Laid, layLabels } from "../src/labels.js";

const runs = (text: string) => AxisRuns.parse(text, AxisRuns.create());
const listed = (r: AxisRuns) => Array.from({ length: r.count }, (_, i) => [r.start[i], r.stop[i]]);
const block = (b: Partial<AxisBlock> & { block: number; offset: number; count: number }): AxisBlock => ({
  first: 0,
  ligand: false,
  shown: null,
  reversed: false,
  layout: null,
  unit: 0,
  members: null,
  alignment: -1,
  ...b,
});
const meta = (...blocks: AxisBlock[]): AxisMeta => ({
  section: 0,
  scrollGroup: null,
  length: blocks.reduce((n, b) => n + b.count, 0),
  blocks,
});
/** What `sink` puts on a block's axis for the block-system run [s, t). */
const through = (b: AxisBlock, ...ranges: [number, number][]) => {
  const sink = new AxisSink();
  sink.block = b;
  for (const [s, t] of ranges) sink.add(s, t);
  return listed(AxisRuns.finish(sink.runs));
};

describe("axis runs", () => {
  it("parse Nightingale's ranges (1-based, inclusive) and print them back; a range without a colon is skipped", () => {
    expect(listed(runs("3:5,9:9,12:20"))).toEqual([
      [2, 5],
      [8, 9],
      [11, 20],
    ]);
    expect(AxisRuns.text(runs("3:5,9:9,12:20"))).toBe("3:5,9:9,12:20");
    expect(listed(runs("7,2:3,0:4,5:4"))).toEqual([[1, 3]]);
    expect(AxisRuns.text(runs(""))).toBe("");
  });

  it("sort when added out of order and join runs that touch or overlap", () => {
    const r = AxisRuns.create(2);
    for (const [s, t] of [
      [10, 12],
      [0, 3],
      [3, 5],
      [11, 15],
      [20, 21],
    ])
      AxisRuns.add(r, s, t);
    expect(r.sorted).toBe(false);
    expect(listed(AxisRuns.finish(r))).toEqual([
      [0, 5],
      [10, 15],
      [20, 21],
    ]);
    expect(r.sorted).toBe(true);
  });

  it("sort many runs: the packed keys keep (start, length) exact", () => {
    const r = AxisRuns.create();
    const starts = Array.from({ length: 5_000 }, (_, i) => (i * 7919) % 5_000).map((x) => x * 3);
    for (const s of starts) AxisRuns.add(r, s, s + 2);
    AxisRuns.finish(r);
    expect(listed(r)).toEqual([...starts].sort((a, b) => a - b).map((s) => [s, s + 2]));
  });

  it("compare and copy by value", () => {
    const [a, b, c] = [runs("1:3,6:6"), runs("1:3,6:6"), AxisRuns.create()];
    expect(AxisRuns.equal(a, b)).toBe(true);
    AxisRuns.copy(c, a);
    expect(AxisRuns.equal(a, c)).toBe(true);
    expect(AxisRuns.equal(a, runs("1:3,6:7"))).toBe(false);
    expect(AxisRuns.equal(a, runs("1:3"))).toBe(false);
  });
});

describe("block-system runs onto the axis (AxisSink)", () => {
  it("a block shown whole: first + k, clipped to the block", () => {
    const b = block({ block: 0, offset: 10, count: 8, first: 5 });
    expect(through(b, [6, 9])).toEqual([[11, 14]]);
    expect(through(b, [0, 7], [12, 40])).toEqual([
      [10, 12],
      [17, 18],
    ]);
    expect(through(b, [0, 5], [13, 20])).toEqual([]);
  });

  it("a ligand: every index is its one position", () => {
    const b = block({ block: 1, offset: 4, count: 1, ligand: true });
    expect(through(b, [0, 1], [3, 9])).toEqual([[4, 5]]);
  });

  it("a positions subset: ranks, cut-out stretches abut so a run stays one run", () => {
    const b = block({ block: 0, offset: 2, count: 5, shown: OrderedSet.ofSortedArray(Int32Array.of(2, 3, 7, 8, 9)) });
    expect(through(b, [0, 4])).toEqual([[2, 4]]);
    expect(through(b, [3, 10])).toEqual([[3, 7]]);
    expect(through(b, [4, 7])).toEqual([]);
  });

  it("a reversed piece and one across a cyclic origin, through the layout map", () => {
    const system = { kind: 1, id: 0, row: -1 };
    const reversed = layoutOf(system, 2, 12, -1, 20).map;
    const r = block({ block: 0, offset: 0, count: 10, first: 11, layout: reversed });
    expect(through(r, [9, 12])).toEqual([[0, 3]]);
    expect(through(r, [2, 4])).toEqual([[8, 10]]);
    const wrapped = layoutOf(system, 16, 4, 1, 20).map;
    const w = block({ block: 0, offset: 5, count: 8, first: 16, layout: wrapped });
    expect(through(w, [18, 20])).toEqual([[7, 9]]);
    expect(through(w, [0, 2])).toEqual([[9, 11]]);
    expect(through(w, [19, 20], [0, 1])).toEqual([[8, 10]]);
  });
});

describe("the axis back to blocks and indices", () => {
  const first = block({ block: 3, offset: 0, count: 5, first: 10 });
  const second = block({ block: 4, offset: 5, count: 4, shown: OrderedSet.ofSortedArray(Int32Array.of(1, 4, 5, 9)) });
  const m = meta(first, second);
  const out: Located = { block: -1, index: -1 };

  it("locates by binary search over the blocks' offsets; none beyond the axis", () => {
    expect([0, 4, 5, 8].map((at) => (locate(m, at, out) ? [out.block, out.index] : null))).toEqual([
      [3, 10],
      [3, 14],
      [4, 1],
      [4, 9],
    ]);
    expect(locate(m, 9, out)).toBe(false);
    expect(locate(m, -1, out)).toBe(false);
  });

  it("locates in a reversed piece through its layout map", () => {
    const layout = layoutOf({ kind: 1, id: 0, row: -1 }, 2, 12, -1, 20).map;
    const r = meta(block({ block: 0, offset: 0, count: 10, first: 11, layout }));
    expect([0, 3, 9].map((at) => (locate(r, at, out) ? out.index : null))).toEqual([11, 8, 2]);
  });

  it("turns axis runs into block runs, across a block edge and splitting a positions subset where its indices jump", () => {
    const x = BlockRuns.create();
    expect(BlockRuns.toArray(extentOf(m, runs("4:7"), x))).toEqual([
      { block: 3, start: 13, stop: 15 },
      { block: 4, start: 1, stop: 2 },
      { block: 4, start: 4, stop: 5 },
    ]);
    // The same buffer again: it is reused.
    expect(extentOf(m, runs("1:1"), x)).toBe(x);
    expect(BlockRuns.toArray(x)).toEqual([{ block: 3, start: 10, stop: 11 }]);
    expect(extentOf(m, AxisRuns.create(), x).count).toBe(0);
  });

  it("turns axis runs into block runs through a layout map, every copy a position has", () => {
    const layout = layoutOf({ kind: 1, id: 0, row: -1 }, 2, 12, -1, 20).map;
    const r = meta(block({ block: 0, offset: 0, count: 10, first: 11, layout }));
    expect(BlockRuns.toArray(extentOf(r, runs("1:3"), BlockRuns.create()))).toEqual([{ block: 0, start: 9, stop: 12 }]);
  });
});

describe("what an event names", () => {
  it("positionOf: its highlight, a residue, a feature's start, or the linegraph's point by its series' name", () => {
    expect(positionOf({ highlight: "12:14,20:21" })).toBe(12);
    expect(positionOf({ type: "highlight", value: "7:9" })).toBe(7);
    expect(positionOf({ feature: { position: 4 } })).toBe(4);
    expect(positionOf({ feature: { start: 9 } })).toBe(9);
    expect(positionOf({ feature: { conservation: { position: 31, value: 0.5 } } }, "conservation")).toBe(31);
    expect(positionOf({ feature: { conservation: { position: 31, value: 0.5 } } })).toBe(31);
    expect(positionOf({ feature: { other: { position: 5 } } }, "conservation")).toBe(-1);
    expect(positionOf({})).toBe(-1);
  });

  it("rangesOf: our feature's fragments as numbers, else Nightingale's own highlight, else the one position", () => {
    const out = AxisRuns.create();
    const svs = { layer: 0 };
    const feature = {
      svs,
      locations: [
        {
          fragments: [
            { start: 8, end: 9 },
            { start: 12, end: 12 },
          ],
        },
      ],
    };
    expect(listed(rangesOf({ feature, highlight: "8:12" }, 8, out))).toEqual([
      [7, 9],
      [11, 12],
    ]);
    expect(listed(rangesOf({ highlight: "3:3,9:11" }, 3, out))).toEqual([
      [2, 3],
      [8, 11],
    ]);
    expect(listed(rangesOf({ feature: { locations: [] } }, 5, out))).toEqual([[4, 5]]);
    expect(rangesOf({}, -1, out).count).toBe(0);
    // A feature of Nightingale's own (no `svs`) names its highlight, not its fragments.
    expect(listed(rangesOf({ feature: { locations: feature.locations }, highlight: "1:2" }, 1, out))).toEqual([[0, 2]]);
  });
});

describe("the arbiter", () => {
  it("a component's hover wins only for an item; its leave only after one; a click changes neither", () => {
    const a = { onItem: false };
    expect(arbitrate(a, "hover", false)).toBe(false);
    expect(arbitrate(a, "leave", true)).toBe(false);
    expect(arbitrate(a, "hover", true)).toBe(true);
    expect(a.onItem).toBe(true);
    // A click on the hovered item must not hand hover back to the pointer.
    expect(arbitrate(a, "click", true)).toBe(true);
    expect(a.onItem).toBe(true);
    expect(arbitrate(a, "leave", true)).toBe(true);
    expect(a.onItem).toBe(false);
    expect(arbitrate(a, "click", true)).toBe(true);
    expect(a.onItem).toBe(false);
  });
});

describe("the pointer's hover", () => {
  const view = axisView({ start: 0, end: 0, left: 0, right: 0 }, 1, -1, 10, 120, 10, 10); // 10 px a position from x = 10
  const m = meta(block({ block: 2, offset: 0, count: 10, first: 3 }));
  const [manager, layer] = [{}, {}];
  const state = (): HoverTarget => ({ manager: null, position: -1, layer: null });
  const out: Located = { block: -1, index: -1 };

  it("reports a change once, then nothing for the same position, manager and layer element", () => {
    const s = state();
    expect(hoverAt(s, manager, m, view, 25, layer, out)).toBe(HoverResult.moved);
    expect([out.block, out.index]).toEqual([2, 4]);
    out.block = out.index = -9;
    expect(hoverAt(s, manager, m, view, 29, layer, out)).toBe(HoverResult.same);
    expect([out.block, out.index]).toEqual([-9, -9]); // nothing was located
    expect(hoverAt(s, manager, m, view, 35, layer, out)).toBe(HoverResult.moved);
    expect(out.index).toBe(5);
    expect(hoverAt(s, {}, m, view, 35, layer, out)).toBe(HoverResult.moved);
    expect(hoverAt(s, s.manager, m, view, 35, {}, out)).toBe(HoverResult.moved);
  });

  it("is off in a gutter and past the shown positions, leaving the last report alone", () => {
    const s = state();
    hoverAt(s, manager, m, view, 25, layer, out);
    expect(hoverAt(s, manager, m, view, 5, layer, out)).toBe(HoverResult.off);
    expect(hoverAt(s, manager, m, view, 111, layer, out)).toBe(HoverResult.off);
    expect(hoverAt(s, manager, m, view, 25, layer, out)).toBe(HoverResult.same);
    // A position the blocks do not hold (the axis is shorter than the track shows).
    const short = meta(block({ block: 2, offset: 0, count: 3, first: 3 }));
    expect(hoverAt(state(), manager, short, view, 75, layer, out)).toBe(HoverResult.off);
  });
});

describe("a manager's highlight (lit)", () => {
  it("builds a string only for a changed manager, and nothing for the same runs again", () => {
    const l = newLight();
    expect(nextHighlight(l, runs(""))).toBeUndefined();
    expect(nextHighlight(l, runs("3:5"))).toBe("3:5");
    expect(nextHighlight(l, runs("3:5"))).toBeUndefined();
    expect(nextHighlight(l, runs("3:5,9:9"))).toBe("3:5,9:9");
    expect(nextHighlight(l, runs(""))).toBeNull();
    expect(nextHighlight(l, runs(""))).toBeUndefined();
  });

  it("sends again after Nightingale showed its own hover, unless it shows the string already", () => {
    const l = newLight();
    nextHighlight(l, runs("3:5"));
    noteLit(l, "8:8"); // the manager's own hover replaced ours
    expect(nextHighlight(l, runs("3:5"))).toBe("3:5");
    noteLit(l, "3:5"); // ... and now shows exactly what we would send
    expect(nextHighlight(l, runs("3:5"))).toBeUndefined();
    expect(l.lit).toBe("3:5");
    noteLit(l, undefined);
    expect(nextHighlight(l, runs(""))).toBeUndefined();
  });
});

describe("fetched values", () => {
  it("fill only the positions with no value: inline values survive, a progressive block's positions are filled", () => {
    const series = [
      { position: 1, value: 0.5 }, // an inline block
      { position: 2, value: 0 },
      { position: 3, value: null }, // a progressive block
      { position: 4, value: null },
    ];
    const fetched = [0.9, 0.9, 0.3, null];
    expect(mergeSeries(series, fetched)).toEqual([
      { position: 1, value: 0.5 },
      { position: 2, value: 0 },
      { position: 3, value: 0.3 },
      { position: 4, value: null },
    ]);
  });
});

describe("labels laid at a view", () => {
  const style = () => ({ display: "none", left: "", width: "", top: "", height: "", lineHeight: "" });
  const laid = (start: number, end: number, text = "Domain"): Laid & { writes: number } => {
    const s = style();
    const box: LabelBox = { mode: "none", left: 0, width: 0 };
    return {
      feature: {},
      text,
      fragments: [{ start, end }],
      style: s,
      width: -1,
      cut: -1,
      box,
      top: 0,
      height: 0,
      visible: false,
      writes: 0,
    };
  };
  const layout = { getFeatureYPos: () => 3, getFeatureHeight: () => 12 };
  const view = axisView({ start: 0, end: 0, left: 0, right: 0 }, 1, -1, 100, 200, 0, 0); // 2 px a position

  it("measures a label once, writes its style when its box changes, and not again for the same view", () => {
    const labels = [laid(10, 40), laid(60, 61, "Hi")];
    const measured: string[] = [];
    const measure = (t: string) => measured.push(t) * 0 + t.length * 6;
    expect(layLabels(labels, view, layout, measure)).toBe(1); // the second does not fit: still hidden, nothing written
    expect(measured).toEqual(["Domain", "Dom…", "Hi"]);
    expect(labels[0].style).toMatchObject({ display: "", left: "18px", width: "62px", top: "3px", height: "12px" });
    expect(labels[1].style.display).toBe("none");
    // The manager's attribute storm on each highlight: the view is unchanged, so nothing is measured or written.
    expect([layLabels(labels, view, layout, measure), layLabels(labels, view, layout, measure)]).toEqual([0, 0]);
    expect(measured.length).toBe(3);
    // Zoomed to 55..70: the first leaves, the second fits.
    const zoomed = axisView({ start: 0, end: 0, left: 0, right: 0 }, 55, 70, 100, 200, 0, 0);
    expect(layLabels(labels, zoomed, layout, measure)).toBe(2);
    expect([labels[0].style.display, labels[1].style.display]).toEqual(["none", ""]);
    // A row that moved (a redrawn track) writes that label again.
    expect(layLabels(labels, zoomed, { ...layout, getFeatureYPos: () => 5 }, measure)).toBe(1);
    expect(labels[1].style.top).toBe("5px");
  });
});
