/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout (LY1; seqspec:model#pieces-and-layout-maps; seqspec:viewer#state): entries side by side, stacked and
// wrapped, lazy rows, row sizing, height modes, groups and clades, `row_subset`, lanes, a filter view through
// `positions`, rebuilds with no surface mounted. LY1's own: each stacked block's entries, a wrapped line's band from
// the features on its block coordinate system runs, no lines before a width, the rows' version kept across a resize, an
// unknown data column type, a layout kept across a removal and return.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { case3, case6 } from "../../core/_test/format/cases.ts";
import type { ShownIndex } from "../src/axis.ts";
import { DefaultRenderProps, EntryKind, KernelRenderers, type LayerRenderer, type RenderProps } from "../src/index.ts";
import { SectionLayouts } from "../src/layout.ts";
import { TextColumn } from "../src/renderers/columns.ts";
import { set } from "./axis-fixtures.ts";
import { familyWithTree } from "./families.ts";
import { genotypeDoc } from "./genotype-fixtures.ts";
import { laneRanges, lazyLoaders, opened, renderers } from "./layout-fixtures.ts";

const [T, G, B] = [EntryKind.track, EntryKind.group, EntryKind.title];
const k = (n: number) => n as ShownIndex;

describe("entries", () => {
  it("side by side (case 6): the ruler pinned, a row per track, the ligand one position over `ligandCells`", async () => {
    const t = await opened(case6());
    const l = t.layout("chains.main");
    const ax = t.vp("chains.main").axis;
    expect(t.entries("chains.main").map(([kind, , , top, h]) => [kind, top, h])).toEqual([
      [T, 20, 30],
      [T, 0, 22],
      [T, 22, 48],
    ]);
    expect([l.pinnedCount, l.navHeight, l.scrolledTop, l.dataLeft, l.dataRight]).toEqual([1, 20, 50, 158, 992]);
    expect([Array.from(ax?.shown ?? []), Array.from(ax?.span ?? []), Array.from(ax?.offset ?? [])]).toEqual([
      [274, 1],
      [1, 3],
      [0, 274],
    ]);
    expect([ax?.gap, ax?.cells, ax?.insetX, ax?.pinnedCount]).toEqual([28, 277, 166, 1]);
    // Fitted: the data width less the gap over 277 cells.
    expect(t.vp("chains.main").live.current.cell).toBeCloseTo((1000 - 166 - 28) / 277);
  });

  it("stacked, unequal blocks: a title per block, every track per block, one scale; each block's own entries", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEKMKVLATGHEK" });
    b.data.sequence({ ref: "b", value: "MKVLA" });
    const v = b.view({ ref: "v", label: "V", arrangement: "stacked" });
    v.track({ ref: "v.ruler" }).ruler();
    v.track({ ref: "v.letters" }).letters();
    const t = await opened(b.getDocument());
    const rows = (block: number, y: number) => [
      [B, `block ${block}`, -1, y, 22],
      [T, "v.ruler", -1, y + 22, 30],
      [T, "v.letters", -1, y + 52, 22],
    ];
    expect(t.entries("v.main")).toEqual([...rows(0, 0), ...rows(1, 74 + 14)]);
    const l = t.layout("v.main");
    expect([l.pinnedCount, Array.from(l.entryBlock.subarray(0, 6))]).toEqual([0, [0, 0, 0, 1, 1, 1]]);
    // A block's marks stay within its title and lane, never the other block's.
    expect(Array.from(l.blockEntries.subarray(0, 4))).toEqual([0, 3, 3, 6]);
    const ax = t.vp("v.main").axis;
    expect([ax?.cells, Array.from(ax?.offset ?? []), ax?.gap]).toEqual([20, [0, 0], 0]);
    expect(t.vp("v.main").live.current.cell).toBe(28);
  });

  it("side by side, every block draws across every scrolled entry", async () => {
    const t = await opened(case6());
    const l = t.layout("chains.main");
    expect(Array.from(l.blockEntries.subarray(0, 4))).toEqual([1, 3, 1, 3]);
  });

  it("wrapped: three blocks each from a new line, a filter's break cuts a group, a ligand spans its cells", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEKMKVLATGHEKMKVLA" });
    const s = b.data.sequence({ ref: "b", value: "MKVLATGHEKMK" });
    const g = b.data.ligand({ ref: "lig", value: ["C", "C", "N", "O"], chem_comp_ids: ["XYZ"] });
    b.data.annotation({ ref: "dom", geometry: "ranges", target: a, value: { start: [0, 10], stop: [5, 25] } });
    b.data.annotation({ ref: "all-b", geometry: "ranges", target: s, value: { start: [0], stop: [12] } });
    b.data.annotation({ ref: "all-l", geometry: "ranges", target: g, value: { start: [0], stop: [4] } });
    const v = b.view({ ref: "v", label: "V", layout: "text", positions: ["dom", "all-b", "all-l"] });
    v.track({ ref: "v.letters" }).letters();
    v.track({ ref: "v.dom" }).ranges({ data: "dom" });
    // Without group labels (their band and placements: text.test.ts).
    const t = await opened(b.getDocument(), {
      size: { "v.main": [200, 400] },
      render: { gutterRight: 8, wrappedGroupLabels: "none", wrappedFontPx: 15 },
    });
    const l = t.layout("v.main");
    const lines = l.lines;
    // Mono 15 px: 9 px cells. a: groups of 5 (the break after block index 4), 10, 5; b: 10, then its short last group
    // (2) on a line of its own, where a full group would not fit; the ligand 3 cells. Lines `wrappedGap` (4) apart, a
    // block's first `stackGap` (14).
    expect(t.vp("v.main").axis?.fixedCell).toBe(9);
    expect([
      lines.count,
      ...[0, 1, 2, 3, 4, 5].map((n) => [lines.block[n], lines.start[n], lines.stop[n], lines.top[n]]),
    ]).toEqual([6, [0, 0, 5, 8], [0, 5, 15, 58], [0, 15, 20, 108], [1, 0, 10, 168], [1, 10, 12, 218], [2, 0, 1, 278]]);
    // The text's inset (8 px) above the first line, left of the gutter and right of the text.
    expect(t.entries("v.main").slice(0, 4)).toEqual([
      [T, "v.letters", 0, 8, 24],
      [T, "v.dom", 0, 32, 22],
      [T, "v.letters", 1, 58, 24],
      [T, "v.dom", 1, 82, 22],
    ]);
    expect([l.entryCount, l.dataLeft, l.dataRight]).toEqual([12, 64, 192]);
    // Each block's lines' entries.
    expect(Array.from(l.blockEntries.subarray(0, 6))).toEqual([0, 6, 6, 10, 10, 12]);
    const runs = t.layouts.stats;
    expect([runs.prepareRuns, runs.layoutRuns]).toEqual([3, 1]);
    // Wider: one line per block; the runs are not prepared again.
    t.vp("v.main").setSize(600, 400, 1);
    expect(t.layout("v.main").lines.count).toBe(3);
    expect([runs.prepareRuns, runs.layoutRuns]).toEqual([3, 2]);
    // The ligand is one position, 3 cells wide, where the axis view draws it.
    const view = t.ctx.layout.axis("v.main");
    expect(view && [view.x(2, k(0)), view.x(2, k(1))]).toEqual([64, 64 + 27]);
  });

  it("wrapped: the gutter takes its width only where it draws the lines' numbers (no group labels)", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEKMKVLATGHEKMKVLA" });
    b.view({ ref: "v", label: "V", layout: "text" }).track({ ref: "v.letters" }).letters();
    const left = async (render: Partial<RenderProps>, style?: string) => {
      const t = await opened(b.getDocument(), { size: { "v.main": [200, 400] }, render });
      if (style) t.ctx.ui.setStyle("v.main", style);
      return t.layout("v.main").dataLeft;
    };
    // The text's inset (8 px), and the gutter inside it.
    expect(await left({})).toBe(8);
    expect(await left({ wrappedGroupLabels: "inlineStart" })).toBe(8);
    expect(await left({ wrappedGroupLabels: "none" })).toBe(64);
    expect(await left({}, "compact")).toBe(8);
  });

  it("a wrapped layout before the viewport has a width has no lines (not one per position)", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(3500) });
    b.view({ ref: "v", label: "V", layout: "text" }).track({ ref: "v.letters" }).letters();
    const t = await opened(b.getDocument(), { size: { "v.main": [0, 400] } });
    expect([t.layout("v.main").lines.count, t.layout("v.main").entryCount, t.layouts.stats.layoutRuns]).toEqual([
      0, 0, 0,
    ]);
    t.vp("v.main").setSize(400, 400, 1);
    expect(t.layout("v.main").lines.count).toBeGreaterThan(100);
  });
});

describe("a wrapped line's band: the features on the line's block coordinate system runs", () => {
  // 40 residues filtered to [0, 5) and [20, 35); features at block indices [2, 4), and [22, 26) over [24, 30) (two lanes).
  const doc = () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(4) });
    b.data.annotation({ ref: "shown", geometry: "ranges", target: a, value: { start: [0, 20], stop: [5, 35] } });
    b.data.annotation({ ref: "f", geometry: "ranges", target: a, value: { start: [2, 22, 24], stop: [4, 26, 30] } });
    const v = b.view({ ref: "v", label: "V", layout: "text", positions: "shown" });
    v.track({ ref: "v.f" }).ranges({ data: "f" });
    return b.getDocument();
  };
  const R = () => renderers().map((r) => (r.name === "ranges" ? laneRanges() : r));

  it("lines in shown positions, features in block indices: each line counts its own lanes", async () => {
    // Groups of 5 shown positions, 45 px each: one group a line in 60 px (and the text's inset, 8 px a side).
    const t = await opened(doc(), {
      renderers: R(),
      size: { "v.main": [76, 400] },
      render: { wrappedGroupLabels: "none", wrappedGutter: 0, gutterRight: 0, wrappedGroup: 5 },
    });
    const l = t.layout("v.main");
    const lines = Array.from({ length: l.lines.count }, (_, n) => [l.lines.start[n], l.lines.stop[n]]);
    expect(lines).toEqual([
      [0, 5],
      [5, 10],
      [10, 15],
      [15, 20],
    ]);
    // Shown [5, 10) is block [20, 25) and [10, 15) is [25, 30): both lanes on each; [15, 20) none. The old comparison
    // of shown bounds with block indices found nothing on lines 1 and 2. A lane is 5 px (`wrappedLaneHeight`), and 2
    // above the first (`wrappedTrackGap`).
    expect(t.entries("v.main").map(([, , line, , h]) => [line, h])).toEqual([
      [0, 7],
      [1, 12],
      [2, 12],
      [3, 2],
    ]);
  });
});

describe("rows", () => {
  // Nothing in the layout is per row.
  it("lazy rows are one entry; the visible rows follow the scroll", async () => {
    const [uri, n] = ["https://example.org/huge", 1_000_000];
    const b = createSeqViewSpec();
    const aln = b.data.alignment({
      ref: "huge",
      row_count: n,
      column_count: 200,
      source: { uri, protocol: "seqstar-tiles/1" },
    });
    const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
    v.track({ ref: "v.rows" }).rows();
    const t = await opened(b.getDocument(), { loaders: lazyLoaders(uri, n, 200) });
    const l = t.layout("v.main");
    expect([l.entryCount, l.entryRun[0], l.entryHeight[0], l.rows.lazy, l.rows.count]).toEqual([1, n, 18 * n, true, n]);
    const vp = t.vp("v.main");
    vp.scrollTo(0, 18 * 50_000, "api");
    const r = t.ctx.layout.visibleRows("v.main")?.current;
    const shown = Math.ceil((400 - l.scrolledTop) / 18);
    expect(r && [r.count, r.displayRow[0], r.dataRow[0], r.entry[0], r.top[0]]).toEqual([
      shown,
      50_000,
      50_000,
      0,
      18 * 50_000,
    ]);
    expect([vp.live.current.firstEntry, vp.live.current.lastEntry]).toEqual([0, 0]);
  });

  it("collapsed clades fold into a wedge row and survive an unrelated patch; row_subset; lanes on the members", async () => {
    const t = await opened(case3(), { size: { "family.main": [1000, 400] } });
    const l = () => t.layout("family.main");
    const rows = () => {
      const r = t.ctx.layout.visibleRows("family.main")?.current;
      return r ? Array.from(r.dataRow.subarray(0, r.count)) : [];
    };
    // Tree order ((SRC, GRB2), (PTN11, P85A)): rows 0, 2, 1, 3; every track but the members' pinned (column data).
    expect([l().pinnedCount, l().rows.count, rows()]).toEqual([5, 4, [0, 2, 1, 3]]);
    // Lanes: the members' phosphosites (rows 0 and 2) draw on those rows' entries.
    const a = t.ctx.graph.assembled.value;
    const v = a?.view?.resolved;
    if (!v || !a) throw new Error("no view");
    const points = Array.from(v.layers.type).findLastIndex((x) => x >= 0 && a.set.strings.get(x) === "points");
    const lanes = Array.from(
      { length: v.layers.placementCount[points] },
      (_, n) => v.placements.lane[v.layers.placementStart[points] + n],
    );
    const at = lanes.map((lane) => {
      const d = SectionLayouts.displayOfLane(l().rows, lane);
      const e = SectionLayouts.entryOfDisplay(l().rows, d);
      return [v.lanes.row[lane], d, l().entryRow[e], l().entryTrack[e]];
    });
    // Row 1's site is the layer's filter's to drop (`kept`): it is placed all the same.
    expect(at).toEqual([
      [0, 0, 0, l().rows.track],
      [2, 1, 1, l().rows.track],
      [1, 2, 2, l().rows.track],
    ]);
    // Collapse (SRC, GRB2): node 1 in preorder.
    t.ctx.ui.toggleClade("sh2", 1);
    expect([l().rows.count, rows(), Array.from(l().rows.node ?? [])]).toEqual([3, [-1, 1, 3], [1, -1, -1]]);
    const r = l().rows;
    expect([1, 2].map((n) => SectionLayouts.displayOfLane(r, r.laneStart + n))).toEqual([-1, 1]);
    await t.ctx.applyPatch(set("family.main/2", { hidden: true }));
    expect([l().rows.count, rows()]).toEqual([3, [-1, 1, 3]]);
    // row_subset: only rows 0 and 2, in tree order.
    t.ctx.ui.toggleClade("sh2", 1);
    const members = v.layers.ref[l().rows.layer];
    await t.ctx.applyPatch(set(a.set.strings.get(members), { row_subset: [2, 0] }));
    expect([l().rows.count, rows()]).toEqual([2, [0, 2]]);
  });

  it("a resize or a row zoom keeps the rows' version (the tree is not laid out again); a fold bumps it", async () => {
    const t = await opened(case3(), { size: { "family.main": [1000, 400] } });
    const version = () => t.layout("family.main").rows.version;
    const v0 = version();
    t.vp("family.main").setSize(800, 300, 1);
    t.vp("family.main").setRowScale(0.5, 0);
    expect(version()).toBe(v0);
    t.ctx.ui.toggleClade("sh2", 1);
    expect(version()).toBe(v0 + 1);
    t.vp("family.main").setSize(900, 300, 1);
    expect([version(), t.layout("family.main").rows.count]).toEqual([v0 + 1, 3]);
    t.ctx.ui.toggleClade("sh2", 1);
    expect([version(), t.layout("family.main").rows.count]).toEqual([v0 + 2, 4]);
  });

  it("groups: a header entry, its children while expanded, nothing under a hidden one", async () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEK" });
    const v = b.view({ ref: "v", label: "V" });
    v.track({ ref: "v.letters" }).letters();
    const g = v.group({ ref: "v.g", label: "Group", collapsed: true });
    g.track({ ref: "v.g.one" }).letters();
    const inner = g.group({ ref: "v.g.inner", label: "Inner" });
    inner.track({ ref: "v.g.two" }).letters();
    const t = await opened(b.getDocument());
    const names = () => t.entries("v.main").map(([kind, n]) => [kind, n]);
    expect(names()).toEqual([
      [T, "v.letters"],
      [G, "v.g"],
    ]);
    t.ctx.ui.toggleGroup("v.g");
    expect(names()).toEqual([
      [T, "v.letters"],
      [G, "v.g"],
      [T, "v.g.one"],
      [G, "v.g.inner"],
      [T, "v.g.two"],
    ]);
    await t.ctx.applyPatch(set("v.g.inner", { hidden: true }));
    expect(names()).toEqual([
      [T, "v.letters"],
      [G, "v.g"],
      [T, "v.g.one"],
    ]);
  });

  it("an unknown data column type is laid out as `text`, reported once per section and type", async () => {
    const t = await opened(case3());
    const issues: string[] = [];
    t.ctx.events.issue.subscribe((i) => issues.push(`${i.code} ${i.message}`));
    const m = t.ctx.graph.assembled.value?.view;
    const rows = t.layout("family.main").rows.layer;
    if (!m) throw new Error("no view");
    const ref = m.set.strings.get(m.resolved.layers.ref[rows]);
    const cols = m.resolved.layers.params[rows].data_columns as Record<string, unknown>[];
    await t.ctx.applyPatch(set(ref, { data_columns: [{ ...cols[0], type: "org.example:sparkline" }] }));
    const l = t.layout("family.main");
    expect([l.columnCount, l.columnWidth[0]]).toEqual([1, TextColumn.width]);
    t.vp("family.main").setSize(800, 400, 1);
    await t.ctx.applyPatch(set("family.main", { label: "Family!" }));
    t.layout("family.main");
    expect(issues).toEqual([
      'svs.render.fallback section "family.main": no data column type "org.example:sparkline"; drawn as "text"',
    ]);
  });
});

describe("data columns yield to a narrow width", () => {
  /** A resize sweep: per width, [columns shown, data left, data area (with the overview)]. */
  const sweep = async (render?: Partial<RenderProps>) => {
    const t = await opened(familyWithTree(160, 20, { columns: true, sites: true }), { render });
    const s = "fam.main";
    const at = (w: number) => {
      t.vp(s).setSize(w, 400, 1);
      const l = t.layout(s);
      return [l.columnCount, l.dataLeft, l.dataRight - l.dataLeft];
    };
    return Object.assign(at, { layout: () => t.layout(s), axis: () => t.vp(s).axis });
  };
  // The fixture presets: a 150 px header, columns 64, 84, 110 and 100 px, an 8 px gutter each side, a 24 px overview:
  // the data area's left at 158, 222, 306, 416 and 516 px with 0 to 4 columns, 32 px right of it.
  const widths = [1000, 700, 600, 560, 400, 300, 1000];

  it("from the last, while the data area would be under `dataFloor` (64 px); back as the width returns; the header never", async () => {
    expect(DefaultRenderProps.dataFloor).toBe(64);
    const at = await sweep();
    expect(at(1000)).toEqual([4, 516, 452]);
    // 600: all four leave 52 px; without the last, 152.
    expect(at(600)).toEqual([3, 416, 152]);
    // 560: three leave 112 px, over 64 (the old 120 px floor took the third).
    expect(at(560)).toEqual([3, 416, 112]);
    // A shrink past the old data area still lays out (the viewport is not left hidden at the old width).
    expect(at(400)).toEqual([1, 222, 146]);
    // None fits: the header alone, whatever data area is left.
    expect(at(300)).toEqual([0, 158, 110]);
    expect(widths.map((w) => at(w)[0])).toEqual([4, 4, 3, 3, 1, 0, 4]);
    for (const w of [1000, 760, 640, 500, 420]) expect(at(w)[2]).toBeGreaterThanOrEqual(64);
  });

  it("at a custom floor, from the last while the data area would be under it", async () => {
    const at = await sweep({ dataFloor: 200 });
    expect(widths.map((w) => at(w)[0])).toEqual([4, 3, 2, 2, 0, 0, 4]);
    expect(at(700)).toEqual([3, 416, 252]);
  });

  it("at 0, never: a section narrower than its header and columns has no data area, its view waits", async () => {
    const at = await sweep({ dataFloor: 0 });
    const t = () => at.layout();
    expect(widths.map((w) => at(w)[0])).toEqual([4, 4, 4, 4, 4, 4, 4]);
    expect(at(600)).toEqual([4, 516, 52]);
    // At 400 the columns end past the width: no data area, so the view is hidden and the layout waits at 600.
    expect([...at(400), t().width]).toEqual([4, 516, 52, 600]);
    at(300);
    expect([at.axis()?.columnsShown, at.axis()?.left, at.axis()?.lefts?.length]).toEqual([4, 516, 5]);
    expect([...at(1000), t().width]).toEqual([4, 516, 452, 1000]);
  });
});

describe("the data columns' heads band", () => {
  // Case 3's family: five pinned column tracks of 18 px under a 20 px navigator, the species column.
  const band = async (render?: Partial<RenderProps>) => {
    const t = await opened(case3(), { size: { "family.main": [1000, 300] }, render });
    const l = t.layout("family.main");
    const last = l.pinnedCount - 1;
    return [
      l.pinnedCount,
      l.columnCount,
      l.entryTop[last] + l.entryHeight[last],
      l.headsTop,
      l.headsHeight,
      l.scrolledTop,
    ];
  };

  it("`shared` (the default): the last pinned entry's line, no row of its own", async () => {
    expect(await band()).toEqual([5, 1, 110, 92, 18, 110]);
  });

  it("`ownRow`: a row of its own directly above the first scrolled entry, pinned entries or none", async () => {
    expect(await band({ dataColumnHeads: "ownRow" })).toEqual([5, 1, 110, 110, 18, 128]);
  });
});

describe("sizes", () => {
  const doc = (params: Record<string, unknown> = {}) => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEKMKVLATGHEK" });
    b.data.annotation({ ref: "dom", geometry: "ranges", target: a, value: { start: [0], stop: [4] } });
    b.data.annotation({ ref: "val", geometry: "values", target: a, value: { index: [0, 1], value: [1, 2] } });
    const v = b.view({ ref: "v", label: "V", ...params });
    v.track({ ref: "v.ruler" }).ruler();
    v.track({ ref: "v.letters" }).letters();
    v.track({ ref: "v.dom" }).ranges({ data: "dom" });
    v.track({ ref: "v.val" }).values({ data: "val" });
    v.track({ ref: "v.tall", height: 2 }).ranges({ data: "dom" });
    return b.getDocument();
  };
  const heights = (t: Awaited<ReturnType<typeof opened>>) => t.entries("v.main").map(([, n, , , h]) => [n, h]);

  it("uniform rows share the row height; natural gives each type its own; a track's height wins in both", async () => {
    const t = await opened(doc());
    expect(heights(t)).toEqual([
      ["v.ruler", 30],
      ["v.letters", 22],
      ["v.dom", 22],
      ["v.val", 22],
      ["v.tall", 44],
    ]);
    await t.ctx.applyPatch(set("v.main", { row_sizing: "natural" }));
    expect(heights(t)).toEqual([
      ["v.ruler", 30],
      ["v.letters", 26],
      ["v.dom", 14],
      ["v.val", 40],
      ["v.tall", 44],
    ]);
  });

  it("`content` is as tall as its rows plus 8 px; `box` prefers `height_rows` rows (no margin: its rows scroll to their last)", async () => {
    const t = await opened(doc({ height: "content" }));
    const l = t.layout("v.main");
    expect(t.ctx.layout.contentHeight("v.main")?.value).toBe(20 + 30 + 22 * 3 + 44 + 8);
    expect([l.heightMode, l.minRowScale, l.scrolledHeight]).toEqual([0, 1, 22 * 3 + 44 + 8]);
    const box = await opened(doc({ height: "box", height_rows: 2 }));
    expect(box.ctx.layout.contentHeight("v.main")?.value).toBe(20 + 30 + 2 * 22);
  });
});

describe("navigator and overview geometry", () => {
  it("the navigator's window across the data area, the overview's down its column, for the gestures", async () => {
    const t = await opened(case3(), { size: { "family.main": [1000, 150] } });
    const vp = t.vp("family.main");
    vp.zoomTo(0, 0, 11, { origin: "api" });
    vp.scrollTo(0, 18, "api");
    const v = t.ctx.layout.axis("family.main");
    const l = t.layout("family.main");
    const s = vp.live.current;
    if (!v) throw new Error("no view");
    const w = l.dataRight - l.dataLeft;
    expect([l.navHeight, l.vnavWidth]).toEqual([20, 24]);
    expect(v.navX0).toBeCloseTo(l.dataLeft + (s.scrollX / s.contentWidth) * w);
    expect(v.navX1).toBeCloseTo(l.dataLeft + ((s.scrollX + w) / s.contentWidth) * w);
    expect(v.navX1 - v.navX0).toBeCloseTo((w * w) / s.contentWidth);
    const rh = l.height - l.scrolledTop;
    expect([v.vnavY0, v.vnavY1]).toEqual([
      l.scrolledTop + (18 / l.scrolledHeight) * rh,
      l.scrolledTop + ((18 + rh) / l.scrolledHeight) * rh,
    ]);
  });
});

describe("the layout's lifecycle", () => {
  it("rebuilds with no surface mounted: on a scheduler frame, once per change", async () => {
    const t = await opened(case6());
    const stats = t.layouts.stats;
    // Nobody read it: a frame builds it.
    expect(stats.builds).toBe(0);
    t.scheduler.advance(16);
    const built = stats.builds;
    expect(built).toBeGreaterThan(0);
    expect(t.layout("chains.main").entryCount).toBe(3);
    t.scheduler.advance(16);
    expect(stats.builds).toBe(built);
    // A change: rebuilt on the next frame, or the first read before it; once.
    t.vp("chains.main").setSize(800, 300, 1);
    t.scheduler.advance(16);
    expect([t.layout("chains.main").width, stats.builds]).toEqual([800, built + 1]);
    t.ctx.ui.toggleGroup("nothing");
    expect(t.layout("chains.main").width).toBe(800);
    t.scheduler.advance(16);
    expect(stats.builds).toBe(built * 2 + 1);
  });

  it("a row zoom rescales the rows without a rebuild: the entries a rebuild at that scale makes (KS2)", async () => {
    // Rows, then a track below them that moves with their height.
    const b = createSeqViewSpec({ title: "Zoom" });
    const aln = b.data.alignment({ ref: "aln", value: Array.from({ length: 12 }, () => "ACDEFGHIKL") });
    const f = b.data.annotation({
      ref: "f",
      geometry: "ranges",
      target: aln,
      value: { start: [1], stop: [4] },
    } as never);
    const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
    v.track({ label: "Rows" }).rows();
    v.group({ ref: "v.below", label: "Below", collapsed: false })
      .track({ label: "Below" })
      .ranges({ data: f } as never);
    const S = "v.main";
    const t = await opened(b.getDocument());
    t.vp(S).setSize(1000, 150, 1);
    const stats = t.layouts.stats;
    const snap = () => {
      const l = t.layout(S);
      const n = l.entryCount;
      return {
        n,
        top: Array.from(l.entryTop.subarray(0, n)),
        height: Array.from(l.entryHeight.subarray(0, n)),
        scrolled: l.scrolledHeight,
        rows: l.rows.height,
        vscroll: l.vscroll,
        // What a host sizes a section's frame from (a box without `height_rows`: its content).
        content: t.ctx.layout.contentHeight(S)?.value,
      };
    };
    const before = snap();
    expect(t.layout(S).rows.count).toBeGreaterThan(1);
    expect(t.layout(S).rows.entry + t.layout(S).rows.count).toBeLessThan(before.n);
    let builds = stats.builds;
    const rescales = stats.rescales;
    const r0 = t.layout(S).rows.height;
    for (const scale of [0.8, 0.55, 0.7]) {
      t.vp(S).setRowView(scale, 0, "api");
      const after = snap();
      expect(after.rows).toBeCloseTo(r0 * scale, 9);
      expect(stats.builds).toBe(builds);
      // A rebuild at the same scale: the same entries.
      t.ctx.ui.toggleGroup("nothing");
      const rebuilt = snap();
      expect(stats.builds).toBe(builds + 1);
      builds++;
      expect(rebuilt.n).toBe(after.n);
      for (let k = 0; k < after.n; k++) {
        expect(after.top[k]).toBeCloseTo(rebuilt.top[k], 6);
        expect(after.height[k]).toBeCloseTo(rebuilt.height[k], 6);
      }
      expect([after.scrolled, after.vscroll]).toEqual([expect.closeTo(rebuilt.scrolled, 6), rebuilt.vscroll]);
      expect(after.content).toBeCloseTo(rebuilt.content as number, 6);
      expect(after.content).not.toBe(before.content);
    }
    expect(stats.rescales - rescales).toBe(3);
  });

  it("a row zoom rescales a dense matrix's track and the tracks below it: the entries a rebuild makes (KS2)", async () => {
    const S = "v.main";
    const doc = genotypeDoc({ samples: 3000, positions: 40, below: 6 });
    const t = await opened(doc, { renderers: [...KernelRenderers] as LayerRenderer[] });
    t.vp(S).setSize(1000, 200, 1);
    const stats = t.layouts.stats;
    const snap = () => {
      const l = t.layout(S);
      const n = l.entryCount;
      return {
        n,
        top: Array.from(l.entryTop.subarray(0, n)),
        height: Array.from(l.entryHeight.subarray(0, n)),
        scrolled: l.scrolledHeight,
        vscroll: l.vscroll,
        content: t.ctx.layout.contentHeight(S)?.value,
      };
    };
    const matrix = () => {
      const l = t.layout(S);
      for (let k = 0; k < l.entryCount; k++) if (l.entryLayer[k] === 0 && l.entryKind[k] === EntryKind.track) return k;
      return -1;
    };
    const before = snap();
    const m = matrix();
    expect(before.height[m]).toBe(160);
    expect(m).toBeLessThan(before.n - 1);
    let builds = stats.builds;
    const rescales = stats.rescales;
    for (const [scale, height] of [
      [1.5, 240],
      [2.6, 416],
      [1.2, 192],
    ]) {
      t.vp(S).setRowView(scale, 0, "api");
      const after = snap();
      expect(after.height[m]).toBe(height);
      expect(stats.builds).toBe(builds);
      t.ctx.ui.toggleGroup("nothing");
      const rebuilt = snap();
      expect(stats.builds).toBe(builds + 1);
      builds++;
      expect(rebuilt.n).toBe(after.n);
      for (let k = 0; k < after.n; k++) {
        expect(after.top[k]).toBeCloseTo(rebuilt.top[k], 6);
        expect(after.height[k]).toBeCloseTo(rebuilt.height[k], 6);
      }
      expect([after.scrolled, after.vscroll]).toEqual([expect.closeTo(rebuilt.scrolled, 6), rebuilt.vscroll]);
      expect(after.content).toBeCloseTo(rebuilt.content as number, 6);
    }
    expect(stats.rescales - rescales).toBe(3);
  });

  it("a section that leaves and returns within one update keeps its layout; one that stays out is closed", async () => {
    const t = await opened(case6());
    const live = t.ctx.layout.get("chains.main");
    const refs = t.ctx.viewports.sections.value;
    t.layouts.sync([]);
    t.layouts.sync(refs);
    await Promise.resolve();
    expect(t.ctx.layout.get("chains.main")).toBe(live);
    let closed = false;
    live?.version.subscribe({
      complete: () => {
        closed = true;
      },
    });
    t.layouts.sync(refs.filter((r) => r !== "chains.main"));
    expect(closed).toBe(false);
    await Promise.resolve();
    expect([closed, t.ctx.layout.get("chains.main")]).toEqual([true, null]);
  });
});
