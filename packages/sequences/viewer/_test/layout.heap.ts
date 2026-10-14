/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for LY1's hot paths (seqspec:viewer#performance), on `heapScaling` (a step is one move or one frame's
// reads): hit testing per pointer move over every region of cases 3, 6 and 9, a wrapped section with inline labels, a
// reversed and a cyclic piece, split rows and the tree; layout lookups per frame with the view moving; fits, numbers
// and cached labels; `layoutRuns` on resize; a wrapped line's lanes over 5,000 features on 35,000 residues. The
// positive control: hit testing's moves allocating an object each (`leaking`), through the same measurement, are caught.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { beforeAll, describe, expect, it } from "vitest";
import { case1, case3, case6, case9, case13 } from "../../core/_test/format/cases.js";
import type { ShownIndex } from "../src/axis.js";
import { Hit, type LayerRenderer } from "../src/index.js";
import { type Layouts, SectionLayouts } from "../src/layout.js";
import { type LaneSegments, laneLineHeights, layoutRuns, lineTable } from "../src/runs.js";
import { fitCount } from "../src/text.js";
import { hitTree, TREE_INSET, type TreeRows, treeOf, treeSpan } from "../src/tree.js";
import { setPiece } from "./axis-fixtures.js";
import { FAMILY_SECTION as F, familyWithTree } from "./families.js";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.js";
import { laneRanges, type Opened, opened, renderers } from "./layout-fixtures.js";

const N = 50_000;
let sink = 0;

/** Points over a section's rows region and chrome, by section number. */
type Point = readonly [Opened, number, number, number];
const points: Point[] = [];

const rowsWithTree: LayerRenderer<TreeRows> = {
  name: "rows",
  geometries: ["linear"],
  prepare: (input) => {
    const a = input.section.blocks[0]?.alignment ?? -1;
    const aln = input.set.alignments[a];
    const track = input.section.tracks[input.layer.track];
    const rowOf = Int32Array.from({ length: track.laneCount }, (_, k) => input.section.lanes.row[track.laneStart + k]);
    return {
      set: input.set,
      alignment: a,
      rowOf,
      order: null,
      source: null,
      tree: aln?.tree ? treeOf(aln.tree) : null,
    };
  },
  header: { inset: () => TREE_INSET, hit: (v, s, _p, out) => hitTree(v, s, out), span: treeSpan },
};

function grid(t: Opened, section: string, dx: number, dy: number) {
  const l = t.layout(section);
  const s = t.ctx.layout.sectionIndex(section);
  for (let y = 1.5; y < l.height; y += dy) for (let x = 2.25; x < l.width; x += dx) points.push([t, s, x, y]);
}

beforeAll(async () => {
  const [c3, c6, c9, c13] = await Promise.all([
    opened(case3(), { size: { "family.main": [1000, 150] } }),
    opened(case6()),
    opened(case9()),
    opened(case13(), { size: { "features.main": [1000, 900] } }),
  ]);
  c9.vp("msa.main").scrollTo(0, 4711 * 18, "api");
  grid(c3, "family.main", 23, 7);
  grid(c6, "chains.main", 17, 9);
  grid(c9, "msa.main", 29, 11);
  grid(c13, "features.main", 31, 5);
  // Wrapped with inline labels: positions under a fractional x.
  const b = createSeqViewSpec();
  b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(30) });
  b.view({ ref: "v", label: "V", layout: "text" }).track({ ref: "v.letters" }).letters();
  const w = await opened(b.getDocument(), {
    size: { "v.main": [400, 900] },
    render: { wrappedGroupLabels: "inlineStart" },
  });
  grid(w, "v.main", 3.7, 13);
  // A reversed and a cyclic piece.
  for (const [start, stop, dir] of [
    [0, 300, -1],
    [150, 150, 1],
  ]) {
    const t = await opened(b.getDocument(), { size: { "v.main": [400, 300] } });
    const m = t.ctx.graph.assembled.value?.view;
    if (!m) throw new Error("no view");
    setPiece(m.resolved, 0, start, stop, dir, 300);
    (t.layouts as Layouts).update(m);
    grid(t, "v.main", 7.3, 11);
  }
  // The tree's header.
  const fam = await opened(familyWithTree(160, 20), {
    renderers: [...renderers().filter((r) => r.name !== "rows"), rowsWithTree as LayerRenderer],
    size: { [F]: [1000, 500] },
  });
  grid(fam, F, 3.1, 4);
}, 120_000);

describe("LY1's hot paths allocate nothing", () => {
  it("hit testing per pointer move: cases 3, 6, 9 and 13 (split rows), wrapped with inline labels, reversed and cyclic pieces, the tree", () => {
    const hit = Hit.create();
    const kinds = new Set<number>();
    let nodes = 0;
    for (const [t, s, x, y] of points) {
      t.ctx.layout.hitTest(s, x, y, hit);
      kinds.add(hit.kind * 16 + hit.region);
      if (hit.treeNode >= 0) nodes++;
    }
    expect([kinds.size > 10, nodes > 10]).toEqual([true, true]);
    const step = (i: number) => {
      const p = points[(i * 7919) % points.length];
      p[0].ctx.layout.hitTest(p[1], p[2], p[3], hit);
      sink = (sink + hit.kind + hit.index + hit.system.kind) & 0xffff;
    };
    const r = heapScaling(N, step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it("layout lookups per frame: the layout, the axis view and its x, the visible rows and entries, with the view moving", async () => {
    const t = await opened(case3(), { size: { "family.main": [1000, 150] } });
    const vp = t.vp("family.main");
    vp.zoomTo(0, 0, 10);
    const layout = t.ctx.layout;
    const r = heapScaling(N, (i) => {
      if (i & 1) vp.panBy(i & 2 ? 4 : -4, i & 4 ? 3 : -3, "pointer");
      const l = layout.get("family.main")?.current;
      const v = layout.axis("family.main");
      const rows = layout.visibleRows("family.main")?.current;
      if (!l || !v || !rows) return;
      sink = (sink + v.visLo[0] + v.x(0, (i & 15) as ShownIndex) + rows.count + SectionLayouts.entryAt(l, 9)) & 0xffff;
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("fits, numbers and cached labels", async () => {
    const t = await opened(case1());
    const f = t.ctx.text.font("ui", 12);
    const p = t.ctx.text.label("Homo sapiens (human) hemoglobin subunit alpha", f);
    const e = t.ctx.text.label("…", f).width;
    const r = heapScaling(N, (i) => {
      sink = (sink + fitCount(p, 40 + (i & 127), e)) & 0xffff;
      sink = (sink + t.ctx.text.numberWidth(f, (i & 4095) - 2048)) & 0xffff;
      sink = (sink + t.ctx.text.label("Homo sapiens (human) hemoglobin subunit alpha", f).ends.length) & 0xffff;
    });
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("layoutRuns on resize (case 1's text view, and 35,000 residues)", async () => {
    const t = await opened(case1(), { size: { "text.main": [900, 400] } });
    const small = t.ctx.layout.axis("text.main")?.runs ?? [];
    const b = createSeqViewSpec();
    const letters = "ACDEFGHIKLMNPQRSTVWY";
    b.data.sequence({ ref: "big", value: Array.from({ length: 35_000 }, (_, i) => letters[(i * 7) % 20]).join("") });
    b.view({ ref: "v", label: "V", layout: "text" }).track({ ref: "v.letters" }).letters();
    const u = await opened(b.getDocument(), { size: { "v.main": [1320, 640] } });
    const big = u.ctx.layout.axis("v.main")?.runs ?? [];
    expect([small.length, big[0]?.count]).toEqual([1, 35_000]);
    const out = lineTable(8192);
    const r = heapScaling(2_000, (i) => {
      sink = (sink + layoutRuns(small, 300 + (i & 511), 8, out)) & 0xffff;
      sink = (sink + layoutRuns(big, 600 + (i & 511), 8, out)) & 0xffff;
    });
    expect(isFlat(r), describeScaling(r, 2_000)).toBe(true);
  });

  it("a wrapped line's lanes: 5,000 features over 35,000 residues, every line's band", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "big", value: "MKVLATGHEK".repeat(3_500) });
    const start = Array.from({ length: 5_000 }, (_, i) => i * 7);
    b.data.annotation({ ref: "f", geometry: "ranges", target: a, value: { start, stop: start.map((s) => s + 4) } });
    const v = b.view({ ref: "v", label: "V", layout: "text" });
    v.track({ ref: "v.f" }).ranges({ data: "f" });
    const t = await opened(b.getDocument(), {
      renderers: renderers().map((x) => (x.name === "ranges" ? laneRanges() : x)),
      size: { "v.main": [1320, 640] },
    });
    const l = t.layout("v.main");
    const view = t.ctx.layout.axis("v.main");
    const m = t.ctx.graph.assembled.value?.view;
    const f = m?.prepared.find((x) => x.renderer === "ranges")?.state as LaneSegments;
    if (!view || !f) throw new Error("not opened");
    expect([l.lines.count > 200, f.seg.length > 1]).toEqual([true, true]);
    const out = new Float32Array(l.lines.count);
    const step = () => {
      laneLineHeights(f, null, view, l.lines, out);
      sink = (sink + out[0]) & 0xffff;
    };
    // 600 steps: at 200 the template's control leaked under the scaled limit, blind (LC-19).
    const r = heapScaling(600, step);
    expect(isFlat(r), describeScaling(r, 600)).toBe(true);
    expect(controlCaught(() => heapScaling(600, leaking(step)))).toBe(true);
  });
});
