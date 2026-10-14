/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hit testing (LY1): every region of cases 3, 6 and 9 fills the right locations and coordinate
// systems; features hit in logarithmic time; a section with no blocks hits nothing; a reversed and a cyclic piece hit
// their block indices through the layout map. Allocation-free: layout.heap.ts.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { CoordinateSystemKind } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { case3, case6, case9 } from "../../core/_test/format/cases.js";
import type { ShownIndex } from "../src/axis.js";
import { Hit, HoverKind as K, HitRegion as R } from "../src/index.js";
import type { Layouts } from "../src/layout.js";
import { setPiece } from "./axis-fixtures.js";
import { opened } from "./layout-fixtures.js";

const hitAt = (t: Awaited<ReturnType<typeof opened>>, section: string, x: number, y: number) => {
  const h = Hit.create();
  t.ctx.layout.hitTest(section, x, y, h);
  return h;
};
const keyOf = (t: Awaited<ReturnType<typeof opened>>, id: number) => {
  const set = t.ctx.graph.assembled.value?.set;
  const a = set?.annotations[id];
  return set && a ? set.strings.get(a.key) : "";
};

describe("CASE_3 (an alignment axis: navigator, overview, header, a data column, pinned columns, rows)", async () => {
  const t = await opened(case3(), { size: { "family.main": [1000, 150] } });
  const l = t.layout("family.main");
  const sh2 = t.ctx.graph.assembled.value?.set.alignments.find((a) => a)?.id;
  const row = (d: number) => l.scrolledTop + d * 18 + 9;

  it("lays out chrome as the msa preset gives it", () => {
    // Five pinned column tracks of 18 px under a 20 px navigator; the species column; the overview at the right.
    expect([l.navHeight, l.pinnedCount, l.scrolledTop, l.columnCount, l.columnX[0], l.dataLeft, l.dataRight]).toEqual([
      20, 5, 110, 1, 150, 268, 968,
    ]);
  });

  it("fills each region", () => {
    const at = (x: number, y: number) => {
      const h = hitAt(t, "family.main", x, y);
      return [h.region, h.kind, h.entry, h.dataColumn, h.displayRow, h.dataRow];
    };
    expect(at(500, 10)).toEqual([R.navigator, K.none, -1, -1, -1, -1]);
    expect(at(990, 60)).toEqual([R.vnav, K.none, -1, -1, -1, -1]);
    expect(at(50, 60)).toEqual([R.header, K.header, 2, -1, -1, -1]);
    expect(at(180, 60)).toEqual([R.dataColumnHeader, K.none, 2, 0, -1, -1]);
    expect(at(180, row(1))).toEqual([R.dataColumn, K.none, 6, 0, 1, 2]);
    expect(at(264, row(1))).toEqual([R.gutter, K.none, 6, -1, 1, 2]);
    // Four rows of 18 px under 110 px of chrome in 150 px: they overflow, and the right gutter is the scrollbar.
    expect(at(972, row(0))).toEqual([R.scrollbar, K.none, 5, -1, 0, 0]);
    expect(at(500, 149.5)[0]).toBe(R.data);
    expect(at(-1, 10)[0]).toBe(R.none);
  });

  it("a pinned entry is a column; its points are items; rows are cells, their residues off a gap", () => {
    const x = t.xOf("family.main", 0, 25);
    const c = hitAt(t, "family.main", x, l.entryTop[1] + 9);
    expect([c.region, c.kind, c.block, c.index, c.column.alignment, c.column.column]).toEqual([
      R.column,
      K.column,
      0,
      25,
      sh2,
      25,
    ]);
    // FLVR arginine (a point at column 25), in the fifth pinned track.
    const p = hitAt(t, "family.main", x, l.entryTop[4] + 9);
    expect([p.region, p.kind, keyOf(t, p.annotation), p.item, p.itemLoc.item]).toEqual([
      R.column,
      K.item,
      "flvr-arg",
      0,
      0,
    ]);
    // SRC (row 0) at column 2 is residue 2; its phosphosite (a point on its lane) is on top.
    const s = hitAt(t, "family.main", t.xOf("family.main", 0, 2), row(0));
    expect([
      s.kind,
      s.lane >= 0,
      keyOf(t, s.annotation),
      s.aligned.row,
      s.aligned.column,
      s.row.row,
      s.row.index,
    ]).toEqual([K.item, true, "sites", 0, 2, 0, 2]);
    // A gap of SRC's (column 10): a cell, no residue.
    const g = hitAt(t, "family.main", t.xOf("family.main", 0, 10), row(0));
    expect([g.kind, g.aligned.row, g.aligned.column, g.row.row]).toEqual([K.aligned, 0, 10, -1]);
    // GRB2 (row 2, display 1) at column 2 holds the other phosphosite; PTN11 (display 2) none there.
    expect(keyOf(t, hitAt(t, "family.main", t.xOf("family.main", 0, 2), row(1)).annotation)).toBe("sites");
    // Scrolled two rows down: PTN11 (row 1, display 2) has no site there.
    t.vp("family.main").scrollTo(0, 36, "api");
    const n = hitAt(t, "family.main", t.xOf("family.main", 0, 2), row(0));
    expect([n.kind, n.displayRow, n.dataRow, n.row.index]).toEqual([K.aligned, 2, 1, 2]);
    t.vp("family.main").scrollTo(0, 0, "api");
  });

  it("a collapsed clade's wedge row names its tree node", () => {
    t.ctx.ui.toggleClade("sh2", 1);
    const w = hitAt(t, "family.main", 500, row(0));
    expect([w.displayRow, w.dataRow, w.treeNode, w.kind]).toEqual([0, -1, 1, K.none]);
    t.ctx.ui.toggleClade("sh2", 1);
  });
});

describe("CASE_3 with the heads on a row of their own (`dataColumnHeads: own-row`)", async () => {
  const t = await opened(case3(), { size: { "family.main": [1000, 200] }, render: { dataColumnHeads: "ownRow" } });
  const l = t.layout("family.main");
  const at = (x: number, y: number) => {
    const h = hitAt(t, "family.main", x, y);
    return [h.region, h.kind, h.entry, h.dataColumn, h.displayRow];
  };

  it("the reserved row is the heads band, no entry's; the last pinned line and the rows below as before", () => {
    // Five pinned tracks end at 110; the heads take 110 to 128, where the rows start.
    expect([l.pinnedCount, l.headsTop, l.headsHeight, l.scrolledTop]).toEqual([5, 110, 18, 128]);
    expect(at(180, 119)).toEqual([R.dataColumnHeader, K.none, -1, 0, -1]);
    expect(at(500, 119)).toEqual([R.data, K.none, -1, -1, -1]);
    expect(at(50, 119)).toEqual([R.header, K.none, -1, -1, -1]);
    expect(at(180, 101)).toEqual([R.dataColumnHeader, K.none, 4, 0, -1]);
    expect(at(180, 128 + 18 + 9)).toEqual([R.dataColumn, K.none, 6, 0, 1]);
  });
});

describe("CASE_6 (blocks side by side, a ligand over three cells)", async () => {
  const t = await opened(case6());
  const set = t.ctx.graph.assembled.value?.set;
  const ax = t.vp("chains.main").axis;
  if (!set || !ax) throw new Error("not opened");
  const y = t.yOf("chains.main", "chains.residues");

  it("a residue of the kinase, every cell of the ligand, nothing in the gap", () => {
    const k = hitAt(t, "chains.main", t.xOf("chains.main", 0, 100), y);
    expect([k.region, k.kind, k.block, k.shown, k.index, k.seq.unit, k.seq.index, k.track]).toEqual([
      R.data,
      K.seq,
      0,
      100,
      100,
      ax.unit[0],
      100,
      t.trackIndex("chains.residues"),
    ]);
    const view = t.ctx.layout.axis("chains.main");
    if (!view) throw new Error("no view");
    const [x0, x1] = [view.x(1, 0 as ShownIndex), view.x(1, 1 as ShownIndex)];
    expect(x1 - x0).toBeCloseTo(3 * view.cell, 0);
    for (const f of [0.1, 0.5, 0.9]) {
      const h = hitAt(t, "chains.main", x0 + f * (x1 - x0), y);
      expect([h.block, h.shown, h.index, h.seq.unit]).toEqual([1, 0, 0, ax.unit[1]]);
    }
    const gap = hitAt(t, "chains.main", view.x(0, 274 as ShownIndex) + 14, y);
    expect([gap.region, gap.kind, gap.block]).toEqual([R.data, K.none, -1]);
  });

  it("the pinned ruler is a position on a blocks axis; header and gutters by x", () => {
    const r = hitAt(t, "chains.main", t.xOf("chains.main", 0, 5), 35);
    expect([r.region, r.kind, r.entry, r.index]).toEqual([R.data, K.seq, 0, 5]);
    expect(hitAt(t, "chains.main", 20, y).region).toBe(R.header);
    expect(hitAt(t, "chains.main", 154, y).region).toBe(R.gutter);
    expect(hitAt(t, "chains.main", 995, y).region).toBe(R.gutter);
    expect(hitAt(t, "chains.main", 500, 5).region).toBe(R.navigator);
  });
});

describe("CASE_9 (lazy rows bound to a structure; ranges on a linked section)", async () => {
  const t = await opened(case9());
  const set = t.ctx.graph.assembled.value?.set;
  if (!set) throw new Error("not opened");
  const l = t.layout("msa.main");
  const hba = 4711;

  it("a lazy row by y; the bound row's cell, residue and its unit", () => {
    expect([l.entryCount - l.pinnedCount, l.entryRun[l.pinnedCount], l.rows.lazy]).toEqual([1, 18000, true]);
    t.vp("msa.main").scrollTo(0, hba * 18, "api");
    const h = hitAt(t, "msa.main", t.xOf("msa.main", 0, 40), l.scrolledTop + 9);
    expect([h.kind, h.displayRow, h.dataRow, h.aligned.row, h.aligned.column, h.row.index]).toEqual([
      K.aligned,
      hba,
      hba,
      hba,
      40,
      0,
    ]);
    // Bound through the row's binding to 1A3N's α chains: the unit's molecule is chain A's.
    const a = t.vp("msa.chains").axis?.unit[0] ?? -1;
    expect(h.seq.unit >= 0 && set.units.molecule[h.seq.unit] === set.units.molecule[a]).toBe(true);
    const next = hitAt(t, "msa.main", t.xOf("msa.main", 0, 40), l.scrolledTop + 18 + 9);
    expect([next.dataRow, next.kind]).toEqual([hba + 1, K.aligned]);
  });

  it("HBB's row is bound to 1A3N's β chains: its cell's unit is chain B's molecule", () => {
    const hbb = 9133;
    t.vp("msa.main").scrollTo(0, hbb * 18, "api");
    const h = hitAt(t, "msa.main", t.xOf("msa.main", 0, 40), l.scrolledTop + 9);
    expect([h.kind, h.dataRow, h.row.index]).toEqual([K.aligned, hbb, 2]);
    const b = t.vp("msa.chains").axis?.unit[1] ?? -1;
    expect(set.strings.get(set.units.key[b])).toBe("1a3n/B");
    expect(h.seq.unit >= 0 && set.units.molecule[h.seq.unit] === set.units.molecule[b]).toBe(true);
  });

  it("the linked section's ranges are items (the secondary structure, drawn as ranges)", () => {
    const unit = t.vp("msa.chains").axis?.unit[0] ?? -1;
    const ss = set.annotations.find((a) => a && set.strings.get(a.key).includes("ss"));
    if (!ss) throw new Error("no ss");
    const r = Array.from({ length: ss.count }, (_, k) => k).find((k) => ss.target.unit[k] === unit);
    if (r === undefined) throw new Error("no range on A");
    const h = hitAt(
      t,
      "msa.chains",
      t.xOf("msa.chains", 0, ss.start[r]),
      t.yOf("msa.chains", t.entries("msa.chains")[0][1]),
    );
    expect([h.kind, h.annotation, h.item, h.block]).toEqual([K.item, ss.id, ss.item[r], 0]);
  });
});

describe("features", () => {
  const doc = (n: number) => {
    const b = createSeqViewSpec();
    const len = 100_000;
    const a = b.data.sequence({ ref: "a", value: "A".repeat(len) });
    const start = Array.from({ length: n }, (_, i) => Math.floor((i * len) / n));
    b.data.annotation({ ref: "f", geometry: "ranges", target: a, value: { start, stop: start.map((s) => s + 5) } });
    b.view({ ref: "v", label: "V" }).track({ ref: "v.f" }).ranges({ data: "f" });
    return b.getDocument();
  };

  it("hit in logarithmic time: 10,000 features cost about what 100 do", { tags: ["slow"] }, async () => {
    const time = async (n: number) => {
      const t = await opened(doc(n));
      const h = Hit.create();
      const y = t.yOf("v.main", "v.f");
      let items = 0;
      const run = () => {
        for (let k = 0; k < 200_000; k++) {
          t.ctx.layout.hitTest("v.main", 158 + ((k * 7919) % 834), y, h);
          if (h.kind === K.item) items++;
        }
      };
      run();
      const t0 = performance.now();
      run();
      return { ms: performance.now() - t0, items };
    };
    const few = await time(100);
    const many = await time(10_000);
    expect(many.items).toBeGreaterThan(few.items);
    expect(many.ms / few.ms).toBeLessThan(4);
  });
});

describe("LY1", () => {
  const seq = () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEKMKVLATGHEK" });
    b.view({ ref: "v", label: "V" }).track({ ref: "v.letters" }).letters();
    return b.getDocument();
  };

  it("the element's coordinate system: a unit's position, a column, a row's residue off a gap, the columns at one", async () => {
    const t = await opened(seq());
    const h = hitAt(t, "v.main", t.xOf("v.main", 0, 3), t.yOf("v.main", "v.letters"));
    const unit = t.vp("v.main").axis?.unit[0];
    expect([h.kind, h.index, h.system.kind, h.system.id, h.system.row]).toEqual([
      K.seq,
      3,
      CoordinateSystemKind.unit,
      unit,
      -1,
    ]);
    const c3 = await opened(case3(), { size: { "family.main": [1000, 150] } });
    const l = c3.layout("family.main");
    const sh2 = c3.ctx.graph.assembled.value?.set.alignments.find((a) => a)?.id;
    const col = hitAt(c3, "family.main", c3.xOf("family.main", 0, 25), l.entryTop[1] + 9);
    expect([col.kind, col.system.kind, col.system.id]).toEqual([K.column, CoordinateSystemKind.columns, sh2]);
    const row = (d: number) => l.scrolledTop + d * 18 + 9;
    const cell = hitAt(c3, "family.main", c3.xOf("family.main", 0, 3), row(0));
    expect([cell.kind, cell.system.kind, cell.system.id, cell.system.row]).toEqual([
      K.aligned,
      CoordinateSystemKind.row,
      sh2,
      0,
    ]);
    const gap = hitAt(c3, "family.main", c3.xOf("family.main", 0, 10), row(0));
    expect([gap.kind, gap.row.row, gap.system.kind, gap.system.row]).toEqual([
      K.aligned,
      -1,
      CoordinateSystemKind.columns,
      -1,
    ]);
  });

  it("a section with no blocks hits nothing, and does not throw", async () => {
    const t = await opened(seq());
    // The view's one section loses its blocks (a failed source): the layout rebuilt from it.
    const m = t.ctx.graph.assembled.value?.view;
    if (!m) throw new Error("no view");
    (m.resolved.sections.blockCount as Int32Array)[0] = 0;
    (t.layouts as Layouts).update(m);
    const h = Hit.create();
    expect(t.ctx.layout.hitTest("v.main", 500, t.yOf("v.main", "v.letters"), h)).toBe(false);
    expect([h.region, h.kind, h.block]).toEqual([R.none, K.none, -1]);
  });

  it("a reversed piece and one across a cyclic origin hit their block indices through the layout map", async () => {
    const at = async (start: number, stop: number, dir: number) => {
      const t = await opened(seq());
      const m = t.ctx.graph.assembled.value?.view;
      if (!m) throw new Error("no view");
      setPiece(m.resolved, 0, start, stop, dir, 20);
      (t.layouts as Layouts).update(m);
      const v = t.ctx.layout.axis("v.main");
      if (!v) throw new Error("no view");
      const y = t.yOf("v.main", "v.letters");
      return [0, 1, 5, 19].map((s) => {
        const h = hitAt(t, "v.main", (v.x(0, s as ShownIndex) + v.x(0, (s + 1) as ShownIndex)) / 2, y);
        return [h.shown, h.index, h.seq.index];
      });
    };
    expect(await at(0, 20, -1)).toEqual([
      [0, 19, 19],
      [1, 18, 18],
      [5, 14, 14],
      [19, 0, 0],
    ]);
    expect(await at(15, 15, 1)).toEqual([
      [0, 15, 15],
      [1, 16, 16],
      [5, 0, 0],
      [19, 14, 14],
    ]);
  });
});
