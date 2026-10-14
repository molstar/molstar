/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Blocks on coordinate systems (VW4): a block's piece and its layout map (the only conversion between shown positions
// and indices; a ranged unit at every index; `dir` -1 and a cyclic origin through `layoutOf`, which AX1 and CY1 give
// the format), and `ResolvedView.systems`: the coordinate systems on each section's axis and the closure index over the
// span maps of shown placements and links.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type CoordinateSystem, Coords, OrderedSet, type SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  compile,
  layoutOf,
  type ResolvedView,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
  validate,
} from "../../src/format/index.ts";
import { case10 } from "./cases.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";

const reg = registries();
const all = supportedAll(reg);
const codes = (issues: { code: string }[]) => issues.map((i) => i.code);
async function resolve(doc: SeqViewSpecDoc, view: string) {
  const { set, issues } = await compile(doc, reg, fixtureLoaders());
  expect(issues.filter((i) => i.severity === "error")).toEqual([]);
  return { set, ...resolveView(doc, set, view, reg, all) };
}
const S: CoordinateSystem = { kind: 1, id: 0, row: -1 };

/** Index → shown and shown → index over the whole domain, through the map alone. */
function readBack(map: ReturnType<typeof layoutOf>["map"], domain: number) {
  const shown = Array.from({ length: domain }, (_, i) => Coords.toAt(map, i, 0));
  const index: number[] = [];
  const n = shown.filter((k) => k >= 0).length;
  for (let k = 0; k < n; k++) Coords.forEachFrom(map, k, k + 1, (s) => index.push(s), null);
  return { shown, index };
}

describe("pieces and their layout maps", () => {
  it("forward: shown k is index first + k, nothing outside the piece", () => {
    const p = layoutOf(S, 5, 15, 1, 10);
    expect(p.first).toBe(5);
    const { shown, index } = readBack(p.map, 20);
    expect(shown).toEqual([...Array(5).fill(-1), ...Array.from({ length: 10 }, (_, k) => k), ...Array(5).fill(-1)]);
    expect(index).toEqual(Array.from({ length: 10 }, (_, k) => 5 + k));
  });

  it("dir -1: from the piece's highest index down; first is stop - 1", () => {
    const p = layoutOf(S, 3, 8, -1, 10);
    expect(p.first).toBe(7);
    const { shown, index } = readBack(p.map, 10);
    expect(shown).toEqual([-1, -1, -1, 4, 3, 2, 1, 0, -1, -1]);
    expect(index).toEqual([7, 6, 5, 4, 3]);
  });

  it("across a cyclic origin: [start, L) then [0, stop), reversed with dir -1", () => {
    const fwd = layoutOf(S, 8, 3, 1, 10);
    expect([fwd.first, readBack(fwd.map, 10).index]).toEqual([8, [8, 9, 0, 1, 2]]);
    const rev = layoutOf(S, 8, 3, -1, 10);
    expect([rev.first, readBack(rev.map, 10).index]).toEqual([2, [2, 1, 0, 9, 8]]);
    expect(fwd.map.count).toBe(2);
  });

  it("a whole cyclic system from an origin (start === stop), either way", () => {
    const fwd = layoutOf(S, 4, 4, 1, 10);
    expect([fwd.first, readBack(fwd.map, 10).index]).toEqual([4, [4, 5, 6, 7, 8, 9, 0, 1, 2, 3]]);
    const rev = layoutOf(S, 4, 4, -1, 10);
    expect([rev.first, readBack(rev.map, 10).index]).toEqual([3, [3, 2, 1, 0, 9, 8, 7, 6, 5, 4]]);
    // Origin 0 is the identity, one segment.
    expect(layoutOf(S, 0, 0, 1, 10).map.count).toBe(1);
  });

  it("a shown run converts to index runs and back through the map (the viewer's run converters read it)", () => {
    const rev = layoutOf(S, 8, 3, -1, 10);
    const runs: number[][] = [];
    Coords.forEachFrom(rev.map, 1, 4, (s, t) => runs.push([s, t]), null);
    expect(runs.sort((a, b) => a[0] - b[0])).toEqual([
      [0, 2],
      [9, 10],
    ]);
    const back: number[][] = [];
    Coords.forEachTo(rev.map, 0, 2, (s, t) => back.push([s, t]), null);
    expect(back).toEqual([[1, 3]]); // indices 0 and 1 are shown at 2 and 1
  });
});

describe("a ranged unit's block indexes from its rangeStart", () => {
  function ranged() {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQRSTVW" });
    const u = b.data.unit({ ref: "u", sequence: s, range: [5, 15] });
    const every = Array.from({ length: 20 }, (_, i) => i);
    const one = b.data.annotation({ ref: "one", geometry: "ranges", target: s, value: { start: every } });
    const site = b.data.annotation({ ref: "site", geometry: "ranges", target: s, value: { start: [7], stop: [10] } });
    const v = b.view({ ref: "v", label: "V" }).block({ unit: u });
    v.track({ label: "Every position" }).ranges({ data: one });
    v.focus({ target: u });
    b.view({ ref: "f", label: "Filtered", positions: site }).block({ unit: u }).track({ label: "L" }).letters();
    b.view({ ref: "out", label: "Out" }).block({ unit: u }).focus({ target: u, start: 2, stop: 6 });
    return b.getDocument();
  }

  it("[5, 15): base, start and first 5; every index maps to its shown position and back", async () => {
    const r = await resolve(ranged(), "v");
    const b = r.view.blocks;
    expect([b.base[0], b.start[0], b.stop[0], b.first[0], b.dir[0], b.length[0]]).toEqual([5, 5, 15, 5, 1, 10]);
    const { shown, index } = readBack(b.layout[0], 20);
    expect(shown).toEqual(Array.from({ length: 20 }, (_, i) => (i >= 5 && i < 15 ? i - 5 : -1)));
    expect(index).toEqual(Array.from({ length: 10 }, (_, k) => 5 + k));
    // Data on the molecule lands directly, in the same indices: each row at its own index, the viewer clips.
    const l = r.view.layers;
    expect([l.placementCount[0], r.view.placements.path[l.placementStart[0]]]).toEqual([1, -1]);
    expect(r.view.sections.focus[0]).toEqual({ block: 0, start: 5, stop: 15 });
  });

  it("a positions filter in block indices, and a focus outside the unit's domain", async () => {
    const f = await resolve(ranged(), "f");
    const shown = f.view.blocks.positions[0];
    expect(shown && [OrderedSet.size(shown), OrderedSet.min(shown), OrderedSet.max(shown)]).toEqual([3, 7, 9]);
    const out = await resolve(ranged(), "out");
    expect(codes(out.issues)).toEqual(["svs.focus.range"]);
  });

  it("case 10's abl1-kd [241, 493): every index, and its sites through the layout map", async () => {
    const r = await resolve(case10(), "text");
    const b = r.view.blocks;
    expect([b.base[0], b.first[0], b.stop[0]]).toEqual([241, 241, 493]);
    const { shown } = readBack(b.layout[0], 600);
    for (let i = 0; i < 600; i++) expect(shown[i]).toBe(i >= 241 && i < 493 ? i - 241 : -1);
    // The emphasized sites (UniProt 253, 255, 315, 351) are molecule indices placed directly: T315 is shown at 73.
    const sites = r.set.annotations.find((a) => a && r.set.strings.get(a.key) === "kd-sites");
    expect(Array.from(sites?.start ?? [], (i) => Coords.toAt(b.layout[0], i, 0))).toEqual([11, 13, 73, 109]);
  });
});

/**
 * Pieces in the format (AX1; seqspec:format#blocks): 20 residues `s`, its unit `u` and a second copy `u2` of one group;
 * a section `-1` whose blocks are `u`'s [2, 8) reversed, `u`'s [10, 16) forward (the block's direction wins) and `u2`'s
 * [2, 8) reversed (a copy of the first piece, merged); `site` covers [4, 6) and [12, 13).
 */
function pieces(o: { block?: Record<string, unknown>; focus?: Record<string, unknown> } = {}) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAQRSTVW" });
  const u = b.data.unit({ ref: "u", sequence: s });
  const u2 = b.data.unit({ ref: "u2", sequence: s });
  const site = b.data.annotation({
    ref: "site",
    geometry: "ranges",
    target: s,
    value: { start: [4, 12], stop: [6, 13] },
  });
  const v = b.view({ ref: "v", label: "V", direction: -1 });
  v.block({ unit: u, range: [2, 8], ...o.block });
  v.block({ unit: u, range: [10, 16], direction: 1 });
  v.block({ unit: u2, range: [2, 8] });
  v.track({ label: "L" }).letters();
  if (o.focus) v.focus(o.focus as never);
  b.view({ ref: "f", label: "F", positions: site })
    .block({ unit: u, range: [2, 8], direction: -1 })
    .track()
    .letters();
  b.view({ ref: "a", label: "A", direction: -1 }).block({ sequence: s }).track().letters();
  return b.getDocument();
}
const piecesOf = (v: ResolvedView) =>
  Array.from({ length: v.blocks.start.length }, (_, k) => [v.blocks.start[k], v.blocks.stop[k], v.blocks.dir[k]]);

describe("pieces in the format: block range and direction, section direction (AX1)", () => {
  it("range and direction make the piece; the block's direction wins over the section's; copies of one piece merge", async () => {
    const r = await resolve(pieces(), "v");
    const b = r.view.blocks;
    expect(piecesOf(r.view)).toEqual([
      [2, 8, -1],
      [10, 16, 1],
    ]);
    expect([b.first[0], b.first[1], b.base[0], b.length[0]]).toEqual([7, 10, 0, 20]);
    expect(readBack(b.layout[0], 20).index).toEqual([7, 6, 5, 4, 3, 2]);
    // u and u2 over [2, 8) reversed are one block of two members; u's [10, 16) is a block of its own.
    expect([b.members[0]?.length, b.members[1]]).toEqual([2, null]);
    // A section's direction alone reverses every block.
    const a = await resolve(pieces(), "a");
    expect(piecesOf(a.view)).toEqual([[0, 20, -1]]);
    expect(a.view.blocks.first[0]).toBe(19);
  });

  it("a positions filter within the piece: only what it covers there; a piece it misses is left out", async () => {
    const f = await resolve(pieces(), "f");
    const shown = f.view.blocks.positions[0];
    expect(shown && Array.from({ length: OrderedSet.size(shown) }, (_, k) => OrderedSet.getAt(shown, k))).toEqual([
      4, 5,
    ]);
    const doc = pieces();
    const fv = doc.root.children?.find((c) => c.ref === "f");
    const block = fv?.children?.[0].children?.find((c) => c.kind === "block");
    if (block) block.params = { ...block.params, range: [14, 18] };
    expect((await resolve(doc, "f")).view.blocks.start.length).toBe(0);
  });

  it("a focus opens on the first piece that shows it; outside every piece it is svs.focus.range", async () => {
    const r = await resolve(pieces({ focus: { target: "u", start: 11, stop: 14 } }), "v");
    expect(r.view.sections.focus[0]).toEqual({ block: 1, start: 11, stop: 14 });
    const out = await resolve(pieces({ focus: { target: "u", start: 17, stop: 19 } }), "v");
    expect(codes(out.issues)).toEqual(["svs.focus.range"]);
  });

  it("a range outside the unit's domain is svs.position.range at resolve, and its block is left out", async () => {
    const r = await resolve(pieces({ block: { range: [15, 25] } }), "v");
    expect(codes(r.issues)).toEqual(["svs.position.range"]);
    expect(r.issues[0].path).toBe("/root/children/1/children/0/children/0/params/range");
    expect(piecesOf(r.view)).toEqual([
      [10, 16, 1],
      [2, 8, -1],
    ]);
  });

  it("an alignment's columns reversed; equal axes compare pieces (svs.section.scroll_group)", async () => {
    const b = createSeqViewSpec();
    const aln = b.data.alignment({ ref: "aln", value: ["MKV-L", "MK-AL"] });
    const v = b.view({ ref: "v", label: "V" });
    v.section({ ref: "v.one", scroll_group: "g", direction: -1 }).block({ alignment: aln }).track().rows();
    v.section({ ref: "v.two", scroll_group: "g" }).block({ alignment: aln }).track().rows();
    const { set, issues } = await compile(b.getDocument(), reg, fixtureLoaders());
    const r = resolveView(b.getDocument(), set, "v", reg, all);
    expect([...codes(issues), ...codes(r.issues)]).toEqual(["svs.section.scroll_group"]);
    expect(piecesOf(r.view)).toEqual([
      [0, 5, -1],
      [0, 5, 1],
    ]);
    expect(r.view.blocks.first[0]).toBe(4);
  });

  it("the validator: range on a complex or an alignment, an empty or negative range, a direction not 1 or -1", () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "MKVLA" });
    const aln = b.data.alignment({ ref: "aln", value: ["MK", "MK"] });
    b.view({ ref: "v", label: "V" })
      .block({ alignment: aln, range: [0, 1] })
      .track()
      .rows();
    b.view({ ref: "w", label: "W" })
      .block({ sequence: s, range: [3, 3] })
      .track()
      .letters();
    b.view({ ref: "x", label: "X" })
      .block({ sequence: s, range: [-1, 2] })
      .track()
      .letters();
    b.view({ ref: "y", label: "Y", direction: 2 as never })
      .block({ sequence: s, direction: 0 as never })
      .track()
      .letters();
    const issues = validate(b.getDocument(), reg).map((i) => [i.code, i.path, i.message]);
    expect(issues).toEqual([
      [
        "svs.params.exclusive",
        "/root/children/1/children/0/children/0/params/range",
        expect.stringMatching(/columns take positions/),
      ],
      [
        "svs.position.range",
        "/root/children/2/children/0/children/0/params/range",
        expect.stringMatching(/give start < stop/),
      ],
      ["svs.position.range", "/root/children/3/children/0/children/0/params/range", expect.any(String)],
      ["svs.params.type", "/root/children/4/children/0/params/direction", expect.any(String)],
      ["svs.params.type", "/root/children/4/children/0/children/0/params/direction", expect.any(String)],
    ]);
  });
});

/** Case 14's gene with one or two proteins, a letters track of them and a protein section (a link). */
function gene(o: { proteins?: 1 | 2; hideTrack?: boolean; hideSection?: boolean } = {}) {
  const b = createSeqViewSpec();
  const g = b.data.sequence({ ref: "gene", alphabet: "dna", value: "ATGGCTTCCTAA" });
  b.data.sequence({ ref: "p1", alphabet: "protein", value: "MAS" });
  b.data.mapping({ ref: "m1", from: "p1", to: "gene", method: "codon", cds: [[0, 9]] } as never);
  if (o.proteins === 2) {
    b.data.sequence({ ref: "p2", alphabet: "protein", value: "AS" });
    b.data.mapping({ ref: "m2", from: "p2", to: "gene", method: "codon", cds: [[3, 9]] } as never);
  }
  const v = b.view({ ref: "v", label: "V" });
  v.block({ sequence: g })
    .track({ label: "Proteins", hidden: o.hideTrack ?? false })
    .letters({ of: "mapped" });
  v.section({ ref: "v.protein", hidden: o.hideSection ?? false })
    .block({ sequence: "p1" })
    .track({ label: "P" })
    .letters();
  return b.getDocument();
}
const systemName = (set: SeqSet, v: ResolvedView, k: number) => {
  const s = v.systems;
  const key = s.kind[k] === 0 ? set.molecules[s.id[k]]?.key : set.units.key[s.id[k]];
  return `${["molecule", "unit"][s.kind[k]] ?? s.kind[k]} ${set.strings.get(key ?? -1)}`;
};

describe("ResolvedView.systems", () => {
  it("the systems on each section's axis, their domains, and every system a placement or link reads", async () => {
    const { set, view } = await resolve(gene(), "v");
    const s = view.systems;
    const names = Array.from({ length: s.count }, (_, k) => systemName(set, view, k));
    expect(names.slice(0, 2)).toEqual(["unit gene/unit", "unit p1/unit"]);
    expect(names).toContain("molecule p1");
    expect([Array.from(s.sectionStart), Array.from(s.sectionSystems)]).toEqual([
      [0, 1, 2],
      [0, 1],
    ]);
    expect([s.start[0], s.stop[0], s.start[1], s.stop[1]]).toEqual([0, 12, 0, 3]);
  });

  it("the closure index: each span map once (the letters' and the link's), by to system and hull", async () => {
    const { set, view } = await resolve(gene(), "v");
    const c = view.systems.closure;
    expect(c.count).toBe(2);
    for (let k = 0; k < c.count; k++) {
      expect(view.paths[c.path[k]].map.span).toBe(3);
      expect([systemName(set, view, c.to[k]), c.toStart[k], c.toStop[k]]).toEqual(["unit gene/unit", 0, 9]);
    }
    expect(new Set(Array.from(c.path)).size).toBe(2);
  });

  it("hidden tracks and sections add no closure maps", async () => {
    expect((await resolve(gene({ hideTrack: true }), "v")).view.systems.closure.count).toBe(1);
    expect((await resolve(gene({ hideSection: true }), "v")).view.systems.closure.count).toBe(1);
    expect((await resolve(gene({ hideTrack: true, hideSection: true }), "v")).view.systems.closure.count).toBe(0);
  });

  it("sorted by toStart within a to system, toMaxStop its running maximum", async () => {
    const { view } = await resolve(gene({ proteins: 2, hideSection: true }), "v");
    const c = view.systems.closure;
    expect([Array.from(c.toStart), Array.from(c.toStop), Array.from(c.toMaxStop)]).toEqual([
      [0, 3],
      [9, 9],
      [9, 9],
    ]);
  });

  it("span-1 placements and copies add none (case 10's chains, alignment and strip)", async () => {
    for (const v of ["chains", "alignment", "strip"])
      expect((await resolve(case10(), v)).view.systems.closure.count).toBe(0);
  });
});
