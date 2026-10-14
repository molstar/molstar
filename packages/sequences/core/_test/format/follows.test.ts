/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Axis links (LK1; seqspec:format#sections): `section.follow` and `follow_anchor` resolved into `ResolvedView.follows`,
// a route per block of either side through the sides' anchors; the validator's ref checks; `svs.section.follow` for a
// section of another view, itself, a loop (through follows and scroll groups), sections no path joins, and a misplaced
// anchor; a pair kept without routes while a side has nothing loaded; and the builder's `codingPair`.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import { Coords, type SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import {
  compile,
  type ResolvedView,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
  validate,
} from "../../src/format/index.js";
import { case9 } from "./cases.js";
import { fixtureLoaders, PF00042, registries } from "./load-fixtures.js";

const reg = registries();
const all = supportedAll(reg);
const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);
async function resolve(doc: SeqViewSpecDoc, view: string) {
  expect(validate(doc, reg).filter((i) => i.severity === "error")).toEqual([]);
  const { set, issues } = await compile(doc, reg, fixtureLoaders());
  expect(issues.filter((i) => i.severity === "error")).toEqual([]);
  return { set, ...resolveView(doc, set, view, reg, all) };
}

/** Case 14's gene and its translation. */
function gene(b = createSeqViewSpec()) {
  const g = b.data.sequence({ ref: "gene", alphabet: "dna", value: `GAATTCGCCACC${"ATGGCTTCCAAGGAA"}CTGACCCCACGTGAC` });
  const p = translate(b.data, { ref: "protein", from: g, cds: [[12, 42]] });
  return { b, g, p };
}

/** Each pair's routes as [section, block, to, forward, span] from both sides. */
function routesOf(v: ResolvedView, set: SeqSet) {
  const f = v.follows;
  const name = (s: number) => set.strings.get(v.sections.ref[s]);
  const out: (string | number)[][] = [];
  for (let p = 0; p < f.count; p++) {
    let r = f.routeStart[p];
    for (const s of [f.a[p], f.b[p]])
      for (let k = 0; k < v.sections.blockCount[s]; k++, r++) {
        const path = f.routes.path[r];
        out.push([name(s), k, f.routes.to[r], f.routes.forward[r], path < 0 ? 0 : v.paths[path].map.span]);
      }
  }
  return out;
}

describe("follows", () => {
  it("a protein following its gene: one pair, a route each way through the codon step, the protein's cell three", async () => {
    const { b, g, p } = gene();
    const v = b.view({ ref: "v", label: "V" });
    const dna = v.section({ ref: "v.dna" }).block({ sequence: g });
    dna.track().letters();
    v.section({ ref: "v.protein", follow: dna.handle }).block({ sequence: p }).track().letters();
    const r = await resolve(b.getDocument(), "v");
    expect(r.issues).toEqual([]);
    const f = r.view.follows;
    expect([f.count, f.a[0], f.b[0], f.scale[0]]).toEqual([1, 0, 1, 3]);
    // The codon path keeps the protein as `from`: read backward from the gene, forward from the protein.
    expect(routesOf(r.view, r.set)).toEqual([
      ["v.dna", 0, 0, 0, 3],
      ["v.protein", 0, 0, 1, 3],
    ]);
    const m = r.view.paths[f.routes.path[1]].map;
    expect([Coords.toAt(m, 0, 0), Coords.toAt(m, 9, 2)]).toEqual([12, 41]);
  });

  it("pieces and blocks: each block of a multi-block side routes to the first block it reaches", async () => {
    const { b, g, p } = gene();
    const other = b.data.sequence({ ref: "other", value: "MKVLAAG" });
    const v = b.view({ ref: "v", label: "V" });
    const top = v.section({ ref: "v.top" });
    b.preset.multiRegion(
      top,
      [
        [g, 0, 20],
        [g, 20, 42],
      ],
      { direction: -1 },
    );
    top.block({ sequence: other }).track().letters();
    v.section({ ref: "v.p", follow: top.handle }).block({ sequence: p }).track().letters();
    const r = await resolve(b.getDocument(), "v");
    expect(r.issues).toEqual([]);
    expect(r.view.follows.scale[0]).toBe(3);
    // `other` reaches nothing; the protein reaches the first piece of the gene (the viewer picks the piece by index).
    expect(routesOf(r.view, r.set)).toEqual([
      ["v.top", 0, 0, 0, 3],
      ["v.top", 1, 0, 0, 3],
      ["v.top", 2, -1, 0, 0],
      ["v.p", 0, 0, 1, 3],
    ]);
  });

  it("chains are allowed; the same unit in two sections follows through the identity", async () => {
    const { b, g, p } = gene();
    const v = b.view({ ref: "v", label: "V" });
    const one = v.section({ ref: "v.one" }).block({ sequence: g });
    one.track().letters();
    const two = v.section({ ref: "v.two", direction: -1, follow: one.handle }).block({ sequence: g });
    two.track().letters();
    v.section({ ref: "v.three", follow: two.handle }).block({ sequence: p }).track().letters();
    const r = await resolve(b.getDocument(), "v");
    expect(r.issues).toEqual([]);
    expect([...r.view.follows.a, ...r.view.follows.b, ...r.view.follows.scale]).toEqual([0, 1, 1, 2, 1, 3]);
    expect(r.view.paths[r.view.follows.routes.path[0]].map.count).toBe(1);
  });

  it("svs.section.follow: another view's section, itself, a loop by follows or by a scroll group, no path", async () => {
    const { b, g, p } = gene();
    const lone = b.data.sequence({ ref: "lone", value: "MKVLA" });
    const w = b.view({ ref: "w", label: "W" });
    w.section({ ref: "w.s" }).block({ sequence: g }).track().letters();
    const v = b.view({ ref: "v", label: "V" });
    const sec = (ref: string, of: typeof g, params: Record<string, unknown> = {}) =>
      v
        .section({ ref, ...params })
        .block({ sequence: of })
        .track()
        .letters();
    sec("v.a", g, { follow: "w.s" });
    sec("v.b", g, { follow: "v.b" });
    sec("v.c", g, { follow: "v.d" });
    sec("v.d", p, { follow: "v.c" });
    sec("v.e", g, { scroll_group: "s" });
    sec("v.f", g, { scroll_group: "s", follow: "v.e" });
    sec("v.g", lone, { follow: "v.e" });
    const r = await resolve(b.getDocument(), "v");
    expect(r.issues.map((i) => [i.code, i.path.split("/").slice(-2).join("/"), i.message.split(":")[0]])).toEqual([
      ["svs.section.follow", "params/follow", '"w.s" is not a section of view "v"'],
      ["svs.section.follow", "params/follow", "a section cannot follow itself"],
      ["svs.section.follow", "params/follow", 'following "v.c" closes a loop'],
      ["svs.section.follow", "params/follow", 'following "v.e" closes a loop'],
      ["svs.section.follow", "params/follow", 'no mapping path joins "v.g" and "v.e"'],
    ]);
    expect(r.issues[2].path).toContain("/children/3/");
    // The one pair kept: v.d follows v.c (v.c's own follow closed the loop).
    expect([r.view.follows.count, r.view.follows.a[0], r.view.follows.b[0]]).toEqual([1, 3, 2]);
  });

  it("the validator: follow names a section, follow_anchor a row or a unit", () => {
    const { b, g } = gene();
    const v = b.view({ ref: "v", label: "V" });
    v.section({ ref: "v.a", follow: "gene", follow_anchor: "v.a" }).block({ sequence: g });
    expect(codes(validate(b.getDocument(), reg))).toEqual(["svs.ref.wrong_kind", "svs.ref.wrong_kind"]);
  });

  it("a side with nothing loaded yet: the pair kept without routes, and no issue", async () => {
    const { b, g } = gene();
    const later = b.data.sequence({ ref: "later", source: { host_file: "missing.fasta", format: "fasta" } });
    const v = b.view({ ref: "v", label: "V" });
    const a = v.section({ ref: "v.a" }).block({ sequence: g });
    a.track().letters();
    v.section({ ref: "v.b", follow: a.handle }).block({ sequence: later }).track().letters();
    const doc = b.getDocument();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const r = resolveView(doc, set, "v", reg, all);
    expect(codes(r.issues)).toEqual([]);
    expect([r.view.follows.count, r.view.sections.blockCount[1], ...r.view.follows.routes.path]).toEqual([1, 0, -1]);
  });
});

describe("alignment axes: anchored on a row", () => {
  /** Two protein rows and their cDNAs, gapped differently. */
  function pair(o: { anchor?: string; reference?: number } = {}) {
    const b = createSeqViewSpec();
    const prot = b.data.alignment({
      ref: "prot",
      alphabet: "protein",
      value: ["MA-S", "MAKS"],
      ...(o.reference !== undefined ? { reference_row: o.reference } : {}),
    });
    const cdna = b.data.alignment({ ref: "cdna", alphabet: "dna", value: ["ATGGCT---TCC", "ATGGCTAAGTCC"] });
    const v = b.view({ ref: "v", label: "V" });
    const { protein, cdna: c } = b.preset.codingPair(v, prot, cdna);
    protein.track().rows();
    c.track().rows();
    const doc = b.getDocument();
    const top = doc.root.children?.[1].children?.[0];
    if (o.anchor && top) top.params = { ...top.params, follow_anchor: o.anchor };
    return doc;
  }

  it("codingPair: a codon mapping per row; the cDNA follows the protein at a third of its cell", async () => {
    const r = await resolve(pair(), "v");
    expect(r.issues).toEqual([]);
    const f = r.view.follows;
    expect([f.count, f.a[0], f.b[0], f.scale[0]]).toEqual([1, 0, 1, 1 / 3]);
    // Columns → reference row 0 → its cDNA row → the cDNA's columns: protein column 3 (S) is cDNA columns 9 to 11.
    const m = r.view.paths[f.routes.path[0]].map;
    expect([f.routes.forward[0], m.span, Coords.toAt(m, 3, 0), Coords.toAt(m, 1, 2)]).toEqual([1, 3, 9, 5]);
    // Column 2 is a gap in row 0: it maps nowhere (the viewer takes the nearest mapped column).
    expect(Coords.toAt(m, 2, 0)).toBe(-1);
    expect(routesOf(r.view, r.set)).toEqual([
      ["v.protein", 0, 0, 1, 3],
      ["v.coding-dna", 0, 0, 0, 3],
    ]);
  });

  it("follow_anchor on the protein section, or the reference row, picks the row", async () => {
    for (const doc of [pair({ anchor: "prot/1" }), pair({ reference: 1 })]) {
      const r = await resolve(doc, "v");
      expect(r.issues).toEqual([]);
      const m = r.view.paths[r.view.follows.routes.path[0]].map;
      // Through row 1 column 2 is K, cDNA columns 6 to 8.
      expect(Coords.toAt(m, 2, 0)).toBe(6);
    }
  });

  it("follow_anchor on a blocks axis, or a row of another alignment, is svs.section.follow", async () => {
    const bad = pair({ anchor: "cdna/0" });
    const r = await resolve(bad, "v");
    expect(r.issues.map((i) => [i.code, i.message])).toEqual([
      ["svs.section.follow", `"cdna/0" is not a row of this section's alignment, nor a unit`],
    ]);
    // Its side has no anchor, so no route from it; the cDNA leads through its own reference row.
    expect([...r.view.follows.routes.path].map((p) => p >= 0)).toEqual([false, true]);
    const { b, g, p } = gene();
    const v = b.view({ ref: "v", label: "V" });
    const a = v.section({ ref: "v.a", follow_anchor: "protein/unit" }).block({ sequence: g });
    a.track().letters();
    v.section({ ref: "v.b", follow: a.handle }).block({ sequence: p }).track().letters();
    const q = await resolve(b.getDocument(), "v");
    expect(q.issues.map((i) => i.message)).toEqual([
      "follow_anchor anchors an alignment axis; a blocks axis follows the block under its center",
    ]);
    expect(q.view.follows.count).toBe(1);
  });

  it("a lazy alignment (case 9): anchored on a held row, the chains follow through it; on a row not held, pending", async () => {
    for (const [row, routes] of [
      [PF00042.hbaRow, [true, true, true]],
      [1, [false, true, true]],
    ] as const) {
      const doc = case9();
      const v = doc.root.children?.find((c) => c.ref === "msa");
      const [msa, chains] = v?.children ?? [];
      if (!msa || !chains) throw new Error("case 9 has two sections");
      msa.params = { ...msa.params, follow_anchor: `pf00042/${row}` };
      chains.params = { ...chains.params, follow: msa.ref };
      const { set } = await compile(doc, reg, fixtureLoaders());
      const r = resolveView(doc, set, "msa", reg, all);
      expect(codes(r.issues)).toEqual([]);
      expect([...r.view.follows.routes.path].map((p) => p >= 0)).toEqual(routes);
      // The chains' unit A reaches the columns through HBA's row; unit B through HBB's.
      if (row !== 1) expect(r.view.paths[r.view.follows.routes.path[0]].map.count).toBeGreaterThan(0);
    }
  });
});
