/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `follow` behavior (LK1; seqspec:viewer#behaviors, seqspec:format#sections): a section's window set from another's
// through a route at a fraction of the cell, both ways; across a codon step (a protein over its gene, case 14's), a
// reversed piece, an intron (the nearest mapped element), an alignment's anchor row; chains and scroll groups; no echo
// and no extra commit; a resize or a late section set from the last leader, which stays; removable. The route
// arithmetic on its own first.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { type CoordinateSystem, Coords } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { imageBackward, imageForward } from "../src/follow.js";
import { SeqViewerConfig } from "../src/index.js";
import { type Driven, driven } from "./input-fixtures.js";

const SYS: CoordinateSystem = { kind: 1, id: 0, row: -1 };
const TO: CoordinateSystem = { kind: 1, id: 1, row: -1 };
const io = new Float64Array(2);
type Cmap = ReturnType<typeof Coords.create>;
const fwd = (m: Cmap, u: number) => {
  io[0] = u;
  imageForward(m, io);
  return io[0];
};
const back = (m: Cmap, u: number) => {
  io[0] = u;
  imageBackward(m, io);
  return io[0];
};

describe("the route's arithmetic: positions with their fraction, through a map and back", () => {
  it("identity, gapped and minus strand: an element's fraction is kept, read in index order", () => {
    const gapped = Coords.create(SYS, TO, 1, [0, 10, 100, 1, 10, 5, 200, 1]);
    expect([fwd(gapped, 3.25), fwd(gapped, 12.5), back(gapped, 203.75)]).toEqual([103.25, 202.5, 13.75]);
    const minus = Coords.create(SYS, TO, 1, [0, 10, 50, -1]);
    // Element 0 is position 59; its low edge (index order) is 59's high edge.
    expect([fwd(minus, 0), fwd(minus, 0.25), fwd(minus, 9.5), back(minus, 59.75), back(minus, 50.5)]).toEqual([
      60, 59.75, 50.5, 0.25, 9.5,
    ]);
  });

  it("a codon map (span 3) across an intron: a residue's fraction is its codon's; nowhere is the nearest element", () => {
    // Case 14's: 34 coding nucleotides from 12, then 29 from 76.
    const codon = Coords.create(SYS, TO, 3, [0, 34, 12, 1, 34, 29, 76, 1]);
    expect([fwd(codon, 0), fwd(codon, 2.5), fwd(codon, 11 + 2 / 3), fwd(codon, 20.5)]).toEqual([
      12,
      19.5,
      77,
      76 + 61.5 - 34,
    ]);
    expect([back(codon, 12), back(codon, 19.5), back(codon, 77)]).toEqual([0, 2.5, 11 + 2 / 3]);
    // In the intron (46 to 76) the nearest coding nucleotide's middle: 45 below 60, 76 above 62.
    expect([back(codon, 60.2), back(codon, 62)]).toEqual([33.5 / 3, 34.5 / 3]);
    // Past the ends: the first or the last element's middle.
    expect([back(codon, 3), back(codon, 110), fwd(codon, -4), fwd(codon, 40)]).toEqual([
      0.5 / 3,
      62.5 / 3,
      12.5,
      104.5,
    ]);
  });

  it("round trips: every quarter position of a gapped minus-strand map comes back", () => {
    const m = Coords.create(SYS, TO, 1, [0, 7, 30, -1, 7, 4, 10, -1, 11, 9, 50, 1]);
    for (let u = 0; u < 20; u += 0.25) expect(back(m, fwd(m, u))).toBeCloseTo(u, 9);
  });

  it("an empty map gives nothing", () => {
    const m = Coords.create(SYS, TO, 1, []);
    io[0] = 3;
    expect([imageForward(m, io), imageBackward(m, io), io[0]]).toEqual([false, false, 3]);
  });
});

/** Case 14's gene and its translation; sections as asked (ref → follow, direction, scroll group, block). */
function geneDoc(
  sections: [string, { of: "gene" | "protein"; follow?: string; direction?: -1; group?: string }][],
): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const GENE =
    "GAATTCGCCACC" +
    "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC" +
    "GTAAGTCTTTCCTAACTTCTGTTTCCACAG" +
    "ATTGGTACTTCAACAGCCAGGTGGAGTAA" +
    "GCGGCCGCATAAAGT";
  const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: GENE.repeat(1) });
  const protein = translate(b.data, {
    ref: "protein",
    from: gene,
    cds: [
      [12, 46],
      [76, 105],
    ],
  });
  const v = b.view({ ref: "v", label: "V" });
  for (const [ref, o] of sections) {
    v.section({
      ref,
      ...(o.follow ? { follow: o.follow } : {}),
      ...(o.direction ? { direction: o.direction } : {}),
      ...(o.group ? { scroll_group: o.group } : {}),
    })
      .block({ sequence: o.of === "gene" ? gene : protein })
      .track()
      .letters();
  }
  return b.getDocument();
}

/** Every section 360 px wide (the protein's 21 residues at the largest cell are wider), with every behavior. */
function drive(doc: SeqViewSpecDoc) {
  const size: Record<string, readonly [number, number]> = {};
  for (const v of doc.root.children ?? [])
    for (const s of v.children ?? []) if (s.kind === "section" && s.ref) size[s.ref] = [360, 200];
  return driven(doc, { labels: false, size });
}

/** A section's center: its block, its block-system position, and px per shown position. */
function center(t: Driven, ref: string) {
  const out = new Float64Array(2);
  const b = t.vp(ref).centerElement(out);
  return { b, u: out[0], px: out[1] };
}
/** The tolerance of a center: scroll snaps to device px (dpr 1: half a px either way). */
const near = (t: Driven, ref: string, u: number) => {
  const c = center(t, ref);
  expect(Math.abs(c.u - u)).toBeLessThanOrEqual(0.5 / c.px + 1e-9);
};
/** Through the codon map: a residue position to its nucleotide position. */
const toGene = (r: number) => (3 * r < 34 ? 12 + 3 * r : 76 + 3 * r - 34);

describe("follow: a protein and its gene through the codon step", () => {
  const doc = () =>
    geneDoc([
      ["v.dna", { of: "gene" }],
      ["v.protein", { of: "protein", follow: "v.dna" }],
    ]);

  it("the gene leads: the protein centers the residue whose codon is at the gene's center, its cell three times", async () => {
    const t = await drive(doc());
    t.vp("v.dna").zoomTo(0, 14, 34);
    const d = center(t, "v.dna");
    const p = center(t, "v.protein");
    expect(t.vp("v.protein").live.current.origin).toBe("behavior");
    expect(p.px).toBeCloseTo(Math.min(28, 3 * d.px), 9);
    near(t, "v.protein", (d.u - 12) / 3);
  });

  it("the protein leads (a drag, then a Ctrl wheel): the gene follows both ways, at a third of the cell", async () => {
    const t = await drive(doc());
    t.vp("v.protein").zoomTo(0, 5, 18);
    const y = t.layout("v.protein").scrolledTop + 8.5;
    t.drag("v.protein", 300.5, y, 240.5, y);
    const p = center(t, "v.protein");
    expect(t.vp("v.dna").live.current.cell).toBeCloseTo(p.px / 3, 9);
    near(t, "v.dna", toGene(p.u));
    t.wheel("v.protein", 250.5, y, 0, -40, { ctrl: true });
    const q = center(t, "v.protein");
    expect(q.px).not.toBeCloseTo(p.px, 3);
    expect(center(t, "v.dna").px).toBeCloseTo(q.px / 3, 9);
    near(t, "v.dna", toGene(q.u));
  });

  it("the gene's center in the intron: the protein stays on the nearest coding nucleotide's residue", async () => {
    const t = await drive(doc());
    t.vp("v.dna").zoomTo(0, 55, 65);
    expect(center(t, "v.dna").u).toBeCloseTo(60, 0);
    // 45 is nearer than 76: residue 11's first codon nucleotide, its middle.
    near(t, "v.protein", 33.5 / 3);
  });

  it("no echo: a drag moves each section once a move; one commit each, the follower's origin `behavior`", async () => {
    const t = await drive(doc());
    const [dna, protein] = [t.vp("v.dna"), t.vp("v.protein")];
    dna.zoomTo(0, 20, 60);
    t.advance(1000);
    const commits: string[] = [];
    dna.committed.subscribe((c) => commits.push(`dna ${c.origin}`));
    protein.committed.subscribe((c) => commits.push(`protein ${c.origin}`));
    commits.length = 0;
    const [v0, w0] = [dna.live.version.value, protein.live.version.value];
    const y = t.layout("v.dna").scrolledTop + 8.5;
    t.pointer("down", "v.dna", 300.5, y);
    for (let i = 1; i <= 6; i++) t.pointer("move", "v.dna", 300.5 - 7 * i, y);
    // At most one change of the follower per move of the leader (a move inside one device px changes nothing there).
    const moves = dna.live.version.value - v0;
    expect(protein.live.version.value - w0).toBeGreaterThan(moves - 2);
    expect(protein.live.version.value - w0).toBeLessThanOrEqual(moves);
    expect(protein.live.current.interacting).toBe(true);
    t.pointer("up", "v.dna", 258.5, y);
    expect(commits).toEqual(["dna pointer"]);
    t.advance(t.ctx.config.get(SeqViewerConfig.Viewport.SettleMs) + 1);
    expect(commits).toEqual(["dna pointer", "protein behavior"]);
    // An api move commits the leader at once and the follower once it settles.
    commits.length = 0;
    dna.zoomTo(0, 30, 50);
    t.advance(1000);
    expect(commits).toEqual(["dna api", "protein behavior"]);
    // A vertical scroll leads nothing.
    const w1 = protein.live.version.value;
    dna.scrollTo(dna.live.current.scrollX, 3);
    expect(protein.live.version.value).toBe(w1);
  });

  it("a resize of the follower re-centers it on the leader, which stays", async () => {
    const t = await drive(doc());
    t.vp("v.dna").zoomTo(0, 20, 40);
    const d = center(t, "v.dna");
    const x = t.vp("v.dna").live.current.scrollX;
    t.vp("v.protein").setSize(300, 200, 1);
    expect(t.vp("v.dna").live.current.scrollX).toBe(x);
    near(t, "v.protein", (d.u - 12) / 3);
  });

  it("removed, nothing follows", async () => {
    const t = await drive(doc());
    t.ctx.behaviors.get("follow")?.remove();
    const x = t.vp("v.protein").live.current.scrollX;
    t.vp("v.dna").zoomTo(0, 20, 40);
    expect(t.vp("v.protein").live.current.scrollX).toBe(x);
  });
});

describe("follow: pieces, chains, scroll groups", () => {
  it("a reversed piece follows a forward one: the same element at both centers, both ways", async () => {
    const t = await drive(
      geneDoc([
        ["v.fwd", { of: "gene" }],
        ["v.rev", { of: "gene", direction: -1, follow: "v.fwd" }],
      ]),
    );
    t.vp("v.fwd").zoomTo(0, 10, 30);
    near(t, "v.rev", center(t, "v.fwd").u);
    // Dragging the reversed one right shows lower indices on its left: the forward one moves on.
    const before = center(t, "v.fwd").u;
    const y = t.layout("v.rev").scrolledTop + 8.5;
    t.drag("v.rev", 220.5, y, 300.5, y);
    expect(center(t, "v.fwd").u).toBeGreaterThan(before);
    near(t, "v.fwd", center(t, "v.rev").u);
    expect(center(t, "v.fwd").px).toBeCloseTo(center(t, "v.rev").px, 9);
  });

  it("a chain, A ← B ← C: either end leads, the middle passes it on", async () => {
    const t = await drive(
      geneDoc([
        ["v.a", { of: "gene" }],
        ["v.b", { of: "gene", direction: -1, follow: "v.a" }],
        ["v.c", { of: "protein", follow: "v.b" }],
      ]),
    );
    t.vp("v.a").zoomTo(0, 15, 35);
    near(t, "v.b", center(t, "v.a").u);
    near(t, "v.c", (center(t, "v.b").u - 12) / 3);
    t.vp("v.c").zoomTo(0, 3, 9);
    near(t, "v.b", toGene(center(t, "v.c").u));
    near(t, "v.a", center(t, "v.b").u);
  });

  it("a scroll group's member follows, and its group moves with it; a member's move leads the follow", async () => {
    const t = await drive(
      geneDoc([
        ["v.dna", { of: "gene" }],
        ["v.p1", { of: "protein", follow: "v.dna", group: "g" }],
        ["v.p2", { of: "protein", group: "g" }],
      ]),
    );
    t.vp("v.dna").zoomTo(0, 14, 34);
    const p1 = t.vp("v.p1").live.current;
    expect([t.vp("v.p2").live.current.scrollX, t.vp("v.p2").live.current.cell]).toEqual([p1.scrollX, p1.cell]);
    near(t, "v.p1", (center(t, "v.dna").u - 12) / 3);
    const y = t.layout("v.p2").scrolledTop + 8.5;
    t.drag("v.p2", 300.5, y, 240.5, y);
    expect(t.vp("v.p1").live.current.scrollX).toBe(t.vp("v.p2").live.current.scrollX);
    near(t, "v.dna", toGene(center(t, "v.p2").u));
  });
});

describe("scroll groups and `graph` changes (LB-05)", () => {
  it("a member fitted by a data edit is set again from the member a person moved last, which keeps its window", async () => {
    const [A, B] = ["v.a", "v.b"];
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "s", value: "MKVLATGHEKMKVLATGHEK".repeat(20) });
    b.data.sequence({ ref: "s2", value: "MKVLATGHEKMKVLATGHEK".repeat(20) });
    const v = b.view({ ref: "v", label: "V" });
    v.section({ ref: A, scroll_group: "g" } as never)
      .block({ sequence: "s" } as never)
      .track()
      .letters();
    v.section({ ref: B, scroll_group: "g" } as never)
      .block({ sequence: "s2" } as never)
      .track()
      .letters();
    const size = { [A]: [600, 200], [B]: [600, 200] } as Record<string, readonly [number, number]>;
    const t = await driven(b.getDocument(), { labels: false, size });
    t.advance(2000);
    t.vp(B).zoomTo(0, 150, 170, { origin: "pointer" });
    t.advance(2000);
    const cell = t.vp(B).live.current.cell;
    // A's sequence shrinks to 100: its view no longer shows, so it fits (origin `graph`).
    await t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "s", params: { value: "MKVLATGHEK".repeat(10) } }],
    });
    t.advance(2000);
    t.vp(B).centerElement(io);
    expect([Math.round(io[0]), t.vp(B).live.current.cell]).toEqual([160, cell]);
    expect(t.vp(A).live.current.cell).toBe(cell);
    t.ctx.dispose();
  });
});

describe("follow: what loads later, and alignments", () => {
  it("a section added later opens where its leader is; the leader stays", async () => {
    const t = await drive(geneDoc([["v.dna", { of: "gene" }]]));
    t.vp("v.dna").zoomTo(0, 20, 40);
    const x = t.vp("v.dna").live.current.scrollX;
    const d = center(t, "v.dna");
    const next = geneDoc([
      ["v.dna", { of: "gene" }],
      ["v.protein", { of: "protein", follow: "v.dna" }],
    ]);
    const done = t.ctx.setDocument(next);
    t.vp("v.protein").setSize(360, 200, 1);
    expect((await done).ok).toBe(true);
    expect(t.vp("v.dna").live.current.scrollX).toBe(x);
    near(t, "v.protein", (d.u - 12) / 3);
  });

  it("a protein alignment and its cDNA's (codingPair): the cDNA's columns follow through the reference row", async () => {
    const b = createSeqViewSpec();
    const rep = (x: string) => x.repeat(6);
    const prot = b.data.alignment({ ref: "prot", alphabet: "protein", value: [rep("MA-SKL"), rep("MAKSKL")] });
    const cdna = b.data.alignment({
      ref: "cdna",
      alphabet: "dna",
      value: [rep("ATGGCT---TCCAAGCTG"), rep("ATGGCTAAGTCCAAGCTG")],
    });
    const v = b.view({ ref: "v", label: "V" });
    const pair = b.preset.codingPair(v, prot, cdna);
    pair.protein.track().rows();
    pair.cdna.track().rows();
    const t = await drive(b.getDocument());
    const [P, C] = ["v.protein", "v.coding-dna"];
    // Six columns a period, eighteen in the cDNA. Protein column 16 (row 0's K) is cDNA columns 48 to 50.
    t.vp(P).zoomTo(0, 16, 17);
    near(t, C, 3 * (center(t, P).u - 16) + 48);
    expect(center(t, C).px).toBeCloseTo(center(t, P).px / 3, 9);
    // cDNA column 51 is row 0's L codon's first nucleotide: protein column 17.
    t.vp(C).zoomTo(0, 51, 54);
    near(t, P, 17 + (center(t, C).u - 51) / 3);
    // Column 14 is a gap in row 0: the nearest residue's (A, before it) last nucleotide, column 41, its middle.
    t.vp(P).zoomTo(0, 14, 15);
    expect(center(t, C).u).toBeCloseTo(41.5, 0);
  });
});

describe("follow: a document's focus change on a linked section (GG-03)", () => {
  const A = "v.a";
  const B = "v.b";
  function doc(): SeqViewSpecDoc {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: "s", value: "MKVLATGHEKMKVLATGHEK".repeat(20) });
    const v = b.view({ ref: "v", label: "V" });
    const a = v.section({ ref: A });
    a.block({ sequence: s });
    a.track().letters();
    a.focus({ ref: "fa", target: "s", start: 20, stop: 40 } as never);
    const sb = v.section({ ref: B, follow: A } as never);
    sb.block({ sequence: s });
    sb.track().letters();
    sb.focus({ ref: "fb", target: "s", start: 300, stop: 320 } as never);
    return b.getDocument();
  }
  const center = (t: Driven, ref: string) => {
    t.vp(ref).centerElement(io);
    return Math.round(io[0]);
  };
  const focus = (t: Driven, ref: string, start: number) =>
    t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref, params: { target: "s", start, stop: start + 20 } }],
    });

  it("moves its section, which then leads: the follower's own, and the leader's after the follower led", async () => {
    const size = { [A]: [600, 200], [B]: [600, 200] } as Record<string, readonly [number, number]>;
    const t = await driven(doc(), { labels: false, size });
    t.advance(2000);
    expect([center(t, A), center(t, B)]).toEqual([30, 30]);
    // The follower's focus: it moves there, and the leader follows it.
    await focus(t, "fb", 100);
    t.advance(2000);
    expect([center(t, A), center(t, B)]).toEqual([110, 110]);
    // The follower leads once (a host's move), then the leader's own focus changes: the leader moves there.
    t.vp(B).zoomTo(0, 150, 170);
    t.advance(2000);
    expect([center(t, A), center(t, B)]).toEqual([160, 160]);
    await focus(t, "fa", 250);
    t.advance(2000);
    expect([center(t, A), center(t, B)]).toEqual([260, 260]);
    t.ctx.dispose();
  });

  // LB-03: a reassembly that moves the focus's block to another index changes no focus: nothing moves or leads.
  const inserted = (t: Driven, section: string) =>
    t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [
        {
          op: "insert",
          parent: section,
          after: "first",
          node: { kind: "block", ref: `${section}.b2`, params: { sequence: "s2" } },
        },
      ],
    } as never);
  const withShort = () => {
    const d = doc();
    const data = d.root.children?.find((c) => c.kind === "data");
    data?.children?.push({ kind: "sequence", ref: "s2", params: { value: "MKVLATGHEK".repeat(5) } } as never);
    return d;
  };
  it("a block inserted before the follower's focus block moves nothing: the person's position stays", async () => {
    const size = { [A]: [600, 200], [B]: [600, 200] } as Record<string, readonly [number, number]>;
    const t = await driven(withShort(), { labels: false, size });
    t.advance(2000);
    t.vp(A).zoomTo(0, 150, 170, { origin: "pointer" });
    t.advance(2000);
    expect([center(t, A), center(t, B)]).toEqual([160, 160]);
    const moved: string[] = [];
    t.ctx.focus.changed.subscribe((c) => moved.push(c.section));
    const r = await inserted(t, B);
    t.advance(2000);
    expect([r.ok, JSON.stringify(r.issues)]).toEqual([true, "[]"]);
    expect(moved).toEqual([]);
    expect(t.vp(B).centerElement(io)).toBe(1);
    expect([center(t, A), Math.round(io[0])]).toEqual([160, 160]);
    // The entry follows its block to the new index.
    expect(t.ctx.focus.current.value.get(B)?.block).toBe(1);
    t.ctx.dispose();
  });

  it("so does one before an unlinked section's focus block", async () => {
    const d = withShort();
    const view = d.root.children?.find((c) => c.kind === "view");
    const b = view?.children?.find((c) => c.ref === B) as { params: Record<string, unknown> };
    delete b.params.follow;
    const t = await driven(d, { labels: false, size: { [A]: [600, 200], [B]: [600, 200] } });
    t.advance(2000);
    t.vp(B).zoomTo(0, 150, 170, { origin: "pointer" });
    t.advance(2000);
    await inserted(t, B);
    t.advance(2000);
    expect(t.vp(B).centerElement(io)).toBe(1);
    expect(Math.round(io[0])).toBe(160);
    t.ctx.dispose();
  });
});
