/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Labels (seqspec:viewer#labels): the position provider at genome scale (a chromosome reaches hundreds of proteins
// through codon mappings, and composing a path to each is O(the chromosome) on its first hover; paths are composed only
// to molecules with a named numbering), with a positive control; then when a label is computed: once per hovered
// element, again for a new provider or view, never on a move; then the value line for each kind of data.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import type * as ViewSpec from "@molstar/seq-core/format";
import type { KeyTarget, SeqSet } from "@molstar/seq-core/model";
import { describe, expect, it, vi } from "vitest";
import { case10 } from "../../core/_test/format/cases.js";
import { Hit, HitRegion, HoverKind, type LabelProvider } from "../src/index.js";
import { PreparedParts } from "../src/prepare.js";
import { type Driven, driven } from "./input-fixtures.js";
import { renderers } from "./layout-fixtures.js";

const calls = { pathTo: 0 };
vi.mock("@molstar/seq-core/format", async (original) => {
  const m = await original<typeof ViewSpec>();
  return {
    ...m,
    pathTo: (...a: Parameters<typeof m.pathTo>) => {
      calls.pathTo++;
      return m.pathTo(...a);
    },
  };
});

const CDS = "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCTAA";
/** A DNA with `n` coding ranges, each translated into a protein numbered by ordinal only (case 14's form). */
function chromosome(n: number, named?: number) {
  const b = createSeqViewSpec();
  const dna = b.data.sequence({
    ref: "chr",
    alphabet: "dna",
    value: Array.from({ length: n }, () => CDS).join("GG"),
  });
  for (let k = 0; k < n; k++) {
    const at = k * (CDS.length + 2);
    const p = translate(b.data, { ref: `p${k}`, from: dna, cds: [[at, at + CDS.length]] });
    // The positive control: one protein with a named numbering, which a label shows and needs a path to.
    if (k === named)
      p.numbering({ name: "kabat", method: "labels", labels: Array.from({ length: 11 }, (_, i) => `${101 + i}`) });
  }
  b.view({ ref: "v", label: "Chromosome" }).block({ sequence: dna }).track({ label: "DNA" }).letters();
  return b.getDocument();
}

/** The hover at DNA position `index` of the chromosome, as the pointer would set it, and its label. */
async function labelAt(doc: ReturnType<typeof chromosome>, index: number) {
  const o = await driven(doc);
  const set = o.ctx.graph.assembled.value?.set as SeqSet;
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  set.byName("chr", t);
  const unit = [...Array(set.units.count).keys()].find((u) => set.units.molecule[u] === t.id) ?? -1;
  const hit = Hit.create();
  Object.assign(hit, {
    section: 0,
    region: HitRegion.data,
    kind: HoverKind.seq,
    layer: 0,
    block: 0,
    shown: index,
    index,
  });
  Object.assign(hit.seq, { set, unit, index });
  calls.pathTo = 0;
  o.ctx.hover.setFromHit(hit, "pointer");
  return { o, label: o.ctx.labels.compute(o.ctx.hover.live.current) };
}

describe("the position label provider over many mapped proteins", () => {
  it("composes no path to a protein without a named numbering", async () => {
    const { label } = await labelAt(chromosome(40), 4);
    expect(label?.lines.map((l) => l.key)).toEqual(["ordinal"]);
    expect(calls.pathTo).toBe(0);
  });

  it("positive control: a protein with a named numbering gets its path, once, and its line", async () => {
    const { o, label } = await labelAt(chromosome(40, 0), 4);
    expect(calls.pathTo).toBeGreaterThan(0);
    expect(label?.lines).toContainEqual({ key: "kabat", value: "102" });
    // Composed once per unit and set: a second element on the same unit reads the cache.
    const composed = calls.pathTo;
    const hit = Hit.create();
    Object.assign(hit, {
      section: 0,
      region: HitRegion.data,
      kind: HoverKind.seq,
      layer: 0,
      block: 0,
      shown: 7,
      index: 7,
    });
    Object.assign(hit.seq, o.ctx.hover.live.current.seq, { index: 7 });
    o.ctx.hover.setFromHit(hit, "pointer");
    expect(o.ctx.labels.compute(o.ctx.hover.live.current)?.lines).toContainEqual({ key: "kabat", value: "103" });
    expect(calls.pathTo).toBe(composed);
  });
});

describe("when a label is computed", () => {
  const S = "chains.main";
  /** A provider that counts its calls and titles the label with the count. */
  const counting = (priority = 2000) => {
    const seen = { n: 0 };
    const provider: LabelProvider = {
      name: "lab:count",
      priority,
      applies: (h) => h.kind === HoverKind.seq,
      label: (_h, _lc, out) => out.title(`count ${++seen.n}`),
    };
    return { seen, provider };
  };

  it("once per hovered element, never on a move; again for a provider registered under a still pointer", async () => {
    const { seen, provider } = counting();
    const t = await driven(case10(), { registry: [{ labels: [provider] }] });
    const y = t.yOf(S, "chains.residues");
    const x = t.xOf(S, 0, 100);
    t.pointer("move", S, x, y);
    expect(t.ctx.labels.current.value?.title).toBe("count 1");
    for (let i = 1; i <= 5; i++) t.pointer("move", S, x + i * 0.05, y);
    expect(seen.n).toBe(1);
    t.pointer("move", S, t.xOf(S, 0, 101), y);
    expect([seen.n, t.ctx.labels.current.value?.title]).toEqual([2, "count 2"]);
    // A provider registered while the pointer is still: the label is made again, with it.
    const other: LabelProvider = {
      name: "lab:other",
      priority: 3000,
      applies: (h) => h.kind === HoverKind.seq,
      label: (_h, _lc, out) => out.title("other"),
    };
    const undo = t.ctx.register({ labels: [other] });
    expect([t.ctx.labels.current.value?.title, seen.n]).toEqual(["other", 3]);
    undo();
    expect([t.ctx.labels.current.value?.title, seen.n]).toEqual(["count 4", 4]);
  });

  it("a numbering switch (a patch, so a new view) recomputes the label under a still pointer", async () => {
    const t = await driven(case10());
    const y = t.yOf(S, "chains.residues");
    const x = t.xOf(S, 0, 100);
    t.pointer("move", S, x, y);
    const before = t.ctx.labels.current.value?.title;
    expect(before).toBeTruthy();
    await t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: S, params: { numbering: "uniprot:P00519-2" } }],
    });
    const after = t.ctx.labels.current.value?.title;
    expect(after).toBeTruthy();
    expect(after).not.toBe(before);
  });
});

describe("the value line: what the layer draws, read from its prepared state", () => {
  /** The fixture renderers, with the real prepare for the layers the value provider reads (LY2's `PerPosition`). */
  const prepared = () =>
    renderers().map((r) =>
      r.name === "values" || r.name === "color" ? { ...r, prepare: PreparedParts[r.name].prepare } : r,
    );
  const ann = (b: ReturnType<typeof createSeqViewSpec>, ref: string, target: unknown, value: unknown) =>
    b.data.annotation({ ref: `${ref}-data`, label: ref, geometry: "values", target: target as never, value } as never);
  /**
   * A chain of 36 with dense values (i / 10) and one sparse value on its molecule (0.5 at 5), and a 12-residue chain
   * mapped onto it at [10, 22) with values 100 + i.
   */
  function chains() {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", label: "A", value: "MKVLATGHEKLV".repeat(3) });
    const q = b.data.sequence({ ref: "q", label: "Q", value: "MKVLATGHEKLV" });
    b.data.mapping({ ref: "m", from: q, to: a, method: "segments", segments: [[0, 10, 12]] });
    const v = b.view({ ref: "v", label: "V" }).block({ sequence: a });
    v.track({ ref: "dense" }).values({
      data: ann(b, "dense", a, { value: Array.from({ length: 36 }, (_, i) => i / 10) }),
    });
    v.track({ ref: "sparse" }).color({ data: ann(b, "sparse", a, { index: [5], value: [0.5] }) });
    const mapped = { value: Array.from({ length: 12 }, (_, i) => 100 + i) };
    v.track({ ref: "mapped" }).values({ data: ann(b, "mapped", q, mapped) });
    return b.getDocument();
  }
  /** Case 14's gene (prepare-fixtures.ts `gene`) with a value per protein residue (k + 1), drawn through codons. */
  function codons() {
    const b = createSeqViewSpec();
    const GENE =
      "GAATTCGCCACCATGGCTTCCAAGGAACTGACCCCACGTGACGGCCGTAAGTCTTTCCTAACTTCTGTTTCCACAGATTGGTACTTCAACAGCCAGGTGGAGTAAGCGGCCGCATAAAGT";
    const g = b.data.sequence({ ref: "gene", label: "Gene", alphabet: "dna", value: GENE });
    const cds: [number, number][] = [
      [12, 46],
      [76, 105],
    ];
    const p = translate(b.data, { ref: "protein", from: g, cds });
    const v = b.view({ ref: "g", label: "G" }).block({ sequence: g });
    v.track({ ref: "codon" }).values({
      data: ann(b, "codon", p, { value: Array.from({ length: 21 }, (_, k) => k + 1) }),
    });
    return b.getDocument();
  }
  const valuesAt = (t: Driven, section: string, track: string, index: number) => {
    t.pointer("move", section, t.xOf(section, 0, index), t.yOf(section, track));
    const h = t.ctx.hover.live.current;
    expect([h.kind, h.index]).toEqual([HoverKind.seq, index]);
    return (t.ctx.labels.current.value?.lines ?? []).filter((l) => l.key === track).map((l) => l.value);
  };

  it("dense and molecule-target data, data mapped from another chain and through codons, each as drawn", async () => {
    const t = await driven(chains(), { renderers: prepared() });
    expect([3, 17].map((i) => valuesAt(t, "v.main", "dense", i))).toEqual([["0.3"], ["1.7"]]);
    expect([5, 6].map((i) => valuesAt(t, "v.main", "sparse", i))).toEqual([["0.5"], []]);
    expect([5, 10, 15, 21, 22].map((i) => valuesAt(t, "v.main", "mapped", i))).toEqual([
      [],
      ["100"],
      ["105"],
      ["111"],
      [],
    ]);
    const g = await driven(codons(), { renderers: prepared() });
    // The UTR has none; a codon's three nucleotides read its residue's, the codon split by the intron too (45, 76, 77).
    expect([11, 12, 14, 15, 45, 76, 77, 78].map((i) => valuesAt(g, "g.main", "codon", i))).toEqual([
      [],
      ["1"],
      ["1"],
      ["2"],
      ["12"],
      ["12"],
      ["12"],
      ["13"],
    ]);
  });
});
