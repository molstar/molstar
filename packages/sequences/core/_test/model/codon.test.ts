/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Codon mappings (seqspec:model#codon-maps): case 14's gene and protein, statuses, strands, genetic codes, paths with
// one codon step, runs and conversions; NCBI's tables against gc.prt itself.

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import {
  type ComposedPath,
  Conversion,
  ConversionStatus,
  type CoordinateMap,
  type CoordinateSystem,
  Coords,
  GeneticCode,
  type Issue,
  type Mapping,
  type MappingInput,
  MappingStatus,
  Mappings,
  OrderedSet,
  type PathStep,
  Paths,
  type SeqIndex,
  SeqLocation,
  SeqLoci,
  type SeqSet,
  SeqSets,
  threeLetter,
  type UnitId,
} from "../../src/model/index.js";
import { CDS, GENE, PROTEIN } from "./codon-fixtures.js";

const COMPLEMENT: Record<string, string> = { A: "T", T: "A", C: "G", G: "C" };
const revcomp = (s: string) =>
  [...s]
    .reverse()
    .map((c) => COMPLEMENT[c] ?? c)
    .join("");

type Codon = Partial<Extract<MappingInput, { method: "codon" }>>;
function gene(codon: Codon = {}, residues = PROTEIN, dna = GENE) {
  const input: MappingInput = {
    key: "protein.codon",
    from: "protein",
    to: "gene",
    method: "codon",
    cds: CDS,
    ...codon,
  };
  return SeqSets.build({
    molecules: [
      { key: "gene", alphabet: "dna", residues: dna, ownUnit: true },
      { key: "protein", alphabet: "protein", residues, ownUnit: true },
    ],
    mappings: [input],
  });
}
const codonMap = (set: SeqSet) => set.mappings[0] as Mapping;
const unit = (set: SeqSet, key: string) => Array.from(set.units.key).indexOf(set.strings.find(key)) as UnitId;
const space = (set: SeqSet, key: string): CoordinateSystem => {
  const s = { kind: -1, id: -1, row: -1 };
  if (!Mappings.systemOf(set, key, s)) throw new Error(`no space ${key}`);
  return s;
};
function compose(set: SeqSet, from: string, to: string): ComposedPath {
  const steps: PathStep[] = [];
  const n = Mappings.path(set, space(set, from), space(set, to), steps);
  if (n < 0) throw new Error(`no path ${from} → ${to}: ${n}`);
  return Mappings.compose(set, space(set, from), steps);
}
/** A range of the path's `from` on its `to`: through the map forward, or backward when its `from` is the path's `to`. */
const runs = (path: ComposedPath, start: number, stop: number) => {
  const out: [number, number][] = [];
  const walk = Coords.sameSystem(path.map.from, path.from) ? Coords.forEachTo : Coords.forEachFrom;
  walk(path.map, start, stop, (s, t, acc) => acc.push([s, t]), out);
  return out;
};
/** The first element holding a to position and its phase ([-1, -1]: none): the old `backward` and `phase`. */
const held = (m: CoordinateMap, to: number) => {
  const h = { element: -1, phase: -1 };
  return Coords.forEachFromAt(m, to, h, () => {}, null) ? [h.element, h.phase] : [-1, -1];
};
const codes = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);
const statuses = (m: Mapping) => Array.from(m.status);

describe("codon mappings", () => {
  it("pairs case 14's protein with its spliced CDS: spans, backward, phase, every status exact", () => {
    const { set, issues } = gene();
    expect(issues).toEqual([]);
    const m = codonMap(set);
    expect([m.method, m.span, m.codon?.strand, m.codon?.geneticCode]).toEqual(["codon", 3, 1, 1]);
    expect(set.strings.get(m.provenance.method)).toBe("translation");
    expect(Array.from(m.spans?.subarray(0, 3) ?? [])).toEqual([12, 13, 14]);
    expect(Array.from(m.spans?.subarray(33, 36) ?? [])).toEqual([45, 76, 77]); // His12, split by the intron
    expect([m.forward[0], m.forward[11], ...[45, 77, 46, 102].map((j) => held(m.map, j)[0])]).toEqual([
      12, 45, 11, 11, -1, -1,
    ]);
    // Phases; the stop codon (102–104) pairs with no residue, so it has none.
    expect([45, 76, 77, 60, 104].map((j) => held(m.map, j)[1])).toEqual([0, 1, 2, -1, -1]);
    expect([m.map.count, Array.from(m.map.at), Array.from(m.map.toLo), [m.map.toStart, m.map.toStop]]).toEqual([
      2,
      [0, 34],
      [12, 76],
      [12, 102],
    ]);
    expect(statuses(m).every((s) => s === MappingStatus.exact)).toBe(true);
  });

  it("places protein ranges as runs of whole codons on the gene, and nucleotides on their residues", () => {
    const { set } = gene();
    const down = compose(set, "protein", "gene/unit");
    expect([down.span, down.heldSpans?.length, Coords.sameSystem(down.map.from, down.from)]).toEqual([3, 60, true]);
    expect(runs(down, 4, 14)).toEqual([
      [24, 46],
      [76, 84],
    ]); // the domain, residues 5–14
    expect(runs(down, 16, 17)).toEqual([[90, 93]]); // pS17
    expect(runs(down, 11, 12)).toEqual([
      [45, 46],
      [76, 78],
    ]);
    const up = compose(set, "gene/unit", "protein/unit");
    // From the nucleotides the protein side is still the map's `from` (span 3); the path's own arrays are span 1.
    expect([up.span, up.map.span, Paths.at(up, 13), Paths.at(up, 5), Paths.at(up, 103)]).toEqual([1, 3, 0, -1, -1]);
    expect([0, 1, 2].map((p) => Coords.toAt(up.map, 11, p))).toEqual([45, 76, 77]);
    expect(runs(up, 40, 80)).toEqual([[9, 13]]);
  });

  it("reads the reverse complement on strand -: the same protein, placements mirrored", () => {
    const mirrored = CDS.map(([s, t]) => [120 - t, 120 - s] as [number, number]);
    const { set, issues } = gene({ cds: mirrored, strand: "-" }, PROTEIN, revcomp(GENE));
    expect(issues).toEqual([]);
    const m = codonMap(set);
    expect(statuses(m).every((s) => s === MappingStatus.exact)).toBe(true);
    expect([m.codon?.strand, Array.from(m.spans?.subarray(0, 3) ?? [])]).toEqual([-1, [107, 106, 105]]);
    const down = compose(set, "protein", "gene/unit");
    expect(runs(down, 4, 14)).toEqual([
      [36, 44],
      [74, 96],
    ]);
    expect(runs(down, 16, 17)).toEqual([[27, 30]]);
  });

  it("gives statuses: a mutation is a conflict, residues past the codons none, a final * the stop", () => {
    const mutant = `${PROTEIN.slice(0, 7)}A${PROTEIN.slice(8)}`; // P8A
    expect(statuses(codonMap(gene({}, mutant).set))[7]).toBe(MappingStatus.conflict);
    const longer = codonMap(gene({}, `${PROTEIN}KK`).set);
    // Residue 21 (K) against the stop codon is a conflict; residue 22, past the codons, none.
    expect(statuses(longer).slice(19)).toEqual([MappingStatus.exact, MappingStatus.conflict, MappingStatus.none]);
    const star = codonMap(gene({}, `${PROTEIN}*`).set);
    expect([star.status[20], star.spans?.[60]]).toEqual([MappingStatus.exact, 102]);
    // A construct coding residues 3–8 of a longer protein (`fromRange`).
    const part = gene({ cds: undefined, start: 18, fromRange: [3, 9] }, `XXX${PROTEIN.slice(2, 8)}`);
    expect(part.issues).toEqual([]);
    expect(statuses(codonMap(part.set))).toEqual([2, 2, 2, 0, 0, 0, 0, 0, 0]);
    expect(codonMap(part.set).forward[3]).toBe(18);
  });

  it("reads starts, Sec, Pyl and ambiguity codes by the table", () => {
    // GTG is a start in table 11, not in table 1: M only as the first coded residue.
    const dna = "GTGGTGTGATAGNCT";
    const build = (residues: string, geneticCode: number) =>
      statuses(codonMap(gene({ cds: [[0, 15]], geneticCode }, residues, dna).set));
    expect(build("MMUOX", 11)).toEqual([0, 1, 0, 0, 0]);
    expect(build("MVUOA", 1)).toEqual([1, 0, 0, 0, 1]);
    // Table 2 (vertebrate mitochondrial) against table 1: TGA is W, AGA a stop, ATA M.
    const mito = "TGAAGAATA";
    const two = (residues: string, geneticCode: number) =>
      statuses(codonMap(gene({ cds: [[0, 9]], geneticCode }, residues, mito).set));
    expect(two("W*M", 2)).toEqual([0, 0, 0]);
    expect(two("W*M", 1)).toEqual([1, 1, 1]);
    expect(two("*RI", 1)).toEqual([0, 0, 0]);
  });

  it("reports a wrong alphabet, ranges out of order or overlapping, a length not a multiple of 3, an unknown code", () => {
    const issue = (codon: Codon, swap = false) => {
      const input = { key: "c", from: swap ? "gene" : "protein", to: swap ? "protein" : "gene", method: "codon" };
      const { issues, set } = SeqSets.build({
        molecules: [
          { key: "gene", alphabet: "dna", residues: GENE },
          { key: "protein", alphabet: "protein", residues: PROTEIN },
        ],
        mappings: [{ ...input, cds: CDS, ...codon } as MappingInput],
      });
      expect(set.mappings[0] ?? null).toBe(null);
      return codes(issues);
    };
    expect(issue({}, true)).toEqual(["model.mapping.codon /mappings/0"]);
    expect(issue({ cds: [CDS[1], CDS[0]] })).toEqual(["model.mapping.codon /mappings/0/cds/1"]);
    expect(issue({ cds: [CDS[0], [40, 69]] })).toEqual(["model.mapping.codon /mappings/0/cds/1"]);
    expect(issue({ cds: [CDS[0], [76, 104]] })).toEqual(["model.mapping.codon /mappings/0/cds"]);
    expect(issue({ geneticCode: 7 })).toEqual(["model.mapping.codon /mappings/0/geneticCode"]);
    expect(issue({ cds: [CDS[0], [100, 121]] })).toEqual(["model.mapping.range /mappings/0/cds/1"]);
    expect(issue({ cds: undefined, start: 90 })).toEqual(["model.mapping.range /mappings/0/start"]);
  });

  it("takes at most one codon step per path", () => {
    const { set } = SeqSets.build({
      molecules: [
        { key: "gene", alphabet: "dna", residues: GENE, ownUnit: true },
        { key: "copy", alphabet: "dna", residues: GENE, ownUnit: true },
        { key: "protein", alphabet: "protein", residues: PROTEIN, ownUnit: true },
        { key: "other", alphabet: "protein", residues: PROTEIN, ownUnit: true },
      ],
      mappings: [
        { key: "a", from: "protein", to: "gene", method: "codon", cds: CDS },
        { key: "b", from: "other", to: "gene", method: "codon", cds: CDS },
        { key: "c", from: "gene", to: "copy", method: "identity" },
      ],
    });
    const steps: PathStep[] = [];
    expect(Mappings.path(set, space(set, "protein"), space(set, "other"), steps)).toBe(-1);
    expect(Mappings.path(set, space(set, "protein"), space(set, "copy/unit"), steps)).toBe(3);
    const p = Mappings.compose(set, space(set, "protein"), steps);
    expect(runs(p, 11, 12)).toEqual([
      [45, 46],
      [76, 78],
    ]);
  });

  it("converts a residue to its three nucleotides (multiplied) and a UTR nucleotide to nothing (dropped)", () => {
    const { set } = gene();
    const [p, g] = [unit(set, "protein/unit"), unit(set, "gene/unit")];
    const c = Conversion.create();
    const loci = (u: UnitId, s: number, t: number): SeqLoci => ({
      kind: "seq-loci",
      set,
      elements: [{ unit: u, indices: OrderedSet.ofRange(s as SeqIndex, t as SeqIndex) }],
    });
    const nts = SeqLoci.toSystem(loci(p, 11, 13), compose(set, "protein/unit", "gene/unit"), c);
    const list = (s: OrderedSet) => Array.from({ length: OrderedSet.size(s) }, (_, i) => OrderedSet.getAt(s, i));
    expect([c.status, c.dropped, c.multiplied, list(nts.elements[0].indices)]).toEqual([
      ConversionStatus.exact,
      0,
      2,
      [45, 76, 77, 78, 79, 80],
    ]);
    const res = SeqLoci.toSystem(loci(g, 10, 16), compose(set, "gene/unit", "protein/unit"), c);
    expect([c.status, c.dropped, c.multiplied, list(res.elements[0].indices)]).toEqual([
      ConversionStatus.partial,
      2,
      0,
      [0, 1],
    ]);
    const l = SeqLocation.set(SeqLocation.create(set), p, 0);
    const seen: number[] = [];
    const n = SeqLocation.toSystemEach(
      l,
      compose(set, "protein/unit", "gene/unit"),
      SeqLocation.create(set),
      (x) => {
        seen.push(x.index);
      },
      null,
    );
    expect([n, seen]).toEqual([3, [12, 13, 14]]);
  });
});

describe("genetic codes", () => {
  const gcPrt = readFileSync(new URL("../fixtures/gc.prt", import.meta.url), "utf8");
  const tables = [
    ...gcPrt
      .split("\n")
      .filter((l) => !l.trim().startsWith("--"))
      .join("\n")
      .matchAll(/\bid (\d+)\s*,\s*ncbieaa\s+"([^"]+)"\s*,\s*sncbieaa\s+"([^"]+)"/g),
  ];
  const BASES = "TCAG";
  const code = (s: string) => s.charCodeAt(0);

  it("translates every table's 64 codons as gc.prt lists them, T and U alike, starts marked", () => {
    expect(tables.length).toBe(27);
    for (const [, id, aas, starts] of tables) {
      expect(GeneticCode.has(Number(id))).toBe(true);
      for (let k = 0; k < 64; k++) {
        const [a, b, c] = [BASES[k >> 4], BASES[(k >> 2) & 3], BASES[k & 3]];
        expect(String.fromCharCode(GeneticCode.translate(Number(id), code(a), code(b), code(c)))).toBe(aas[k]);
        const u = (x: string) => code(x === "T" ? "U" : x);
        expect(GeneticCode.translate(Number(id), u(a), u(b), u(c))).toBe(code(aas[k]));
        expect(GeneticCode.isStart(Number(id), code(a), code(b), code(c))).toBe(starts[k] === "M");
      }
    }
    expect([GeneticCode.has(7), GeneticCode.has(0), GeneticCode.translate(7, 65, 65, 65)]).toEqual([false, false, -1]);
    expect(GeneticCode.translate(1, code("A"), code("N"), code("G"))).toBe(code("X"));
  });

  it("names residues by their three-letter codes", () => {
    expect([..."MUOBZJX*W?"].map((c) => threeLetter(code(c)))).toEqual(
      "Met Sec Pyl Asx Glx Xle Xaa Ter Trp Xaa".split(" "),
    );
  });
});
