/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import {
  AlignedLocation,
  AlignedLoci,
  type AlignmentId,
  type ColumnIndex,
  ColumnLocation,
  ColumnLoci,
  type ComposedPath,
  Conversion,
  ConversionStatus,
  type CoordinateMap,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  Fasta,
  type Issue,
  type MappingInput,
  MappingStatus,
  Mappings,
  type Molecule,
  OrderedSet,
  type PathStep,
  Paths,
  Providers,
  type RowIndex,
  RowLocation,
  RowLoci,
  type SeqIndex,
  SeqLocation,
  SeqLoci,
  type SeqSet,
  SeqSets,
  type SetInput,
  type UnitId,
} from "../../src/model/index.ts";

/** A composed path's first to positions per from element, dense (the sparse path read through `Paths`). */
const dense = (p: ComposedPath) => Array.from({ length: p.length }, (_, i) => Paths.at(p, i));
const denseSpans = (p: ComposedPath) =>
  Array.from({ length: p.span * p.length }, (_, k) => Paths.spanAt(p, (k / p.span) | 0, k % p.span));

import { pf00042, SH2_A3M } from "./alignment-fixtures.ts";
import { a3n1, brs1 } from "./fixtures.ts";
import { alignment as generated } from "./support/generators.ts";

const spike = (name: string) => JSON.parse(readFileSync(new URL(`../fixtures/${name}.json`, import.meta.url), "utf8"));
const codes = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);
/** The first from element holding a to position, -1: none (the old `backward`). */
const back = (m: CoordinateMap | undefined, to: number) => {
  const h = { element: -1, phase: -1 };
  let first = -1;
  if (m) Coords.forEachFromAt(m, to, h, (x) => (first = first < 0 ? x.element : first), null);
  return first;
};
const space = (set: SeqSet, key: string): CoordinateSystem => {
  const s = { kind: -1, id: -1, row: -1 };
  if (!Mappings.systemOf(set, key, s)) throw new Error(`no space ${key}`);
  return s;
};
const mapping = (set: SeqSet, key: string) => {
  const m = set.mappings.find((x) => x && set.strings.get(x.key) === key);
  if (!m) throw new Error(`no mapping ${key}`);
  return m;
};
const unit = (set: SeqSet, key: string) => Array.from(set.units.key).indexOf(set.strings.find(key)) as UnitId;
const molecule = (set: SeqSet, key: string) =>
  set.molecules.find((m) => m && set.strings.get(m.key) === key) as Molecule;
const list = (s: OrderedSet) => Array.from({ length: OrderedSet.size(s) }, (_, i) => OrderedSet.getAt(s, i));
/** A path's composed arrays, or the path's code (-1 none, -2 ambiguous). */
function along(set: SeqSet, from: string, to: string) {
  const steps: PathStep[] = [];
  const n = Mappings.path(set, space(set, from), space(set, to), steps);
  return { n, steps, path: n >= 0 ? Mappings.compose(set, space(set, from), steps) : null };
}
function seqsOf(set: SeqSet, run: (out: SeqLocation, cb: (l: SeqLocation) => void) => number) {
  const seen: [string, number][] = [];
  const n = run(SeqLocation.create(set), (l) => seen.push([set.strings.get(set.units.key[l.unit]), l.index]));
  expect(n).toBe(seen.length);
  return seen;
}

const barnase = spike("1brs").barnase.sequence as string;
const barstar = spike("1brs").barstar.sequence as string;

describe("mappings: methods, statuses and issues", () => {
  it("builds barstar by alignment, with statuses (C40A, C82A)", () => {
    const P11540 = { key: "P11540", residues: barstar, ownUnit: true };
    const draft = SeqSets.build({ ...brs1(), molecules: [...brs1().molecules, P11540] }).set;
    const cigar = "1I89M"; // what the aligners in @molstar/seq-core/algorithms give for barstar against 1BRS
    const input: MappingInput = { key: "sifts", from: "P11540", to: "1brs/entity:2", method: "pairwise", cigar };
    const { set, issues } = SeqSets.withChanges(draft, { upsert: { mappings: [input] } });
    expect(issues).toEqual([]);
    const map = mapping(set, "sifts");
    expect([map.forward[0], map.forward[1], back(map.map, 88)]).toEqual([-1, 0, 89]);
    expect(map.status[0]).toBe(MappingStatus.none);
    const conflicts = Array.from(map.status.keys()).filter((i) => map.status[i] === MappingStatus.conflict);
    expect(conflicts).toEqual([40, 82]);
    expect(set.strings.get(map.provenance.method)).toBe("alignment");
  });

  it("builds ABL1 1a ↔ 1b by segments with an open length, and the 1IEP chain by numbering", () => {
    const abl = spike("1iep").sequence as string;
    const filler = "G".repeat(224 - 26);
    const set = SeqSets.build({
      molecules: [
        { key: "ABL1-1a", residues: `${"A".repeat(26)}${filler}${abl}`, ownUnit: true },
        { key: "ABL1-1b", residues: `${"C".repeat(45)}${filler}${abl}`, ownUnit: true },
        { key: "1iep/entity:1", residues: abl },
      ],
      units: [
        {
          key: "1iep/A",
          molecule: "1iep/entity:1",
          numberings: [{ name: "auth", method: "offset", first: 225, isDefault: true }],
        },
      ],
      mappings: [
        { key: "isoforms", from: "ABL1-1a", to: "ABL1-1b", method: "segments", segments: [[26, 45, null]] },
        {
          key: "chain",
          from: "1iep/A",
          to: "ABL1-1a",
          method: "by_numbering",
          numbering: "auth",
          toNumbering: "ordinal",
        },
      ],
    }).set;
    const iso = mapping(set, "isoforms");
    expect([iso.forward[25], iso.forward[26], iso.forward[497], back(iso.map, 44), back(iso.map, 45)]).toEqual([
      -1, 45, 516, -1, 26,
    ]);
    expect(iso.status.filter((s) => s === MappingStatus.exact).length).toBe(198 + abl.length);
    expect(set.strings.get(iso.provenance.method)).toBe("curated");
    expect(mapping(set, "chain").forward[90]).toBe(314);

    // T315 (1a numbering) is T334 in 1b: unit → 1a → 1b → 1b's own unit.
    const { n, path } = along(set, "1iep/A", "ABL1-1b/unit");
    expect(n).toBe(3);
    const out = SeqLocation.create(set);
    const l = SeqLocation.set(SeqLocation.create(set), unit(set, "1iep/A"), 90);
    expect(path && SeqLocation.toSystem(l, path, out)).toBe(true);
    expect([out.unit, out.index + 1, String.fromCharCode(molecule(set, "ABL1-1b").code[out.index])]).toEqual([
      unit(set, "ABL1-1b/unit"),
      334,
      "T",
    ]);
  });

  it("reports range, CIGAR, endpoint and numbering errors, and an empty mapping", () => {
    const ms: MappingInput[] = [
      { key: "a", from: "1brs/A", to: "1brs/D", method: "segments", segments: [[100, 0, 20]] },
      { key: "b", from: "1brs/A", to: "1brs/D", method: "pairwise", cigar: "100M" },
      { key: "c", from: "1brs/A", to: "1brs/D", method: "pairwise", cigar: "10Q" },
      { key: "d", from: "1brs/A", to: "nothing", method: "identity" },
      { key: "e", from: "1brs/A", to: "1brs/D", method: "by_numbering", numbering: "kabat" },
      { key: "f", from: "1brs/A", to: "1brs/D", method: "identity", fromRange: [0, 100] },
      {
        key: "g",
        from: "1brs/A",
        to: "1brs/D",
        method: "segments",
        segments: [
          [0, 0, 5],
          [3, 10, 5],
        ],
      },
      { key: "h", from: "1brs/A", to: "1brs", method: "identity" },
      { key: "i", from: "1brs/A", to: "1brs/entity:2", method: "pairwise", cigar: "110I89D" },
    ];
    const { set, issues } = SeqSets.build({ ...brs1(), mappings: ms });
    expect(codes(issues)).toEqual([
      "model.mapping.range /mappings/0/segments/0",
      "model.mapping.cigar /mappings/1/cigar",
      "model.mapping.cigar /mappings/2/cigar",
      "model.mapping.endpoint /mappings/3",
      "model.mapping.numbering /mappings/4/numbering",
      "model.mapping.range /mappings/5/fromRange",
      "model.mapping.range /mappings/6/segments/1",
      "model.mapping.endpoint /mappings/7",
      "model.mapping.empty /mappings/8",
    ]);
    expect(set.mappings.filter(Boolean).length).toBe(1);
    expect(set.byKey(set.strings.find("a"), { kind: 0, id: 0, row: 0 })).toBe(false);
  });

  it("keeps ids, rebuilds on an endpoint's change and drops a mapping that lost one", () => {
    const input: MappingInput = { key: "m", from: "1brs/A", to: "1brs/B", method: "identity" };
    const v1 = SeqSets.build({ ...brs1(), mappings: [input] }).set;
    const m1 = mapping(v1, "m");
    const v2 = SeqSets.withChanges(v1, { upsert: { complexes: [{ key: "x", units: [] }] } }).set;
    expect([v2.mappings === v1.mappings, mapping(v2, "m") === m1]).toEqual([true, true]);
    const mutant = { key: "1brs/entity:1", residues: `G${barnase.slice(48)}`, numberings: [] };
    const v3 = SeqSets.withChanges(v2, { upsert: { molecules: [mutant] } }).set;
    const m3 = mapping(v3, "m");
    expect([m3 === m1, m3.id, m3.status[0], m3.status[1]]).toEqual([false, m1.id, MappingStatus.exact, 0]);
    const r = SeqSets.withChanges(v3, { remove: ["1brs/B"] });
    expect([codes(r.issues), r.set.mappings[m1.id]]).toEqual([["model.mapping.orphaned /remove/0"], null]);
    const v5 = SeqSets.withChanges(r.set, { upsert: { mappings: [{ ...input, to: "1brs/C" }] } }).set;
    // ST1b-a, seqspec:model#identity-and-lineage: a key that went and comes back gets its old id.
    expect(mapping(v5, "m").id).toBe(m1.id);
  });
});

/** PF00042 with HBA_HUMAN (row 0, from residue 27) bound to 1A3N's α entity, and SH2 as A3M. */
function globins(extra: Partial<SetInput> = {}) {
  const rows = pf00042();
  const base = a3n1();
  const HBA = base.molecules[0];
  return SeqSets.build({
    ...base,
    molecules: [
      { ...HBA, numberings: [{ name: "uniprot:P69905", method: "offset", first: 2 }] },
      ...base.molecules.slice(1),
      ...(extra.molecules ?? []),
    ],
    units: [...(base.units ?? []), ...(extra.units ?? [])],
    alignments: [
      {
        key: "PF00042",
        flavor: "aligned",
        gapped: rows.map((r) => r.gapped),
        labels: rows.map((r) => r.label),
        rowNumberingFirst: rows.map((r) => r.start),
      },
      Fasta.alignment(Fasta.parse(SH2_A3M).records, "sh2", "a3m").input,
      ...(extra.alignments ?? []),
    ],
    mappings: [
      { key: "hba", from: "PF00042/0", to: "1a3n/entity:1", method: "segments", segments: [[0, 25, null]] },
      ...(extra.mappings ?? []),
    ],
  });
}

describe("rows, columns and units", () => {
  const { set, issues } = globins();
  const msa = 0 as AlignmentId;
  const a = set.alignments[msa] as NonNullable<SeqSet["alignments"][number]>;
  const [A, C] = [unit(set, "1a3n/A"), unit(set, "1a3n/C")];

  it("binds a row by numbering as by segments, and synthesizes row numberings", () => {
    expect(issues.filter((i) => i.severity === "error").map((i) => i.message)).toEqual([]);
    const ns = a.numberingsOf(0 as RowIndex);
    expect([ns?.count, ns?.byName(set.strings.find("row"))?.isDefault]).toEqual([2, true]);
    const more = SeqSets.withChanges(set, {
      upsert: {
        mappings: [
          {
            key: "sifts",
            from: "PF00042/0",
            to: "1a3n/entity:1",
            method: "by_numbering",
            numbering: "row",
            toNumbering: "uniprot:P69905",
          },
        ],
      },
    }).set;
    expect(Array.from(mapping(more, "sifts").forward)).toEqual(Array.from(mapping(set, "hba").forward));
    expect(mapping(set, "hba").status.every((s) => s === MappingStatus.exact)).toBe(true);
  });

  it("converts columns ↔ positions through the row's binding, with gaps", () => {
    const col = ColumnLocation.set(ColumnLocation.create(set), msa, 20);
    expect(seqsOf(set, (out, cb) => ColumnLocation.toSeqLocations(col, out, cb, null))).toEqual([
      ["1a3n/A", 45],
      ["1a3n/C", 45],
    ]);
    ColumnLocation.set(col, msa, 21);
    expect(seqsOf(set, (out, cb) => ColumnLocation.toSeqLocations(col, out, cb, null))).toEqual([]);
    const cell = AlignedLocation.set(AlignedLocation.create(set), msa, 0, 21);
    const res = RowLocation.create(set);
    expect(AlignedLocation.toRow(cell, res)).toBe(false);
    AlignedLocation.set(cell, msa, 0, 25);
    expect([AlignedLocation.toRow(cell, res), res.index]).toEqual([true, 21]);
    expect(seqsOf(set, (out, cb) => RowLocation.toSeqLocations(res, out, cb, null))).toEqual([
      ["1a3n/A", 46],
      ["1a3n/C", 46],
    ]);
    const back = AlignedLocation.create(set);
    expect([RowLocation.toAligned(res, back), back.column]).toEqual([true, 25]);
    const rows: [number, number][] = [];
    const seq = SeqLocation.set(SeqLocation.create(set), C, 46);
    SeqLocation.toRowLocations(seq, msa, RowLocation.create(set), (l) => rows.push([l.row, l.index]), null);
    expect(rows).toEqual([[0, 21]]);
    SeqLocation.set(seq, A, 3);
    expect(SeqLocation.toRowLocations(seq, msa, RowLocation.create(set), () => {}, null)).toBe(0);
  });

  it("converts loci with dropped gaps, dropped inserts and copies multiplied", () => {
    const c = Conversion.create();
    const rect: AlignedLoci = {
      kind: "aligned-loci",
      set,
      alignment: msa,
      elements: {
        rows: OrderedSet.ofRange(0 as RowIndex, 2 as RowIndex),
        columns: OrderedSet.ofRange(18 as ColumnIndex, 26 as ColumnIndex),
      },
    };
    const residues = AlignedLoci.toRowLoci(rect, c);
    // Row 0 has gaps at 21-24, row 1 at 20-24.
    expect([c.status, c.dropped, residues.elements.map((e) => list(e.indices))]).toEqual([
      ConversionStatus.partial,
      9,
      [
        [18, 19, 20, 21],
        [16, 17, 18],
      ],
    ]);
    const units = RowLoci.toSeqLoci(residues, c);
    expect([c.status, c.dropped, c.multiplied]).toEqual([ConversionStatus.partial, 3, 4]);
    expect(units.elements.map((e) => [e.unit, list(e.indices)])).toEqual([
      [A, [43, 44, 45, 46]],
      [C, [43, 44, 45, 46]],
    ]);
    const rows = SeqLoci.toRowLoci(units, msa, c);
    expect([c.status, c.multiplied, rows.elements.map((e) => [e.row, list(e.indices)])]).toEqual([
      ConversionStatus.exact,
      0,
      [[0, [18, 19, 20, 21]]],
    ]);
    const columns = AlignedLoci.toColumnLoci(rect, c);
    expect(list(columns.columns)).toEqual([18, 19, 20, 21, 22, 23, 24, 25]);
    const band = ColumnLoci.toAlignedLoci(columns, OrderedSet.ofRange(3 as RowIndex, 5 as RowIndex), c);
    expect(Array.isArray(band.elements)).toBe(false);

    // SH2 as A3M: PTN11 (row 1) has the insert "kd" after column 9; its residues 10, 11 have no column.
    const sh2 = 1 as AlignmentId;
    const ptn11: RowLoci = {
      kind: "row-loci",
      set,
      alignment: sh2,
      elements: [{ row: 1 as RowIndex, indices: OrderedSet.ofRange(8 as SeqIndex, 14 as SeqIndex) }],
    };
    const cells = RowLoci.toAlignedLoci(ptn11, c);
    expect([c.status, c.dropped]).toEqual([ConversionStatus.partial, 2]);
    expect(Array.isArray(cells.elements) && list(cells.elements[0].columns)).toEqual([8, 9, 10, 11]);
    const unbound = RowLoci.toSeqLoci(ptn11, c);
    expect([c.status, unbound.elements.length]).toEqual([ConversionStatus.unmapped, 0]);
  });

  it("finds a row ↔ chain ↔ another row path, and nothing past 4 steps", () => {
    const HBA = a3n1().molecules[0].residues as string;
    const two = globins({
      alignments: [{ key: "pair", flavor: "aligned", gapped: [`--${HBA}`, `MM${HBA.replace(/K/g, "R")}`] }],
      mappings: [{ key: "pairA", from: "pair/0", to: "1a3n/A", method: "identity" }],
    }).set;
    const { n, steps, path } = along(two, "PF00042/0", "pair/0");
    expect(n).toBe(3);
    expect(steps.map((s) => [s.mapping >= 0, s.reversed, s.unit])).toEqual([
      [true, false, -1],
      [false, true, unit(two, "1a3n/A")],
      [true, true, -1],
    ]);
    expect([
      path ? Paths.at(path, 0) : undefined,
      path ? Paths.at(path, 21) : undefined,
      back(path?.map, 25),
      back(path?.map, 0),
    ]).toEqual([25, 46, 0, -1]);
    expect(along(two, "PF00042/1", "1a3n/A").n).toBe(4);
    expect(along(two, "PF00042/0", "PF00042/0").n).toBe(0);
    // Unit steps are free: columns → row 0 → entity → A → pair/0 → columns, and the other row of "pair" through
    // them, are 4 counted steps and 5 in all (the old bound refused both).
    expect(along(two, "PF00042", "pair").n).toBe(5);
    expect(along(two, "PF00042/0", "pair/1").n).toBe(5);
    // Nothing past 4 counted steps: PF00042/1 → columns → row 0 → entity → (A) → pair/0 → columns is 5.
    expect(along(two, "PF00042/1", "pair").n).toBe(-1);
    expect(along(two, "PF00042/1", "pair/1").n).toBe(-1);
    // Composed arrays are cached per set and path.
    expect(Mappings.compose(two, space(two, "PF00042/0"), steps)).toBe(path);
  });
});

describe("paths", () => {
  const seq = (key: string) => ({ key, residues: "ACDEFGHIKLMNPQRSTVWY", ownUnit: true });
  const seg = (key: string, from: string, to: string, shift: number): MappingInput => ({
    key,
    from,
    to,
    method: "segments",
    segments: [[0, shift, null]],
  });

  it("goes from UniProt through the entity to a unit", () => {
    const P11540 = { key: "P11540", residues: barstar, ownUnit: true };
    const { set } = SeqSets.build({
      ...brs1(),
      molecules: [...brs1().molecules, P11540],
      mappings: [{ key: "u", from: "P11540", to: "1brs/entity:2", method: "pairwise", cigar: "1I89M" }],
    });
    const { n, path } = along(set, "P11540", "1brs/E");
    expect(n).toBe(2);
    const out = SeqLocation.create(set);
    const l = SeqLocation.set(SeqLocation.create(set), unit(set, "P11540/unit"), 40);
    expect([path && SeqLocation.toSystem(l, path, out), out.unit, out.index]).toEqual([true, unit(set, "1brs/E"), 39]);
    SeqLocation.set(l, unit(set, "P11540/unit"), 0);
    expect(path && SeqLocation.toSystem(l, path, out)).toBe(false);
    const c = Conversion.create();
    const loci: SeqLoci = {
      kind: "seq-loci",
      set,
      elements: [{ unit: unit(set, "P11540/unit"), indices: OrderedSet.ofRange(0 as SeqIndex, 3 as SeqIndex) }],
    };
    const e = path ? SeqLoci.toSystem(loci, path, c) : null;
    expect([c.status, c.dropped, e?.elements.map((x) => list(x.indices))]).toEqual([1, 1, [[0, 1]]]);
  });

  it("reports several shortest paths with different results as ambiguous", () => {
    const molecules = ["P", "Q1", "Q2", "T"].map(seq);
    const equal = SeqSets.build({
      molecules,
      mappings: [seg("a", "P", "Q1", 1), seg("b", "P", "Q2", 1), seg("c", "Q1", "T", 0), seg("d", "Q2", "T", 0)],
    }).set;
    expect(along(equal, "P", "T").n).toBe(2);
    const different = SeqSets.withChanges(equal, { upsert: { mappings: [seg("d", "Q2", "T", 2)] } }).set;
    expect(along(different, "P", "T").n).toBe(-2);
  });

  it("leaves a row not held unmapped", () => {
    const rows = pf00042();
    const { set, issues } = SeqSets.build({
      ...a3n1(),
      alignments: [
        {
          key: "lazy",
          flavor: "aligned",
          gapped: [],
          lazy: {
            rowCount: 10,
            columnCount: rows[0].gapped.length,
            rows: new Map([[3, { label: rows[0].label, gapped: rows[0].gapped }]]),
            labelOf: () => null,
          },
        },
      ],
      mappings: [
        { key: "held", from: "lazy/3", to: "1a3n/A", method: "segments", segments: [[0, 25, null]] },
        { key: "not", from: "lazy/5", to: "1a3n/A", method: "segments", segments: [[0, 25, null]] },
      ],
    });
    expect(codes(issues)).toEqual(["model.mapping.endpoint /mappings/1"]);
    const steps: PathStep[] = [];
    const notHeld = { kind: CoordinateSystemKind.row, id: 0, row: 5 };
    expect(Mappings.path(set, notHeld, space(set, "lazy"), steps)).toBe(-1);
    expect(along(set, "lazy/3", "1a3n/A").n).toBe(1);
    const col = ColumnLocation.set(ColumnLocation.create(set), 0, 20);
    expect(seqsOf(set, (out, cb) => ColumnLocation.toSeqLocations(col, out, cb, null))).toEqual([["1a3n/A", 45]]);
    const c = Conversion.create();
    const loci: RowLoci = {
      kind: "row-loci",
      set,
      alignment: 0 as AlignmentId,
      elements: [{ row: 5 as RowIndex, indices: OrderedSet.ofRange(0 as SeqIndex, 4 as SeqIndex) }],
    };
    expect([RowLoci.toSeqLoci(loci, c).elements.length, c.status, c.dropped]).toEqual([0, 3, 4]);
    expect([RowLoci.toAlignedLoci(loci, c).elements, c.status]).toEqual([[], ConversionStatus.unmapped]);
    const r = RowLocation.set(RowLocation.create(set), 0, 5, 0);
    expect(RowLocation.toSeqLocations(r, SeqLocation.create(set), () => {}, null)).toBe(0);
  });

  it("marks a row with several bindings ambiguous", () => {
    const { set } = globins({
      mappings: [{ key: "again", from: "PF00042/0", to: "1a3n/A", method: "segments", segments: [[0, 25, null]] }],
    });
    const c = Conversion.create();
    const loci: RowLoci = {
      kind: "row-loci",
      set,
      alignment: 0 as AlignmentId,
      elements: [{ row: 0 as RowIndex, indices: OrderedSet.ofSingleton(0 as SeqIndex) }],
    };
    const units = RowLoci.toSeqLoci(loci, c);
    expect([c.status, c.multiplied, units.elements.length]).toEqual([ConversionStatus.ambiguous, 1, 2]);
    expect(Mappings.bindingsOf(set, 0)?.bindings[0]).toBe(2);
  });
});

describe("the path search (M1)", () => {
  const sys = (set: SeqSet, key: string) => space(set, key);
  const n = (set: SeqSet, from: string, to: string) => along(set, from, to).n;

  it("keeps copies apart in the search, and finds the longer legal path through their rows", () => {
    const BARNASE = brs1().molecules[0].residues as string;
    const rows = {
      alignments: [{ key: "aln", flavor: "aligned" as const, gapped: [BARNASE, BARNASE] }],
      mappings: [
        { key: "r1", from: "aln/0", to: "1brs/A", method: "identity" as const },
        { key: "r2", from: "aln/1", to: "1brs/B", method: "identity" as const },
      ],
    };
    const { set } = SeqSets.build({ ...brs1(), ...rows });
    // A → R1 → columns → R2 → B, never A → entity → B.
    const { steps, path } = along(set, "1brs/A", "1brs/B");
    expect(steps.map((x) => [x.mapping, x.row, x.unit])).toEqual([
      [0, -1, -1],
      [-1, 0, -1],
      [-1, 1, -1],
      [1, -1, -1],
    ]);
    expect([path ? Paths.at(path, 50) : undefined, path?.codon]).toEqual([50, -1]);
    expect(n(SeqSets.build(brs1()).set, "1brs/A", "1brs/B")).toBe(-1);
    // Data on the molecule still reaches every unit; a unit's own data none of its copies.
    expect([n(set, "1brs/entity:1", "1brs/C"), n(set, "aln/0", "1brs/C")]).toEqual([1, -1]);

    // Coverage against a reference mapped to chain A only does not leak to B.
    const P00648 = { key: "P00648", residues: spike("1brs").barnase.sequence as string };
    const sifts = {
      key: "sifts",
      from: "P00648",
      to: "1brs/A",
      method: "segments" as const,
      segments: [[47, 0, null]],
    };
    const ref = SeqSets.build({
      ...brs1(),
      molecules: [...brs1().molecules, P00648],
      mappings: [sifts as MappingInput],
    });
    const cov = Providers.coverage(ref.set, { key: "cov", of: ["1brs/A", "1brs/B"], reference: "P00648" });
    expect(codes(cov.issues)).toEqual(["model.provider.unmapped /reference"]);
    expect(cov.issues[0].message).toContain('"1brs/B"');
  });

  it("compares whole maps: a split codon placed differently is ambiguous though the codon starts agree", () => {
    // Two transcripts of one gene whose splice acceptor sits in residue 1's codon at different places.
    const { set } = SeqSets.build({
      molecules: [
        { key: "P", alphabet: "protein", residues: "MK" },
        { key: "T1", alphabet: "dna", residues: "ATGAAA" },
        { key: "T2", alphabet: "dna", residues: "ATGAAA" },
        { key: "G", alphabet: "dna", residues: "ATGAAAAAAAAA" },
      ],
      mappings: [
        { key: "c1", from: "P", to: "T1", method: "codon", cds: [[0, 6]] },
        { key: "c2", from: "P", to: "T2", method: "codon", cds: [[0, 6]] },
        {
          key: "e1",
          from: "T1",
          to: "G",
          method: "segments",
          segments: [
            [0, 0, 4],
            [4, 10, 2],
          ],
        },
        {
          key: "e2",
          from: "T2",
          to: "G",
          method: "segments",
          segments: [
            [0, 0, 5],
            [5, 11, 1],
          ],
        },
      ],
    });
    const a = Mappings.compose(set, sys(set, "P"), [step(0, false), step(2, false)]);
    const b = Mappings.compose(set, sys(set, "P"), [step(1, false), step(3, false)]);
    expect([dense(a), dense(b)]).toEqual([
      [0, 3],
      [0, 3],
    ]);
    expect([denseSpans(a), denseSpans(b)]).toEqual([
      [0, 1, 2, 3, 10, 11],
      [0, 1, 2, 3, 4, 11],
    ]);
    expect([Mappings.sameComposed(a, b), a.codon, n(set, "P", "G")]).toEqual([false, 0, -2]);
    // From the nucleotides too: G 4 is residue 1 through T2 only.
    expect(n(set, "G", "P")).toBe(-2);
  });

  it("drops a path that pairs nothing before comparing (a chimera's parts)", () => {
    const seq = (key: string, residues: string) => ({ key, residues });
    const GFP = "MSKGEELFTG";
    const P53 = "MEEPQSDPSV";
    const { set } = SeqSets.build({
      molecules: [seq("gfp", GFP), seq("p53", P53), seq("fusion", `${GFP}GGG${P53}`), seq("other", P53)],
      mappings: [
        { key: "a", from: "fusion", to: "gfp", method: "segments", segments: [[0, 0, 10]] },
        { key: "b", from: "fusion", to: "p53", method: "segments", segments: [[13, 0, 10]] },
      ],
    });
    // GFP → fusion → p53 composes to nothing: no path.
    expect(n(set, "gfp", "p53")).toBe(-1);
    // With a real path of the same length beside it, that one is not ambiguous.
    const more = SeqSets.withChanges(set, {
      upsert: {
        mappings: [
          { key: "c", from: "gfp", to: "other", method: "identity", fromRange: [0, 3] },
          { key: "d", from: "other", to: "p53", method: "identity" },
        ],
      },
    }).set;
    const { n: steps, path } = along(more, "gfp", "p53");
    expect([steps, path ? Paths.at(path, 2) : undefined, path ? Paths.at(path, 3) : undefined]).toEqual([2, 2, -1]);
  });

  it("searches a 20,000-row alignment through its two bound rows only", () => {
    const g = generated();
    const molecules = g.bound.map((b, i) => ({ key: `s${i}`, residues: b.sequence, ownUnit: true }));
    const { set, issues } = SeqSets.build({
      molecules: [...molecules, { key: "lone", residues: "ACDEF", ownUnit: true }],
      alignments: [{ key: "msa", flavor: "aligned", gapped: [...g.rows] }],
      mappings: g.bound.map((b, i) => ({
        key: `m${i}`,
        from: `msa/${b.row}`,
        to: `s${i}`,
        method: "identity" as const,
      })),
    });
    expect(issues).toEqual([]);
    expect(n(set, "s0/unit", "s1/unit")).toBe(6);
    expect(Mappings.expanded()).toBeLessThan(20);
    expect(n(set, "s0/unit", "lone/unit")).toBe(-1);
    expect(Mappings.expanded()).toBeLessThan(20);
    // An unbound row as either end: its own edge to the columns.
    expect([n(set, "msa/7", "s1"), n(set, "s0", "msa/7")]).toEqual([3, 3]);
    expect(Mappings.expanded()).toBeLessThan(20);
    // Cached per set: the same search again expands nothing.
    expect(n(set, "s0/unit", "s1/unit")).toBe(6);
  });

  it("reads the ambiguity case and free unit steps together", () => {
    const seq = (key: string) => ({ key, residues: "ACDEFGHIKLMNPQRSTVWY", ownUnit: true });
    const seg = (key: string, from: string, to: string, shift: number): MappingInput => ({
      key,
      from,
      to,
      method: "segments",
      segments: [[0, shift, null]],
    });
    const molecules = ["P", "Q1", "Q2", "T"].map(seq);
    const equal = SeqSets.build({
      molecules,
      mappings: [seg("a", "P", "Q1", 1), seg("b", "P", "Q2", 1), seg("c", "Q1", "T", 0), seg("d", "Q2", "T", 0)],
    }).set;
    // Two counted steps and the two unit steps, either way: one result.
    expect(n(equal, "P/unit", "T/unit")).toBe(4);
    const different = SeqSets.withChanges(equal, { upsert: { mappings: [seg("d", "Q2", "T", 2)] } }).set;
    expect([n(different, "P/unit", "T/unit"), n(different, "P", "T/unit")]).toEqual([-2, -2]);
  });
});

function step(mapping: number, reversed: boolean): PathStep {
  return { mapping, reversed, alignment: -1, row: -1, unit: -1 };
}
