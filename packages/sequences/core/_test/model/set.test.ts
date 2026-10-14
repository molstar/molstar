/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
// Molecule tests are molecule.test.ts (MD1); the asserts that go through the set are here: row refs, alignment ids
// across versions, a unit's numbering winning.
import {
  type AlignmentInput,
  KeyKind,
  type KeyTarget,
  type NameId,
  Numbering,
  type SeqIndex,
  type SeqSet,
  SeqSets,
  type SetInput,
  type UnitId,
} from "../../src/model/index.js";
import { a3n1, brs1 } from "./fixtures.js";

const keyOf = (set: SeqSet, id: number) => set.strings.get(set.units.key[id]);
const target = (): KeyTarget => ({ kind: -1, id: -1, row: -1 });
function lookup(set: SeqSet, key: string) {
  const out = target();
  return set.byKey(set.strings.find(key), out) ? out : null;
}
/** Copy groups as lists of unit keys. */
const groups = (set: SeqSet) => set.groups.members.map((m) => Array.from(m, (u) => keyOf(set, u)));
const codes = (issues: { code: string; path: string }[]) => issues.map((i) => `${i.code} ${i.path}`);

describe("units, complexes and copy groups", () => {
  it("groups 1BRS's barnase and barstar chains", () => {
    const { set, issues } = SeqSets.build(brs1());
    expect(issues).toEqual([]);
    expect(groups(set)).toEqual([
      ["1brs/A", "1brs/B", "1brs/C"],
      ["1brs/D", "1brs/E", "1brs/F"],
    ]);
    expect(Array.from(set.groups.first)).toEqual([0, 3]);
    expect(set.units.observed[0]).toEqual(Int32Array.of(2, 110));
    expect(set.units.observed[3]).toEqual(Int32Array.of(0, 63, 65, 89));
    const cx = set.complexes[0];
    expect(cx && Array.from(cx.units)).toEqual([0, 1, 2, 3, 4, 5]);
    expect(Array.from(set.units.complex)).toEqual([0, 0, 0, 0, 0, 0]);
    expect(lookup(set, "1brs")).toEqual({ kind: KeyKind.complex, id: 0, row: -1 });
    expect(lookup(set, "1brs/E")).toEqual({ kind: KeyKind.unit, id: 4, row: -1 });
    expect(lookup(set, "1brs/entity:2")).toEqual({ kind: KeyKind.molecule, id: 1, row: -1 });
    expect(lookup(set, "1brs/Z")).toBeNull();
  });

  it("groups 1A3N's α and β chains and its four hemes", () => {
    const { set, issues } = SeqSets.build(a3n1());
    expect(issues).toEqual([]);
    expect(groups(set)).toEqual([
      ["1a3n/A", "1a3n/C"],
      ["1a3n/B", "1a3n/D"],
      ["1a3n/E", "1a3n/F", "1a3n/G", "1a3n/H"],
    ]);
  });

  it("splits copies by range, and keeps own units apart", () => {
    const { set, issues } = SeqSets.build({
      molecules: [{ key: "p53", residues: "MEEPQSDPSV", ownUnit: true }],
      units: [
        { key: "full", molecule: "p53" },
        { key: "dbd", molecule: "p53", range: [2, 8] },
        { key: "dbd2", molecule: "p53", range: [2, 8], score: 0.5 },
      ],
    });
    expect(issues).toEqual([]);
    expect(groups(set)).toEqual([["p53/unit"], ["full"], ["dbd", "dbd2"]]);
    expect(set.strings.get(set.units.origin[0])).toBe("sequence");
    expect([set.units.rangeStart[1], set.units.rangeStop[1]]).toEqual([0, 10]);
    expect(set.units.score[3]).toBeCloseTo(0.5);
    expect(Number.isNaN(set.units.score[1])).toBe(true);
  });

  it("reports bad units and complexes with their input paths", () => {
    const { issues } = SeqSets.build({
      molecules: [{ key: "m", residues: "ACDEF", ownUnit: true }],
      units: [
        { key: "a", molecule: "nope" },
        { key: "b", molecule: "m", range: [3, 9] },
        {
          key: "c",
          molecule: "m",
          observed: [
            [2, 4],
            [3, 5],
          ],
        },
        { key: "d", molecule: "m", origin: "sequence" },
        { key: "m", molecule: "m" },
        { key: "e", molecule: "m" },
      ],
      complexes: [
        { key: "x", units: ["e", "zz"] },
        { key: "y", units: ["e"] },
      ],
    });
    expect(codes(issues)).toEqual([
      "model.unit.molecule /units/0",
      "model.unit.range /units/1",
      "model.unit.observed /units/2",
      "model.unit.origin /units/3",
      "model.key.duplicate /units/4",
      "model.complex.unit /complexes/0/units/1",
      "model.complex.unit /complexes/1/units/0",
    ]);
  });
});

describe("withChanges", () => {
  const base = (): SetInput => ({
    molecules: [
      { key: "a", residues: "MKV", ownUnit: true },
      { key: "b", residues: "ACDEFG" },
    ],
    units: [
      { key: "b1", molecule: "b" },
      { key: "b2", molecule: "b", numberings: [{ name: "auth", method: "offset", first: 10 }] },
    ],
    complexes: [{ key: "cx", units: ["b1", "b2"] }],
  });

  it("keeps ids, appends new ones, never gives a removed one to another key, and gives it back to its own", () => {
    const v1 = SeqSets.build(base()).set;
    const r2 = SeqSets.withChanges(v1, {
      upsert: { molecules: [{ key: "c", residues: "WW" }], units: [{ key: "c1", molecule: "c" }] },
      remove: ["b1"],
    });
    expect(r2.issues).toEqual([]);
    const v2 = r2.set;
    expect(v2.lineage).toBe(v1.lineage);
    expect(v2.strings).toBe(v1.strings);
    expect(lookup(v2, "c")?.id).toBe(2);
    expect(lookup(v2, "c1")?.id).toBe(3);
    expect(lookup(v2, "b1")).toBeNull();
    expect(v2.units.alive[1]).toBe(0);
    expect(v2.complexes[0] && Array.from(v2.complexes[0].units)).toEqual([2]);
    // seqspec:model#identity-and-lineage: a key removed and upserted again gets its old id back, in no complex; a new
    // key the next id.
    const v3 = SeqSets.withChanges(v2, {
      upsert: {
        units: [
          { key: "b1", molecule: "b" },
          { key: "d1", molecule: "c" },
        ],
      },
    }).set;
    expect(lookup(v3, "b1")?.id).toBe(1);
    expect(lookup(v3, "d1")?.id).toBe(4);
    expect(v3.units.complex[1]).toBe(-1);
    expect(groups(v3)).toEqual([["a/unit"], ["b1", "b2"], ["c1", "d1"]]);
  });

  it("evicts a key with its dependents, quietly, and gives every id back when they return", () => {
    const input = (): SetInput => ({
      molecules: [
        { key: "doc", residues: "MKVLA", ownUnit: true },
        { key: "s/p", residues: "MKV", ownUnit: true },
      ],
      units: [{ key: "s/u", molecule: "s/p" }],
      mappings: [{ key: "s/m", from: "s/p/unit", to: "doc/unit", method: "identity" }],
      annotations: [
        { key: "s/a", geometry: "ranges", columns: { target: ["s/u"], start: [0], stop: [2] } },
        { key: "own", geometry: "ranges", columns: { target: ["doc/unit"], start: [0], stop: [2] } },
      ],
    });
    const v1 = SeqSets.build(input()).set;
    const ids = (set: SeqSet) => ["s/p", "s/p/unit", "s/u", "s/m", "s/a", "own"].map((k) => lookup(set, k)?.id ?? null);
    const before = ids(v1);
    const r2 = SeqSets.withChanges(v1, { evict: ["s/p", "nothing"] });
    expect(r2.issues).toEqual([]);
    expect(ids(r2.set)).toEqual([null, null, null, null, null, before[5]]);
    expect(r2.set.annotations[before[5] as number]).toBe(v1.annotations[before[5] as number]);
    // The same removal is reported: the dependents are orphaned, and an unknown key named.
    expect(codes(SeqSets.withChanges(v1, { remove: ["s/p", "nothing"] }).issues)).toEqual([
      "model.key.unknown /remove/1",
      "model.unit.orphaned /remove/0",
      "model.mapping.orphaned /remove/0",
      "model.annotation.orphaned /remove/0",
    ]);
    const r3 = SeqSets.withChanges(r2.set, { upsert: input() });
    expect(r3.issues).toEqual([]);
    expect(ids(r3.set)).toEqual(before);
    expect(r3.set.lineage).toBe(v1.lineage);
  });

  it("shares what a change does not touch", () => {
    const v1 = SeqSets.build(base()).set;
    const v2 = SeqSets.withChanges(v1, { upsert: { molecules: [{ key: "b", residues: "ACDEFH" }] } }).set;
    expect(v2.molecules[0]).toBe(v1.molecules[0]);
    expect(v2.molecules[1]).not.toBe(v1.molecules[1]);
    expect(v2.molecules[1]?.id).toBe(1);
    expect(v2.units).toBe(v1.units);
    expect(v2.groups).toBe(v1.groups);
    expect(v2.complexes).toBe(v1.complexes);
    const v3 = SeqSets.withChanges(v2, { upsert: { complexes: [{ key: "cx", units: ["b2"] }] } }).set;
    expect(v3.molecules).toBe(v2.molecules);
    expect(v3.units).not.toBe(v2.units);
    expect(v3.units.numberings[2]).toBe(v2.units.numberings[2]);
    expect(Array.from(v3.units.complex)).toEqual([-1, -1, 0]);
  });

  it("removes units over a removed molecule, and units that no longer fit a changed one", () => {
    const v1 = SeqSets.build(base()).set;
    const r2 = SeqSets.withChanges(v1, {
      upsert: { molecules: [{ key: "b", residues: "ACDEFGHI" }] },
      remove: ["x", "a"],
    });
    expect(codes(r2.issues)).toEqual(["model.key.unknown /remove/0", "model.unit.stale /upsert/molecules/0"]);
    const v2 = r2.set;
    expect(lookup(v2, "a")).toBeNull();
    expect(lookup(v2, "a/unit")).toBeNull();
    expect(lookup(v2, "b2")).toBeNull();
    expect(v2.units.rangeStop[1]).toBe(8);
    expect(v2.complexes[0] && Array.from(v2.complexes[0].units)).toEqual([1]);
    const r3 = SeqSets.withChanges(v2, { remove: ["b"] });
    expect(codes(r3.issues)).toEqual(["model.unit.orphaned /remove/0"]);
    expect(r3.set.groups.count).toBe(0);
  });

  it("keys annotations and answers none for unknown keys", () => {
    const { set, issues } = SeqSets.build({
      molecules: [],
      annotations: [{ key: "f", geometry: "ranges", columns: {} }],
    });
    expect(codes(issues)).toEqual([]);
    expect(lookup(set, "f")).toEqual({ kind: KeyKind.annotation, id: 0, row: -1 });
    expect(lookup(SeqSets.build({ molecules: [] }).set, "msa/0")).toBeNull();
    expect(SeqSets.build({ molecules: [] }).set.byKey(-1 as NameId, target())).toBe(false);
  });
});

describe("through the set (MD1's deferred asserts)", () => {
  it("names a lazy alignment's rows by key and by text, interning nothing (alignment.test)", () => {
    const lazy = {
      rowCount: 1_000_000,
      columnCount: 6,
      rows: new Map([[3, { label: "a", gapped: "-VLSPA" }]]),
      labelOf: () => null,
    };
    const { set, issues } = SeqSets.build({
      molecules: [],
      alignments: [{ key: "big", flavor: "aligned", gapped: [], lazy }],
    });
    expect(issues).toEqual([]);
    const t: KeyTarget = { kind: -1, id: -1, row: -1 };
    expect(set.byKey(set.strings.intern("big/999999"), t) && [t.kind, t.id, t.row]).toEqual([KeyKind.row, 0, 999_999]);
    expect(set.byKey(set.strings.intern("big/1000000"), t)).toBe(false);
    // By text: a row ref no one interned.
    const size = set.strings.size;
    expect(set.byName("big/123", t) && [t.kind, t.id, t.row]).toEqual([KeyKind.row, 0, 123]);
    expect([set.byName("big/12x", t), set.byName("bi/1", t), set.byName("big/", t)]).toEqual([false, false, false]);
    expect(set.strings.size).toBe(size);
  });

  it("keeps alignment ids across withChanges and removes by key (alignment.test)", () => {
    const one = (key: string, gapped: string[]): AlignmentInput => ({ key, flavor: "aligned", gapped });
    const v1 = SeqSets.build({ molecules: [], alignments: [one("a", ["AC"]), one("b", ["GT"])] }).set;
    const v2 = SeqSets.withChanges(v1, { upsert: { alignments: [one("b", ["GTT"]), one("c", ["A"])] } }).set;
    expect(v2.alignments.map((a) => a && [v2.strings.get(a.key), a.id, a.columnCount])).toEqual([
      ["a", 0, 2],
      ["b", 1, 3],
      ["c", 2, 1],
    ]);
    expect(v2.alignments[0]).toBe(v1.alignments[0]);
    const v3 = SeqSets.withChanges(v2, { remove: ["a"] }).set;
    expect(v3.alignments[0]).toBeNull();
    expect(v3.byKey(v3.strings.intern("a/0"), { kind: -1, id: -1, row: -1 })).toBe(false);
    const clash = SeqSets.build({ molecules: [{ key: "a", residues: "AC" }], alignments: [one("a", ["AC"])] });
    expect(codes(clash.issues)).toEqual(["model.key.duplicate /alignments/0"]);
  });

  it("keeps a unit's auth numbering apart from the molecule's label numbering, and lets it win (numbering.test)", () => {
    const { set, issues } = SeqSets.build({
      molecules: [
        {
          key: "abl1",
          residues: "GQYGEVYEGVW",
          numberings: [
            { name: "label", method: "offset", first: 1 },
            { name: "auth", method: "offset", first: 1 },
          ],
        },
      ],
      units: [{ key: "1iep/A", molecule: "abl1", numberings: [{ name: "auth", method: "offset", first: 248 }] }],
    });
    expect(issues).toEqual([]);
    const unit = 0 as UnitId;
    const auth = SeqSets.numberingOf(set, unit, set.strings.find("auth"));
    const label = SeqSets.numberingOf(set, unit, set.strings.find("label"));
    const ordinal = SeqSets.numberingOf(set, unit, set.strings.find("ordinal"));
    expect(auth && Numbering.format(auth, 0 as SeqIndex, set.strings)).toBe("248");
    expect(label && Numbering.format(label, 0 as SeqIndex, set.strings)).toBe("1");
    expect(ordinal && Numbering.indexOf(ordinal, 11, -1 as NameId)).toBe(10);
    expect(auth && Numbering.indexOf(auth, 250, -1 as NameId)).toBe(2);
  });
});

describe("lineage after a failed re-upsert (CM-07)", () => {
  it("a removed key whose upsert fails keeps its old id for the next one", () => {
    const idOf = (set: ReturnType<typeof SeqSets.build>["set"], key: string) =>
      set.molecules.findIndex((m) => m && set.strings.get(m.key) === key);
    const a = SeqSets.build({
      molecules: [
        { key: "x", residues: "ACDE" },
        { key: "p", residues: "MKTA" },
      ],
    }).set;
    const b = SeqSets.withChanges(a, { remove: ["p"] }).set;
    const failed = SeqSets.withChanges(b, { upsert: { molecules: [{ key: "p", residues: "MKé" }] } });
    expect(failed.issues.map((i) => i.code)).toEqual(["model.residues.invalid"]);
    const after = SeqSets.withChanges(failed.set, { upsert: { molecules: [{ key: "p", residues: "MKTA" }] } }).set;
    expect(idOf(after, "p")).toBe(idOf(a, "p"));
    // Given back once, the id is the key's again: a removal and return after it get it back again.
    const again = SeqSets.withChanges(SeqSets.withChanges(after, { remove: ["p"] }).set, {
      upsert: { molecules: [{ key: "p", residues: "MKTA" }] },
    }).set;
    expect(idOf(again, "p")).toBe(idOf(a, "p"));
  });
});
