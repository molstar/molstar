/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Dependents after a change (seqspec:model#identity-and-lineage; M2): a kept entry is checked against this version, a
// consequence is a warning that points at its cause, and the unit table is re-frozen only when units change. Also
// parent rows.

import { describe, expect, it } from "vitest";
import {
  type Annotation,
  type AnnotationInput,
  type Issue,
  type Mapping,
  type MappingInput,
  Mappings,
  type SeqSet,
  SeqSets,
  type SetInput,
} from "../../src/model/index.ts";
import { contacts, iep1 } from "./annotation-fixtures.ts";

const codes = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);
const byKey = <T extends { key: number } | null>(set: SeqSet, list: readonly T[], key: string) =>
  list.find((x) => x && set.strings.get(x.key) === key) ?? null;

const SEQ = "ACDEFGHIKL";
const base = (): SetInput => ({
  molecules: [
    { key: "a", residues: SEQ },
    { key: "b", residues: SEQ },
  ],
  units: [
    { key: "b1", molecule: "b" },
    { key: "b2", molecule: "b", numberings: [{ name: "auth", method: "offset", first: 10 }] },
  ],
  mappings: [{ key: "m", from: "a", to: "b", method: "pairwise", cigar: "10M" }],
  annotations: [
    { key: "f", geometry: "ranges", target: "b1", columns: { start: [2, 8], stop: [4, 10] } },
    { key: "o", geometry: "ranges", target: "b2", columns: { start: [2] } },
  ],
});
const typo: MappingInput = { key: "m", from: "a", to: "b", method: "pairwise", cigar: "10Q" };

describe("dependents of a change", () => {
  it("revalidate a mapping kept after a failed upsert: the endpoint lengthened in the same change", () => {
    const v1 = SeqSets.build(base()).set;
    const r = SeqSets.withChanges(v1, {
      upsert: { molecules: [{ key: "a", residues: `${SEQ}MN` }], mappings: [typo] },
    });
    const m = byKey(r.set, r.set.mappings, "m") as Mapping | null;
    // The previous version is never kept unchecked: rebuilt against the new `a`, or removed with an issue.
    expect(m === null || m.forward.length === 12).toBe(true);
    expect(codes(r.issues)).toContain("model.mapping.cigar /upsert/mappings/0/cigar");
    expect(r.issues.filter((i) => !i.path)).toEqual([]);
    if (!m) expect(codes(r.issues).some((c) => c.endsWith(" /upsert/molecules/0"))).toBe(true);
  });

  it("orphan a mapping kept after a failed upsert when the same change removes its endpoint", () => {
    const v1 = SeqSets.build(base()).set;
    const r = SeqSets.withChanges(v1, { upsert: { mappings: [typo] }, remove: ["b"] });
    expect(codes(r.issues)).toContain("model.mapping.orphaned /remove/0");
    expect([byKey(r.set, r.set.mappings, "m"), r.set.byName("m", { kind: -1, id: -1, row: -1 })]).toEqual([
      null,
      false,
    ]);
  });

  it("orphan an annotation kept after a failed upsert when its target goes in the same change", () => {
    const v1 = SeqSets.build(base()).set;
    const bad: AnnotationInput = { key: "f", geometry: "nope" as never, target: "b1", columns: {} };
    const r = SeqSets.withChanges(v1, { upsert: { annotations: [bad] }, remove: ["x", "b1"] });
    expect(codes(r.issues)).toEqual([
      "model.key.unknown /remove/0",
      "model.annotation.geometry /upsert/annotations/0/geometry",
      "model.annotation.orphaned /remove/1",
    ]);
    expect(byKey(r.set, r.set.annotations, "f")).toBeNull();
  });

  it("point rebuild and orphan issues at the change that caused them, naming the dependent", () => {
    const v1 = SeqSets.build(base()).set;
    // b shortened: b2 is stale (a warning), b1 follows; the mapping no longer fits, f's second row is out of range, o
    // loses its target. Every issue names its dependent and points at the change.
    const r = SeqSets.withChanges(v1, { upsert: { molecules: [{ key: "b", residues: "ACDEFGHI" }] } });
    expect(r.issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ["model.unit.stale", "warning", "/upsert/molecules/0"],
      ["model.mapping.cigar", "error", "/upsert/molecules/0"],
      ["model.position.range", "error", "/upsert/molecules/0"],
      ["model.annotation.orphaned", "warning", "/upsert/molecules/0"],
    ]);
    expect(r.issues.map((i) => i.message.match(/^(Unit|Mapping|Annotation) "(\w+)"/)?.[2])).toEqual([
      "b2",
      "m",
      "f",
      "o",
    ]);
    expect((byKey(r.set, r.set.annotations, "f") as Annotation).count).toBe(1);
    // A removal: the orphan points at it.
    const r2 = SeqSets.withChanges(v1, { remove: ["a"] });
    expect(codes(r2.issues)).toEqual(["model.mapping.orphaned /remove/0"]);
  });

  it("re-freeze the unit table only when units change: removing an annotation shares it", () => {
    const input = iep1();
    const v1 = SeqSets.build({
      ...input,
      alignments: [{ key: "aln", flavor: "aligned", gapped: ["MKV", "MKV"] }],
      annotations: [contacts(), { key: "h", geometry: "ranges", target: "1iep/A", columns: { start: [3] } }],
    }).set;
    const pairs = byKey(v1, v1.annotations, "contacts") as Annotation;
    const r = SeqSets.withChanges(v1, { remove: ["h"] });
    expect(r.issues).toEqual([]);
    const v2 = r.set;
    expect([v2.units === v1.units, v2.groups === v1.groups, byKey(v2, v2.annotations, "contacts") === pairs]).toEqual([
      true,
      true,
      true,
    ]);
    expect(Mappings.bindingsOf(v2, 0) === Mappings.bindingsOf(v1, 0)).toBe(true);
    // Removing a unit still does.
    expect(SeqSets.withChanges(v1, { remove: ["1iep/C"] }).set.units).not.toBe(v1.units);
  });
});

describe("parent rows", () => {
  it("warn on a bad parent and on the parent that closes a cycle; the alignment is built, that row has none", () => {
    const { set, issues } = SeqSets.build({
      molecules: [],
      alignments: [
        { key: "two", flavor: "aligned", gapped: ["MK", "MK"], parentRows: [1, 0] },
        { key: "four", flavor: "aligned", gapped: ["M", "M", "M", "M"], parentRows: [1, 2, 0, 9] },
      ],
    });
    expect(issues.map((i) => [i.code, i.severity, i.path])).toEqual([
      ["model.alignment.parent", "warning", "/alignments/0/parentRows/1"],
      ["model.alignment.parent", "warning", "/alignments/1/parentRows/3"],
      ["model.alignment.parent", "warning", "/alignments/1/parentRows/2"],
    ]);
    expect(set.alignments.map((a) => a && Array.from(a.parent))).toEqual([
      [1, -1],
      [1, 2, -1, -1],
    ]);
  });
});
