/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `Extents.addLoci` (NG1; IN1 shares it): a host's loci as entries in each element's own coordinate system.

import { describe, expect, it } from "vitest";
import {
  type AlignmentId,
  type AnnotationId,
  Extents,
  type ItemLoci,
  CoordinateSystemKind as K,
  type Loci,
  OrderedSet,
  type RowIndex,
  type SeqSet,
  SeqSets,
  type UnitId,
} from "../../src/model/index.js";
import { contacts, iep1 } from "./annotation-fixtures.js";
import { brs1 } from "./fixtures.js";

const range = <T extends number = never>(a: number, b: number) => OrderedSet.ofRange(a as T, b as T);
const sorted = (...xs: number[]) => OrderedSet.ofSortedArray<never>(Int32Array.from(xs));
const listed = (e: Extents) =>
  Array.from({ length: e.count }, (_, i) => [e.kind[i], e.id[i], e.row[i], e.start[i], e.stop[i]]);
const from = (loci: Loci) => {
  const e = Extents.create(1);
  return [Extents.addLoci(e, loci), listed(e)] as const;
};

const stubs = (base: SeqSet): SeqSet =>
  SeqSets.withChanges(base, {
    upsert: {
      alignments: [{ key: "msa", gapped: ["AC-DE", "ACGD-", "AC-DE"], flavor: "aligned" }],
      annotations: [
        { key: "f", geometry: "ranges", target: "1brs/A", columns: { start: [3, 10], stop: [5, 12], id: ["a", "b"] } },
      ],
    },
  }).set;

describe("Extents.addLoci", () => {
  const set = stubs(SeqSets.build(brs1()).set);

  it("a seq-loci interval is one entry, a sorted array one per run of consecutive indices; boundaries are (b, b)", () => {
    expect(from({ kind: "seq-loci", set, elements: [{ unit: 2 as UnitId, indices: range(4, 9) }] })).toEqual([
      1,
      [[K.unit, 2, -1, 4, 9]],
    ]);
    const mixed = from({
      kind: "seq-loci",
      set,
      elements: [{ unit: 0 as UnitId, indices: sorted(1, 2, 3, 7, 9, 10) }],
    });
    expect(mixed).toEqual([
      3,
      [
        [K.unit, 0, -1, 1, 4],
        [K.unit, 0, -1, 7, 8],
        [K.unit, 0, -1, 9, 11],
      ],
    ]);
    expect(from({ kind: "boundary-loci", set, elements: [{ unit: 1 as UnitId, indices: sorted(0, 5) }] })[1]).toEqual([
      [K.unit, 1, -1, 0, 0],
      [K.unit, 1, -1, 5, 5],
    ]);
  });

  it("a row's residues are in the row's system, columns in the alignment's", () => {
    const aln = 0 as AlignmentId;
    expect(
      from({ kind: "row-loci", set, alignment: aln, elements: [{ row: 1 as RowIndex, indices: range(2, 4) }] })[1],
    ).toEqual([[K.row, 0, 1, 2, 4]]);
    expect(from({ kind: "column-loci", set, alignment: aln, columns: sorted(0, 1, 4) })[1]).toEqual([
      [K.columns, 0, -1, 0, 2],
      [K.columns, 0, -1, 4, 5],
    ]);
  });

  it("an aligned cell is its residue in the row's system, a gap its column; both forms, the rectangle row by row", () => {
    const aln = 0 as AlignmentId;
    // Row 0 is AC-DE: columns 0, 1 are residues 0, 1; column 2 is a gap; columns 3, 4 are residues 2, 3.
    const one = from({
      kind: "aligned-loci",
      set,
      alignment: aln,
      elements: [{ row: 0 as RowIndex, columns: range(0, 5) }],
    });
    expect(one[1]).toEqual([
      [K.row, 0, 0, 0, 2],
      [K.columns, 0, -1, 2, 3],
      [K.row, 0, 0, 2, 4],
    ]);
    const rect = from({
      kind: "aligned-loci",
      set,
      alignment: aln,
      elements: { rows: sorted(0, 1), columns: range(3, 5) },
    });
    expect(rect[1]).toEqual([
      [K.row, 0, 0, 2, 4],
      [K.row, 0, 1, 3, 4],
      [K.columns, 0, -1, 4, 5],
    ]);
  });

  it("an item is its rows' ranges in their targets' systems", () => {
    const loci: ItemLoci = {
      kind: "item-loci",
      set,
      elements: [{ annotation: 0 as AnnotationId, items: range(0, 2) }],
    };
    const [n, entries] = from(loci);
    expect(n).toBe(2);
    expect(entries.map((e) => [e[0], e[3], e[4]])).toEqual([
      [K.unit, 3, 5],
      [K.unit, 10, 12],
    ]);
  });

  it("an item of a pairs table is both its sides", () => {
    const iep = SeqSets.build({ ...iep1(), annotations: [contacts()] }).set;
    const [n, entries] = from({
      kind: "item-loci",
      set: iep,
      elements: [{ annotation: 0 as AnnotationId, items: range(0, 1) }],
    });
    expect(n).toBe(2);
    expect(entries.map((e) => e[0])).toEqual([K.unit, K.unit]);
    expect(entries[0][1]).not.toBe(entries[1][1]);
  });

  it("appends: the entries already in `out` stay, and the count is what was added", () => {
    const e = Extents.create(1);
    const loci: Loci = { kind: "seq-loci", set, elements: [{ unit: 0 as UnitId, indices: range(0, 2) }] };
    expect(Extents.addLoci(e, loci)).toBe(1);
    expect(Extents.addLoci(e, loci)).toBe(1);
    expect(e.count).toBe(2);
    expect(Extents.addLoci(e, { kind: "seq-loci", set, elements: [] })).toBe(0);
  });
});
