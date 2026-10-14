/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import {
  AlignedLocation,
  AlignedLoci,
  type AlignmentId,
  type AnnotationId,
  BoundaryLocation,
  BoundaryLoci,
  ColumnLocation,
  ColumnLoci,
  Interval,
  ItemLocation,
  ItemLoci,
  Loci,
  OrderedSet,
  type RowIndex,
  RowLocation,
  RowLoci,
  SeqLocation,
  SeqLoci,
  type SeqSet,
  SeqSets,
  type UnitId,
} from "../../src/model/index.js";
import { brs1 } from "./fixtures.js";

/** Any brand: `never` elements assign to every OrderedSet<T>. */
const range = <T extends number = never>(a: number, b: number) => OrderedSet.ofRange(a as T, b as T);
const list = (s: OrderedSet) => Array.from({ length: OrderedSet.size(s) }, (_, i) => OrderedSet.getAt(s, i));

/** A set with an alignment "msa" and an annotation "features" whose rows are items with these ids ("": none). */
function withStubs(base: SeqSet, columns: number, rows: number, items: string[]): SeqSet {
  const gapped = Array.from({ length: rows }, () => "A".repeat(columns));
  const columnsOf = { start: items.map((_, i) => i), id: items.map((k) => k || null) };
  const annotations = [{ key: "features", geometry: "ranges" as const, target: "1brs/A", columns: columnsOf }];
  return SeqSets.withChanges(base, { upsert: { alignments: [{ key: "msa", gapped, flavor: "aligned" }], annotations } })
    .set;
}

describe("locations", () => {
  it("are set in place and compared by lineage and fields", () => {
    const { set } = SeqSets.build(brs1());
    const a = SeqLocation.create(set);
    expect([a.kind, a.unit, a.index]).toEqual(["seq", -1, -1]);
    SeqLocation.set(a, 2, 15);
    const b = SeqLocation.set(SeqLocation.create(set), 2, 15);
    expect(SeqLocation.areEqual(a, b)).toBe(true);
    expect(SeqLocation.areEqual(a, SeqLocation.set(SeqLocation.create(set), 2, 16))).toBe(false);
    expect(SeqLocation.areEqual(a, SeqLocation.set(SeqLocation.create(SeqSets.build(brs1()).set), 2, 15))).toBe(false);
    const cell = AlignedLocation.set(AlignedLocation.create(set), 0, 3, 7);
    expect([cell.alignment, cell.row, cell.column]).toEqual([0, 3, 7]);
  });
});

describe("loci", () => {
  const { set } = SeqSets.build(brs1());

  it("hold positions per unit: has, forEachLocation, union, isEmpty, areEqual", () => {
    const a = SeqLoci.fromLocation(SeqLocation.set(SeqLocation.create(set), 3, 10));
    const b: SeqLoci = { kind: "seq-loci", set, elements: [{ unit: 0 as UnitId, indices: range(5, 8) }] };
    const u = SeqLoci.union(a, b);
    expect(u.elements.map((e) => [e.unit, list(e.indices)])).toEqual([
      [0, [5, 6, 7]],
      [3, [10]],
    ]);
    const l = SeqLocation.create(set);
    expect(SeqLoci.has(u, SeqLocation.set(l, 0, 6))).toBe(true);
    expect(SeqLoci.has(u, SeqLocation.set(l, 0, 8))).toBe(false);
    expect(SeqLoci.has(u, SeqLocation.set(l, 1, 6))).toBe(false);
    expect(Loci.has(u, BoundaryLocation.set(BoundaryLocation.create(set), 0, 6))).toBe(false);
    const seen: string[] = [];
    const n = SeqLoci.forEachLocation(u, l, (x, ctx) => ctx.push(`${x.unit}:${x.index}`), seen);
    expect([n, seen]).toEqual([4, ["0:5", "0:6", "0:7", "3:10"]]);
    expect(Loci.areEqual(u, SeqLoci.union(b, a))).toBe(true);
    expect(Loci.areEqual(u, a)).toBe(false);
    expect(Loci.isEmpty(u)).toBe(false);
    expect(Loci.isEmpty({ kind: "seq-loci", set, elements: [{ unit: 0 as UnitId, indices: OrderedSet.Empty }] })).toBe(
      true,
    );
  });

  it("hold boundaries, row residues, columns and items", () => {
    const at = BoundaryLocation.set(BoundaryLocation.create(set), 1, 110);
    expect(BoundaryLoci.has(BoundaryLoci.fromLocation(at), at)).toBe(true);
    const r = RowLocation.set(RowLocation.create(set), 0, 4, 9);
    const rows = RowLoci.union(
      RowLoci.fromLocation(r),
      RowLoci.fromLocation(RowLocation.set(RowLocation.create(set), 0, 4, 10)),
    );
    expect(rows.elements.map((e) => list(e.indices))).toEqual([[9, 10]]);
    expect(RowLoci.has(rows, RowLocation.set(RowLocation.create(set), 1, 4, 9))).toBe(false);
    const c = ColumnLoci.union(
      ColumnLoci.fromLocation(ColumnLocation.set(ColumnLocation.create(set), 0, 3)),
      ColumnLoci.fromLocation(ColumnLocation.set(ColumnLocation.create(set), 0, 4)),
    );
    expect(Interval.is(c.columns) && list(c.columns)).toEqual([3, 4]);
    const item = ItemLocation.set(ItemLocation.create(set), 0, 2);
    expect(ItemLoci.has(ItemLoci.fromLocation(item), item)).toBe(true);
    expect(() => Loci.union(c, { ...c, alignment: 1 as AlignmentId })).toThrow();
  });

  it("keep a rectangle of aligned cells, and a union that is one", () => {
    const rect: AlignedLoci = {
      kind: "aligned-loci",
      set,
      alignment: 0 as AlignmentId,
      elements: { rows: range<RowIndex>(0, 3), columns: range(5, 10) },
    };
    const l = AlignedLocation.create(set);
    expect(AlignedLoci.has(rect, AlignedLocation.set(l, 0, 2, 9))).toBe(true);
    expect(AlignedLoci.has(rect, AlignedLocation.set(l, 0, 3, 9))).toBe(false);
    expect(AlignedLoci.forEachLocation(rect, l, () => {}, null)).toBe(15);
    const nextRow = AlignedLoci.fromLocation(AlignedLocation.set(l, 0, 3, 5));
    const grown = AlignedLoci.union(rect, { ...nextRow, elements: [{ row: 3 as RowIndex, columns: range(5, 10) }] });
    expect("rows" in grown.elements && list(grown.elements.rows)).toEqual([0, 1, 2, 3]);
    const ragged = AlignedLoci.union(rect, nextRow);
    expect(Array.isArray(ragged.elements) && ragged.elements.length).toBe(4);
    expect(AlignedLoci.has(ragged, AlignedLocation.set(l, 0, 3, 5))).toBe(true);
    expect(AlignedLoci.has(ragged, AlignedLocation.set(l, 0, 3, 6))).toBe(false);
    const wider = AlignedLoci.union(rect, {
      ...rect,
      elements: { rows: range<RowIndex>(0, 3), columns: range(10, 12) },
    });
    expect("rows" in wider.elements && list(wider.elements.columns)).toEqual([5, 6, 7, 8, 9, 10, 11]);
    expect(
      Loci.areEqual(rect, { ...rect, elements: [0, 1, 2].map((r) => ({ row: r as RowIndex, columns: range(5, 10) })) }),
    ).toBe(true);
  });
});

describe("Loci.remap", () => {
  it("goes by key across lineages and drops what is gone", () => {
    const from = SeqSets.build(brs1()).set;
    const input = brs1();
    input.units = input.units?.filter((u) => u.key !== "1brs/A").reverse();
    input.complexes = [];
    const to = SeqSets.build(input).set;
    expect(to.lineage).not.toBe(from.lineage);
    const loci: SeqLoci = {
      kind: "seq-loci",
      set: from,
      elements: [
        { unit: 0 as UnitId, indices: range(0, 5) },
        { unit: 1 as UnitId, indices: range(100, 120) },
      ],
    };
    const out = Loci.remap(loci, from, to);
    expect(out.set).toBe(to);
    expect(out.elements.map((e) => [to.strings.get(to.units.key[e.unit]), list(e.indices)])).toEqual([
      ["1brs/B", [100, 101, 102, 103, 104, 105, 106, 107, 108, 109]],
    ]);
  });

  it("remaps alignments by key with bounds, and items by their ids", () => {
    const from = withStubs(SeqSets.build(brs1()).set, 20, 5, ["x", "y", "z", ""]);
    const to = withStubs(SeqSets.build(brs1()).set, 8, 3, ["z", "x", "w"]);
    const cells: AlignedLoci = {
      kind: "aligned-loci",
      set: from,
      alignment: 0 as AlignmentId,
      elements: { rows: range<RowIndex>(1, 5), columns: range(4, 12) },
    };
    const r = Loci.remap(cells, from, to);
    expect("rows" in r.elements && [list(r.elements.rows), list(r.elements.columns)]).toEqual([
      [1, 2],
      [4, 5, 6, 7],
    ]);
    const items: ItemLoci = {
      kind: "item-loci",
      set: from,
      elements: [{ annotation: 0 as AnnotationId, items: OrderedSet.ofSortedArray(Int32Array.of(0, 2, 3)) }],
    };
    expect(Loci.remap(items, from, to).elements.map((e) => list(e.items))).toEqual([[0, 1]]);
    const cols = ColumnLoci.fromLocation(ColumnLocation.set(ColumnLocation.create(from), 0, 15));
    expect(Loci.isEmpty(Loci.remap(cols, from, to))).toBe(true);
    expect(Loci.remap(cols, from, SeqSets.build(brs1()).set).alignment).toBe(-1);
  });
});
