/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tables (V3; seqspec:format#positions, seqspec:format#annotation-tables): one column-role table per geometry,
// the model's `ColumnRoles`, shared by inline validation, file typing and the model; typed columns and codes; repeated
// failures as one issue.

import { ColumnRoles } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { GeometryColumns, normalize, validate } from "../../src/format/index.ts";
import { renameFields, typeColumns } from "../../src/format/tables/table.ts";
import { dataOf, fullDoc, type J, kernel, sectionOf } from "./fixtures.ts";

const reg = kernel();
const AT = "/root/children/0/children/8/params";
/** The issues of the full document with annotation 8 given these params. */
function issues(params: J): string[] {
  const d = fullDoc();
  dataOf(d)[8].params = { target: "p53", ...params };
  // Only the annotation's own: other nodes read it as a ranges annotation.
  return validate(d, reg)
    .map((i) => `${i.code} ${i.path}`)
    .filter((i) => i.includes(AT));
}

describe("column roles", () => {
  it("the validator reads the model's table", () => {
    for (const [g, cols] of Object.entries(GeometryColumns))
      expect(cols).toEqual(Object.keys(ColumnRoles[g as "rows"]));
  });

  it("a rows table carries label, category and value as fields", () => {
    const value = { row: [0, 1], label: ["a", "b"], category: ["x", "y"], value: [0.5, 1] };
    expect(issues({ geometry: "rows", target: "aln", value })).toEqual([]);
  });

  it("a target or position column of another geometry is svs.columns.unknown, in either form", () => {
    expect(issues({ geometry: "rows", target: "aln", value: { row: [0], start: [1] } })).toEqual([
      `svs.columns.unknown ${AT}/value/start`,
    ]);
    expect(issues({ geometry: "ranges", value: [{ start: 0, row: 1 }] })).toEqual([`svs.columns.unknown ${AT}/value`]);
  });

  it("rows may name a state or a node instead of a row", () => {
    expect(issues({ geometry: "rows", target: "aln", value: { state: ["A"], label: ["a"] } })).toEqual([]);
  });

  it("inline values over a range validate, as the same table from a file does", () => {
    const inline = { start: [0, 5], stop: [4, 9], value: [0.5, 1] };
    expect(issues({ geometry: "values", value: inline })).toEqual([]);
    const typed = typeColumns({ start: ["0", "5"], stop: ["4", "9"], value: ["0.5", "1"] });
    expect(issues({ geometry: "values", value: typed })).toEqual([]);
    expect(Array.from(typed.start as Int32Array)).toEqual(inline.start);
    expect(Array.from(typed.value as Float64Array)).toEqual(inline.value);
  });
});

describe("empty tables", () => {
  it("[] and {} need no columns", () => {
    for (const geometry of ["ranges", "values", "pairs", "rows"]) {
      expect(issues({ geometry, value: [] })).toEqual([]);
      expect(issues({ geometry, value: {} })).toEqual([]);
    }
    expect(issues({ geometry: "ranges", value: { label: ["a"] } })).toEqual([`svs.params.missing ${AT}/value`]);
  });
});

describe("fields", () => {
  const file = { name: ["n1"], label: ["l1"], category: ["c1"], score: ["1"] };

  it("renames simultaneously: a swap, a fan-out and a chain", () => {
    expect(renameFields(file, { category: "label", label: "category" })).toEqual({
      name: ["n1"],
      score: ["1"],
      category: ["l1"],
      label: ["c1"],
    });
    expect(renameFields(file, { label: "name", category: "name" })).toEqual({
      score: ["1"],
      label: ["n1"],
      category: ["n1"],
    });
    expect(renameFields(file, { id: "label", label: "category", category: "name" })).toEqual({
      score: ["1"],
      id: ["l1"],
      label: ["c1"],
      category: ["n1"],
    });
    expect(renameFields(file, { score: null })).toEqual({ name: ["n1"], label: ["l1"], category: ["c1"] });
  });

  it("with an inline value is svs.params.exclusive", () => {
    expect(issues({ geometry: "ranges", value: { start: [0] }, fields: { label: "name" } })).toEqual([
      `svs.params.exclusive ${AT}/fields`,
    ]);
  });
});

describe("file cells typed by role", () => {
  it("01 stays text in an id or a category; indices and values are numbers; a numeric field too", () => {
    const t = typeColumns({
      id: ["01", "1"],
      category: ["01", "10"],
      target: ["1", "2"],
      beg: ["27", "28"],
      start: ["0", "7"],
      stop: ["4", null],
      value: ["1E5", "0.5"],
      score: ["007", "2"],
      note: ["a", "1"],
    });
    expect(t.id).toEqual(["01", "1"]);
    expect(t.category).toEqual(["01", "10"]);
    expect(t.target).toEqual(["1", "2"]);
    expect(t.beg).toEqual(["27", "28"]);
    expect(t.start).toBeInstanceOf(Int32Array);
    expect(t.stop).toBeInstanceOf(Float64Array);
    expect(Array.from(t.stop as Float64Array)).toEqual([4, Number.NaN]);
    expect(Array.from(t.value as Float64Array)).toEqual([100000, 0.5]);
    expect(Array.from(t.score as Float64Array)).toEqual([7, 2]);
    expect(t.note).toEqual(["a", "1"]);
  });
});

describe("typed columns and codes", () => {
  it("typed numbers and typed codes validate; bad ones are one issue each, with their count", () => {
    const value = {
      start: Int32Array.of(0, 5, 9),
      stop: Float64Array.of(4, Number.NaN, 12),
      value: Float32Array.of(0.15, 1, 2),
      label: { dict: ["TA", "DBD"], codes: Uint8Array.of(0, 1, 0) },
    };
    expect(issues({ geometry: "ranges", value })).toEqual([]);
    const bad = { ...value, start: Int32Array.of(0, -1, -2), label: { dict: ["TA"], codes: Uint8Array.of(0, 3, 3) } };
    const d = fullDoc();
    dataOf(d)[8].params = { geometry: "ranges", target: "p53", value: bad };
    const found = validate(d, reg);
    expect(found.map((i) => `${i.code} ${i.path}`)).toEqual([
      `svs.columns.type ${AT}/value/start/1`,
      `svs.columns.type ${AT}/value/label/codes/1`,
    ]);
    expect(found.map((i) => i.message)).toEqual([
      "expected an index (int >= 0), got -1 (2 rows)",
      "expected a code in 0..0 or null, got 3 (2 rows)",
    ]);
  });

  it("a coded column's entry is checked once for all its rows", () => {
    const value = { start: { dict: ["x", "y"], codes: [0, 1, 1, null] } };
    expect(issues({ geometry: "ranges", value })).toEqual([`svs.columns.type ${AT}/value/start/dict/0`]);
  });
});

describe("normalize", () => {
  it("a column named json is a column of numbers, not of JSON", () => {
    const d = fullDoc();
    dataOf(d)[8].params.value = { start: [0, 1], json: [1, 2] };
    const n: J = normalize(d, reg);
    expect(dataOf(n)[8].params.value.json).toEqual([1, 2]);
    dataOf(d)[8].params.value = { start: [0, 1], extra: { json: [{ a: 1 }, null] } };
    expect(dataOf(normalize(d, reg))[8].params.value.extra).toEqual({ json: [{ a: 1 }, null] });
  });

  it("a split by a column an inline table has passes, in either form", () => {
    const d = fullDoc();
    const layer = sectionOf(d)[3].children[0].children[0];
    dataOf(d)[9].params.value = { beg: [1, "1A"], label: ["a", "b"] };
    layer.params.split = { by: "label" };
    expect(validate(d, reg).filter((i) => i.code === "svs.layer.split")).toEqual([]);
    layer.params.split = { by: "nothing" };
    expect(validate(d, reg).map((i) => i.code)).toContain("svs.layer.split");
  });
});
