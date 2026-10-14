/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CF-02: a column named like an `Object.prototype` member (`constructor`, `toString`) or `__proto__` is a column like
// any other, in every table form and through `fields` (seqspec:format#annotation-tables).
import { Annotations } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { compile } from "../../src/format/index.js";
import { KernelFormats, tableColumns } from "../../src/format/load/formats.js";
import { renameFields, typeColumns } from "../../src/format/tables/table.js";
import { CASES } from "./cases.js";
import { fixtureLoaders, registries } from "./load-fixtures.js";

/** The key, read as data (a literal `.__proto__` would read the prototype). */
const PROTO = "__proto__";
const keysOf = (o: object) => [
  Object.keys(o),
  Object.getPrototypeOf(o) === Object.prototype || Object.getPrototypeOf(o) === null,
];

describe("data-keyed objects", () => {
  it("keeps constructor, toString and __proto__ columns in the rows and columns forms", () => {
    const rows = tableColumns(
      JSON.parse('[{"start": 0, "stop": 2, "constructor": "pET", "toString": "a", "__proto__": "x"}]'),
    );
    expect(rows && keysOf(rows.columns)).toEqual([["start", "stop", "constructor", "toString", "__proto__"], true]);
    expect(rows?.columns.constructor).toEqual(["pET"]);
    expect(rows?.columns[PROTO]).toEqual(["x"]);
    const cols = tableColumns(JSON.parse('{"start": [0], "stop": [2], "__proto__": ["x"]}'));
    expect(cols && keysOf(cols.columns)).toEqual([["start", "stop", "__proto__"], true]);
  });

  it("keeps them in a TSV, typed, and through fields", async () => {
    const tsv = KernelFormats.find((f) => f.name === "tsv");
    const text = "start\tstop\tconstructor\t__proto__\n0\t2\tpET\t7\n";
    const r = await tsv?.read(new TextEncoder().encode(text), {} as never, { report: () => {} } as never);
    const columns = (r as unknown as { columns: Record<string, unknown[]> }).columns;
    expect(keysOf(columns)).toEqual([["start", "stop", "constructor", "__proto__"], true]);
    const typed = typeColumns(columns as never);
    expect(keysOf(typed)).toEqual([["start", "stop", "constructor", "__proto__"], true]);
    expect(Array.from(typed[PROTO] as Float64Array)).toEqual([7]);
    const renamed = renameFields(columns, JSON.parse('{"__proto__": "constructor", "constructor": null}'));
    expect(keysOf(renamed)).toEqual([["start", "stop", "__proto__"], true]);
    expect(renamed[PROTO]).toEqual(["pET"]);
  });

  it("reads a constructor column's values end to end", async () => {
    const doc = CASES.case1();
    const data = (doc.root.children ?? []).find((c) => c.kind === "data") as {
      children: { ref?: string; params: Record<string, unknown> }[];
    };
    const regions = data.children.find((c) => c.ref === "regions");
    if (!regions) throw new Error("no regions");
    regions.params.value = [
      { start: 0, stop: 44, constructor: "pET-28" },
      { start: 101, stop: 292, constructor: "pGEX" },
    ];
    const { set, issues } = await compile(doc, registries(), fixtureLoaders());
    expect(issues.filter((i) => i.severity === "error")).toEqual([]);
    const a = set.annotations.find((x) => x && set.strings.get(x.key) === "regions");
    if (!a) throw new Error("no annotation");
    expect([0, 1].map((r) => Annotations.rowValue(a, set, r, "constructor"))).toEqual(["pET-28", "pGEX"]);
  });

  it("reads a __proto__ column's values end to end, in the rows and columns forms (normalized params included)", async () => {
    for (const value of [
      '[{"start": 0, "stop": 44, "__proto__": "x"}, {"start": 101, "stop": 292, "__proto__": "y"}]',
      '{"start": [0, 101], "stop": [44, 292], "__proto__": ["x", "y"]}',
    ]) {
      const doc = CASES.case1();
      const data = (doc.root.children ?? []).find((c) => c.kind === "data") as {
        children: { ref?: string; params: Record<string, unknown> }[];
      };
      const regions = data.children.find((c) => c.ref === "regions");
      if (!regions) throw new Error("no regions");
      regions.params.value = JSON.parse(value);
      const { set, issues } = await compile(doc, registries(), fixtureLoaders());
      expect(issues.filter((i) => i.severity === "error")).toEqual([]);
      const a = set.annotations.find((x) => x && set.strings.get(x.key) === "regions");
      if (!a) throw new Error("no annotation");
      expect(
        [0, 1].map((r) => Annotations.rowValue(a, set, r, PROTO)),
        value,
      ).toEqual(["x", "y"]);
    }
  });
});
