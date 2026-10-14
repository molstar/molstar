/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import { fromJSON, normalize, printReference, toJSON, validate } from "../../src/format/index.ts";
import { dataOf, fullDoc, type J, kernel, sectionOf } from "./fixtures.ts";

describe("JSON", () => {
  it("fromJSON(toJSON(d)) deep-equals d, compact and pretty", () => {
    const d = fullDoc();
    for (const pretty of [false, true]) {
      const { doc, issues } = fromJSON(toJSON(d, { pretty }));
      expect(issues).toEqual([]);
      expect(doc).toStrictEqual(d);
    }
  });

  it("text that is not JSON, or not a SeqViewSpec, is svs.format.invalid", () => {
    expect(fromJSON("{").issues.map((i) => i.code)).toEqual(["svs.format.invalid"]);
    expect(fromJSON("{").doc).toBeNull();
    const other = fromJSON(JSON.stringify({ format: "mvs", version: "1.0.0" }));
    expect([other.doc, other.issues.map((i) => i.path)]).toEqual([null, ["/format"]]);
  });
});

describe("normalize", () => {
  const reg = kernel();

  it("is idempotent and keeps the document valid", () => {
    const n = normalize(fullDoc(), reg);
    expect(normalize(n, reg)).toStrictEqual(n);
    expect(validate(n, reg)).toEqual([]);
    expect(fromJSON(toJSON(n)).doc).toStrictEqual(n);
  });

  it("fills every default, nested ones included", () => {
    const n: J = normalize(fullDoc(), reg);
    expect(n.requires).toEqual([]);
    expect(n.snapshots[0]).toMatchObject({ description_format: "markdown", duration_ms: 1000, transition_ms: 0 });
    expect(n.root.children[0].ref).toBe("data");
    expect(n.root.children[0].params).toEqual({});
    const ranges = sectionOf(n)[3].children[0].children[0].params;
    expect(ranges).toMatchObject({ opacity: 1, filter: null, fallback: null, display: "lanes", labels: false });
    expect(ranges.color).toEqual({
      theme: "field_palette",
      params: {
        field: null,
        palette: {
          kind: "categorical",
          colors: "Set2",
          repeat_color_list: true,
          sort: "none",
          sort_direction: "ascending",
          case_insensitive: false,
          missing_color: null,
        },
      },
      fallback: null,
    });
    const rows = n.root.children[2].children[0].children[1].children[0].params;
    expect(rows.color.params.base).toEqual({ theme: "letter", params: { colors: null }, fallback: null });
    expect(dataOf(n)[4].custom).toEqual({
      mvs: { label_asym_id: "A", structure_ref: null, auth_asym_id: null, instance_id: null },
      "org.example": { anything: 1 },
    });
    expect(dataOf(n)[10].params).toMatchObject({ row_subset: null, gaps: "count", reference_row: null });
  });

  it("ranges' labels: a boolean hint, false unless the document asks", () => {
    const d: J = fullDoc();
    const ranges = sectionOf(d)[3].children[0].children[0].params;
    ranges.labels = true;
    expect(validate(d, reg)).toEqual([]);
    expect(sectionOf(normalize(d, reg))[3].children[0].children[0].params.labels).toBe(true);
    // Documents are immutable once validated: a change is a new params object.
    sectionOf(d)[3].children[0].children[0].params = { ...ranges, labels: "yes" };
    expect(validate(d, reg).map((i) => i.severity)).toContain("error");
  });

  it("writes tables as columns; rows, columns and dictionary codes normalize equal", () => {
    const rows = fullDoc();
    dataOf(rows)[8].params.value = [
      { start: 0, stop: 4, label: "TA" },
      { start: 5, label: "DBD", category: "domain" },
      { start: 9, stop: 12, label: "TA" },
    ];
    const columns = fullDoc();
    dataOf(columns)[8].params.value = {
      start: [0, 5, 9],
      stop: [4, null, 12],
      label: { dict: ["TA", "DBD"], codes: [0, 1, 0] },
      category: [null, "domain", null],
    };
    expect(validate(rows, reg)).toEqual([]);
    expect(validate(columns, reg)).toEqual([]);
    const a: J = normalize(rows, reg);
    expect(a).toStrictEqual(normalize(columns, reg));
    expect(dataOf(a)[8].params.value).toEqual({
      start: [0, 5, 9],
      stop: [4, null, 12],
      label: ["TA", "DBD", "TA"],
      category: [null, "domain", null],
    });
  });

  it("keeps JSON columns and what it does not know", () => {
    const d = fullDoc();
    dataOf(d)[8].params.value = { start: [0], extra: { json: [{ a: 1 }] } };
    dataOf(d).push({ kind: "computed", ref: "x", params: { provider: "org.example:x", options: { k: 1 } } });
    const n: J = normalize(d, reg);
    expect(dataOf(n)[8].params.value.extra).toEqual({ json: [{ a: 1 }] });
    expect(dataOf(n)[12].params).toEqual({ provider: "org.example:x", options: { k: 1 } });
  });
});

describe("the reference", () => {
  it("is generated from the tables and the registries", () => {
    const text = printReference(kernel());
    for (const kind of ["root", "annotation", "layer", "focus"]) expect(text).toContain(`### \`${kind}\``);
    expect(text).toContain('`geometry: "pairs"`');
    expect(text).toContain(
      "| `letters` | residues (element symbols for atoms), one alignment row's, or another polymer's through a mapping | none |",
    );
    expect(text).toContain("### `conservation_shading`");
    expect(text).toContain("### `mvs` on `unit`");
    expect(text).toContain("| `svs.layer.unplaced` | warning |");
    expect(text).toContain("| `opacity` | `float` | `1` |");
  });
});
