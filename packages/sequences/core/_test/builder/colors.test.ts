/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `b.colors`: the schemes as `scheme` theme refs in any layer's `color`, value schemes with their options, custom and
// derived schemes, the default per subject and the list per subject. The table's twin check against the format's
// `KernelColorSchemes` is in the format's schemes test (the builder has no runtime dependency).

import { describe, expect, expectTypeOf, it } from "vitest";
import { Colors, createSeqViewSpec, type KernelSchemeName, type SchemeRef } from "../../src/builder/index.ts";

describe("b.colors", () => {
  it("refs go into a layer's color as written; value schemes take a domain, a field and a direction", () => {
    const b = createSeqViewSpec({ title: "colors" });
    const aln = b.data.alignment({ ref: "a", value: ["MKV", "MRV"] });
    const cons = b.data.annotation({
      ref: "c",
      geometry: "values",
      target: aln,
      value: { index: [0, 1, 2], value: [1, 0.5, 1] },
    });
    const v = b.view({ ref: "v", label: "V" }).block({ alignment: aln });
    v.track({ ref: "t.rows", label: "Rows" }).rows({ ref: "rows", color: b.colors.clustalX });
    v.track({ ref: "t.cons", label: "Conservation" }).values({
      ref: "cons",
      data: cons,
      color: b.colors.viridis({ domain: [0, 1] }),
    });
    const layers = JSON.stringify(b.getDocument());
    expect(layers).toContain('"color":{"theme":"scheme","params":{"name":"clustal_x"}}');
    expect(layers).toContain('"color":{"theme":"scheme","params":{"name":"viridis","domain":[0,1]}}');
    expect(b.colors.rdBu({ field: "score", reverse: true })).toEqual({
      theme: "scheme",
      params: { name: "rd_bu", field: "score", reverse: true },
    });
  });

  it("custom, a registered name, and one derived from another with its overrides merged", () => {
    const own = Colors.custom({ kind: "letters", alphabets: ["dna"], colors: { A: "#0a0" }, default: "#ccc" });
    expect(own).toEqual({
      theme: "scheme",
      params: { custom: { kind: "letters", alphabets: ["dna"], colors: { A: "#0a0" }, default: "#ccc" } },
    });
    expect(Colors.scheme("org.example:brand")).toEqual({ theme: "scheme", params: { name: "org.example:brand" } });
    // CB-11: no overrides on either side, none written (format.md#builders: defaults are left out).
    expect(Colors.from(Colors.viridis(), { reverse: true }).params).toEqual({ name: "viridis", reverse: true });
    const once = Colors.from(Colors.zappo, { overrides: { C: "#f0c" } });
    expect(Colors.from(once, { overrides: { W: "#000" } }).params).toEqual({
      name: "zappo",
      overrides: { C: "#f0c", W: "#000" },
    });
  });

  it("auto picks by subject; list gives the schemes that fit, column tables only on an alignment", () => {
    expect(Colors.auto({ alphabet: "protein", alignment: true })).toBe(Colors.clustalX);
    expect(Colors.auto({ alphabet: "protein" })).toBe(Colors.zappo);
    expect(Colors.auto({ alphabet: "rna" })).toBe(Colors.nucleotide);
    expect(Colors.auto({ alphabet: "generic" })).toBeNull();
    expect(Colors.auto({ values: "diverging" })?.params.name).toBe("rd_bu");
    expect(Colors.auto({ categories: true })).toBe(Colors.okabeIto);
    expect(Colors.list({ alphabet: "protein" })).toEqual(["clustal_x", "zappo", "taylor", "hydrophobicity"]);
    expect(Colors.list({ alphabet: "dna", alignment: true })).toEqual([
      "percent_identity",
      "nucleotide",
      "purine_pyrimidine",
      "nucleotide_cb",
    ]);
    expect(Colors.list({ values: true })).toEqual(["viridis", "cividis", "magma", "blues", "rd_bu"]);
    expect(Colors.list({ categories: true })).toEqual(["category10", "okabe_ito", "tableau10"]);
  });

  it("types: names are the kernel's; a custom scheme's kind decides its fields", () => {
    expectTypeOf(Colors.list({ values: true })).toEqualTypeOf<KernelSchemeName[]>();
    expectTypeOf(Colors.zappo).toEqualTypeOf<SchemeRef>();
    // @ts-expect-error a ramp has no alphabets
    Colors.custom({ kind: "sequential", colors: ["#fff", "#000"], alphabets: ["dna"] });
    // @ts-expect-error letters need their alphabets
    Colors.custom({ kind: "letters", colors: { A: "#0a0" } });
    // @ts-expect-error not an alphabet
    Colors.list({ alphabet: "amino" });
  });
});
