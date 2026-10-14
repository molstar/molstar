/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for LY3a's per-element path (seqspec:viewer#performance; seqspec:viewer#extension-points): `theme.value`
// (and `alpha`) per position for every built-in and every named rule table compiled into `rules`, on `heapScaling` (a
// step one element). The positive control: `letter`'s step allocating an object too (`leaking`), measured alike, is
// caught.

import { SeqSets } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { type Theme, ThemeLocation } from "../src/index.js";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.js";
import { context } from "./theme-fixtures.js";

const N = 50_000;
let sink = 0;

const set = SeqSets.build({
  molecules: [{ key: "p", residues: "ACDEFGHIKLMNPQRSTVWY" }],
  units: [{ key: "A", molecule: "p" }],
  alignments: [
    {
      key: "aln",
      gapped: ["ACDEFGHIKL", "ACDKFGHIKL", "WCD-FGHVKL", "ACKKEGHIKL"],
      flavor: "aligned",
      alphabet: "protein",
    },
  ],
  annotations: [
    {
      key: "f",
      geometry: "ranges",
      target: "p",
      columns: {
        start: [0, 2, 4, 6],
        stop: [2, 4, 6, 8],
        category: ["x", "y", "x", "z"],
        glyph: ["diamond", "★", "square", "pin"],
        color: ["#ff0000", "", "#0000ff", "not a color"],
        score: [1, 5, 9, 2],
      },
    },
    {
      key: "cons",
      geometry: "values",
      target: "aln",
      columns: { start: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], value: [1, 0.2, 0.75, 0.5, 1, 0.9, 0.4, 0.8, 1, 0.3] },
    },
    { key: "d", geometry: "values", target: "p", columns: { index: [0, 1, 2, 3], value: [1, 2, 2, 3] } },
  ],
}).set;

describe("theme values per element allocate nothing", () => {
  const onSeq = context(set, { data: ["f"] });
  const onAln = context(set, { data: ["f"], blocks: [{ alignment: "aln" }] });
  const onValues = context(set, { data: ["d"], domain: [0, 3] });
  const cases: [string, Theme][] = [
    ["letter", onSeq.color("letter")],
    ["field_palette (rows)", onSeq.color("field_palette")],
    ["field_palette (continuous, by value)", onValues.color("field_palette", { palette: { kind: "continuous" } })],
    ["field_palette (categorical, by value)", onValues.color("field_palette", { field: "value" })],
    [
      "field_palette (discrete, by value)",
      onValues.color("field_palette", { palette: { kind: "discrete", colors: ["#ff0000", "#00ff00", "#0000ff"] } }),
    ],
    ["uniform (color)", onSeq.color("uniform", { color: "#123456" })],
    ["field_color", onSeq.color("field_color")],
    ["by_unit", onSeq.color("by_unit")],
    ["by_copy_group", onSeq.color("by_copy_group")],
    ["clustal_x (rules)", onAln.color("clustal_x")],
    ["clustal_x on a sequence (rules)", onSeq.color("clustal_x")],
    ["identity_to_row over clustal_x (rules)", onAln.color("identity_to_row", { base: { theme: "clustal_x" } })],
    [
      "rules on a share, no alpha",
      onAln.color("rules", { rules: [{ when: [{ measure: { share: "letter" }, min: 0.5 }], color: "#ff0000" }] }),
    ],
    ["conservation_shading (rules, alpha)", onAln.color("conservation_shading", { data: "cons", threshold: 0.3 })],
    [
      "rules over fields, rank and a default",
      onAln.color("rules", {
        rules: [
          { when: [{ measure: { field: "score" }, min: 5 }], color: "#ff0000" },
          { when: [{ measure: { field: "category" }, in: ["z"] }], color: "#00ff00" },
          { when: [{ measure: { rank: "letter" }, max: 0 }], color: "#0000ff" },
        ],
        default: { theme: "letter" },
      }),
    ],
    // Color schemes (seqspec:format#color-schemes): the theme a scheme stands for, nothing of its own per element.
    ["scheme: zappo (letter)", context(set, { data: ["f"], type: "letters" }).color("scheme", { name: "zappo" })],
    [
      "scheme: a custom one with a default and overrides (letter)",
      context(set, { data: ["f"], type: "letters" }).color("scheme", {
        custom: { kind: "letters", alphabets: ["protein"], groups: [{ letters: "DE", color: "red" }], default: "#ccc" },
        overrides: { K: "#00f" },
      }),
    ],
    [
      "scheme: percent identity (rules)",
      context(set, { data: ["f"], blocks: [{ alignment: "aln" }], type: "rows" }).color("scheme", {
        name: "percent_identity",
      }),
    ],
    [
      "scheme: blues, dark (continuous, by value)",
      context(set, { data: ["d"], domain: [0, 3], type: "values", dark: true }).color("scheme", { name: "blues" }),
    ],
    ["field_glyph", onSeq.glyph("field_glyph")],
    ["field_map (glyph)", onSeq.glyph("field_map", { map: { x: "diamond", y: { char: "✶" } }, default: "cross" })],
    ["uniform (size)", onSeq.size("uniform", { size: 3 })],
    ["field_scale (rows, sqrt)", onSeq.size("field_scale", { field: "score", range: [2, 6], scale: "sqrt" })],
    ["field_scale (a string field's map)", onSeq.size("field_scale", { field: "category", map: { x: 3 }, default: 1 })],
    ["field_scale (by value)", onValues.size("field_scale", { range: [1, 5], scale: "log" })],
    ["rules (size)", onSeq.size("rules", { rules: [{ when: [{ measure: { field: "score" }, min: 5 }], size: 4 }] })],
    ["rules (glyph)", onSeq.glyph("rules", { rules: [{ when: [{ measure: "value", min: 0.5 }], glyph: "star" }] })],
  ];
  const f = set.annotations.findIndex((a) => a && set.strings.get(a.key) === "f");
  /** One element a step: block, unit, annotation, row, column, letter and value cycling. */
  const reading = (theme: Theme) => {
    const loc = ThemeLocation.create();
    return (k: number) => {
      loc.block = 0;
      loc.unit = 0;
      loc.annotation = k & 4 ? f : -1;
      loc.row = k & 3;
      loc.column = k % 10;
      loc.letter = 65 + (k % 26);
      loc.value = (k & 3) * (k & 8 ? 1 : 0.37);
      sink = (sink + theme.value(loc) + (theme.alpha ? theme.alpha(loc) : 0)) & 0xffff;
    };
  };
  for (const [name, theme] of cases) {
    it(name, () => {
      expect(theme).not.toBeNull();
      const r = heapScaling(N, reading(theme));
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
    });
  }

  it("the same elements allocating an object each are caught (positive control)", () => {
    expect(controlCaught(() => heapScaling(N, leaking(reading(cases[0][1]))))).toBe(true);
  });
});
