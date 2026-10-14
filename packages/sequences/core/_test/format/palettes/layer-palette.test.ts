/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// One palette per layer (seqspec:format#themes; PL): built over all rows of all the layer's data, a normalized
// palette's null ends from the layer's domain, the legend's keys in the palette's order.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type Annotation, SeqSets } from "@molstar/seq-core/model";
import { describe, expect, it } from "vitest";
import { Color, ColorLists, compile, layerPalette, type Palette } from "../../../src/format/index.js";
import { UNNAMED_ITEM } from "../../../src/format/palettes/fields.js";
import { kernel } from "../fixtures.js";
import { loaders } from "../load-support.js";

const reg = kernel();
async function data() {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: "s", value: "MKVLAAGIVGLLLAG" });
  b.data.annotation({
    ref: "a",
    geometry: "ranges",
    target: s,
    value: { start: [0, 1, 2], category: ["SH3", "sh2", "Kinase"] },
  });
  b.data.annotation({ ref: "b", geometry: "ranges", target: s, value: { start: [3, 4], category: ["SH2", "Linker"] } });
  b.data.annotation({ ref: "n", geometry: "ranges", target: s, value: { start: [5, 6], id: ["x", null] } });
  b.data.annotation({ ref: "v", geometry: "values", target: s, value: { index: [0, 1, 2], value: [2, 10, 6] } });
  b.data.annotation({
    ref: "m",
    geometry: "matrix",
    target: s,
    states: ["x", "y"],
    value: {
      index: [0, 1],
      value: [
        [0, 20],
        [5, 15],
      ],
    },
  });
  b.view({ ref: "end", label: "V" }).track({ label: "L" }).letters();
  const { set, issues } = await compile(b.getDocument(), reg, loaders({}));
  expect(issues).toEqual([]);
  const ann = (k: string) => set.annotations.find((a) => a && set.strings.get(a.key) === k) as Annotation;
  return { set, ann };
}

describe("layerPalette", () => {
  it("first occurrence over every row of the data, annotations in data order", async () => {
    const { set, ann } = await data();
    const f = layerPalette(set, [ann("a"), ann("b")], { kind: "categorical" }, null, [Number.NaN, Number.NaN]);
    const c10 = ColorLists.Category10;
    expect(["SH3", "sh2", "Kinase", "SH2", "Linker"].map((v) => f(v))).toEqual([0, 1, 2, 3, 4].map((k) => c10[k]));
    // Without "a", "b"'s values take the first colors: the palette is the layer's data's, the same in every viewer.
    expect(layerPalette(set, [ann("b")], { kind: "categorical" }, null, [0, 0])("SH2")).toBe(c10[0]);
  });

  it("legend keys in the palette's order, case-folded, unnamed items left out", async () => {
    const { set, ann } = await data();
    const p: Palette = { kind: "categorical", sort: "lexical", case_insensitive: true };
    expect(layerPalette(set, [ann("a"), ann("b")], p, null, [0, 0]).keys).toEqual(["KINASE", "LINKER", "SH2", "SH3"]);
    const named = layerPalette(set, [ann("n")], { kind: "categorical" }, null, [0, 0]);
    expect(named.keys).toEqual(["x"]);
    expect(named(UNNAMED_ITEM)).toBe(ColorLists.Category10[1]); // the unnamed items' one value still has a color
    const dict = layerPalette(
      set,
      [ann("a")],
      { kind: "categorical", colors: { Kinase: "red", SH3: "blue" } },
      null,
      [0, 0],
    );
    expect(dict.keys).toEqual(["SH3", "Kinase"]);
  });

  it("a normalized palette's null ends from the layer's domain when the field is the value", async () => {
    const { set, ann } = await data();
    const viridis: Palette = { kind: "continuous", colors: "Viridis" };
    const v = ColorLists.Viridis;
    const at = (t: number) => {
      const x = t * (v.length - 1);
      const i = Math.floor(x);
      return i + 1 < v.length ? Color.interpolate(v[i], v[i + 1], x - i) : v[i];
    };
    // A 0..20 matrix: 15 is 0.75 of the way (the old viewer read [0, 1] and left it uncolored).
    expect(layerPalette(set, [ann("m")], viridis, null, [0, 20])(15)).toBe(at(0.75));
    // Without a resolved domain, the data's: the matrix's cells and the values' rows.
    expect(layerPalette(set, [ann("m")], viridis, null, [Number.NaN, Number.NaN])(15)).toBe(at(0.75));
    expect(layerPalette(set, [ann("v")], viridis, null, [Number.NaN, Number.NaN])(6)).toBe(at(0.5));
    // A served or given domain wins over the rows; a given end over the domain.
    expect(layerPalette(set, [ann("v")], viridis, null, [0, 12])(6)).toBe(at(0.5));
    const half = { ...viridis, value_domain: [null, 8] } as Palette;
    expect(layerPalette(set, [ann("v")], half, null, [0, 12])(4)).toBe(at(0.5));
    // Another field than the value: the domain is not its.
    expect(layerPalette(set, [ann("v")], viridis, "start", [0, 100])(1)).toBe(at(0.5));
  });

  it("decodeFirst: a split color's first half; undefined when invalid or not a string", () => {
    expect([Color.decodeFirst("#f00/#00f"), Color.decodeFirst("nope"), Color.decodeFirst(3)]).toEqual([
      0xff0000,
      undefined,
      undefined,
    ]);
  });
});

describe("a streamed annotation's values outside the palette's keys (LB-08)", () => {
  const NONE: [number, number] = [Number.NaN, Number.NaN];
  const c10 = ColorLists.Category10;
  /** A streamed annotation (no rows), its manifest's vocabulary as given. */
  const streamed = (vocabulary?: Record<string, string[]>) => {
    const set = SeqSets.build({
      molecules: [{ key: "p", residues: "ACDEFGHIKL" }],
      units: [{ key: "A", molecule: "p" }],
      annotations: [{ key: "g", geometry: "ranges", streamed: ["p"], columns: {} }],
    } as never).set;
    const a = set.annotations[0] as Annotation;
    if (vocabulary) a.custom.props.set("svs:stream", { manifest: { vocabulary } });
    return { set, a };
  };

  it("a list palette colors each by FNV-1a of its UTF-8 bytes into the list, the same in every palette made", () => {
    const { set, a } = streamed();
    const f = layerPalette(set, [a], { kind: "categorical" }, null, NONE);
    // FNV-1a (32-bit): "a" 0xe40c292c, "b" 0xe70c2de5, "c" 0xe60c2c52, "rs1042522" 0x23f12a94, "é" 0x1e9de8c1; mod 10.
    expect(["a", "b", "c", "rs1042522", "é"].map((v) => f(v))).toEqual([0, 7, 8, 6, 7].map((i) => c10[i]));
    // The legend cannot list them, and says so.
    expect([f.keys, f.unlisted]).toEqual([[], true]);
    // Another set and palette (a reload): the same colors; a missing value is still uncolored.
    const again = streamed();
    const g = layerPalette(again.set, [again.a], { kind: "categorical" }, null, NONE);
    expect(["a", "b", "c"].map((v) => g(v))).toEqual(["a", "b", "c"].map((v) => f(v)));
    expect([g(null), g(".")]).toEqual([undefined, undefined]);
    // `case_insensitive` hashes the value folded: "a" as "A" (0xc40bf6cc).
    expect(layerPalette(set, [a], { kind: "categorical", case_insensitive: true }, null, NONE)("a")).toBe(c10[2]);
  });

  it("a vocabulary's values keep their order's colors; one it lacks is hashed, the legend lists the vocabulary", () => {
    const { set, a } = streamed({ id: ["zz", "a"] });
    const f = layerPalette(set, [a], { kind: "categorical" }, null, NONE);
    expect([f("zz"), f("a"), f("b")]).toEqual([c10[0], c10[1], c10[7]]);
    expect([f.keys, f.unlisted]).toEqual([["zz", "a"], undefined]);
  });

  it("a dictionary palette and held rows hash nothing: a value they lack is the missing color", async () => {
    const { set, a } = streamed();
    const dict: Palette = { kind: "categorical", colors: { x: "red" }, missing_color: "gray" };
    expect(layerPalette(set, [a], dict, null, NONE)("a")).toBe(Color.decodePrimary("gray"));
    const held = await data();
    const f = layerPalette(held.set, [held.ann("a")], { kind: "categorical" }, null, NONE);
    expect([f("unknown"), f.unlisted]).toEqual([undefined, undefined]);
  });
});
