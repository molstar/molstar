/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import {
  Color,
  ColorDictionaries,
  ColorDictionaryNames,
  ColorListNames,
  ColorLists,
  ColorNames,
  ElementSymbolColors,
  makePaletteFunction,
  type Palette,
  type PaletteValue,
} from "../../../src/format/index.js";

/**
 * Outputs of MVS's own code (vendored Mol* 5.11.0): its named lists and dictionaries, Mol*'s color names, and every
 * palette kind and mode evaluated by `palettePropsFromMVSPalette` + `makePaletteFunction` (recorded once by a
 * scratch script importing the vendored source; this test imports no Mol*).
 */
const recorded = JSON.parse(readFileSync(new URL("./mvs-recorded.json", import.meta.url), "utf8")) as {
  lists: Record<string, string[]>;
  dicts: Record<string, Record<string, string>>;
  colorNames: Record<string, string>;
  elementSymbolWithCarbon: Record<string, string>;
  cases: {
    name: string;
    palette: Palette;
    values: (string | number)[];
    queries: (string | number)[];
    expected: (string | null)[];
  }[];
};

const hex = (c: Color | undefined) => (c === undefined ? null : Color.toHex(c));
const primary = (c: string) => Color.toHex(Color.decodePrimary(c));

describe("named lists and dictionaries are MVS's", () => {
  it("lists", () => {
    expect(ColorListNames).toEqual(Object.keys(recorded.lists));
    for (const name of ColorListNames) expect(ColorLists[name].map(Color.toHex), name).toEqual(recorded.lists[name]);
  });

  it("dictionaries (split colors by their first half)", () => {
    expect(ColorDictionaryNames).toEqual(Object.keys(recorded.dicts));
    for (const name of ColorDictionaryNames) {
      const expected = Object.fromEntries(Object.entries(recorded.dicts[name]).map(([k, c]) => [k, primary(c)]));
      const actual = Object.fromEntries(Object.entries(ColorDictionaries[name]).map(([k, c]) => [k, Color.toHex(c)]));
      expect(actual, name).toEqual(expected);
    }
    expect(Object.keys(ColorDictionaries.CarbohydrateSymbol)).toHaveLength(1229);
  });

  it("element symbols with carbon, and color names", () => {
    const elements = Object.fromEntries(Object.entries(ElementSymbolColors).map(([k, c]) => [k, Color.toHex(c)]));
    expect(elements).toEqual(recorded.elementSymbolWithCarbon);
    expect(elements.C).toBe("#909090");
    expect(Object.fromEntries([...ColorNames].map(([k, c]) => [k, Color.toHex(c)]))).toEqual(recorded.colorNames);
  });
});

describe("palettes evaluate as MVS's", () => {
  for (const c of recorded.cases) {
    it(c.name, () => {
      // An MVS palette in SVS: MVS's default `repeat_color_list` (false) written out, SVS's being true
      // (seqspec:format#themes).
      const mvs = c.palette.kind === "categorical" ? { repeat_color_list: false, ...c.palette } : c.palette;
      const f = makePaletteFunction(mvs, c.values);
      expect(c.queries.map((q) => hex(f(q)))).toEqual(c.expected);
    });
  }

  it("covers every kind and mode", () => {
    const kinds = new Set(recorded.cases.map((c) => `${c.palette.kind}:${"mode" in c.palette ? c.palette.mode : ""}`));
    expect([...kinds].sort()).toEqual([
      "categorical:",
      "continuous:",
      "continuous:absolute",
      "discrete:",
      "discrete:absolute",
    ]);
  });
});

describe("palette rules", () => {
  const colorsOf = (palette: Palette, values: PaletteValue[], queries: PaletteValue[]) => {
    const f = makePaletteFunction(palette, values);
    return queries.map((q) => hex(f(q)));
  };

  it("first occurrence order, skipping missing values, `.` and `?`; they take `missing_color` and no list color", () => {
    const palette: Palette = {
      kind: "categorical",
      colors: ["red", "green", "blue"],
      repeat_color_list: false,
      missing_color: "gray",
    };
    const values = [null, ".", "b", undefined, "?", "a", "b", "c", "d"];
    expect(colorsOf(palette, values, ["b", "a", "c", "d", ".", "?", null, undefined])).toEqual([
      "#ff0000",
      "#008000",
      "#0000ff",
      "#808080",
      "#808080",
      "#808080",
      "#808080",
      "#808080",
    ]);
    // Without a `missing_color` they stay uncolored (MVS's), and the list's colors go to the real values alone.
    const none: Palette = { kind: "categorical", colors: ["red", "green", "blue"], repeat_color_list: false };
    expect(colorsOf(none, values, ["b", "a", ".", null, undefined])).toEqual(["#ff0000", "#008000", null, null, null]);
  });

  it("a list repeats by default (MVS's does not): the 11th value of Category10 takes its first color", () => {
    const values = Array.from({ length: 12 }, (_, i) => `Ig ${i + 1}`);
    const ten = ColorLists.Category10.map(Color.toHex);
    expect(colorsOf({ kind: "categorical" }, values, values)).toEqual([...ten, ten[0], ten[1]]);
    // Off, the values past the list take `missing_color` (unmapped: by default not colored).
    const off = colorsOf({ kind: "categorical", repeat_color_list: false }, values, values.slice(9));
    expect(off).toEqual([ten[9], null, null]);
    // A custom list repeats too; a dictionary's unlisted keys stay unmapped.
    expect(colorsOf({ kind: "categorical", colors: ["red", "blue"] }, ["a", "b", "c"], ["a", "b", "c"])).toEqual([
      "#ff0000",
      "#0000ff",
      "#ff0000",
    ]);
    expect(colorsOf({ kind: "categorical", colors: { a: "red" } }, ["a", "b"], ["a", "b"])).toEqual(["#ff0000", null]);
  });

  it("numbers and strings compare as strings", () => {
    const f = makePaletteFunction({ kind: "categorical" }, [1, "1", 2.5, "2.50"]);
    expect([f(1), f("1"), f(2.5), f("2.50")].map(hex)).toEqual(["#1f77b4", "#1f77b4", "#ff7f0e", "#2ca02c"]);
  });

  it("reverse flips a list, not explicit checkpoints", () => {
    const values = [0, 1];
    expect(colorsOf({ kind: "continuous", colors: ["red", "blue"], reverse: true }, values, [0, 1])).toEqual([
      "#0000ff",
      "#ff0000",
    ]);
    const pairs: Palette = {
      kind: "continuous",
      colors: [
        ["red", 0],
        ["blue", 1],
      ],
      reverse: true,
    };
    expect(colorsOf(pairs, values, [0, 1])).toEqual(["#ff0000", "#0000ff"]);
    expect(colorsOf({ kind: "discrete", colors: ["red", "blue"], reverse: true }, values, [0.2, 0.8])).toEqual([
      "#0000ff",
      "#ff0000",
    ]);
  });

  it("value_domain fixes the normalized domain; absolute mode ignores it", () => {
    const values = [0, 100];
    const palette: Palette = { kind: "continuous", colors: ["#000000", "#ffffff"], value_domain: [0, 10] };
    expect(colorsOf(palette, values, [0, 5, 10, 11])).toEqual(["#000000", "#7f7f7f", "#ffffff", null]);
    expect(colorsOf({ ...palette, value_domain: [null, 10] }, [4, 100], [4, 7])).toEqual(["#000000", "#7f7f7f"]);
    expect(colorsOf({ ...palette, mode: "absolute" }, values, [0.5, 2])).toEqual(["#7f7f7f", null]);
  });

  it('"auto" underflow and overflow take the end colors; a color or null otherwise', () => {
    const base: Palette = { kind: "continuous", colors: ["red", "blue"], value_domain: [0, 10] };
    const queries = [-1, 11];
    expect(colorsOf({ ...base, underflow_color: "auto", overflow_color: "auto" }, [], queries)).toEqual([
      "#ff0000",
      "#0000ff",
    ]);
    expect(colorsOf({ ...base, underflow_color: "white", overflow_color: null }, [], queries)).toEqual([
      "#ffffff",
      null,
    ]);
    expect(colorsOf(base, [], queries)).toEqual([null, null]);
  });

  it("values outside every bin of a discrete palette are not colored", () => {
    const palette: Palette = {
      kind: "discrete",
      mode: "absolute",
      colors: [
        ["red", 0, 1],
        ["blue", 2, 3],
      ],
    };
    expect(colorsOf(palette, [], [0.5, 1.5, 2.5, "x"])).toEqual(["#ff0000", null, "#0000ff", null]);
  });

  it("Color strings as MVS's ColorT", () => {
    for (const ok of ["#abc", "#A1B2C3", "red", "cornflower", "red/#00f"]) expect(Color.isValid(ok), ok).toBe(true);
    for (const bad of ["#abcd", "Red", "rgb(1,2,3)", "red/blue/green", "", 3])
      expect(Color.isValid(bad), `${bad}`).toBe(false);
    expect(Color.decode("#f0a")).toBe(0xff00aa);
    expect(Color.decode("nope")).toBeUndefined();
    expect(Color.decodePrimary("nope")).toBe(0);
  });
});
