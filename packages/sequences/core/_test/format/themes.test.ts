/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Channels and themes per kind (seqspec:format#extensions; seqspec:format#layer-types, seqspec:format#themes; S4
// decision 6, VW4a): each layer type's channels, `ThemeRef`s checked by the theme's table for the channel's kind, the
// `rules` params, and the named rule tables' `expand`, which resolution applies (resolve.test.ts) so a viewer
// implements `rules` alone.

import { describe, expect, it } from "vitest";
import {
  type ChannelKind,
  Color,
  countColumn,
  numericDefaultColor,
  type ThemeRef,
  themeRef,
  validate,
} from "../../src/format/index.ts";
import { clustalXColor, createColumnComposition } from "../../src/format/palettes/clustal-x.ts";
import { quiet } from "../../src/format/schema/field-schema.ts";
import { KernelLayerTypes } from "../../src/format/tables/kernel.ts";
import { type J, kernel } from "./fixtures.ts";

const reg = kernel();
const issuesOf = (kind: ChannelKind, ref: unknown) => {
  const ctx = quiet(reg);
  Object.assign(ctx, { lenient: false });
  themeRef(kind).check(ref, "/t", ctx);
  return ctx.issues.map((i) => `${i.code} ${i.path}`);
};

describe("channels", () => {
  it("every type has color first; ranges, points, variants (ST4b) and pairs a glyph; the defaults are the table's", () => {
    const of = (name: string) => KernelLayerTypes.find((t) => t.name === name)?.channels ?? [];
    for (const t of KernelLayerTypes) expect(t.channels[0]?.name).toBe("color");
    const glyphs = KernelLayerTypes.filter((t) => t.channels.some((c) => c.kind === "glyph")).map((t) => t.name);
    expect(glyphs).toEqual(["ranges", "points", "variants", "pairs"]);
    expect(of("letters")[0].default).toEqual({ theme: "letter" });
    expect(of("values")[0].default).toBeNull();
    expect(of("ranges")[1].default).toEqual({ theme: "uniform", params: { glyph: "tick" } });
  });
});

describe("ThemeRef per kind", () => {
  it("uniform serves color and glyph, each with its own table; a theme of another kind is svs.theme.channel", () => {
    expect(issuesOf("color", { theme: "uniform", params: { color: "red" } })).toEqual([]);
    expect(issuesOf("glyph", { theme: "uniform", params: { glyph: "star" } })).toEqual([]);
    expect(issuesOf("glyph", { theme: "uniform", params: { color: "red" } })).toContain(
      "svs.params.missing /t/params/glyph",
    );
    expect(issuesOf("color", { theme: "field_map", params: { map: {} } })).toEqual(["svs.theme.channel /t/theme"]);
    expect(issuesOf("glyph", { theme: "letter" })).toEqual(["svs.theme.channel /t/theme"]);
  });

  it("rules: rules, clauses and measures checked; a nested theme of the same kind", () => {
    const rules = (r: unknown[], extra = {}) => ({ theme: "rules", params: { rules: r, ...extra } });
    const ok = rules(
      [
        { when: [{ measure: { share: "letter" }, min: 0.5 }], color: "#ff0000", label: "conserved" },
        { when: [{ measure: { equals_row: null }, min: 1, in: ["A", "G"] }], theme: { theme: "letter" } },
        { when: [{ measure: { column: "cons" }, min: 0.2 }], theme: { theme: "letter" }, alpha: { column: "cons" } },
        { when: [{ measure: "value", max: 3 }, { measure: { field: "score" } }, { measure: { rank: "letter" } }] },
      ],
      { default: { theme: "uniform", params: { color: "gray" } } },
    );
    expect(issuesOf("color", ok)).toEqual([]);
    expect(issuesOf("glyph", rules([{ when: [], glyph: "star" }]))).toEqual([]);
    expect(issuesOf("color", rules([{ when: [{ measure: { share: 3 } }] }])).length).toBeGreaterThan(0);
    expect(issuesOf("color", rules([{ when: [{ measure: { nope: 1 } }] }])).length).toBeGreaterThan(0);
    expect(issuesOf("color", rules([{ when: [], color: "nope" }]))).toEqual([
      "svs.color.invalid /t/params/rules/0/color",
    ]);
    expect(issuesOf("color", rules([{ when: [], theme: { theme: "field_map", params: { map: {} } } }]))).toEqual([
      "svs.theme.channel /t/params/rules/0/theme/theme",
    ]);
    expect(issuesOf("color", { theme: "rules", params: {} })).toEqual(["svs.params.missing /t/params/rules"]);
  });
});

describe("named rule tables", () => {
  const expand = (theme: string, params: Record<string, unknown> = {}): ThemeRef => {
    const e = reg.get("themes", theme);
    const normalized = themeRef("color").normalize?.({ theme, params }, reg) as ThemeRef;
    const out = e?.expand?.(normalized.params ?? {}, "color") as ThemeRef;
    expect(issuesOf("color", out)).toEqual([]);
    return out;
  };

  it("identity_to_row and conservation_shading: one rule on their measure, the base theme, nothing else", () => {
    expect(expand("identity_to_row", { row: 3 })).toEqual({
      theme: "rules",
      params: {
        rules: [
          { when: [{ measure: { equals_row: 3 }, min: 1 }], theme: expect.objectContaining({ theme: "letter" }) },
        ],
      },
    });
    expect(expand("identity_to_row").params?.rules).toMatchObject([{ when: [{ measure: { equals_row: null } }] }]);
    const base = { theme: "uniform", params: { color: "red" } };
    expect(expand("conservation_shading", { data: "cons", threshold: 0.3, base })).toMatchObject({
      theme: "rules",
      params: {
        rules: [{ when: [{ measure: { column: "cons" }, min: 0.3 }], theme: base, alpha: { column: "cons" } }],
      },
    });
    expect(["uniform", "letter", "rules", "field_palette"].filter((t) => reg.get("themes", t)?.expand)).toEqual([]);
  });

  it("clustal_x's rules color every column as the Clustal X table does", () => {
    type Rule = { when: { measure: { share: string }; min: number; in: string[] }[]; color: string };
    const rules = expand("clustal_x").params?.rules as Rule[];
    expect(rules.length).toBeGreaterThan(30);
    // A plain evaluator of these rules (the viewer's is LY3's): share among the column's letters, gaps left out.
    const byRules = (column: string, letter: string) => {
      const letters = [...column.toUpperCase()].filter((c) => c >= "A" && c <= "Z");
      const share = (set: string) => letters.filter((c) => set.includes(c)).length / letters.length;
      const L = letter.toUpperCase();
      const r = rules.find((x) => x.when.every((c) => c.in.includes(L) && share(c.measure.share) >= c.min));
      return r ? Color.decodePrimary(r.color) : undefined;
    };
    const c = createColumnComposition();
    let seed = 7;
    const rnd = (n: number) => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed % n;
    };
    const ALPHA = "ACDEFGHIKLMNPQRSTVWY-";
    for (let k = 0; k < 2000; k++) {
      // Columns biased to a few letters, so the thresholds are crossed often; exact 50% and 60% shares included.
      const pool = Array.from({ length: 1 + rnd(4) }, () => ALPHA[rnd(ALPHA.length)]).join("");
      const column = Array.from({ length: 5 + rnd(16) }, () => pool[rnd(pool.length)]).join("");
      countColumn(column, c);
      for (const l of new Set(column.replaceAll("-", ""))) expect(byRules(column, l), column).toBe(clustalXColor(c, l));
    }
    for (const column of ["KKKRRRDEST", "KKKKRRRDES", "LLLLIV----", "NNNNNDDDDD", "CCCCCCCCCA"])
      for (const l of new Set(column.replaceAll("-", "")))
        expect(byRules(column, l)).toBe(clustalXColor(countColumn(column, c), l));
  });
});

describe("layer params (VW4a)", () => {
  const withLayer = (layer: J) => ({
    format: "seq_view_spec",
    version: "1.0.0",
    root: {
      kind: "root",
      children: [
        {
          kind: "data",
          children: [
            { kind: "sequence", ref: "s", params: { value: "MKVLA" } },
            { kind: "annotation", ref: "a", params: { geometry: "ranges", target: "s", value: { start: [0] } } },
            {
              kind: "annotation",
              ref: "m",
              params: { geometry: "matrix", target: "s", states: ["x"], value: { index: [0], value: [[1]] } },
            },
          ],
        },
        {
          kind: "view",
          ref: "v",
          params: { label: "V" },
          children: [{ kind: "section", ref: "sec", children: [{ kind: "track", ref: "t", children: [layer] }] }],
        },
      ],
    },
  });
  it("density on ranges, points and variants; ruler's of and numbering lists; matrix's domain", () => {
    for (const type of ["ranges", "points", "variants"])
      expect(
        validate(withLayer({ kind: "layer", ref: "l", params: { type, data: "a", density: "majority" } }), reg),
      ).toEqual([]);
    expect(
      validate(withLayer({ kind: "layer", ref: "l", params: { type: "letters", density: "majority" } }), reg),
    ).not.toEqual([]);
    const ruler = { kind: "layer", ref: "l", params: { type: "ruler", of: ["s"], numbering: ["uniprot:P1", "auth"] } };
    expect(validate(withLayer(ruler), reg)).toEqual([]);
    const matrix = { kind: "layer", ref: "l", params: { type: "matrix", data: "m", domain: [0, 1] } };
    expect(validate(withLayer(matrix), reg)).toEqual([]);
  });

  it("the numeric default: Blues from 0 without negatives, RdBu centered with them", () => {
    const blues = (d: number[]) => ({
      theme: "field_palette",
      params: { palette: { kind: "continuous", colors: "Blues", value_domain: d } },
    });
    expect(numericDefaultColor(2, 9, false)).toEqual(blues([0, 9]));
    expect(numericDefaultColor(2, 9, true)).toEqual(blues([2, 9]));
    expect(numericDefaultColor(-4, 1, false).params?.palette).toEqual({
      kind: "continuous",
      colors: "RdBu",
      value_domain: [-4, 4],
    });
  });
});
