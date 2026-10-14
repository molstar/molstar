/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shared param types (seqspec:format#positions, seqspec:format#themes): labels, colors, palettes, theme references,
// glyphs, tables, provenance, filters, URIs.

import type { Glyph, GlyphShape, ThemeRef } from "../doc.js";
import { Color } from "../palettes/color.js";
import { ColorDictionaryNames } from "../palettes/dictionaries.js";
import { isMapLetter } from "../palettes/letter-maps.js";
import { LetterDictionaryNames } from "../palettes/letters.js";
import { ColorListNames } from "../palettes/lists.js";
import { PaletteDefaults } from "../palettes/palette.js";
import type { ChannelKind } from "../registry.js";
import {
  bool,
  type CheckContext,
  dict,
  type FieldType,
  float,
  int,
  isNamespaced,
  isPlainObject,
  list,
  literal,
  nullable,
  ptr,
  ref,
  refined,
  report,
  reportUnknown,
  show,
  str,
  tuple,
  typeError,
  union,
} from "../schema/field-schema.js";
import {
  checkParams,
  normalizeParams,
  OptionalField,
  object,
  paramsField,
  RequiredField,
  SimpleParamsSchema,
  UnionParamsSchema,
} from "../schema/params-schema.js";
import { resolveArchivePath } from "../uri-path.js";

const LABEL_INT = /^-?[0-9]+([A-Za-z]+|\.[0-9]+)?$/;

/** A position label: an int, a string in an integer numbering (`"27A"`, `"111.1"`), or a name (`"N13"`). */
export const label: FieldType<number | string> = {
  name: "Label",
  check(v, path, ctx) {
    if (Number.isInteger(v) || (typeof v === "string" && /^\S+$/.test(v))) return;
    typeError(ctx, path, 'Label (315, "27A", "111.1", or a name such as "N13")', v, "svs.position.syntax");
  },
};
/** A label in an integer numbering: an int or `^-?[0-9]+([A-Za-z]+|\.[0-9]+)?$`. */
export const integerLabel: FieldType<number | string> = {
  name: "Label",
  check(v, path, ctx) {
    if (Number.isInteger(v) || (typeof v === "string" && LABEL_INT.test(v))) return;
    typeError(ctx, path, 'an integer label (315, "27A", "111.1")', v, "svs.position.syntax");
  },
};

export const color = refined(
  "Color",
  Color.isValid,
  "svs.color.invalid",
  "a Color: #rgb, #rrggbb, a color name, or a/b",
);
export const httpsUrl = refined("str", (s) => s.startsWith("https://"), "svs.source.uri", "an absolute https:// URL");

/**
 * A source's `uri`: `https://`, or inside an archive (`ctx.base` set) a path within it (seqspec:sources#archives,
 * seqspec:sources#data-by-reference).
 */
export const uri: FieldType<string> = {
  name: "str",
  check(v, path, ctx) {
    if (typeof v !== "string") return typeError(ctx, path, "str", v);
    if (v.startsWith("https://") || (ctx.base !== undefined && resolveArchivePath(v, ctx.base) !== null)) return;
    const what = ctx.base === undefined ? "an absolute https:// URL" : "an https:// URL or a path inside the archive";
    report(ctx, "svs.source.uri", path, `${show(v)} is not ${what}`);
  },
};

const named = (names: readonly string[], what: string) =>
  refined("str", (s) => names.includes(s), "svs.color.invalid", `${what} (${names.slice(0, 4).join(", ")}, ...)`);
const listNames = named(ColorListNames, "a color list name");
const flowColor = OptionalField(nullable(union(literal("auto"), color)), null, '"auto", a color, or null: not colored');
const domain = OptionalField(tuple(nullable(float), nullable(float)), [null, null], "[min, max]; null: from the data");
const mode = OptionalField(literal("normalized", "absolute"), "normalized", "normalized: over value_domain");
const reverse = OptionalField(bool, false, "reverse the list");
const P = PaletteDefaults;

/** MVS's palettes, discriminated by `kind`. */
export const PaletteSchema = UnionParamsSchema("kind", "the palette's kind", {
  categorical: SimpleParamsSchema({
    colors: OptionalField(
      union(named([...ColorListNames, ...ColorDictionaryNames], "a list or dictionary name"), list(color), dict(color)),
      P.categorical.colors as unknown,
      "a list name, a dictionary name, colors, or {value: Color}",
    ),
    repeat_color_list: OptionalField(bool, true, "start the list over when it runs out (MVS: false)"),
    sort: OptionalField(literal("none", "lexical", "numeric"), "none", "none: first occurrence"),
    sort_direction: OptionalField(literal("ascending", "descending"), "ascending", ""),
    case_insensitive: OptionalField(bool, false, ""),
    missing_color: OptionalField(nullable(color), null, "null: not colored by this layer"),
  }),
  discrete: SimpleParamsSchema({
    colors: OptionalField(
      union(
        listNames,
        list(color),
        list(tuple(color, float)),
        list(tuple(nullable(color), nullable(float), nullable(float))),
      ),
      P.discrete.colors as unknown,
      "a list name, colors, [Color, checkpoint] pairs, or [Color | null, from | null, to | null] ranges",
    ),
    reverse,
    mode,
    value_domain: domain,
  }),
  continuous: SimpleParamsSchema({
    colors: OptionalField(
      union(listNames, list(color), list(tuple(color, float))),
      P.continuous.colors as unknown,
      "a list name, colors, or [Color, checkpoint] pairs",
    ),
    reverse,
    mode,
    value_domain: domain,
    underflow_color: flowColor,
    overflow_color: flowColor,
  }),
});
export const palette = paramsField(PaletteSchema, "Palette");

/** A letter dictionary or one of MVS's dictionaries, by name. */
export const letterDictionary = union(
  named([...LetterDictionaryNames, ...ColorDictionaryNames], "a dictionary name"),
  dict(color),
);

/**
 * A letter map (seqspec:format#layer-types): a registered map's name, or a dict from single stored letters to the
 * letter shown (printable ASCII; unlisted letters are kept). Errors name the registered maps and the entry at fault.
 */
export const letterMap: FieldType<string | Record<string, string>> = {
  name: "str | dict<str, str>",
  check(v, path, ctx) {
    const names = ctx.registries.list("letterMaps").map((e) => e.name);
    if (typeof v === "string") {
      if (!names.includes(v))
        typeError(ctx, path, `a letter map name (${names.join(", ")}) or a {letter: letter} dict`, v);
      return;
    }
    if (!isPlainObject(v))
      return typeError(ctx, path, `a letter map name (${names.join(", ")}) or a {letter: letter} dict`, v);
    const one = "one printable character";
    for (const [k, x] of Object.entries(v)) {
      if (!isMapLetter(k)) typeError(ctx, ptr(path, k), `a key of ${one} (the stored letter)`, k);
      else if (typeof x !== "string" || !isMapLetter(x)) typeError(ctx, ptr(path, k), `${one} (the letter shown)`, x);
    }
  },
};

/**
 * `{theme, params, fallback}`: params checked by the registered theme's table for the channel's kind, at most 4 deep; a
 * theme that does not serve the kind is `svs.theme.channel`.
 */
export function themeRef(kind: ChannelKind): FieldType<ThemeRef> {
  const code = kind === "glyph" ? "svs.glyph.invalid" : "svs.color.invalid";
  const self: FieldType<ThemeRef> = {
    name: `ThemeRef<${kind}>`,
    check(v, path, ctx: CheckContext) {
      if (!isPlainObject(v) || typeof v.theme !== "string") return typeError(ctx, path, `a ${kind} ThemeRef`, v);
      for (const k of Object.keys(v)) {
        if (k !== "theme" && k !== "params" && k !== "fallback") {
          reportUnknown(
            ctx,
            "svs.params.unknown",
            ptr(path, k),
            `unknown key "${k}"; a ThemeRef has theme, params, fallback`,
          );
        }
      }
      if (ctx.depth >= 4) return report(ctx, code, path, "themes nest at most 4 deep");
      const entry = ctx.registries.get("themes", v.theme);
      const params = entry?.params[kind];
      ctx.depth++;
      if (params) {
        const n = ctx.issues.length;
        checkParams(params, v.params ?? {}, ptr(path, "params"), ctx);
        if (ctx.issues.length === n && entry?.check && isPlainObject(v.params ?? {}))
          entry.check((v.params ?? {}) as Record<string, unknown>, ptr(path, "params"), ctx);
      } else if (entry) {
        report(ctx, "svs.theme.channel", ptr(path, "theme"), `"${v.theme}" is not a ${kind} theme`);
      } else if (!isNamespaced(v.theme)) {
        const known = ctx.registries.list("themes").flatMap((e) => (e.params[kind] ? [e.name] : []));
        reportUnknown(ctx, code, ptr(path, "theme"), `unknown ${kind} theme "${v.theme}"; known: ${known.join(", ")}`);
      } else if (v.params !== undefined && !isPlainObject(v.params)) {
        typeError(ctx, ptr(path, "params"), "object", v.params);
      }
      if (isNamespaced(v.theme) && v.fallback == null) {
        report(
          ctx,
          "svs.layer.fallback",
          ptr(path, "fallback"),
          `extension theme "${v.theme}" needs a kernel fallback`,
        );
      }
      if (v.fallback != null) self.check(v.fallback, ptr(path, "fallback"), ctx);
      ctx.depth--;
    },
    normalize(v, reg) {
      if (!isPlainObject(v) || typeof v.theme !== "string") return v;
      const schema = reg.get("themes", v.theme)?.params[kind];
      const params = schema ? normalizeParams(schema, v.params ?? {}, reg) : (v.params ?? {});
      return { ...v, theme: v.theme, params, fallback: v.fallback ? self.normalize?.(v.fallback, reg) : null };
    },
  };
  return self;
}

/** A rule's measure per element: NaN where it does not apply, and a clause on NaN fails. */
export const measure = union(
  literal("value"),
  object({ field: RequiredField(str, "an item or row field") }),
  object({ column: RequiredField(ref("annotation", "computed"), "values over the alignment's columns") }),
  object({ share: RequiredField(str, '"letter", or a set of letters written as a string, among the column\'s') }),
  object({ rank: RequiredField(literal("letter"), "the cell's letter's rank in its column, 0 the most frequent") }),
  object({ equals_row: RequiredField(nullable(int), "1 where the letter equals that row's; null: the reference row") }),
);
const clause = object({
  measure: OptionalField(nullable(measure), null, "lies in [min, max]"),
  min: OptionalField(nullable(float), null, ""),
  max: OptionalField(nullable(float), null, ""),
  in: OptionalField(nullable(list(str)), null, "the element's letter (a string field's value) is one of these"),
});

/** The `rules` theme's params for one channel kind: the first rule that holds gives its value. */
export function rulesSchema(kind: ChannelKind, value: FieldType<unknown>) {
  const rule = object({
    when: OptionalField(list(clause), [], "every clause holds; []: always"),
    [kind]: OptionalField(nullable(value), null, "its value"),
    theme: OptionalField(nullable(themeRef(kind)), null, "instead of a value, a nested theme's"),
    ...(kind === "color" ? { alpha: OptionalField(nullable(measure), null, "fades it, clamped to 0..1") } : {}),
    label: OptionalField(nullable(str), null, "names it in the legend"),
  });
  return SimpleParamsSchema({
    rules: RequiredField(list(rule), "the first that holds wins"),
    default: OptionalField(nullable(themeRef(kind)), null, "null: no value"),
    profile: OptionalField(nullable(ref("computed")), null, "what share and rank read; null: the alignment's columns"),
  });
}

export const GlyphShapes: readonly GlyphShape[] = [
  "circle",
  "diamond",
  "triangle",
  "square",
  "pin",
  "star",
  "cross",
  "tick",
];
const GLYPH_KEYS = ["shape", "char", "path", "box", "label", "fallback"];
const graphemes = new Intl.Segmenter(undefined, { granularity: "grapheme" });

/** A shape name, or `{shape | char | path, box?, label?, fallback?}`. */
export const glyph: FieldType<Glyph> = {
  name: "Glyph",
  check(v, path, ctx) {
    const bad = (at: string, why: string) => report(ctx, "svs.glyph.invalid", at, why);
    const shapes = GlyphShapes.join(", ");
    if (typeof v === "string") {
      if (!GlyphShapes.includes(v as GlyphShape)) bad(path, `"${v}" is not a shape: ${shapes}`);
      return;
    }
    if (!isPlainObject(v)) return typeError(ctx, path, "Glyph", v);
    for (const k of Object.keys(v)) {
      if (!GLYPH_KEYS.includes(k)) reportUnknown(ctx, "svs.glyph.invalid", ptr(path, k), `unknown key "${k}"`);
    }
    const given = ["shape", "char", "path"].filter((k) => v[k] !== undefined);
    if (given.length !== 1) bad(path, `give exactly one of shape, char, path (given: ${given.join(", ") || "none"})`);
    for (const k of ["shape", "fallback"]) {
      if (v[k] !== undefined && !GlyphShapes.includes(v[k] as GlyphShape)) bad(ptr(path, k), `not a shape: ${shapes}`);
    }
    if (v.char !== undefined && (typeof v.char !== "string" || [...graphemes.segment(v.char)].length !== 1)) {
      bad(ptr(path, "char"), "char is one grapheme");
    }
    if (v.path !== undefined) {
      if (typeof v.path !== "string" || v.path.length > 4096)
        bad(ptr(path, "path"), "path is SVG path data, at most 4 KB");
      if (v.box === undefined) bad(ptr(path, "box"), "a path needs box: [w, h]");
    } else if (v.box !== undefined) bad(ptr(path, "box"), "box goes with path");
    const box = v.box;
    const size = (x: unknown) => typeof x === "number" && Number.isFinite(x) && x > 0;
    if (box !== undefined && !(Array.isArray(box) && box.length === 2 && box.every(size))) {
      bad(ptr(path, "box"), "box is [w, h]: two numbers above 0");
    }
    if (v.label !== undefined && typeof v.label !== "string") bad(ptr(path, "label"), "label is a string");
  },
};

export const provenance = object({
  source: RequiredField(str, "who or what made it"),
  method: OptionalField(nullable(str), null, "how, in an open vocabulary: sifts, alignment, identity, curated, ..."),
  note: OptionalField(nullable(str), null, ""),
  confidence: OptionalField(nullable(float), null, "0..1"),
  references: OptionalField(list(str), [], "DOIs or URLs"),
  parameters: OptionalField(nullable(object()), null, "a program's parameters as given"),
});

export const filter = object({
  field: OptionalField(nullable(str), null, "null: category"),
  values: RequiredField(list(union(str, float)), "rows whose field holds one of these are kept"),
});

const keyTable = RequiredField(ref("annotation", "computed"), "a key table: a rows annotation or a computed one");
/** One key of an order (keyed sets): the member's data order, the tree's, or a key table's field. */
export const orderKey = union(
  literal("data", "tree"),
  object({
    data: keyTable,
    field: RequiredField(str, "a field of the key table"),
    direction: OptionalField(literal("asc", "desc"), "asc", ""),
  }),
);
/** `order`: one key or several, lexicographic; ties in data order, missing values last. */
export const orderField = nullable(union(orderKey, list(orderKey)));
/** A keyed set's filter: a member is kept when its field holds one of `values` and lies in `[min, max]`. */
export const keyFilter = object({
  data: keyTable,
  field: RequiredField(str, "a field of the key table"),
  values: OptionalField(nullable(list(union(str, float, bool))), null, "kept: the field holds one of these"),
  min: OptionalField(nullable(float), null, "kept: the field is at least this"),
  max: OptionalField(nullable(float), null, "kept: the field is at most this"),
});
