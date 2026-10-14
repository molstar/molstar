/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's registry entries: the fourteen layer types (seqspec:format#layer-types), the color and glyph themes
// (seqspec:format#themes) and the `svs.*` issue codes (seqspec:format#validation). Formats, services, providers and
// aligners are in load/.

import type { Severity } from "@molstar/seq-core/model";
import type { Geometry, ThemeRef } from "../doc.ts";
import { clustalXRules } from "../palettes/clustal-x.ts";
import { checkSchemeDef, checkSchemeParams, percentIdentityRules } from "../palettes/schemes.ts";
import type { ChannelKind, LayerTypeEntry, SvsRegistryEntry, ThemeEntry } from "../registry.ts";
import { bool, dict, float, int, list, literal, nullable, ref, str, tuple, union } from "../schema/field-schema.ts";
import {
  type Fields,
  OptionalField,
  object,
  RequiredField,
  SimpleParamsSchema,
  type ValuesFor,
} from "../schema/params-schema.ts";
import {
  color,
  glyph,
  keyFilter,
  letterDictionary,
  letterMap,
  orderField,
  palette,
  rulesSchema,
  themeRef,
} from "./param-types.ts";

const fieldPalette = (params: Record<string, unknown> = {}): ThemeRef => ({ theme: "field_palette", params });
const shape = (g: string): ThemeRef => ({ theme: "uniform", params: { glyph: g } });
/** The default themes the layer type rows name. */
const DEFAULTS: Readonly<Record<string, ThemeRef>> = {
  letter: { theme: "letter" },
  palette: fieldPalette(),
  alt: fieldPalette({ field: "alt" }),
  secondary: fieldPalette({ palette: { kind: "categorical", colors: "SecondaryStructure" } }),
  blues: fieldPalette({ palette: { kind: "continuous", colors: "Blues" } }),
  circle: shape("circle"),
  tick: shape("tick"),
};

/**
 * The numeric default color over a layer's domain (`matrix` and `color` over values; its `domain`, else the data's):
 * with no negative values, sequential `Blues` (light low) from 0 (from `lo` when the domain is given) to `hi`; with
 * negative values, diverging `RdBu` centered at 0. No domain (no data): `Blues` over the values.
 */
export function numericDefaultColor(lo: number, hi: number, given: boolean): ThemeRef {
  if (!(lo <= hi)) return DEFAULTS.blues;
  if (lo >= 0)
    return fieldPalette({ palette: { kind: "continuous", colors: "Blues", value_domain: [given ? lo : 0, hi] } });
  const m = Math.max(-lo, Math.abs(hi));
  return fieldPalette({ palette: { kind: "continuous", colors: "RdBu", value_domain: [-m, m] } });
}

/** The layer types' table, a row each: name | shows | data | layouts | axes | sizing | fallback | color | glyph. */
const LAYER_TYPES = `
letters | residues (element symbols for atoms), one alignment row's, or another polymer's through a mapping | none | tracks text | blocks columns | row | - | letter | -
rows | an alignment's rows in tree order, with the tree and data columns | none | tracks | columns | free | letters | letter | -
ranges | labeled ranges, discontinuous items, boundaries | ranges | tracks text | blocks columns | row | - | palette | tick
points | point features | ranges | tracks text | blocks columns | row | ranges | palette | circle
emphasis | positions emphasized in the letters or rows of its track | ranges | tracks text | blocks columns | row | points | - | -
variants | substitutions per position | ranges | tracks | blocks columns | row | points | alt | circle
secondary | helices, strands, turns | ranges | tracks text | blocks columns | row | ranges | secondary | -
values | a number per position or column | values | tracks | blocks columns | row | color | - | -
color | a color per position, from categories or values | ranges values | tracks text | blocks columns | row | ranges | palette | -
ruler | positions labeled in a numbering | none | tracks | blocks columns | free | omit | - | -
pairs | pairs of positions, by category | pairs | tracks | blocks columns | free | points | palette | circle
logo | letter frequencies per column | matrix | tracks | columns | free | values | letter | -
matrix | a value per position and state | matrix | tracks | blocks columns | free | values | blues | -
dot_plot | positions against positions on two axes: pairs as segments, a grid's values as cells | pairs grid | tracks | blocks | free | omit | palette | -`;

const letterMapField = OptionalField(
  nullable(letterMap),
  null,
  "the letter shown for each stored one: a registered map by name (complement, complement_rna, lower) or single letters mapped; null: as stored",
);
const density = OptionalField(
  literal("stacked", "majority"),
  "stacked",
  "a device column with more items than it draws one by one: each color's share stacked, or the most frequent; a hint",
);
const display = <V extends string>(d: V, ...values: V[]) => OptionalField(literal(...values), d, "a hint");
const split = OptionalField(
  nullable(
    object({
      by: RequiredField(str, "a field of the data's rows: category (a feature type), alt (a substitution)"),
      order: OptionalField(nullable(list(str)), null, "these rows first, every one; the others in data order"),
      empty: OptionalField(
        literal("show", "hide"),
        "show",
        "a value of order no row has: a row (a fixed grid) or none",
      ),
    }),
  ),
  null,
  "one labeled row per value of the field; the track's expanded opens them",
);
/** Each type's size channels: name and unit; none has a default, the viewer's size. */
const SIZE_CHANNELS: Readonly<Record<string, readonly (readonly [string, "px" | "band"])[]>> = {
  ranges: [
    ["thickness", "band"],
    ["radius", "px"],
  ],
  points: [["radius", "px"]],
  variants: [["thickness", "band"]],
  secondary: [["thickness", "band"]],
  values: [["thickness", "band"]],
  pairs: [
    ["width", "px"],
    ["height", "band"],
    ["radius", "px"],
  ],
  matrix: [["height", "band"]],
};
const sizeParam = OptionalField(nullable(themeRef("size")), null, "a size theme; null: the viewer's size");
/** Each type's params (top-level layer params; the list under the table). Typed as written: `KernelLayerParams`. */
const LAYER_PARAMS = {
  letters: {
    row: OptionalField(nullable(int), null, "an alignment row's letters; null: the block's (the reference row's)"),
    of: OptionalField(
      union(literal("mapped"), nullable(ref("sequence", "unit")), list(ref("sequence", "unit"))),
      null,
      "other polymers' letters where mapping paths place them: one, a list, or mapped (every codon mapping's protein)",
    ),
    residue_form: OptionalField(literal("auto", "one", "three"), "auto", "a protein's residues as M or Met; a hint"),
    letter_map: letterMapField,
  },
  rows: {
    order: OptionalField(
      orderField,
      null,
      "tree, data, a key table's field ({data, field, direction}) or a list of them, lexicographic; null: tree when there is one",
    ),
    filter: OptionalField(
      list(keyFilter),
      [],
      "rows kept when every filter keeps them: {data, field, values, min, max} over a key table",
    ),
    show_tree: OptionalField(bool, true, "drawn while the order is the tree"),
    row_subset: OptionalField(nullable(list(int)), null, "only these rows, applied before the filter"),
    data_columns: OptionalField(
      list(
        object({
          data: RequiredField(ref("annotation", "computed"), "a rows annotation"),
          field: RequiredField(str, "its field"),
          label: OptionalField(nullable(str), null, ""),
          type: OptionalField(str, "text", "number, fraction, text, category, <ns>:<name>; unknown: text"),
        }),
      ),
      [],
      "columns beside the rows",
    ),
    letter_map: letterMapField,
  },
  ranges: {
    thickness: sizeParam,
    radius: sizeParam,
    display: display("lanes", "lanes", "brackets", "underline", "union"),
    labels: OptionalField(
      bool,
      false,
      "each item's label inside its feature where it fits (on rings: outside the outermost ring, with a leader); a hint",
    ),
    direction_field: OptionalField(
      nullable(str),
      null,
      "a field of the rows giving each item's direction (+, 1, forward; -, -1, reverse): an arrowhead at its end",
    ),
    collapse: OptionalField(
      nullable(str),
      null,
      "a field of the rows: where the track's lanes do not fit, items sharing its value draw as one item in one lane (a gene's transcripts), until the track is expanded; streamed items; a hint",
    ),
    split,
    density,
  },
  points: { split, density, radius: sizeParam },
  variants: { split, density, thickness: sizeParam },
  secondary: { thickness: sizeParam },
  values: {
    thickness: sizeParam,
    encoding: display("bars", "bars", "line", "area"),
    domain: OptionalField(nullable(tuple(float, float)), null, "null: over all the layer's data; outside: clamped"),
  },
  ruler: {
    of: OptionalField(
      union(nullable(ref("sequence", "unit", "row")), list(ref("sequence", "unit", "row"))),
      null,
      "the coordinate system whose numbering it shows, placed as data is; a list: the first that maps a position",
    ),
    numbering: OptionalField(
      union(nullable(str), list(str)),
      null,
      "null: the section's or block's, or of's default; a list: the first that labels a position",
    ),
  },
  pairs: {
    width: sizeParam,
    height: sizeParam,
    radius: sizeParam,
    display: display("both", "arcs", "ends", "both"),
    side: OptionalField(
      nullable(literal("above", "below")),
      null,
      "null: toward the letters (the side facing the nearest letters track; above when there is none); a hint",
    ),
    arcs: display("always", "always", "hover"),
    labels: OptionalField(
      bool,
      false,
      "each arc's label, else its width theme's field's value, inside it under its apex where it fits; a hint",
    ),
  },
  logo: { scale: display("bits", "bits", "probability") },
  matrix: {
    height: sizeParam,
    density: OptionalField(
      literal("majority", "maximum"),
      "majority",
      "states thinner than a device pixel, per device row: majority (the most common color, faded by the share of states with one) or maximum (the state with the largest value, so a rare signal shows); a hint",
    ),
    domain: OptionalField(nullable(tuple(float, float)), null, "null: over all the layer's data; outside: clamped"),
    order: OptionalField(
      orderField,
      null,
      "the states' order: data, a key table's field ({data, field, direction}) over a table keyed by state, or a list of them; null: data",
    ),
    filter: OptionalField(
      list(keyFilter),
      [],
      "states kept when every filter keeps them: {data, field, values, min, max} over a table keyed by state",
    ),
  },
  dot_plot: {
    domain: OptionalField(nullable(tuple(float, float)), null, "null: over all the layer's data; outside: clamped"),
    aspect_ratio: OptionalField(
      nullable(float),
      1,
      "the plot's width ÷ height (1: square), its height following the data area's width within what the section shows, never flatter than a viewer's stretch of it; null: the track's own height; a hint",
    ),
  },
} as const satisfies Readonly<Record<string, Fields>>;

/** The fourteen kernel layer types, in the table's order (a test keeps it equal to the table and the builder's). */
export const KERNEL_LAYER_TYPE_NAMES = [
  "letters",
  "rows",
  "ranges",
  "points",
  "emphasis",
  "variants",
  "secondary",
  "values",
  "color",
  "ruler",
  "pairs",
  "logo",
  "matrix",
  "dot_plot",
] as const;
export type KernelLayerTypeName = (typeof KERNEL_LAYER_TYPE_NAMES)[number];
/** Each kernel type's own params (the layer's common ones are the tree table's), for the builder's typed methods. */
export type KernelLayerParams = {
  [T in KernelLayerTypeName]: T extends keyof typeof LAYER_PARAMS
    ? ValuesFor<SimpleParamsSchema<(typeof LAYER_PARAMS)[T]>>
    : Record<never, never>;
};

const rows = (table: string) =>
  table
    .trim()
    .split("\n")
    .map((line) => line.split(" | ").map((c) => c.trim()));
const orNull = <T>(c: string, value: T) => (c === "-" ? null : value);
const channel = (kind: ChannelKind, d: ThemeRef | null) => ({ name: kind, kind, unit: null, default: d });

export const KernelLayerTypes: readonly LayerTypeEntry[] = rows(LAYER_TYPES).map(
  ([name, description, accepts, layouts, axes, sizing, fallback, color, glyph]) => ({
    name,
    description,
    options: SimpleParamsSchema((LAYER_PARAMS as Readonly<Record<string, Fields>>)[name] ?? {}),
    accepts: accepts === "none" ? "none" : (accepts.split(" ") as Geometry[]),
    layouts: layouts.split(" ") as LayerTypeEntry["layouts"],
    axes: axes.split(" ") as LayerTypeEntry["axes"],
    sizing: sizing as LayerTypeEntry["sizing"],
    fallback: orNull(fallback, fallback),
    channels: [
      channel("color", orNull(color, DEFAULTS[color])),
      ...(glyph === "-" ? [] : [channel("glyph", DEFAULTS[glyph])]),
      ...(SIZE_CHANNELS[name] ?? []).map(([n, unit]) => ({ name: n, kind: "size" as const, unit, default: null })),
    ],
  }),
);

type P = Readonly<Record<string, unknown>>;
/** A `rules` theme over `rules`, nothing by default (the named tables). */
const ruled = (rules: unknown[]): ThemeRef => ({ theme: "rules", params: { rules } });
/** The named rule tables: each the `rules` theme it stands for, its params filling the table. */
const NAMED_TABLES: Readonly<Record<string, (p: P) => ThemeRef>> = {
  identity_to_row: (p) => ruled([{ when: [{ measure: { equals_row: p.row ?? null }, min: 1 }], theme: p.base }]),
  conservation_shading: (p) => {
    const measure = { column: p.data };
    return ruled([{ when: [{ measure, min: p.threshold ?? 0 }], theme: p.base, alpha: measure }]);
  },
  clustal_x: () => ruled(clustalXRules()),
  percent_identity: () => ruled(percentIdentityRules()),
};

/** Themes from rows (name | kinds | what it gives) and their params per kind. */
const themes = (table: string, params: Readonly<Record<string, Fields>>): ThemeEntry[] =>
  rows(table).map(([name, kinds, description]) => ({
    name,
    description,
    params: Object.fromEntries(
      kinds.split(" ").map((k) => [k, SimpleParamsSchema(params[`${name}/${k}`] ?? params[name] ?? {})]),
    ),
    fallback: null,
    expand: NAMED_TABLES[name],
    check: name === "scheme" ? checkSchemeParams : undefined,
  }));
/** An inline color scheme (seqspec:format#color-schemes), checked by `checkSchemeDef`. */
const schemeDef = { name: "SchemeDef", check: checkSchemeDef };
const categorical = () => OptionalField(palette, { kind: "categorical", colors: "Category10" }, "");
const base = OptionalField(themeRef("color"), DEFAULTS.letter, "the theme it shades");

const KernelThemeList: readonly ThemeEntry[] = themes(
  `
uniform | color glyph size | one value
rules | color glyph size | the first rule that holds
field_scale | size | the row's field (a number through a linear, sqrt or log scale; a string through a map) into a range
field_palette | color | the row's field through the palette
field_color | color | the field's values are Colors; others not colored
letter | color | the residue letter, or a ligand atom's element, through a letter dictionary
clustal_x | color | Clustal X by the column's composition (protein); a named rule table
percent_identity | color | the consensus letter where more than 80, 60, 40% of the column agree (Jalview's); a named rule table
scheme | color | a color scheme by name, or an inline one: the letter, field_palette or rules theme it stands for, its dark variant on a dark skin
by_unit | color | one color per unit, by unit ref in data order
by_copy_group | color | one color per copy group, by its first unit
identity_to_row | color | the base color where the letter equals that row's in the column; none elsewhere; a named rule table
conservation_shading | color | the base color faded by the column's value where it is at least threshold; none below; a named rule table
field_map | glyph | the field's value through a map
field_glyph | glyph | the field's values are shape names or one character`,
  {
    "uniform/color": { color: RequiredField(color, "the color") },
    "uniform/glyph": { glyph: RequiredField(glyph, "the glyph") },
    "uniform/size": { size: RequiredField(float, "the size, in the channel's unit") },
    field_scale: {
      field: OptionalField(nullable(str), null, "null: value"),
      domain: OptionalField(nullable(tuple(float, float)), null, "null: the layer's data's, or its domain_group's"),
      range: RequiredField(tuple(float, float), "the sizes the domain's ends map to, in the channel's unit"),
      scale: OptionalField(
        literal("linear", "sqrt", "log"),
        "linear",
        "log: of one plus the distance from the domain's start",
      ),
      map: OptionalField(nullable(dict(float)), null, "a string field: value -> size, in place of domain"),
      default: OptionalField(
        nullable(float),
        null,
        "a string field: for values the map lacks; null: the viewer's size",
      ),
    },
    field_palette: {
      field: OptionalField(
        nullable(str),
        null,
        "null: value for values and matrix; else category, or the item (id, label) when there is no category column",
      ),
      palette: categorical(),
    },
    field_color: { field: OptionalField(str, "color", "") },
    scheme: {
      name: OptionalField(nullable(str), null, "a registered color scheme; exactly one of name and custom"),
      custom: OptionalField(
        nullable(schemeDef),
        null,
        "an inline scheme: letters, sequential, diverging or categorical",
      ),
      overrides: OptionalField(dict(color), {}, "a letters scheme's letters recolored"),
      field: OptionalField(nullable(str), null, "values and categories: the field, as field_palette's"),
      domain: OptionalField(
        nullable(tuple(nullable(float), nullable(float))),
        null,
        "values: [min, max]; null: the layer's",
      ),
      reverse: OptionalField(bool, false, "values: the ramp reversed"),
    },
    letter: {
      colors: OptionalField(
        nullable(letterDictionary),
        null,
        "a dictionary name or {letter: Color}; null: the alphabet's",
      ),
    },
    by_unit: { palette: categorical() },
    by_copy_group: { palette: categorical() },
    identity_to_row: { row: OptionalField(nullable(int), null, "null: the reference row"), base },
    conservation_shading: {
      data: RequiredField(ref("annotation", "computed"), "values over the alignment's columns"),
      threshold: OptionalField(float, 0, ""),
      base,
    },
    field_map: {
      field: OptionalField(nullable(str), null, "null: category"),
      map: RequiredField(dict(glyph), "value -> glyph"),
      default: OptionalField(nullable(glyph), null, "for values the map lacks; null: none"),
    },
    field_glyph: { field: OptionalField(str, "glyph", "") },
  },
).map((t) =>
  t.name === "rules"
    ? {
        ...t,
        params: {
          color: rulesSchema("color", color),
          glyph: rulesSchema("glyph", glyph),
          size: rulesSchema("size", float),
        },
      }
    : t,
);

const ERRORS = `svs.format.invalid svs.version.missing svs.version.unsupported_major svs.node.kind svs.node.parent
  svs.node.count svs.params.missing svs.params.unknown svs.params.type svs.params.exclusive svs.ref.syntax svs.ref.missing
  svs.ref.duplicate svs.ref.unresolved svs.ref.wrong_kind svs.ref.cycle svs.ref.derived_unknown svs.custom.namespace
  svs.custom.invalid svs.requires.unknown svs.position.syntax svs.position.form svs.position.range svs.columns.length
  svs.columns.unknown svs.columns.type svs.color.invalid svs.glyph.invalid svs.source.uri svs.layer.data
  svs.layer.fallback svs.layer.annotation svs.layer.split svs.section.axis svs.section.positions svs.section.scroll_group
  svs.tree.leaves svs.patch.base svs.patch.ref svs.patch.anchor svs.patch.cycle svs.snapshot.ref
  svs.source.load svs.source.format svs.source.record svs.source.category svs.source.size_bound svs.sequence.alphabet
  svs.ligand.unsupported svs.numbering.length svs.numbering.unknown svs.alignment.ragged svs.alignment.shape
  svs.mapping.cigar svs.mapping.range svs.mapping.codon svs.mapping.ambiguous svs.focus.range
  svs.theme.channel svs.tree.streamed svs.numbering.invalid svs.section.follow
  svs.scheme.unknown svs.scheme.letter svs.query.invalid svs.query.version svs.query.search
  svs.query.input`;
const WARNINGS = `svs.version.newer_minor svs.snapshot.invalid svs.position.unknown svs.numbering.descending
  svs.mapping.empty svs.numbering.fallback svs.provider.unavailable svs.layer.unplaced svs.layer.insert svs.render.fallback svs.render.omitted
  svs.section.style svs.render.theme svs.render.glyph svs.render.preset svs.render.direction svs.alphabet.detected
  svs.alignment.parent svs.annotation.rows svs.unit.orphaned svs.unit.stale svs.mapping.orphaned svs.annotation.orphaned
  svs.provider.unmapped svs.provider.size svs.provider.lazy
  svs.provider.alphabet svs.coverage.ambiguous svs.layer.rows svs.stream.ref svs.stream.conflict svs.stream.depth
  svs.stream.fallback svs.query.stale svs.query.warning svs.layer.streamed_merged
  svs.stream.item_map svs.mapping.align svs.source.id svs.stream.load`;
/** Quiet notes on how valid input was read (seqspec:format#validation): nothing is wrong, nothing is lost. */
const INFOS = `svs.alignment.inserts svs.source.empty_record svs.source.skipped svs.source.unterminated svs.layer.empty
  svs.ref.pending svs.query.note svs.legend.unlisted`;
const codes = (s: string, severity: Severity) =>
  s
    .trim()
    .split(/\s+/)
    .map((code) => ({ code, severity }));

/** The kernel's tables: layer types, themes, issue codes (kernel-registry.ts adds the rest). */
export const KernelTables: SvsRegistryEntry = {
  layerTypes: KernelLayerTypes,
  themes: KernelThemeList,
  issueCodes: [...codes(ERRORS, "error"), ...codes(WARNINGS, "warning"), ...codes(INFOS, "info")],
};

const mvsField = (description: string) => OptionalField(nullable(str), null, description);
const structureRef = mvsField("the MVS ref of the structure node (the application writes it)");

/** The `mvs` extension (seqspec:format#extensions): links units to a MolViewSpec state; nothing in Seq* reads it. */
export const MvsRegistry: SvsRegistryEntry = {
  custom: [
    { namespace: "mvs", kind: "complex", schema: SimpleParamsSchema({ structure_ref: structureRef }) },
    {
      namespace: "mvs",
      kind: "unit",
      schema: SimpleParamsSchema({
        structure_ref: structureRef,
        label_asym_id: mvsField("MVS's label_asym_id"),
        auth_asym_id: mvsField("MVS's auth_asym_id"),
        instance_id: mvsField("MVS's instance_id: the operator"),
      }),
    },
  ],
};
