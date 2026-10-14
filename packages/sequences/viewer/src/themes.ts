/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Themes (seqspec:viewer#extension-points; seqspec:format#themes): one provider kind for every channel (`color`,
// `size`, `glyph`), the drawing half of the format's theme entries. A provider creates a theme per layer and channel in
// the graph's theme stage from a `ThemeContext`, precomputing tables, so `value` reads typed arrays and the location's
// numbers and allocates nothing. Nesting (a base) resolves through `tc.resolve`; a theme with no implementation, or not
// applicable here, falls back (the reference's, then its entry's; a color theme implemented but not applicable here:
// `letter`) and is reported; with no fallback left nothing is drawn. Built-ins:
// `uniform` and `rules` (rules.ts) for every kind; color `field_palette`, `field_color`, `letter`, `by_unit`,
// `by_copy_group`; glyph `field_map`, `field_glyph`. The named rule tables arrive expanded (VW4); a reference that
// names one is expanded here too. Every value goes through the format's palette code (`layerPalette`, `letterColors`).

import {
  type ChannelEntry,
  type ChannelKind,
  Color,
  canonical,
  layerPalette,
  letterColors,
  makePaletteFunction,
  numberColumn,
  type Palette,
  type PaletteFunction,
  paletteValue,
  type ResolvedLayer,
  type ResolvedSection,
  type SvsRegistries,
  type ThemeEntry,
  type ThemeRef,
  themeRef,
} from "@molstar/seq-core/format";
import {
  type Annotation,
  Annotations,
  CoordinateSystemKind,
  Issue,
  Molecules,
  type SeqSet,
  type StringTable,
} from "@molstar/seq-core/model";
import type { Axis } from "./axis.ts";
import { type GlyphRecord, GlyphShape, NO_COLOR, SHAPE_COUNT } from "./draw/backend.ts";
import type { SeqViewerRegistry } from "./registry.ts";

/** What a theme's value varies with: a renderer hoists `value` out of loops it does not vary in. */
export const Granularity = {
  uniform: 0,
  unit: 1,
  letter: 2,
  column: 3,
  row: 4,
  item: 5,
  value: 6,
  position: 7,
} as const;
const G = Granularity;
/**
 * A theme over its bases (a `rules` theme's rule themes and default): equal ones stay, `uniform` takes the other, a
 * letter under a column is the column (a column theme varies with the letter too); any two other kinds, `position`.
 */
export function coverGranularity(own: number, base: number): number {
  if (own === base || base === G.uniform) return own;
  if (own === G.uniform) return base;
  if ((own === G.column && base === G.letter) || (own === G.letter && base === G.column)) return G.column;
  return G.position;
}
/** This layer gives the element no glyph. */
export const NO_GLYPH = -1;
/**
 * A size theme's value is an integer, in 1/`SIZE_ONE` of its channel's unit (px or a band share), so no double crosses
 * a call on the draw's path; `NO_SIZE`: the viewer's own size.
 */
export const SIZE_ONE = 256;
export const NO_SIZE = -1;
/** A size in a channel's unit as a theme gives it; a negative or non-finite one is none. Edge-rate. */
export const sizeOf = (x: number) =>
  Number.isFinite(x) && x >= 0 ? Math.min(Math.round(x * SIZE_ONE), 0x3fffffff) : NO_SIZE;
/** A theme's alpha as the batch's byte: opaque. */
export const OPAQUE = 255;

/**
 * One per pass, reused; -1 (NaN for `value`): not applicable. An item's `annotation` and `row` (its data row) read its
 * fields; `letter` is an ASCII residue code, or a ligand atom's atomic number.
 */
export interface ThemeLocation {
  layer: number;
  block: number;
  unit: number;
  group: number;
  sequence: number;
  index: number;
  alignment: number;
  row: number;
  column: number;
  annotation: number;
  item: number;
  letter: number;
  value: number;
  /** A renderer-defined category index, read by its own extension's themes (the fingerprint's kind); -1: none. */
  category: number;
  state: number;
}
export const ThemeLocation = {
  create(): ThemeLocation {
    const l = { layer: 0, block: 0, unit: 0, group: 0, sequence: 0, index: 0, alignment: 0, row: 0, column: 0 };
    const loc = { ...l, annotation: 0, item: 0, letter: 0, value: 0, category: 0, state: 0 };
    ThemeLocation.reset(loc);
    return loc;
  },
  reset(loc: ThemeLocation) {
    loc.layer = loc.block = loc.unit = loc.group = loc.sequence = loc.index = -1;
    loc.alignment = loc.row = loc.column = loc.annotation = loc.item = loc.letter = loc.category = loc.state = -1;
    loc.value = Number.NaN;
  },
  /**
   * What every renderer fills for block b's index i: the layer, block, its unit (merged copies: set the member drawn
   * after), the index, and on an alignment block its alignment and column. Allocation-free.
   */
  at(loc: ThemeLocation, axis: Axis, layer: number, b: number, i: number): ThemeLocation {
    loc.layer = layer;
    loc.block = b;
    loc.unit = axis.unit[b];
    loc.index = i;
    loc.alignment = axis.alignment[b];
    loc.column = loc.alignment >= 0 ? i : -1;
    return loc;
  },
};

/** A size sample's unit in `LegendSink.size`. */
export const SizeUnit = { px: 0, band: 1 } as const;
/**
 * What a size sample is of (`LegendSink.size`): its channel's name (`width`, `radius`, `thickness`, `height`: how it is
 * drawn), and for a scale over a field the field and the sample's value of it (a color ramp over the same field in the
 * same layer colors the sample: the legend merges them).
 */
export interface LegendSizeOf {
  readonly channel: string;
  readonly field?: string;
  readonly at?: number;
}
/** What a legend writes into (seqspec:interaction#the-legend). */
export interface LegendSink {
  /** `glyph`: NO_GLYPH for a swatch; `color`: NO_COLOR for a glyph in the ink. */
  item(label: string, color: number, glyph: number): void;
  /**
   * A continuous palette over [lo, hi]: `stops` (t, 0xRRGGBB) pairs, t in 0..1; `caption` names what it measures;
   * `field` the field it is over (a size scale over it is merged with it).
   */
  ramp(stops: Float64Array, lo: number, hi: number, caption?: string, field?: string): void;
  /** A sample size (a width, a radius) drawn as its channel draws it; `unit` a `SizeUnit`. */
  size(label: string, value: number, unit: number, of?: LegendSizeOf): void;
}

export interface ThemeProvider<C extends ChannelKind = ChannelKind, P = unknown> {
  readonly name: string;
  readonly label: string;
  readonly category: string;
  /** The kinds it serves, as its format entry's `params`. */
  readonly kinds: readonly C[];
  /** The provider's default `Granularity`. */
  readonly granularity: number;
  /** False: `create` never reads residues or alignments (seqspec:viewer#updates-and-status). */
  readonly content?: boolean;
  /** False: draw its fallback (a measure this layer cannot compute: `share` on a lazy alignment). */
  isApplicable(tc: ThemeContext, kind: C, params: P): boolean;
  /** In the graph's theme stage; precomputes tables. */
  create(tc: ThemeContext, params: P, channel: ChannelEntry): Theme<C>;
}
export interface Theme<C extends ChannelKind = ChannelKind> {
  readonly kind: C;
  /** This instance's, after its params and bases. */
  readonly granularity: number;
  /**
   * How zoomed-out cells reduce it (seqspec:sources#levels-and-the-zoom-gate): the plurality class, or `palette(mean)`
   * of a value.
   */
  readonly reduction: "classes" | "mean";
  /** Color: 0xRRGGBB or NO_COLOR; size: in 1/`SIZE_ONE` of the channel's unit, NO_SIZE the default; glyph: a GlyphId or NO_GLYPH. */
  value(loc: ThemeLocation): number;
  /** Color only, as the batch's byte 0..255 (a double returned from a call is boxed): themes that fade. */
  readonly alpha?: ((loc: ThemeLocation) => number) | undefined;
  /** A `mean` theme over a row field: a data row's value (NaN: none), what `rowValues` reads for its pyramid. */
  readonly rowValue?: ((annotation: number, row: number) => number) | undefined;
  /**
   * A field theme's value for an item with no row in the set (a streamed region's, ST3b): its fields read from its
   * patch's `table` at `row`. Edge rate: a draw keeps the values per region and theme (`featureKey`).
   */
  readonly raw?: ((table: FieldTable, row: number) => number) | undefined;
  /** Its memo key with name, params and channel: unique per theme made in a viewer context. */
  readonly contextHash: number;
  readonly legend?: ((out: LegendSink) => void) | undefined;
}
export interface ThemeContext {
  readonly set: SeqSet;
  readonly section: ResolvedSection;
  readonly layer: ResolvedLayer;
  readonly annotations: readonly Annotation[];
  /** The layer's resolved domain; null: none. */
  readonly domain: readonly [number, number] | null;
  readonly dark: boolean;
  readonly strings: StringTable;
  /** The viewer context's glyph table (`DrawResources.glyph` reads it). */
  readonly glyphs: GlyphTable;
  /** A param typed ThemeRef (a base), for `channel`, created once here; null past 4 deep (reported) or none drawn. */
  resolve(ref: ThemeRef, channel: ChannelEntry): Theme | null;
  /** A new `contextHash`, counted per viewer context. */
  stamp(): number;
  /** The format's registries: field providers' names read as fields (seqspec:format#extensions). */
  readonly svs?: SvsRegistries;
  /** Takes a theme's notes on what it draws (a palette the legend cannot list), on the layer, at the theme stage. */
  readonly report?: (i: Issue) => void;
}

// --- glyphs -------------------------------------------------------------------------------------------------------

const SHAPES = Object.keys(GlyphShape);
/** Shapes are GlyphIds 0..7 in every table; characters and paths interned from SHAPE_COUNT on. */
export class GlyphTable {
  private readonly records: GlyphRecord[] = [];
  private readonly keys = new Map<string, number>();
  /** A document `Glyph` (a shape name or an object); NO_GLYPH when it is none of them. Edge-rate. */
  intern(g: unknown): number {
    if (typeof g === "string") return SHAPES.indexOf(g);
    if (g === null || typeof g !== "object") return NO_GLYPH;
    const o = g as { shape?: string; char?: string; path?: string; box?: [number, number]; fallback?: string };
    if (typeof o.shape === "string") return this.intern(o.shape);
    if (typeof o.char !== "string" && typeof o.path !== "string") return NO_GLYPH;
    // Its keys in any order: one glyph.
    const key = canonical(g);
    let id = this.keys.get(key);
    if (id !== undefined) return id;
    const fallback = Math.max(0, this.intern(o.fallback ?? "circle"));
    this.records.push({ char: o.char ?? null, path: o.path ?? null, box: o.box ?? [1, 1], fallback });
    id = SHAPE_COUNT + this.records.length - 1;
    this.keys.set(key, id);
    return id;
  }
  /** A shape's index; for a character or path, its fallback shape. */
  shape(id: number): number {
    return id < SHAPE_COUNT ? id : (this.records[id - SHAPE_COUNT]?.fallback ?? 0);
  }
  /** A character or path glyph's record; null for a shape. */
  glyph(id: number): GlyphRecord | null {
    return id < SHAPE_COUNT ? null : (this.records[id - SHAPE_COUNT] ?? null);
  }
}

// --- creation -----------------------------------------------------------------------------------------------------

const MAX_DEPTH = 4;
/** Per viewer context (its registry): the glyph table and the theme counter; never per module. */
const tables = new WeakMap<SeqViewerRegistry, { glyphs: GlyphTable; hashes: number }>();
export function glyphTableOf(registry: SeqViewerRegistry): GlyphTable {
  return tablesOf(registry).glyphs;
}
function tablesOf(registry: SeqViewerRegistry) {
  const held = tables.get(registry);
  if (held) return held;
  const t = { glyphs: new GlyphTable(), hashes: 0 };
  tables.set(registry, t);
  return t;
}

export interface ThemeInput {
  readonly set: SeqSet;
  readonly section: ResolvedSection;
  readonly layer: ResolvedLayer;
  readonly annotations: readonly Annotation[];
}
/** The layer's theme context; `report` takes the fallbacks and nesting issues (the theme stage's). Edge-rate. */
export function themeContext(
  input: ThemeInput,
  dark: boolean,
  registry: SeqViewerRegistry,
  report: (i: Issue) => void,
  depth = 0,
): ThemeContext {
  const l = input.layer;
  const t = tablesOf(registry);
  return {
    ...input,
    domain: l.domainMin <= l.domainMax ? [l.domainMin, l.domainMax] : null,
    dark,
    strings: input.set.strings,
    glyphs: t.glyphs,
    stamp: () => ++t.hashes,
    svs: registry.svs,
    report,
    resolve(ref, channel) {
      if (depth + 1 < MAX_DEPTH)
        return createTheme(themeContext(input, dark, registry, report, depth + 1), ref, channel, registry, report);
      const code = channel.kind === "glyph" ? "svs.glyph.invalid" : "svs.color.invalid";
      report(Issue.error(code, l.ref, `themes nest at most ${MAX_DEPTH} deep`));
      return null;
    },
  };
}

/**
 * The theme a reference names for a channel: a named table expanded, its provider when one serves the kind and
 * applies here, else down the fallbacks (the reference's, its entry's; a color theme implemented but not applicable
 * here, with neither: `letter`), each step reported (`svs.render.theme`); null when none is drawn.
 */
export function createTheme(
  tc: ThemeContext,
  given: ThemeRef,
  channel: ChannelEntry,
  registry: SeqViewerRegistry,
  report: (i: Issue) => void,
): Theme | null {
  const kind = channel.kind;
  let next: ThemeRef | null = given;
  for (let k = 0; next && k < 8; k++) {
    const norm = (r: ThemeRef) => (themeRef(kind).normalize?.(r, registry.svs) ?? r) as ThemeRef;
    let ref = norm(next);
    const entry: ThemeEntry | undefined = registry.svs.get("themes", ref.theme);
    if (entry?.expand) ref = norm({ ...entry.expand(ref.params ?? {}, kind), table: ref.theme });
    const p: ThemeProvider | undefined = registry.get("themes", ref.theme);
    const params = ref.params ?? {};
    if (p?.kinds.includes(kind) && p.isApplicable(tc, kind, params)) {
      if (k > 0) report(Issue.warning("svs.render.theme", tc.layer.ref, `"${given.theme}" is drawn as "${ref.theme}"`));
      return p.create(tc, params, channel);
    }
    const letter = p && kind === "color" && ref.theme !== "letter" ? { theme: "letter" } : null;
    next = ref.fallback ?? entry?.fallback ?? letter;
  }
  report(Issue.warning("svs.render.theme", tc.layer.ref, `"${given.theme}" has no implementation and no fallback`));
  return null;
}

/** A theme's value function, or its parts: the provider fills the rest, one shape for every theme. */
export type ThemeValueFn = (loc: ThemeLocation) => number;
export type ThemeParts = Pick<Theme, "value"> &
  Partial<Pick<Theme, "granularity" | "reduction" | "alpha" | "rowValue" | "legend" | "raw">>;
/** A patch's columns by name (a streamed region's items, `StreamItems.table`). */
export type FieldTable = Readonly<Record<string, ArrayLike<unknown>>>;
/** A streamed item's field as a palette reads it (`Theme.raw`): a number or a string; null: none. */
function cellOf(t: FieldTable, field: string, row: number): string | number | null {
  const c = t[field];
  const v = c && row < c.length ? c[row] : null;
  return v == null ? null : typeof v === "number" ? v : String(v);
}
export type ThemeParams = Readonly<Record<string, unknown>>;
export function themeProvider(
  name: string,
  kinds: readonly ChannelKind[],
  granularity: number,
  create: (tc: ThemeContext, p: ThemeParams, channel: ChannelEntry) => ThemeParts | ThemeValueFn,
  o: { content?: boolean; isApplicable?: (tc: ThemeContext, kind: ChannelKind, p: ThemeParams) => boolean } = {},
): ThemeProvider<ChannelKind, ThemeParams> {
  return {
    name,
    label: name,
    category: "kernel",
    kinds,
    granularity,
    content: o.content ?? false,
    isApplicable: o.isApplicable ?? (() => true),
    create(tc, p, channel) {
      const made = create(tc, p, channel);
      const t = typeof made === "function" ? { value: made } : made;
      return {
        kind: channel.kind,
        granularity: t.granularity ?? granularity,
        reduction: t.reduction ?? "classes",
        value: t.value,
        alpha: t.alpha,
        rowValue: t.rowValue,
        raw: t.raw,
        contextHash: tc.stamp(),
        legend: t.legend,
      };
    },
  };
}

/** The index of an annotation id in the layer's data (a few entries), -1: not its. */
export function dataIndex(tc: { readonly layer: ResolvedLayer }, annotation: number): number {
  const d = tc.layer.data;
  for (let k = 0; k < d.length; k++) if (d[k] === annotation) return k;
  return -1;
}
/** The block's alphabet: its alignment's or its unit's molecule's; null: a ligand's atoms. */
function blockAlphabet(set: SeqSet, section: ResolvedSection, b: number) {
  const block = section.blocks[b];
  if (!block) return null;
  if (block.alignment >= 0) return set.alignments[block.alignment]?.alphabet ?? null;
  return set.molecules[set.units.molecule[block.unit]]?.alphabet ?? null;
}
/** Per data annotation, a value per row; read at the location's annotation and data row (`other` with no row). */
export function perRow(tc: ThemeContext, of: (a: Annotation, r: number, k: number) => number): Int32Array[] {
  return tc.annotations.map((a, k) => Int32Array.from({ length: a?.count ?? 0 }, (_, r) => of(a, r, k)));
}
function rowReader(tc: ThemeContext, rows: Int32Array[], other: ThemeValueFn = () => -1): ThemeValueFn {
  return (loc) => {
    if (loc.annotation < 0 || loc.row < 0) return other(loc);
    const k = dataIndex(tc, loc.annotation);
    return k >= 0 && loc.row < rows[k].length ? rows[k][loc.row] : -1;
  };
}

// --- the built-ins ------------------------------------------------------------------------------------------------

const CATEGORY10: Palette = { kind: "categorical", colors: "Category10" };
const orNone = (c: number | undefined) => (c === undefined ? NO_COLOR : c);
const numeric = (a: Annotation | undefined) =>
  a?.geometry === "values" || a?.geometry === "matrix" || a?.geometry === "grid";
const RAMP = 16;

/**
 * An absolute palette's span: its first and last checkpoints (`[color, at]` pairs; a discrete one's finite range
 * ends); null for a normalized palette or a plain list.
 */
function absoluteSpan(palette: Palette): [number, number] | null {
  if (palette.kind === "categorical" || palette.mode !== "absolute" || typeof palette.colors === "string") return null;
  const at: number[] = [];
  for (const c of palette.colors ?? [])
    if (Array.isArray(c)) for (const x of c.slice(1)) if (typeof x === "number" && Number.isFinite(x)) at.push(x);
  return at.length >= 2 ? [Math.min(...at), Math.max(...at)] : null;
}

/** A categorical palette at a value (a row's, dense, a matrix's cells): the values' colors by a binary search. */
function byValue(tc: ThemeContext, f: PaletteFunction): ThemeValueFn {
  const seen = new Set<number>();
  for (const a of tc.annotations) {
    for (const x of (a && numberColumn(a, tc.set, "value")) ?? []) seen.add(x);
    for (const d of a?.dense ?? []) for (let i = 0; i < d.length; i++) seen.add(d[i]);
    for (const x of a?.matrix?.values ?? []) seen.add(x);
    for (const row of a?.grid?.rows ?? []) for (let i = 0; i < row.length; i++) seen.add(row[i]);
  }
  const keys = Float64Array.from([...seen].filter((x) => !Number.isNaN(x))).sort();
  const colors = Int32Array.from(keys, (x) => orNone(f(x)));
  return (loc) => {
    let lo = 0;
    for (let hi = keys.length - 1; lo < hi; ) {
      const mid = (lo + hi) >> 1;
      if (keys[mid] < loc.value) lo = mid + 1;
      else hi = mid;
    }
    return lo < keys.length && keys[lo] === loc.value ? colors[lo] : NO_COLOR;
  };
}
/** A continuous palette at the location's value (`palette(mean)`): through `at`, so no double crosses a call. */
function atValue(f: PaletteFunction): ThemeValueFn {
  const at = f.at;
  if (at) return (loc) => orNone(at(loc));
  return (loc) => (Number.isNaN(loc.value) ? -1 : orNone(f(loc.value)));
}
/**
 * A palette's ramp over its domain. Continuous: RAMP samples. Discrete: its bands as hard steps, a pair of stops at
 * each edge (the color before and after, at one `t`; the edge found to 1/1024 of the domain). Edge-rate.
 */
function ramp(tc: ThemeContext, f: PaletteFunction, palette: Palette, out: LegendSink, field: string) {
  // An absolute palette ignores `value_domain`: the layer's domain, else its checkpoints' span.
  const span = absoluteSpan(palette);
  const given = palette.kind === "categorical" || span ? null : palette.value_domain;
  const lo = f.domain?.[0] ?? given?.[0] ?? tc.domain?.[0] ?? span?.[0] ?? 0;
  const hi = f.domain?.[1] ?? given?.[1] ?? tc.domain?.[1] ?? span?.[1] ?? 1;
  const at = (t: number) => orNone(f(lo + t * (hi - lo)));
  if (palette.kind === "discrete") {
    const n = 1024;
    const s: number[] = [0, at(0)];
    for (let k = 1, prev = at(0); k <= n; k++) {
      const c = at(k / n);
      if (c === prev) continue;
      s.push((k - 0.5) / n, prev, (k - 0.5) / n, c);
      prev = c;
    }
    s.push(1, at(1));
    out.ramp(Float64Array.from(s), lo, hi, undefined, field);
    return;
  }
  const stops = new Float64Array(2 * RAMP);
  for (let k = 0; k < RAMP; k++) {
    stops[2 * k] = k / (RAMP - 1);
    stops[2 * k + 1] = at(stops[2 * k]);
  }
  out.ramp(stops, lo, hi, undefined, field);
}

/** The row's field through the layer's one palette; dense values by value; a ramp or keys. */
const fieldPalette = themeProvider("field_palette", ["color"], G.item, (tc, p) => {
  const field = typeof p.field === "string" ? p.field : null;
  const palette = (p.palette ?? CATEGORY10) as Palette;
  const cat = palette.kind === "categorical";
  const f = layerPalette(tc.set, tc.annotations, palette, field, tc.domain ?? [Number.NaN, Number.NaN], tc.svs);
  // Streamed values with no vocabulary are colored by their hash (seqspec:format#palettes), and an item that sends no
  // value of the field has the missing color: one note per layer, true of both.
  if (f.unlisted) {
    const what = field ? `"${field}"` : "its items' category, id or label";
    const message =
      `the stream lists no values of ${what}: an item's value is colored by its hash, an item without one has the ` +
      "missing color, and the legend cannot list them; list the values in the manifest's vocabulary for a legend";
    tc.report?.(Issue.info("svs.legend.unlisted", tc.layer.ref, message));
  }
  const columns = tc.annotations.map((a) => (a && !cat ? numberColumn(a, tc.set, field ?? "value", tc.svs) : null));
  const rows = perRow(tc, (a, r, k) => {
    const column = columns[k];
    return orNone(f(column ? column[r] : paletteValue(a, tc.set, r, field, tc.svs)));
  });
  const value = field === "value" || (field === null && numeric(tc.annotations[0]));
  return {
    value: rowReader(tc, rows, cat ? byValue(tc, f) : atValue(f)),
    // A streamed item's own value as a held row's (`paletteValue`): its category, else its id, else its label.
    raw: (t, r) => orNone(f(cellOf(t, field ?? (t.category ? "category" : t.id ? "id" : "label"), r))),
    granularity: value ? G.value : G.item,
    reduction: cat ? "classes" : "mean",
    rowValue: cat
      ? undefined
      : (annotation, row) => {
          const c = columns[dataIndex(tc, annotation)];
          return c && row >= 0 && row < c.length ? c[row] : Number.NaN;
        },
    legend: cat
      ? (out) => {
          for (const key of f.keys ?? []) out.item(key, orNone(f(key)), NO_GLYPH);
        }
      : (out) => ramp(tc, f, palette, out, field ?? "value"),
  };
});

const fieldColor = themeProvider("field_color", ["color"], G.item, (tc, p) => {
  const field = typeof p.field === "string" ? p.field : "color";
  return {
    value: rowReader(
      tc,
      perRow(tc, (a, r) => orNone(Color.decodeFirst(Annotations.rowValue(a, tc.set, r, field)))),
    ),
    raw: (t, r) => orNone(Color.decodeFirst(cellOf(t, field, r))),
  };
});

/** Every kind: one value. */
const uniform = themeProvider("uniform", ["color", "glyph", "size"], G.uniform, (tc, p, { kind }) => {
  const v =
    kind === "color"
      ? Color.decodePrimary(p.color ?? "black")
      : kind === "glyph"
        ? tc.glyphs.intern(p.glyph)
        : sizeOf(Number(p.size));
  return () => v;
});

/** A number's place in [0, 1] through `scale` (linear, `sqrt`, `log`) over [lo, hi], clamped. */
function scaleAt(x: number, lo: number, hi: number, scale: unknown): number {
  const d = hi - lo;
  const u = d > 0 ? Math.min(1, Math.max(0, (x - lo) / d)) : 0;
  return scale === "sqrt" ? Math.sqrt(u) : scale === "log" ? (d > 0 ? Math.log1p(u * d) / Math.log1p(d) : 0) : u;
}
/**
 * A number through `scale` over [lo, hi] (clamped) into `range`, as a size (`sizeOf`); NaN: none. Per row at the theme
 * stage; per element on the by-value path (no allocation).
 */
function scaled(x: number, lo: number, hi: number, scale: unknown, range: readonly [number, number]): number {
  if (Number.isNaN(x)) return NO_SIZE;
  return sizeOf(range[0] + scaleAt(x, lo, hi, scale) * (range[1] - range[0]));
}
/**
 * `field_scale`: a row's number field through a linear, sqrt or log scale over the domain (given, else the layer's for
 * the value, else the field's over the layer's data) into `range`, clamped; a string field through `map` (`default` for
 * a value it lacks); a dense position's value (no row) by the value.
 */
const fieldScale = themeProvider("field_scale", ["size"], G.item, (tc, p, channel) => {
  const field = typeof p.field === "string" ? p.field : null;
  const name = field ?? "value";
  const range = p.range as [number, number];
  const map = (p.map ?? null) as Record<string, number> | null;
  const columns = tc.annotations.map((a) => (a && !map ? numberColumn(a, tc.set, name, tc.svs) : null));
  let [lo, hi] = Array.isArray(p.domain) ? (p.domain as [number, number]) : [Number.NaN, Number.NaN];
  if (!(lo <= hi) && name === "value" && tc.domain) [lo, hi] = tc.domain;
  if (!(lo <= hi)) {
    [lo, hi] = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
    for (const c of columns)
      for (let r = 0; c && r < c.length; r++) {
        lo = Math.min(lo, c[r]);
        hi = Math.max(hi, c[r]);
      }
  }
  const none = p.default == null ? NO_SIZE : sizeOf(Number(p.default));
  const rows = perRow(tc, (a, r, k) => {
    if (map) {
      const v = paletteValue(a, tc.set, r, name, tc.svs);
      return v == null ? none : Object.hasOwn(map, String(v)) ? sizeOf(map[String(v)]) : none;
    }
    const c = columns[k];
    return scaled(c ? c[r] : Number(paletteValue(a, tc.set, r, name, tc.svs)), lo, hi, p.scale, range);
  });
  return {
    value: rowReader(tc, rows, (loc) => (map ? none : scaled(loc.value, lo, hi, p.scale, range))),
    legend: (out) => scaleLegend(out, channel, name, map, lo, hi, p.scale, range),
  };
});
/** A number as a legend shows it (a sample, a ramp's end): integers whole, others to three significant digits. */
export const legendNumber = (n: number) => (Number.isInteger(n) ? String(n) : String(Number(n.toPrecision(3))));
/** Most sizes a map lists in the legend. */
const MAP_SAMPLES = 8;
/**
 * `field_scale`'s legend: a number field's domain ends and its middle, each with its size through the scale; a string
 * field's map, entry by entry. Edge-rate.
 */
function scaleLegend(
  out: LegendSink,
  channel: ChannelEntry,
  field: string,
  map: Record<string, number> | null,
  lo: number,
  hi: number,
  scale: unknown,
  range: readonly [number, number],
) {
  const unit = channel.unit === "band" ? SizeUnit.band : SizeUnit.px;
  if (map) {
    for (const [key, size] of Object.entries(map).slice(0, MAP_SAMPLES))
      out.size(key, size, unit, { channel: channel.name, field });
    return;
  }
  if (!(lo <= hi) || !Number.isFinite(lo) || !Number.isFinite(hi)) return;
  const d = hi - lo;
  for (const x of d > 0 ? [lo, lo + d / 2, hi] : [lo]) {
    const size = range[0] + scaleAt(x, lo, hi, scale) * (range[1] - range[0]);
    out.size(legendNumber(x), size, unit, { channel: channel.name, field, at: x });
  }
}

/**
 * A letter table: ASCII code (either case) → color; for atoms (alphabet null), atomic number → color, the symbols from
 * the model's one element table matched case-blind, whether the dictionary is the default or given.
 */
function letterTable(colors: unknown, alphabet: string | null): Int32Array {
  const out = new Int32Array(128).fill(NO_COLOR);
  const dict = letterColors(colors as never, alphabet as never);
  if (!dict) return out;
  if (alphabet === null) {
    const upper = new Map(Object.entries(dict).map(([k, c]) => [k.toUpperCase(), c]));
    for (let z = 1; z < 128; z++) out[z] = upper.get(Molecules.elementSymbol(z).toUpperCase()) ?? NO_COLOR;
    return out;
  }
  for (const [key, c] of Object.entries(dict)) {
    if (key.length !== 1) continue;
    const code = key.toUpperCase().charCodeAt(0);
    if (code < 128) out[code] = out[code | 32] = c;
  }
  return out;
}
/** Letters `of` another polymer are in its alphabet (a protein's over the DNA it is coded in); undefined: no `of`. */
export function ofAlphabet(tc: ThemeContext): string | null | undefined {
  const l = tc.layer;
  if (l.params?.of == null) return undefined;
  const from = l.placements.length > 0 && l.placements[0].path >= 0 ? l.paths[l.placements[0].path].from : null;
  if (!from) return undefined;
  const K = CoordinateSystemKind;
  const m = from.kind === K.unit ? tc.set.units.molecule[from.id] : from.kind === K.molecule ? from.id : -1;
  return tc.set.molecules[m]?.alphabet;
}
const letter = themeProvider("letter", ["color"], G.letter, (tc, p) => {
  const of = ofAlphabet(tc);
  const made = new Map<string | null, Int32Array>();
  const table = (alphabet: string | null) => {
    const t = made.get(alphabet) ?? letterTable(p.colors ?? null, alphabet);
    made.set(alphabet, t);
    return t;
  };
  const tables = tc.section.blocks.map((_, b) => table(of !== undefined ? of : blockAlphabet(tc.set, tc.section, b)));
  const one = tables[0] ?? table(null);
  return (loc) => {
    const t = loc.block >= 0 && loc.block < tables.length ? tables[loc.block] : one;
    return loc.letter >= 0 && loc.letter < 128 ? t[loc.letter] : NO_COLOR;
  };
});

/** One color per unit (by unit ref, data order) or per copy group (by its first unit), the same in every block. */
function byUnit(name: string, group: boolean) {
  return themeProvider(name, ["color"], G.unit, (tc, p) => {
    const u = tc.set.units;
    const key = (id: number) => tc.set.strings.get(u.key[id]);
    const first = (id: number) => (group ? tc.set.groups.first[u.group[id]] : id);
    const keys: string[] = [];
    for (let id = 0; id < u.count; id++) if (u.alive[id] && first(id) === id) keys.push(key(id));
    const f = makePaletteFunction((p.palette ?? CATEGORY10) as Palette, keys);
    const table = Int32Array.from({ length: u.count }, (_, id) => (u.alive[id] ? orNone(f(key(first(id)))) : NO_COLOR));
    return (loc) => (loc.unit >= 0 && loc.unit < table.length ? table[loc.unit] : NO_COLOR);
  });
}

/** Per row, the field's value through `pick` into a GlyphId; a streamed item's from its region's table (ST4b). */
function glyphRows(tc: ThemeContext, field: string, pick: (value: unknown) => unknown): ThemeParts {
  return {
    value: rowReader(
      tc,
      perRow(tc, (a, r) => tc.glyphs.intern(pick(Annotations.rowValue(a, tc.set, r, field)))),
    ),
    raw: (t, r) => tc.glyphs.intern(pick(cellOf(t, field, r))),
  };
}
const fieldMap = themeProvider("field_map", ["glyph"], G.item, (tc, p) => {
  const map = (p.map ?? {}) as ThemeParams;
  const fallback = p.default ?? null;
  return glyphRows(tc, typeof p.field === "string" ? p.field : "category", (v) =>
    v != null && Object.hasOwn(map, String(v)) ? map[String(v)] : fallback,
  );
});
/** One grapheme (as the format checks `char`): a plain character without the segmenter, else segmented (edge rate). */
const graphemes = new Intl.Segmenter(undefined, { granularity: "grapheme" });
const oneGrapheme = (v: string) => v.length === 1 || (v.length > 1 && [...graphemes.segment(v)].length === 1);
const fieldGlyph = themeProvider("field_glyph", ["glyph"], G.item, (tc, p) =>
  glyphRows(tc, typeof p.field === "string" ? p.field : "glyph", (v) =>
    typeof v === "string" && !SHAPES.includes(v) && oneGrapheme(v) ? { char: v } : v,
  ),
);

/** The built-in themes but `rules` (rules.ts's `DefaultThemes` adds it). */
export const KernelThemes: readonly ThemeProvider[] = [
  uniform,
  fieldScale,
  fieldPalette,
  fieldColor,
  letter,
  byUnit("by_unit", false),
  byUnit("by_copy_group", true),
  fieldMap,
  fieldGlyph,
] as ThemeProvider[];

// --- ink on a cell ------------------------------------------------------------------------------------------------

/** The ink that reads on a cell: dark on light fills, light on dark ones, by luminance. */
export const CellInk = { dark: 0x181a1e, light: 0xffffff } as const;
const inkKey = new Int32Array(1024).fill(-1);
const inkOf = new Int32Array(1024);
const linear = (c: number) => (c <= 10 ? c / 3294.6 : ((c / 255 + 0.055) / 1.055) ** 2.4);
/**
 * `CellInk.dark` or `.light` for a fill (0xRRGGBB), cached per color in a direct-mapped table: multiplicative hashing,
 * since a palette's colors often share bytes (an XOR of them kept two letters' colors evicting each other per frame).
 */
export function cellInk(fill: number): number {
  const h = Math.imul(fill, 0x9e3779b1) >>> 22;
  if (inkKey[h] === fill) return inkOf[h];
  const l = 0.2126 * linear((fill >> 16) & 255) + 0.7152 * linear((fill >> 8) & 255) + 0.0722 * linear(fill & 255);
  // Contrast against the dark ink (luminance about 0.01) or white, whichever is higher.
  const ink = (l + 0.05) / 0.06 >= 1.05 / (l + 0.05) ? CellInk.dark : CellInk.light;
  inkKey[h] = fill;
  inkOf[h] = ink;
  return ink;
}
