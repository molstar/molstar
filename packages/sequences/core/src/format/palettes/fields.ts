/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a `field_palette` theme reads per row (seqspec:format#themes), decided once for every viewer, and the layer's
// one palette over it (`layerPalette`): its `field`, else the default, `value` for values and matrix data, `category`
// for the others; and where the data has no `category` column at all, the item's name (its `id`, else its `label`), so
// a track is never left uncolored by default. An item with neither has nothing to key on: every such item takes one
// value (`UNNAMED_ITEM`), so one color, and legends leave it out, since it names nothing. A streamed value the palette's
// keys lack (no vocabulary of the field) takes a list color by its hash (`paletteHash`), which no legend lists.

import { type Annotation, Annotations, ColumnRoles, type SeqSet } from "@molstar/seq-core/model";
import type { FieldProviderEntry, SvsRegistries } from "../registry.ts";
import {
  categoricalList,
  makePaletteFunction,
  type Palette,
  type PaletteFunction,
  type PaletteValue,
  paletteHash,
} from "./palette.ts";

/** The default's value for an item with no category, id or label: one for all of them; never listed in a legend. */
export const UNNAMED_ITEM = "\u0000unnamed";

const provided = new WeakMap<Annotation, Map<FieldProviderEntry, Float64Array | Int32Array>>();
/**
 * A field provider's values on an annotation (seqspec:format#extensions), by its name: read once per annotation and
 * provider, then kept; numbers, or string codes into the set's strings. Null: none registered under the name, or it
 * does not apply.
 */
function providedField(reg: SvsRegistries | undefined, a: Annotation, set: SeqSet, name: string) {
  const e = reg && name.includes(":") ? reg.get("fields", name) : undefined;
  if (!e?.applies(a)) return null;
  let held = provided.get(a);
  if (!held) {
    held = new Map();
    provided.set(a, held);
  }
  let out = held.get(e);
  if (!out) {
    out = e.type === "number" ? new Float64Array(a.count).fill(Number.NaN) : new Int32Array(a.count).fill(-1);
    e.read(a, out, set.strings);
    held.set(e, out);
  }
  return out;
}

/** The value a row takes through a `field_palette` theme with this `field` (null: the default); `reg`: its providers. */
export function paletteValue(
  ann: Annotation,
  set: SeqSet,
  row: number,
  field: string | null,
  reg?: SvsRegistries,
): PaletteValue {
  const value = Annotations.rowValue as (a: Annotation, s: SeqSet, r: number, f: string) => PaletteValue;
  const p = field !== null ? providedField(reg, ann, set, field) : null;
  if (p) return p instanceof Float64Array ? p[row] : p[row] >= 0 ? set.strings.get(p[row]) : null;
  if (field !== null) return value(ann, set, row, field);
  if (ann.geometry === "values" || ann.geometry === "matrix") return value(ann, set, row, "value");
  if (ann.category !== null) return value(ann, set, row, "category");
  return value(ann, set, row, "id") ?? value(ann, set, row, "label") ?? UNNAMED_ITEM;
}

/** A number column by name, read as `rowValue` reads it (`value`, else a number field, else a provider's); null: none. */
export function numberColumn(a: Annotation, set: SeqSet, name: string, reg?: SvsRegistries): Float64Array | null {
  if (name === "value" && ColumnRoles[a.geometry].value !== undefined) return a.sparse;
  const id = set.strings.find(name);
  for (const f of a.fields) if (f.name === id) return f.kind === "number" ? f.data : null;
  const p = providedField(reg, a, set, name);
  return p instanceof Float64Array ? p : null;
}

type Vocab = { manifest?: { vocabulary?: Record<string, string[]> } } | undefined;
const numeric = (a: Annotation) => a.geometry === "values" || a.geometry === "matrix" || a.geometry === "grid";
const num = (v: unknown) => (typeof v === "number" ? v : typeof v === "string" ? Number.parseFloat(v) : Number.NaN);

/**
 * A layer's one palette: first occurrence and a normalized palette's domain over all rows of all its data, whatever its
 * filter keeps, annotations in `data` order; a `null` end of `value_domain` from the layer's resolved `domain` when the
 * colored field is the value. A categorical one carries its legend `keys`, the unnamed items' value left out. Both
 * viewers color through it, once per layer.
 */
export function layerPalette(
  set: SeqSet,
  data: readonly Annotation[],
  palette: Palette,
  field: string | null,
  domain: readonly [number, number],
  reg?: SvsRegistries,
): PaletteFunction {
  const value = field === "value" || (field === null && data.length > 0 && numeric(data[0]));
  const dense = field === null || field === "value";
  if (palette.kind !== "categorical") {
    const given = palette.value_domain ?? [null, null];
    if (palette.mode === "absolute" || (given[0] !== null && given[1] !== null))
      return makePaletteFunction(palette, []);
    let [lo, hi] = value ? domain : [Number.NaN, Number.NaN];
    if (!(lo <= hi)) {
      [lo, hi] = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
      const see = (x: number) => {
        if (x < lo) lo = x;
        if (x > hi) hi = x;
      };
      for (const a of data) {
        const column = numberColumn(a, set, field ?? "value", reg);
        for (let r = 0; r < a.count; r++) see(column ? column[r] : num(paletteValue(a, set, r, field, reg)));
        for (const d of dense ? (a.dense ?? []) : []) for (let i = 0; i < d.length; i++) see(d[i]);
        if (dense && a.matrix) for (let i = 0; i < a.matrix.values.length; i++) see(a.matrix.values[i]);
        for (const row of dense ? (a.grid?.rows ?? []) : []) for (let i = 0; i < row.length; i++) see(row[i]);
      }
    }
    const ends: [number | null, number | null] = lo <= hi ? [given[0] ?? lo, given[1] ?? hi] : [given[0], given[1]];
    return Object.assign(makePaletteFunction({ ...palette, value_domain: ends }, []), { domain: ends });
  }
  const values: PaletteValue[] = [];
  let [streamed, listed] = [false, false];
  for (const a of data) {
    for (let r = 0; r < a.count; r++) values.push(paletteValue(a, set, r, field, reg));
    // A streamed annotation has no rows: its manifest's vocabulary for the field gives the order
    // (seqspec:sources#streamed-sources; `svs:stream`, load/compile.ts' `StreamProp`, read by name: no import of the
    // compiler here). With no field, the one `paletteValue` reads: the category, else the id, else the label.
    const vocab = a.streamed ? (a.custom.props.get("svs:stream") as Vocab)?.manifest?.vocabulary : undefined;
    const given = vocab ? (field ? vocab[field] : (vocab.category ?? vocab.id ?? vocab.label)) : undefined;
    if (given) values.push(...given);
    streamed ||= !!a.streamed;
    listed ||= !!given?.length;
    // Values with no rows (dense, seqspec:model#annotations; a matrix's cells): each distinct one once, in order.
    const seen = new Set<number>();
    const add = (x: number) => {
      if (Number.isNaN(x) || seen.has(x)) return;
      seen.add(x);
      values.push(x);
    };
    for (const d of dense ? (a.dense ?? []) : []) for (let i = 0; i < d.length; i++) add(d[i]);
    if (dense && a.matrix) for (let i = 0; i < a.matrix.values.length; i++) add(a.matrix.values[i]);
    for (const row of dense ? (a.grid?.rows ?? []) : []) for (let i = 0; i < row.length; i++) add(row[i]);
  }
  const f = makePaletteFunction(palette, values);
  const fold = palette.case_insensitive === true;
  const unnamed = fold ? UNNAMED_ITEM.toUpperCase() : UNNAMED_ITEM;
  const keys = f.keys?.filter((k) => k !== unnamed) ?? [];
  // A streamed value outside the keys (no vocabulary, or one it lacks): a list palette's color at its hash
  // (seqspec:format#palettes), the same in every region and reload. Read at a region's first draw (`Theme.raw`), whose
  // colors the draw keeps: no cache here, so a stream of unique ids holds nothing.
  const list = streamed ? categoricalList(palette) : null;
  if (!list?.length) return Object.assign((v: PaletteValue) => f(v), { keys });
  const known = new Set(f.keys);
  const n = list.length;
  const hashed = (v: PaletteValue) => {
    if (v === null || v === undefined || v === "." || v === "?") return f(v);
    const s = fold ? `${v}`.toUpperCase() : `${v}`;
    return known.has(s) ? f(v) : list[paletteHash(s) % n];
  };
  return Object.assign(hashed, { keys, unlisted: listed ? undefined : true });
}
