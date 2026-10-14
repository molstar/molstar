/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Keyed sets (seqspec:format#layer-types, seqspec:format#keyed-sets): the members of a set (an alignment's rows, a
// matrix's states) drawn in an order and kept by a filter, both over key tables: `rows` annotations that give members
// fields (a document's, a provider's `row_stats` or `redundancy`, a field provider's). `keyed` resolves one set to a
// permutation of its members and a kept mask, rebuilt whenever it is called (the viewer's graph calls it when a key
// table changes). Resolve-time code: it allocates per call, never per frame. What a field means is the table's, never
// read here.

import { type Annotation, ColumnRoles, KeyKind, type KeyTarget, type SeqSet } from "@molstar/seq-core/model";
import { paletteValue } from "../palettes/fields.js";
import type { SvsRegistries } from "../registry.js";

/** One key of an order, as the document writes it (defaults filled by normalizing). */
export type OrderKeyValue = "data" | "tree" | { data: string; field: string; direction?: "asc" | "desc" };
/** One filter, as the document writes it. */
export interface KeyFilterValue {
  data: string;
  field: string;
  values?: readonly (string | number | boolean)[] | null;
  min?: number | null;
  max?: number | null;
}

/** What a set's members are: how many, and which member a key table's row names (-1: none). */
export interface Members {
  readonly count: number;
  /** The member annotation row `i` of `ann` stands for. */
  of(ann: Annotation, i: number): number;
  /** Members in tree order (an alignment's leaf order); null: the set has no tree. */
  readonly tree: ArrayLike<number> | null;
}

/** The result: all members in order (`permutation`), and each member's `kept` (1) or not (0). */
export interface Keyed {
  readonly permutation: Int32Array;
  readonly kept: Uint8Array;
}

/** A report of a problem in the spec: the issue code is `svs.layer.data`; `where` is `order/0` or `filter/1`. */
export type KeyedSay = (where: string, message: string) => void;

/** Whether a rows layer's rows follow the tree: its order starts with `tree`, or is `null` and the alignment has one. */
export const followsTree = (order: unknown, hasTree: boolean): boolean =>
  order == null ? hasTree : (Array.isArray(order) ? order[0] : order) === "tree";

/** The members of an alignment's key tables: a table row stands for the row in its `target.row`, of this alignment. */
export function alignmentRows(alignment: number, count: number, tree: ArrayLike<number> | null): Members {
  return { count, tree, of: (ann, i) => (ann.target.alignment[i] === alignment ? ann.target.row[i] : -1) };
}

const kind: KeyTarget = { kind: -1, id: -1, row: -1 };
const MISSING = 2;
const type = (x: unknown) =>
  x === undefined || x === null || (typeof x === "number" && Number.isNaN(x)) ? MISSING : typeof x === "string" ? 1 : 0;
const num = (x: unknown) => (typeof x === "number" ? x : x === true ? 1 : 0);

/** One key table's field read for every member: its value, or undefined where the table has no row for it. */
function column(
  set: SeqSet,
  reg: SvsRegistries,
  members: Members,
  ref: string,
  field: string,
  say: (message: string) => void,
): unknown[] | null {
  const ok = set.byName(ref, kind) && kind.kind === KeyKind.annotation;
  const ann = ok ? set.annotations[kind.id] : undefined;
  if (!ann) {
    say(`"${ref}" names no annotation or computed node here: give a rows annotation (a key table)`);
    return null;
  }
  if (ann.geometry !== "rows") {
    say(`"${ref}" is a ${ann.geometry} annotation: a key table is a rows annotation, keyed by row`);
    return null;
  }
  const out = new Array<unknown>(members.count).fill(undefined);
  let named = 0;
  let found = 0;
  for (let i = 0; i < ann.count; i++) {
    const m = members.of(ann, i);
    if (m < 0 || m >= members.count) continue;
    named++;
    const x = paletteValue(ann, set, i, field, reg);
    if (x !== undefined && x !== null) found++;
    out[m] = x ?? undefined;
  }
  if (named === 0) {
    say(`"${ref}" has no row for any member of this set: its rows must name this set's members`);
    return null;
  }
  if (found === 0) {
    const fields = [
      ...Object.keys(ColumnRoles.rows).filter((k) => k !== "row" && k !== "state" && k !== "node"),
      ...ann.fields.map((f) => set.strings.get(f.name)),
    ];
    // A field the table has but gives no member a value is a field with nothing: the filter keeps none, as for any missing.
    if (fields.includes(field) || reg.get("fields", field)) {
      say(`"${ref}" has no values for "${field}" (every member's is missing)`);
      return out;
    }
    say(`"${ref}" has no field "${field}"${fields.length ? `; its fields: ${[...new Set(fields)].join(", ")}` : ""}`);
    return null;
  }
  return out;
}

/** A comparison of two members' values on one key: numbers (booleans as 0 and 1) before strings, missing last. */
function compare(a: unknown, b: unknown, sign: number): number {
  const [ka, kb] = [type(a), type(b)];
  if (ka === MISSING || kb === MISSING) return ka === kb ? 0 : ka === MISSING ? 1 : -1;
  if (ka !== kb) return ka - kb;
  const c = ka === 0 ? num(a) - num(b) : (a as string) < (b as string) ? -1 : (a as string) > (b as string) ? 1 : 0;
  return sign * (c < 0 ? -1 : c > 0 ? 1 : 0);
}

const holds = (x: unknown, values: readonly (string | number | boolean)[]) =>
  values.some((v) => String(v) === String(x));

/**
 * The order and kept mask of `members`. `subset` (`row_subset`) leaves the others out first; a key table that names
 * nothing or lacks its field is reported through `say` and the key or filter is skipped. `order` `null` is the tree's
 * order when the set has one, else data order; later keys break ties, and the last tie is data order.
 */
export function keyed(
  set: SeqSet,
  reg: SvsRegistries,
  members: Members,
  spec: {
    order?: OrderKeyValue | readonly OrderKeyValue[] | null;
    filter?: readonly KeyFilterValue[] | null;
    subset?: ReadonlySet<number> | null;
  },
  say: KeyedSay,
): Keyed {
  const n = members.count;
  const keys =
    spec.order == null ? [members.tree ? "tree" : "data"] : Array.isArray(spec.order) ? spec.order : [spec.order];
  const rank = new Int32Array(n).fill(n);
  const tree = members.tree;
  for (let k = 0; tree && k < tree.length; k++) if (tree[k] >= 0 && tree[k] < n) rank[tree[k]] = k;
  const sorts: ((a: number, b: number) => number)[] = [];
  keys.forEach((key, k) => {
    const where = `order${Array.isArray(spec.order) ? `/${k}` : ""}`;
    if (key === "data") sorts.push((a, b) => a - b);
    else if (key === "tree") sorts.push(tree ? (a, b) => rank[a] - rank[b] : (a, b) => a - b);
    else if (key && typeof key === "object") {
      const values = column(set, reg, members, key.data, key.field, (m) => say(where, m));
      const sign = key.direction === "desc" ? -1 : 1;
      if (values) sorts.push((a, b) => compare(values[a], values[b], sign));
    }
  });
  const kept = new Uint8Array(n);
  const pool: number[] = [];
  for (let m = 0; m < n; m++) {
    if (spec.subset && !spec.subset.has(m)) continue;
    kept[m] = 1;
    pool.push(m);
  }
  (spec.filter ?? []).forEach((f, k) => {
    const values = column(set, reg, members, f.data, f.field, (m) => say(`filter/${k}`, m));
    if (!values) return;
    for (const m of pool) {
      const x = values[m];
      const pass =
        type(x) !== MISSING &&
        (f.values == null || holds(x, f.values)) &&
        (f.min == null || (typeof x === "number" && x >= f.min)) &&
        (f.max == null || (typeof x === "number" && x <= f.max));
      if (!pass) kept[m] = 0;
    }
  });
  pool.sort((a, b) => {
    for (const s of sorts) {
      const c = s(a, b);
      if (c !== 0) return c;
    }
    return a - b;
  });
  // Members a subset left out come last, unkept: the permutation lists every member once.
  const left: number[] = [];
  for (let m = 0; m < n; m++) if (spec.subset && !spec.subset.has(m)) left.push(m);
  return { permutation: Int32Array.from([...pool, ...left]), kept };
}
