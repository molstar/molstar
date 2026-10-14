/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QR2b's properties: the position slot's (`index`, `letter`, `label`, `value-of`), the system slot's (`ref`, `kind`,
// `label`, `length`, `origin`, `score`, `alphabet`, `field`) and the group slot's (`length`, `runs`, `score`, and
// `overlaps`, `within`, `count`, whose selection is evaluated once per evaluation of their `pick`, before its groups
// are walked: seqspec:seqql#elements-and-groups, seqspec:seqql#performance). Each reads `ctx.slot`; a value not there
// is missing. Tests run per element: no destructuring and no closure per call.

import { type Annotation, Annotations } from "../annotation.ts";
import type { NameId, SeqIndex } from "../common.ts";
import { CoordinateSystemKind as K } from "../coordinates.ts";
import { Numbering, type PositionLabel } from "../numbering.ts";
import type { SeqSet } from "../set.ts";
import type { CompileContext, QueryFn, SymbolRuntime } from "./compile.ts";
import { F64, type QueryContext } from "./context.ts";
import { firstAfter } from "./kernels.ts";
import { annotationOf, rt, unionOf, Wants, wrongKind } from "./runtime.ts";
import type { SeqSelection } from "./selection.ts";
import { domainOf, labelAt, nameOf, numberingOf } from "./systems.ts";

const dom = new Int32Array(2);
const KINDS = ["molecule", "unit", "row", "columns"];
const str = (set: SeqSet, id: number) => (id >= 0 ? set.strings.get(id) : undefined);

/** The letter code at a position; 0 where letters are not held (the columns, a lazy molecule, a row not held). */
function codeAt(set: SeqSet, kind: number, id: number, row: number, index: number): number {
  if (kind === K.unit || kind === K.molecule) {
    const m = set.molecules[kind === K.unit ? set.units.molecule[id] : id];
    return m && index < m.code.length ? m.code[index] : 0;
  }
  const rows = kind === K.row ? set.alignments[id]?.rows : undefined;
  const s = rows ? rows.slotOf(row as never) : -1;
  return rows && s >= 0 ? rows.residues[rows.residueStart[s] + index] : 0;
}

const label: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };
const DIGITS = /^(-?[0-9]+)(.*)$/s;

/**
 * `(= (label) "27A")` and `(in (label) …)` against constants (CM-04): each constant read once per numbering into
 * (num, suffix) or a name, then compared with the position's label as the numbering holds it, so no label is made
 * per position. Equal exactly where the label's text, as `(label)` gives it, is the constant.
 */
function labelEquals(c: CompileContext, values: readonly string[]): QueryFn<number> {
  const name = c.raw<string>("numbering") ?? null;
  const at = { set: null as SeqSet | null, kind: -1, id: -1, row: -1, n: null as Numbering | null | undefined };
  const nums = new Float64Array(values.length);
  const ids = new Int32Array(values.length);
  const empty = values.includes("") ? 1 : 0;
  const parsed = { set: null as SeqSet | null, n: null as Numbering | null };
  /** The constants as labels of `n`'s kind: `ids` a suffix (-1 none, -2 none can match) or a name. */
  const parse = (set: SeqSet, n: Numbering | null) => {
    for (let k = 0; k < values.length; k++) {
      const v = values[k] as string;
      if (n?.method === "names") {
        ids[k] = set.strings.find(v);
        continue;
      }
      const m = DIGITS.exec(v);
      const ok = m !== null && String(Number(m[1])) === m[1];
      nums[k] = ok ? Number(m[1]) : Number.NaN;
      ids[k] = !ok
        ? -2
        : m[2] === ""
          ? -1
          : n === null
            ? -2
            : set.strings.find(m[2]) >= 0
              ? set.strings.find(m[2])
              : -2;
    }
  };
  return (ctx) => {
    const s = ctx.slot;
    if (ctx.set !== at.set || s.kind !== at.kind || s.id !== at.id || s.row !== at.row) {
      const n = numberingOf(ctx.set, s.kind, s.id, s.row, name);
      if (n !== undefined && (n !== parsed.n || ctx.set !== parsed.set)) {
        parse(ctx.set, n);
        parsed.set = ctx.set;
        parsed.n = n;
      }
      at.n = n;
      at.set = ctx.set;
      at.kind = s.kind;
      at.id = s.id;
      at.row = s.row;
    }
    const n = at.n;
    // No such numbering, or an unlabeled position: the text is "".
    if (n === undefined) return empty;
    if (n === null) {
      for (let k = 0; k < values.length; k++) if (ids[k] === -1 && nums[k] === s.index + 1) return 1;
      return 0;
    }
    if (!Numbering.labelOf(n, s.index as SeqIndex, label)) return empty;
    if (n.method === "names") {
      for (let k = 0; k < values.length; k++) if (ids[k] === label.name) return 1;
      return 0;
    }
    for (let k = 0; k < values.length; k++) if (nums[k] === label.num && ids[k] === label.suffix) return 1;
    return 0;
  };
}

/** `(label :numbering N)`: the numbering looked up once per system, the label per position. */
function positionLabel(c: CompileContext): QueryFn<unknown> {
  const name = c.raw<string>("numbering") ?? null;
  const at = { set: null as SeqSet | null, kind: -1, id: -1, row: -1, n: null as Numbering | null | undefined };
  return (ctx) => {
    const s = ctx.slot;
    if (ctx.set !== at.set || s.kind !== at.kind || s.id !== at.id || s.row !== at.row) {
      at.n = numberingOf(ctx.set, s.kind, s.id, s.row, name);
      at.set = ctx.set;
      at.kind = s.kind;
      at.id = s.id;
      at.row = s.row;
    }
    return at.n === undefined ? "" : labelAt(ctx.set, at.n, s.index);
  };
}

const space = { kind: 0, id: 0, row: 0 };
/** A values annotation's slot over a system: its own, or, for a unit, its molecule's. */
function slotFor(set: SeqSet, a: Annotation, kind: number, id: number, row: number): number {
  space.kind = kind;
  space.id = id;
  space.row = row;
  const s = Annotations.slotOf(a, space);
  if (s >= 0 || kind !== K.unit) return s;
  space.kind = K.molecule;
  space.id = set.units.molecule[id];
  return Annotations.slotOf(a, space);
}

/** A number property from its writer of `F64[0]`: the relations call the writer, the rest read it back. */
const numeric = (name: string, writer: (c: CompileContext) => QueryFn<void>): SymbolRuntime => ({
  ...rt(name, (c) => {
    const w = writer(c);
    return (ctx) => {
      w(ctx);
      return F64[0];
    };
  }),
  number: writer,
});

/** `(value-of A)`: A resolved once per set, its slot once per system, its value per position (NaN: none); a dense
 * array is read inline (a double returned from a call is boxed). */
function valueOfAnnotation(c: CompileContext): QueryFn<void> {
  const [ref, path] = [String(c.raw("0")), `${c.path}/args/0`];
  const at = { set: null as SeqSet | null, a: null as Annotation | null, kind: -1, id: -1, row: -1, slot: -1 };
  return (ctx) => {
    const s = ctx.slot;
    if (ctx.set !== at.set) {
      const a = ctx.set.annotations[annotationOf(ctx, ref, path, Wants.values)] as Annotation;
      if (a.geometry !== "values") wrongKind(ctx, ref, path, "", `${a.geometry} annotation`);
      at.a = a;
      at.set = ctx.set;
      at.kind = -1;
    }
    if (s.kind !== at.kind || s.id !== at.id || s.row !== at.row) {
      at.slot = slotFor(ctx.set, at.a as Annotation, s.kind, s.id, s.row);
      at.kind = s.kind;
      at.id = s.id;
      at.row = s.row;
    }
    const a = at.a as Annotation;
    const d = at.slot >= 0 && a.dense ? a.dense[at.slot] : null;
    if (d) F64[0] = s.index < d.length ? d[s.index] : Number.NaN;
    else F64[0] = at.slot < 0 ? Number.NaN : Annotations.valueAt(a, at.slot, s.index);
  };
}

// The system slot: missing when the slot has no one system (a group over several, in `pick`).

const system = (name: string, read: (set: SeqSet, kind: number, id: number, row: number) => unknown) =>
  rt(name, () => (ctx) => (ctx.slot.kind < 0 ? undefined : read(ctx.set, ctx.slot.kind, ctx.slot.id, ctx.slot.row)));
const moleculeOf = (set: SeqSet, kind: number, id: number) =>
  kind === K.unit ? set.molecules[set.units.molecule[id]] : kind === K.molecule ? set.molecules[id] : null;

/** A rows annotation's record per row of one alignment (-1: none), built once per annotation and alignment. */
const RECORDS = new WeakMap<Annotation, Map<number, Int32Array>>();
function recordsOf(set: SeqSet, a: Annotation, aln: number): Int32Array {
  let by = RECORDS.get(a);
  if (!by) {
    by = new Map();
    RECORDS.set(a, by);
  }
  let x = by.get(aln);
  if (!x) {
    x = new Int32Array(set.alignments[aln]?.rowCount ?? 0).fill(-1);
    for (let r = a.count - 1; r >= 0; r--)
      if (a.target.alignment[r] === aln && a.target.row[r] >= 0) x[a.target.row[r]] = r;
    by.set(aln, x);
  }
  return x;
}
const recordOf = (x: Int32Array, row: number) => (row >= 0 && row < x.length ? x[row] : -1);

/** `(field name :of K :ns N)`: a built-in, a rows key table's column for the row, or an app's custom field. */
function systemField(c: CompileContext): QueryFn<unknown> {
  const [name, of, ns, path] = [String(c.raw("0")), c.raw<string>("of"), c.raw<string>("ns"), `${c.path}/args/of`];
  const at = { set: null as SeqSet | null, a: null as Annotation | null };
  return (ctx) => {
    const set = ctx.set;
    const kind = ctx.slot.kind;
    const id = ctx.slot.id;
    const row = ctx.slot.row;
    if (kind < 0) return undefined;
    if (of !== undefined) {
      if (at.set !== set) at.a = keyTable(ctx, of, name, path, c.path);
      at.set = set;
      const r = kind === K.row ? recordOf(recordsOf(set, at.a as Annotation, id), row) : -1;
      return r < 0 ? undefined : Annotations.rowValue(at.a as Annotation, set, r, name);
    }
    const m = moleculeOf(set, kind, id);
    if (ns !== undefined) {
      const own = kind === K.unit ? set.units.custom[id] : kind === K.molecule ? m?.custom : set.alignments[id]?.custom;
      const x = (own?.doc[ns] as Record<string, unknown> | undefined)?.[name];
      return x !== undefined || kind !== K.unit
        ? x
        : (m?.custom.doc[ns] as Record<string, unknown> | undefined)?.[name];
    }
    if (name === "molecule") return m ? set.strings.get(m.key) : undefined;
    if (name === "complex")
      return kind === K.unit ? str(set, set.complexes[set.units.complex[id]]?.key ?? -1) : undefined;
    if (name === "alignment") return kind >= K.row ? str(set, set.alignments[id]?.key ?? -1) : undefined;
    return name === "row" && kind === K.row ? row : undefined;
  };
}
/** `(field name :of K)` as a number: a number column read inline; anything else through the field, NaN unless a
 * number. */
function fieldNumber(c: CompileContext): QueryFn<void> {
  const [field, of, name] = [systemField(c), c.raw<string>("of"), String(c.raw("0"))];
  // Mutable bindings, not an object's fields: code that read a field written once embedded what it held, and was
  // thrown away when that was collected (QR2b's heap guards).
  let set: SeqSet | null = null;
  let a: Annotation | null = null;
  let data: Float64Array | null = null;
  let aln = -1;
  let records = new Int32Array(0) as Int32Array;
  return (ctx) => {
    if (of !== undefined && ctx.slot.kind === K.row) {
      if (set !== ctx.set) {
        a = keyTable(ctx, of, name, `${c.path}/args/of`, c.path);
        const f = a.fields.find((x) => ctx.set.strings.get(x.name) === name);
        data = f?.kind === "number" ? f.data : null;
        set = ctx.set;
        aln = -1;
      }
      if (data) {
        // The records looked up once per alignment, then read by row: per row, two typed-array reads.
        if (ctx.slot.id !== aln) {
          records = recordsOf(ctx.set, a as Annotation, ctx.slot.id);
          aln = ctx.slot.id;
        }
        const r = recordOf(records, ctx.slot.row);
        F64[0] = r < 0 ? Number.NaN : data[r];
        return;
      }
    }
    const v = field(ctx);
    F64[0] = typeof v === "number" ? v : Number.NaN;
  };
}

/** `:of`'s annotation: a rows key table with the column; else an error naming what it is or what it has. */
function keyTable(ctx: QueryContext, of: string, name: string, path: string, at: string): Annotation {
  const a = ctx.set.annotations[annotationOf(ctx, of, path, Wants.annotation)] as Annotation;
  if (a.geometry !== "rows") wrongKind(ctx, of, path, ":of takes a rows key table", `${a.geometry} annotation`);
  const fields = a.fields.map((f) => ctx.set.strings.get(f.name));
  if (!fields.includes(name))
    ctx.fail("seqql.field.unknown", at, `${of} has no field "${name}"; its fields: ${fields.join(", ") || "none"}`);
  return a;
}

// The group slot: the group `pick` tests.

/** The test of a `pick`'s :where and the group tests' selections in it, evaluated before the groups are walked. */
let preparing: QueryFn<void>[] | null = null;
export function compileWhere(c: CompileContext): { where: QueryFn<unknown>; prepare: readonly QueryFn<void>[] } {
  const [outer, prepare] = [preparing, [] as QueryFn<void>[]];
  preparing = prepare;
  try {
    return { where: c.value("where") as QueryFn<unknown>, prepare };
  } finally {
    preparing = outer;
  }
}

const group = (name: string, read: (s: SeqSelection, g: number) => unknown) =>
  rt(name, () => (ctx) => {
    const s = ctx.slot.selection;
    return s && ctx.slot.group >= 0 ? read(s, ctx.slot.group) : undefined;
  });

/** Group g's positions within `n` of `u`'s (0: shared), in the same system; or, with `count`, how many it shares. */
function near(s: SeqSelection, g: number, u: SeqSelection, n: number, count: boolean): number {
  let shared = 0;
  for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
    const kind = s.kind[k];
    const id = s.id[k];
    const row = s.row[k];
    const st = s.start[k];
    const sp = s.stop[k];
    for (let j = firstAfter(u, kind, id, row, st - n); j < u.runCount; j++) {
      if (u.kind[j] !== kind || u.id[j] !== id || u.row[j] !== row || u.start[j] >= sp + n) break;
      if (!count) return 1;
      shared += Math.min(sp, u.stop[j]) - Math.max(st, u.start[j]);
    }
  }
  return shared;
}
/** `overlaps`, `within`, `count`: Q unioned once per evaluation of the `pick`, each group by binary search. */
function nearTo(c: CompileContext, mode: "overlaps" | "within" | "count"): QueryFn<unknown> {
  const [q, n] = [c.query("0") as QueryFn<SeqSelection>, mode === "within" ? (c.raw<number>("1") ?? 0) : 0];
  const at = { u: null as SeqSelection | null };
  preparing?.push((ctx) => {
    at.u = unionOf(ctx, q(ctx));
  });
  return (ctx) => {
    const s = ctx.slot.selection;
    const g = ctx.slot.group;
    const u = at.u;
    if (!s || g < 0 || !u) return undefined;
    return mode === "count" ? near(s, g, u, 0, true) : near(s, g, u, n, false) > 0;
  };
}

export const PROPERTY_RUNTIMES = [
  rt("seq.position.index", () => (ctx) => ctx.slot.index),
  rt("seq.position.letter", () => (ctx) => {
    const s = ctx.slot;
    const code = codeAt(ctx.set, s.kind, s.id, s.row, s.index);
    return code ? String.fromCharCode(code) : "";
  }),
  { ...rt("seq.position.label", positionLabel), equals: labelEquals },
  numeric("seq.position.value_of", valueOfAnnotation),
  system("seq.system.ref", (set, kind, id, row) =>
    kind === K.columns ? str(set, set.alignments[id]?.key ?? -1) : nameOf(set, kind, id, row),
  ),
  system("seq.system.kind", (_, kind) => KINDS[kind]),
  system("seq.system.label", (set, kind, id, row) =>
    kind === K.unit
      ? str(set, set.units.label[id])
      : kind === K.molecule
        ? str(set, set.molecules[id]?.label ?? -1)
        : kind === K.row
          ? (set.alignments[id]?.labelText(row as never) ?? undefined)
          : undefined,
  ),
  system("seq.system.length", (set, kind, id, row) => (domainOf(set, kind, id, row, dom) ? dom[1] - dom[0] : 0)),
  system("seq.system.origin", (set, kind, id) => (kind === K.unit ? str(set, set.units.origin[id]) : undefined)),
  numeric("seq.system.score", () => (ctx) => {
    F64[0] = ctx.slot.kind === K.unit ? ctx.set.units.score[ctx.slot.id] : Number.NaN;
  }),
  system(
    "seq.system.alphabet",
    (set, kind, id) => (moleculeOf(set, kind, id)?.alphabet ?? set.alignments[id]?.alphabet) || undefined,
  ),
  { ...rt("seq.system.field", systemField), number: fieldNumber },
  group("seq.group.length", (s, g) => {
    let n = 0;
    for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) n += s.stop[k] - s.start[k];
    return n;
  }),
  group("seq.group.runs", (s, g) => s.groupStart[g + 1] - s.groupStart[g]),
  numeric("seq.group.score", () => (ctx) => {
    const s = ctx.slot.selection;
    F64[0] = s && ctx.slot.group >= 0 ? s.groupScore[ctx.slot.group] : Number.NaN;
  }),
  rt("seq.group.overlaps", (c) => nearTo(c, "overlaps")),
  rt("seq.group.within", (c) => nearTo(c, "within")),
  rt("seq.group.count", (c) => nearTo(c, "count")),
];
