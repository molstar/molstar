/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The base runtimes of QR1: the core's logic and relations; refs and the generators `all`, `empty`, `systems`, `rows`,
// `items`, `values`, `input`; the combinators `and`, `or`, `not` over selections; `range` and `slice`. Each emits runs
// in canonical order; nothing is visited per position but the values scan.

import { BASE_SYMBOLS, type Expression, isApply, type SymbolDef, suggest } from "@molstar/seq-core/query";
import { type Annotation, Annotations, type Targets } from "../annotation.js";
import { CoordinateSystemKind as K } from "../coordinates.js";
import { NO_LABEL, Numbering, type PositionLabel } from "../numbering.js";
import { KeyKind, type KeyTarget, type SeqSet } from "../set.js";
import type { CompileContext, QueryFn, SymbolRuntime } from "./compile.js";
import { F64, type QueryContext } from "./context.js";
import {
  addRuns,
  eachRun,
  firstAfter,
  intersectRun,
  liftMolecule,
  type SelectionBuilder,
  sameGroup,
  sameSystem,
  subtractRun,
} from "./kernels.js";
import { SeqSelection } from "./selection.js";
import {
  circleOf,
  domainOf,
  labelAt,
  nameOf,
  numberingNames,
  numberingNamesUnder,
  numberingOf,
  toLabel,
  unitsOf,
} from "./systems.js";

const SYMBOLS = new Map<string, SymbolDef>(BASE_SYMBOLS.map((s) => [s.name, s]));
export const rt = (name: string, compile: (c: CompileContext) => QueryFn<unknown>): SymbolRuntime => ({
  symbol: SYMBOLS.get(name) as SymbolDef,
  compile,
});
export const KIND_NAMES = ["molecule", "unit", "alignment", "row", "mapping", "annotation", "complex"];
/** The ref `resolve` found last. */
export const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const dom = new Int32Array(2);
const pair = new Int32Array(2);
const la: PositionLabel = { num: NO_LABEL, suffix: -1 as never, name: -1 as never };
const lb: PositionLabel = { ...la };

// The core: two-valued, a missing value fails every comparison (seqspec:seqql#labels-and-indices).
const missing = (v: unknown) => v == null || (typeof v === "number" && Number.isNaN(v));
const nums = (a: unknown, b: unknown) => typeof a === "number" && typeof b === "number";
/**
 * A relation; a missing side is noted for the empty-result diagnosis. Two sides that are numbers (constants, or
 * properties with a `number` runtime) are compared as doubles read from `F64`, inline: nothing boxed per element
 * (seqspec:seqql#performance).
 */
const relation = (name: string, test: (a: unknown, b: unknown) => boolean, op: number) =>
  rt(name, (c) => {
    if (op <= 1) {
      const probe = textEquals(c, 0, [c.raw("1")]) ?? textEquals(c, 1, [c.raw("0")]);
      if (probe) return probe.f(op === 0);
    }
    const [na, nb, pa, pb] = [c.number("0"), c.number("1"), at(c, 0), at(c, 1)];
    if (na && nb)
      return (ctx) => {
        na(ctx);
        const x = F64[0];
        nb(ctx);
        const y = F64[0];
        const nx = Number.isNaN(x);
        const ny = Number.isNaN(y);
        if (nx) ctx.miss(pa);
        else if (ny) ctx.miss(pb);
        return op === 0
          ? x === y
          : op === 1
            ? !nx && !ny && x !== y
            : op === 2
              ? x < y
              : op === 3
                ? x <= y
                : op === 4
                  ? x > y
                  : x >= y;
      };
    const [a, b] = [c.value("0") as QueryFn<unknown>, c.value("1") as QueryFn<unknown>];
    return (ctx) => {
      const x = a(ctx);
      const y = b(ctx);
      if (missing(x)) ctx.miss(pa);
      else if (missing(y)) ctx.miss(pb);
      return test(x, y);
    };
  });
const at = (c: CompileContext, k: number) => `${c.path}/args/${k}`;
/**
 * `=`, `!=` and `in` of a text property against constant strings, through its runtime's `equals` (CM-04: `(label)`
 * made a string per position): `f(true)` is true where the text is one of the values, `f(false)` where it is none.
 */
function textEquals(c: CompileContext, k: number, values: readonly unknown[]) {
  if (values.length === 0 || !values.every((v) => typeof v === "string")) return undefined;
  const probe = c.equals(String(k), values as string[]);
  if (!probe) return undefined;
  const p = at(c, k);
  return {
    f:
      (equal: boolean): QueryFn<boolean> =>
      (ctx) => {
        const r = probe(ctx);
        if (r < 0) {
          ctx.miss(p);
          return false;
        }
        return (r === 1) === equal;
      },
  };
}
const logic = (name: string, all: boolean) =>
  rt(name, (c) => {
    const tests = [c.value("0") as QueryFn<unknown>, ...c.values("1")];
    return (ctx) => {
      for (let i = 0; i < tests.length; i++) if ((tests[i](ctx) === true) !== all) return !all;
      return all;
    };
  });
const CORE = [
  logic("core.logic.and", true),
  logic("core.logic.or", false),
  rt("core.logic.not", (c) => {
    const f = c.value("0") as QueryFn<unknown>;
    return (ctx) => f(ctx) !== true;
  }),
  relation("core.rel.eq", (a, b) => !missing(a) && a === b, 0),
  relation("core.rel.neq", (a, b) => !missing(a) && !missing(b) && a !== b, 1),
  relation("core.rel.lt", (a, b) => nums(a, b) && (a as number) < (b as number), 2),
  relation("core.rel.lte", (a, b) => nums(a, b) && (a as number) <= (b as number), 3),
  relation("core.rel.gr", (a, b) => nums(a, b) && (a as number) > (b as number), 4),
  relation("core.rel.gre", (a, b) => nums(a, b) && (a as number) >= (b as number), 5),
  rt("core.rel.in", (c) => {
    const probe = textEquals(c, 0, c.raw<unknown[]>("1") ?? []);
    if (probe) return probe.f(true);
    const [x, vs, p0] = [c.value("0") as QueryFn<unknown>, c.values("1"), at(c, 0)];
    return (ctx) => {
      const v = x(ctx);
      if (missing(v)) {
        ctx.miss(p0);
        return false;
      }
      for (let i = 0; i < vs.length; i++) if (vs[i](ctx) === v) return true;
      return false;
    };
  }),
  rt("core.str.matches", (c) => {
    const [s, p, flags, p0] = [
      c.value("0") as QueryFn<unknown>,
      c.value("1") as QueryFn<unknown>,
      c.raw("case") ? "" : "i",
      at(c, 0),
    ];
    let re = /(?:)/;
    let source: unknown = null;
    return (ctx) => {
      // Tests run per element: no destructuring, which allocates its array where the call is not inlined.
      const text = s(ctx);
      const pattern = p(ctx);
      if (pattern !== source) {
        source = pattern;
        try {
          re = new RegExp(String(pattern), flags);
        } catch (e) {
          ctx.fail("seqql.type", c.path, `matches takes a regular expression: ${(e as Error).message}`);
        }
      }
      if (typeof text === "string") return re.test(text);
      ctx.miss(p0);
      return false;
    };
  }),
];

// Refs and the generators.

/** The refs a slot wants, for an unknown ref's message (`unknown ref "stats"; values annotations: cons, plddt`). */
export const Wants = { any: 0, annotation: 1, values: 2, alignment: 3, system: 4 } as const;
const WANTED = ["refs", "annotations", "values annotations", "alignments", "sequences"];
function refsOf(set: SeqSet, want: number): string[] {
  const key = (x: { key: number } | null) => (x ? [set.strings.get(x.key)] : []);
  const anns = set.annotations.filter((a) => a && (want !== Wants.values || a.geometry === "values"));
  const out = want === Wants.alignment ? [] : anns.flatMap(key);
  if (want === Wants.annotation || want === Wants.values) return out;
  out.push(...set.alignments.flatMap(key));
  if (want === Wants.alignment) return out;
  out.push(...set.molecules.flatMap(key), ...set.complexes.flatMap(key));
  for (let u = 0; u < set.units.count; u++) if (set.units.alive[u]) out.push(set.strings.get(set.units.key[u]));
  return out;
}
function unknownRef(ctx: QueryContext, ref: string, path: string, want: number): never {
  const refs = refsOf(ctx.set, want);
  const near = suggest(ref, refs)[0];
  const list = refs.length > 8 ? `${refs.slice(0, 8).join(", ")} (and ${refs.length - 8} more)` : refs.join(", ");
  const say = near
    ? `did you mean "${near}"?`
    : refs.length
      ? `${WANTED[want]}: ${list}`
      : `there are no ${WANTED[want]}`;
  ctx.fail("seqql.ref.unknown", path, `unknown ref "${ref}"; ${say}`);
}
/** Resolves a ref into `t`: the host's names first (a selection), else the set's keys; an error when unknown. */
export function resolve(ctx: QueryContext, ref: string, path: string, want: number = Wants.any): SeqSelection | null {
  const host = ctx.options.resolveRef?.(ref, ctx.set);
  if (host) return host;
  if (!ctx.set.byName(ref, t)) unknownRef(ctx, ref, path, want);
  return null;
}
/** `seqql.ref.kind`: what the ref is, then what was wanted (else the symbol's hint, added at the run's end). */
export const wrongKind = (ctx: QueryContext, ref: string, path: string, wants = "", what = KIND_NAMES[t.kind]): never =>
  ctx.fail("seqql.ref.kind", path, `"${ref}" is ${article(what)}${wants && `; ${wants}`}`);
const article = (w: string) => `${/^[aeio]/.test(w) ? "an" : "a"} ${w}`;

/** A whole system as one group, when it has positions held. */
function whole(b: SelectionBuilder, set: SeqSet, kind: number, id: number, row: number) {
  if (domainOf(set, kind, id, row, dom)) b.group(kind, id, row, dom[0], dom[1]);
}
/** A molecule's systems: its alive units, else itself. */
function moleculeGroups(b: SelectionBuilder, set: SeqSet, m: number) {
  const units = unitsOf(set);
  for (let k = units.start[m]; k < units.start[m + 1]; k++) whole(b, set, K.unit, units.units[k], -1);
  if (units.start[m + 1] === units.start[m]) whole(b, set, K.molecule, m, -1);
}
const sortedUnits = (set: SeqSet, complex: number) => Int32Array.from(set.complexes[complex]?.units ?? []).sort();

/** What a word names where a selection goes: its systems, one group each; an annotation its items. */
export function refQuery(ref: string, path: string): QueryFn<SeqSelection> {
  return (ctx) => {
    const host = resolve(ctx, ref, path);
    if (host) return host;
    if (t.kind === KeyKind.annotation) return items(ctx, ref, t.id, undefined, path);
    if (t.kind === KeyKind.mapping) wrongKind(ctx, ref, path, "a mapping relates sequences: name one of them");
    const [set, b, kind, id, row] = [ctx.set, ctx.builder(), t.kind, t.id, t.row];
    if (kind === KeyKind.unit) whole(b, set, K.unit, id, -1);
    else if (kind === KeyKind.molecule) moleculeGroups(b, set, id);
    else if (kind === KeyKind.complex) for (const u of sortedUnits(set, id)) whole(b, set, K.unit, u, -1);
    else if (kind === KeyKind.row) whole(b, set, K.row, id, row);
    else whole(b, set, K.columns, id, -1);
    return ctx.finish(b);
  };
}

/** Every alive unit, then every held row. */
function all(ctx: QueryContext): SeqSelection {
  const [set, b] = [ctx.set, ctx.builder()];
  for (let u = 0; u < set.units.count; u++) whole(b, set, K.unit, u, -1);
  for (const a of set.alignments) for (let r = 0; a && r < a.rowCount; r++) whole(b, set, K.row, a.id, r);
  return ctx.finish(b);
}

/**
 * Each system of `kind` within `:of` (its key kind and id; -1: everywhere), in canonical order; false when `:of`
 * cannot hold one of that kind.
 */
function eachSystem(
  ctx: QueryContext,
  kind: number,
  of: KeyTarget,
  emit: (ctx: QueryContext, id: number, row: number) => void,
): boolean {
  const set = ctx.set;
  const [ofKind, ofId] = [of.kind, of.id];
  if (kind === K.unit) {
    const units = unitsOf(set);
    if (ofKind < 0) for (let u = 0; u < set.units.count; u++) emit(ctx, u, -1);
    else if (ofKind === KeyKind.unit) emit(ctx, ofId, -1);
    else if (ofKind === KeyKind.complex) for (const u of sortedUnits(set, ofId)) emit(ctx, u, -1);
    else if (ofKind !== KeyKind.molecule) return false;
    else for (let k = units.start[ofId]; k < units.start[ofId + 1]; k++) emit(ctx, units.units[k], -1);
  } else if (kind === K.molecule) {
    if (ofKind < 0) for (let m = 0; m < set.molecules.length; m++) emit(ctx, m, -1);
    else if (ofKind === KeyKind.molecule) emit(ctx, ofId, -1);
    else if (ofKind === KeyKind.unit) emit(ctx, set.units.molecule[ofId], -1);
    else if (ofKind !== KeyKind.complex) return false;
    // In canonical order, by id (CM-12: the units' order put a later unit's lower molecule after a higher one).
    else
      for (const m of Int32Array.from(new Set(Array.from(sortedUnits(set, ofId), (u) => set.units.molecule[u]))).sort())
        emit(ctx, m, -1);
  } else {
    if (ofKind >= 0 && ofKind !== KeyKind.alignment && ofKind !== KeyKind.row) return false;
    for (let a = 0; a < set.alignments.length; a++) {
      const rows = kind === K.row ? (set.alignments[a]?.rowCount ?? 0) : 0;
      if (ofKind >= 0 && a !== ofId) continue;
      if (kind === K.columns) emit(ctx, a, -1);
      for (let r = 0; r < rows; r++) if (ofKind !== KeyKind.row || r === of.row) emit(ctx, a, r);
    }
  }
  return true;
}

const SYSTEM_KIND: Record<string, number> = { unit: K.unit, molecule: K.molecule, row: K.row, columns: K.columns };
const OF: Record<string, string> = { unit: "a molecule or a complex", molecule: "a unit or a complex" };
const EVERYWHERE: KeyTarget = { kind: -1, id: -1, row: -1 };
/** `:of`'s target, copied: a test that resolves a ref (`value-of`, `field :of`) rewrites `t`. */
const OF_T: KeyTarget = { kind: -1, id: -1, row: -1 };
/** The `systems` run in progress: one top-level `emit`, so each run calls the same function (a closure made per run
 * kept the caller from staying optimized, and its tests boxed their doubles: QR2b's heap guards). */
const E = { b: null as SelectionBuilder | null, where: null as QueryFn<unknown> | null };
const S = { kind: 0, tested: 0 };
function emit(ctx: QueryContext, id: number, row: number): void {
  const kind = S.kind;
  if (!domainOf(ctx.set, kind, id, row, dom)) return;
  ctx.slot.kind = kind;
  ctx.slot.id = id;
  ctx.slot.row = row;
  S.tested++;
  if (!E.where || E.where(ctx) === true) (E.b as SelectionBuilder).group(kind, id, row, dom[0], dom[1]);
}

/** `systems` and `rows`: the systems of a kind within `:of`, passing `:where` (the system slot), one group each. */
function systems(c: CompileContext, kindName: string): QueryFn<SeqSelection> {
  const [of, where, kind, path] = [c.raw<string>("of"), c.value("where"), SYSTEM_KIND[kindName], `${c.path}/args/of`];
  const wants =
    c.call.head === "seq.generator.rows"
      ? ""
      : `systems of kind ${kindName} take :of ${OF[kindName] ?? "an alignment"}`;
  return (ctx) => {
    const [b, missed] = [ctx.builder(), ctx.missed];
    const host = of === undefined ? null : resolve(ctx, of, path, kind >= K.row ? Wants.alignment : Wants.system);
    const target = host || of === undefined ? EVERYWHERE : OF_T;
    OF_T.kind = t.kind;
    OF_T.id = t.id;
    OF_T.row = t.row;
    const [b0, where0, kind0, tested0] = [E.b, E.where, S.kind, S.tested];
    E.b = b;
    E.where = where ?? null;
    S.kind = kind;
    S.tested = 0;
    try {
      for (let k = 0; host && k < host.runCount; k++) {
        const repeat = k > 0 && host.id[k] === host.id[k - 1] && host.row[k] === host.row[k - 1];
        if (host.kind[k] === kind && !repeat) emit(ctx, host.id[k], host.row[k]);
      }
      if (!host && !eachSystem(ctx, kind, target, emit)) wrongKind(ctx, of as string, path, wants);
      if (where && b.groupCount === 0) ctx.tested(c.path, S.tested, kind === K.row ? "row" : "system", missed);
    } finally {
      E.b = b0;
      E.where = where0;
      S.kind = kind0;
      S.tested = tested0;
    }
    return ctx.finish(b);
  };
}

/** A row's target (a slot of `tg`, or a pair's side b) [s, t) as runs: a unit's clipped to its range, a molecule's on
 * every alive unit over it (seqspec:seqql#elements-and-groups), its own system when it has none; an alignment's row or
 * columns. */
export function addTarget(b: SelectionBuilder, set: SeqSet, tg: Targets, i: number, s: number, e: number) {
  const unit = tg.unit[i];
  const mol = tg.molecule[i];
  const row = tg.row[i];
  if (unit >= 0) {
    b.add(K.unit, unit, -1, Math.max(s, set.units.rangeStart[unit]), Math.min(e, set.units.rangeStop[unit]));
    return;
  }
  if (mol >= 0 && liftMolecule(b, set, mol, s, e, unitsOf(set))) return;
  const kind = mol >= 0 ? K.molecule : row >= 0 ? K.row : K.columns;
  const id = mol >= 0 ? mol : tg.alignment[i];
  if (id >= 0 && domainOf(set, kind, id, row, dom)) b.add(kind, id, row, Math.max(s, dom[0]), Math.min(e, dom[1]));
}

export const hasItems = (a: Annotation) => !a.dense && !a.streamed && a.geometry !== "matrix";
/** One annotation's items passing `where` (the item slot), one tagged group each, with their extents. */
function itemsInto(ctx: QueryContext, b: SelectionBuilder, a: Annotation, where: QueryFn<unknown> | undefined) {
  const [set, slot, byItem, p] = [ctx.set, ctx.slot, Annotations.itemRows(a), a.pairs];
  for (let item = 0; item < a.itemCount; item++) {
    slot.annotation = a.id;
    slot.item = item;
    if (where && where(ctx) !== true) continue;
    b.begin(a.id, item);
    itemExtent(b, set, a, item, byItem, p);
    b.end();
  }
  return a.itemCount;
}
/** An item's extent into `b`'s open group: its rows' targets, a pair's both ends. */
export function itemExtent(
  b: SelectionBuilder,
  set: SeqSet,
  a: Annotation,
  item: number,
  byItem = Annotations.itemRows(a),
  p = a.pairs,
) {
  for (let k = byItem.start[item]; k < byItem.start[item + 1]; k++) {
    const r = byItem.rows[k];
    if (a.geometry === "rows") addTarget(b, set, a.target, r, 0, 0x7fffffff);
    else addTarget(b, set, a.slots, a.slot[r], a.start[r], a.stop[r]);
    if (p)
      addTarget(b, set, p.bTarget, r, p.bStart[r] < 0 ? 0 : p.bStart[r], p.bStart[r] < 0 ? 0x7fffffff : p.bStop[r]);
  }
}
function items(
  ctx: QueryContext,
  of: string | undefined,
  id: number,
  where: QueryFn<unknown> | undefined,
  path: string,
  node = path,
) {
  const [b, missed] = [ctx.builder(), ctx.missed];
  let n = 0;
  if (of === undefined) {
    for (const a of ctx.set.annotations) if (a && hasItems(a)) n += itemsInto(ctx, b, a, where);
  } else {
    const a = ctx.set.annotations[id] as Annotation;
    if (!hasItems(a)) wrongKind(ctx, of, path, "", `${a.dense ? "values" : a.geometry} annotation`);
    n = itemsInto(ctx, b, a, where);
  }
  if (where && b.groupCount === 0) ctx.tested(node, n, "item", missed);
  return ctx.finish(b);
}

/** A ref that must name an annotation: its id. */
export function annotationOf(ctx: QueryContext, ref: string, path: string, want: number, wants = ""): number {
  if (resolve(ctx, ref, path, want) || t.kind !== KeyKind.annotation) wrongKind(ctx, ref, path, wants);
  return t.id;
}

/** `values`' bounds: min, max, above, below. */
const LIMITS = new Float64Array(4);
/** The runs of `d` passing LIMITS, onto slot s's target: one pass, the comparison inline (no value boxed). */
function scan(b: SelectionBuilder, set: SeqSet, a: Annotation, s: number, d: Float64Array) {
  const lo = LIMITS[0];
  const hi = LIMITS[1];
  const above = LIMITS[2];
  const below = LIMITS[3];
  const n = d.length;
  for (let i = 0; i < n; ) {
    for (; i < n; i++) {
      const v = d[i];
      if (v >= lo && v <= hi && v > above && v < below) break;
    }
    const start = i;
    for (; i < n; i++) {
      const v = d[i];
      if (!(v >= lo && v <= hi && v > above && v < below)) break;
    }
    if (i > start) addTarget(b, set, a.slots, s, start, i);
  }
}

/** seqspec:seqql#empty-results: what a values annotation holds, for a `values` that selected nothing. */
function spanOf(a: Annotation, of: string): string {
  let [lo, hi] = [Infinity, -Infinity];
  for (const d of [...(a.dense ?? []), a.sparse ?? []])
    for (let i = 0; i < d.length; i++) {
      if (d[i] < lo) lo = d[i];
      if (d[i] > hi) hi = d[i];
    }
  const num = (v: number) => String(Number(v.toPrecision(3)));
  return lo > hi ? `${of} holds no values` : `${of} runs ${num(lo)} to ${num(hi)}`;
}

/** Positions of a values annotation passing the bounds, one group per run. */
function values(c: CompileContext): QueryFn<SeqSelection> {
  const [of, path] = [c.raw<string>("of") as string, `${c.path}/args/of`];
  const bounds = ["min", "max", "above", "below"].map((k) => c.value(k));
  return (ctx) => {
    const a = ctx.set.annotations[annotationOf(ctx, of, path, Wants.values)] as Annotation;
    if (a.geometry !== "values") wrongKind(ctx, of, path, "", `${a.geometry} annotation`);
    for (let i = 0; i < 4; i++) {
      const v = bounds[i]?.(ctx);
      LIMITS[i] = typeof v === "number" ? v : i % 2 ? Infinity : -Infinity;
    }
    const [set, b] = [ctx.set, ctx.builder()];
    b.begin();
    if (a.streamed) ctx.warn("seqql.lazy", path, `the values of ${of} are streamed and were not read`);
    for (let s = 0; a.dense && s < a.dense.length; s++) {
      ctx.visit(a.dense[s].length, path, of);
      ctx.tick();
      scan(b, set, a, s, a.dense[s]);
    }
    if (a.sparse) ctx.visit(a.count, path, of);
    for (let r = 0; a.sparse && r < a.count; r++) {
      const v = a.sparse[r];
      if (v >= LIMITS[0] && v <= LIMITS[1] && v > LIMITS[2] && v < LIMITS[3])
        addTarget(b, set, a.target, r, a.start[r], a.stop[r]);
    }
    b.end();
    const u = ctx.finish(b);
    if (u.runCount > ctx.limits.groups) ctx.over("groups", u.runCount);
    if (u.runCount === 0) ctx.because(c.path, spanOf(a, of));
    return eachRun(u);
  };
}

// The combinators: `(and A (not B))` removes B from A's groups, never building the complement; a scan that is
// a direct operand with no `:in` reads only the selection so far.

const SCANS = new Set(["seq.generator.letters", "seq.generator.positions", "seq.generator.search"]);

export function unionOf(ctx: QueryContext, s: SeqSelection): SeqSelection {
  if (s.groupCount === 1 && s.groupAnnotation[0] < 0) return s;
  const b = ctx.builder();
  b.begin();
  addRuns(b, s);
  b.end();
  return ctx.finish(b);
}

/** Whether a system of `a`'s runs has runs in `u` (one canonical group). */
function sharesSystem(a: SeqSelection, u: SeqSelection): boolean {
  for (let k = 0; k < a.runCount; k++) {
    const j = firstAfter(u, a.kind[k], a.id[k], a.row[k], -1);
    if (j < u.runCount && u.kind[j] === a.kind[k] && u.id[j] === a.id[k] && u.row[j] === a.row[k]) return true;
  }
  return false;
}

function intersect(c: CompileContext): QueryFn<SeqSelection> {
  const first = c.query("0") as QueryFn<SeqSelection>;
  const rest = (c.raw<Expression[]>("1") ?? []).map((e, i) => {
    const minus = isApply(e) && e.head === "seq.modifier.complement" && e.args?.in === undefined;
    const scan = isApply(e) && SCANS.has(e.head) && e.args?.in === undefined;
    const key = minus ? `1/${i}/args/0` : `1/${i}`;
    return { minus, scan, key: `${c.path}/args/${key}`, q: c.query(key) as QueryFn<SeqSelection> };
  });
  return (ctx) => {
    let a = first(ctx);
    for (const r of rest) {
      if (r.scan) ctx.scope = unionOf(ctx, a);
      const [u, b] = [unionOf(ctx, r.q(ctx)), ctx.builder()];
      for (let g = 0; g < a.groupCount; g++) {
        b.begin();
        for (let k = a.groupStart[g]; k < a.groupStart[g + 1]; k++) (r.minus ? subtractRun : intersectRun)(b, a, k, u);
        b.end();
      }
      const out = ctx.finish(b);
      if (out.runCount === 0 && a.runCount > 0 && u.runCount > 0) {
        const left = SeqSelection.summary(a, { maxParts: 3 });
        if (r.minus) ctx.because(c.path, `${ctx.say(r.key)} covers all of what it is taken from (${left})`);
        else if (sharesSystem(a, u))
          ctx.because(c.path, `${ctx.say(r.key)} meets none of what came before it (${left})`);
        else {
          const [x, y] = [nameOf(ctx.set, a.kind[0], a.id[0], a.row[0]), nameOf(ctx.set, u.kind[0], u.id[0], u.row[0])];
          const say = `and meets nothing: its operands share no sequence (${x}, ${y}); carry one onto the other's first: (onto Q R)`;
          ctx.warn("seqql.empty", c.path, say);
        }
      }
      a = out;
    }
    return a;
  };
}

/** Every group of every operand, in order; a group equal to an earlier one (runs and item) dropped. */
function merge(c: CompileContext): QueryFn<SeqSelection> {
  const qs = [c.query("0") as QueryFn<SeqSelection>, ...c.queries("1")];
  return (ctx) => {
    const parts = qs.map((q) => q(ctx));
    const total = parts.reduce((n, p) => n + p.groupCount, 0);
    const size = 1 << Math.ceil(Math.log2(2 * total + 2));
    const table = new Int32Array(size);
    const kept = { part: new Int32Array(total), group: new Int32Array(total), count: 0 };
    const b = ctx.builder();
    for (let p = 0; p < parts.length; p++) mergeInto(b, parts, p, table, kept);
    return ctx.finish(b);
  };
}
/** Part p's groups into `b`, each unless an equal one is kept: an open-addressing table of kept groups by hash. */
function mergeInto(b: SelectionBuilder, parts: SeqSelection[], p: number, table: Int32Array, kept: Kept) {
  const s = parts[p];
  const mask = table.length - 1;
  for (let g = 0; g < s.groupCount; g++) {
    let h = s.groupItem[g] | 0;
    for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) h = Math.imul(h ^ s.start[k], 0x9e3779b1) ^ s.id[k];
    let slot = h & mask;
    while (table[slot] && !sameGroup(parts[kept.part[table[slot] - 1]], kept.group[table[slot] - 1], s, g)) {
      slot = (slot + 1) & mask;
    }
    if (table[slot]) continue;
    kept.part[kept.count] = p;
    kept.group[kept.count] = g;
    table[slot] = ++kept.count;
    b.beginAs(s, g);
    addRuns(b, s, g);
    b.end();
  }
}
type Kept = { part: Int32Array; group: Int32Array; count: number };

/** One group: U's positions (default all) not in Q. */
function complement(c: CompileContext): QueryFn<SeqSelection> {
  const [q, within] = [c.query("0") as QueryFn<SeqSelection>, c.query("in") ?? all];
  return (ctx) => {
    const [minus, u, b] = [unionOf(ctx, q(ctx)), unionOf(ctx, within(ctx)), ctx.builder()];
    b.begin();
    for (let k = 0; k < u.runCount; k++) subtractRun(b, u, k, minus);
    b.end();
    return ctx.finish(b);
  };
}

// `range` and `slice`: labels converted once per system, slices clipped to each domain.

const PART =
  /^(-?[0-9]+(?:[A-Za-z]+|\.[0-9]+)?|[A-Za-z][A-Za-z0-9_]*)(?:-(-?[0-9]+(?:[A-Za-z]+|\.[0-9]+)?|[A-Za-z][A-Za-z0-9_]*))?$/;

/** The parts' index ranges in a numbering (null: labels are index + 1), sorted by start, into `out`; marks the parts
 * found. Once per system. */
/**
 * The label ranges as index runs in start order into `out` (two per range: one that wraps a circle's origin is its two
 * runs, CM-10); how many. `circle`: the system's (`circleOf`), 0 for none.
 */
function convert(
  set: SeqSet,
  n: Numbering | null,
  ends: readonly string[],
  out: Int32Array,
  hit: Uint8Array,
  circle: number,
) {
  let m = 0;
  for (let i = 0; i < hit.length; i++) {
    if (!toLabel(set, n, ends[2 * i], la) || !toLabel(set, n, ends[2 * i + 1], lb)) continue;
    if (n) {
      if (!Numbering.rangeOf(n, la, lb, pair, circle > 0)) continue;
    } else if (la.suffix < 0 && lb.suffix < 0 && (lb.num >= la.num || circle > 0)) {
      pair[0] = la.num - 1;
      pair[1] = lb.num;
    } else continue;
    hit[i] = 1;
    if (pair[0] < pair[1]) m = insertRun(out, m, pair[0], pair[1]);
    else {
      m = insertRun(out, m, pair[0], circle);
      m = insertRun(out, m, 0, pair[1]);
    }
  }
  return m;
}

/** Run [s, t) into the first `m` runs of `out`, kept in start order; m + 1. */
function insertRun(out: Int32Array, m: number, s: number, t: number): number {
  let j = m;
  for (; j > 0 && out[2 * j - 2] > s; j--) {
    out[2 * j] = out[2 * j - 2];
    out[2 * j + 1] = out[2 * j - 1];
  }
  out[2 * j] = s;
  out[2 * j + 1] = t;
  return m + 1;
}

/** Run k's system's labels from first to last ("1-108"); sets `dom` to its domain. */
function labelSpan(set: SeqSet, s: SeqSelection, k: number, numbering: string | null): string {
  const n = numberingOf(set, s.kind[k], s.id[k], s.row[k], numbering) ?? null;
  domainOf(set, s.kind[k], s.id[k], s.row[k], dom);
  return `${labelAt(set, n, dom[0])}-${labelAt(set, n, dom[1] - 1)}`;
}

/** "1brs/A has 108 positions, indices 0 to 108 (labels 1-108)": run k's system, for slice's messages. */
function describe(set: SeqSet, s: SeqSelection, k: number): string {
  const labels = labelSpan(set, s, k, null);
  const has =
    dom[0] > 0
      ? `holds indices ${dom[0]} to ${dom[1]} of its molecule`
      : `has ${dom[1]} positions, indices 0 to ${dom[1]}`;
  return `${nameOf(set, s.kind[k], s.id[k], s.row[k])} ${has} (labels ${labels})`;
}

function range(c: CompileContext): QueryFn<SeqSelection> {
  const q = c.query("0") as QueryFn<SeqSelection>;
  const parts = (c.raw<(string | number)[]>("1") ?? []).map(String);
  const ends = parts.flatMap((p) => {
    const m = PART.exec(p);
    return m ? [m[1], m[2] ?? m[1]] : [p, p];
  });
  const name = c.raw<string>("numbering") ?? null;
  const out = new Int32Array(4 * parts.length);
  return (ctx) => {
    const s = q(ctx);
    const set = ctx.set;
    const b = ctx.builder();
    const hit = new Uint8Array(parts.length);
    let m = -1;
    let last = -1;
    let found = 0;
    let lacking = -1;
    for (let g = 0; g < s.groupCount; g++) {
      b.begin();
      for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
        if (last < 0 || !sameSystem(s, last, k)) {
          const n = numberingOf(set, s.kind[k], s.id[k], s.row[k], name);
          m = n === undefined ? -1 : convert(set, n, ends, out, hit, circleOf(set, s.kind[k], s.id[k]));
          if (n === undefined) lacking = k;
          else found++;
          last = k;
        }
        for (let i = 0; i < m; i++) {
          b.add(s.kind[k], s.id[k], s.row[k], Math.max(s.start[k], out[2 * i]), Math.min(s.stop[k], out[2 * i + 1]));
        }
      }
      b.end();
    }
    if (found === 0 && lacking >= 0) {
      const [kind, id, row] = [s.kind[lacking], s.id[lacking], s.row[lacking]];
      const under = numberingNamesUnder(set, kind, id, row, name ?? "");
      if (under.length > 1) {
        const msg = `${nameOf(set, kind, id, row)} has ${under.length} numberings under "${name}": ${under.join(", ")}; name one, as :numbering ${under[0]}`;
        ctx.fail("seqql.numbering.ambiguous", `${c.path}/args/numbering`, msg);
      }
      const has = numberingNames(set, kind, id, row).join(", ") || "none";
      const msg = `${nameOf(set, kind, id, row)} has no numbering "${name}"; it has ${has}`;
      ctx.fail("seqql.numbering.unknown", `${c.path}/args/numbering`, msg);
    }
    for (let i = 0; last >= 0 && i < parts.length; i++) {
      if (hit[i]) continue;
      const span = labelSpan(set, s, last, name);
      const msg = `${nameOf(set, s.kind[last], s.id[last], s.row[last])} has no labels ${parts[i]}${name ? ` in ${name}` : ""}; its labels run ${span}`;
      ctx.warn("seqql.label.unknown", `${c.path}/args/1/${i}`, msg);
    }
    return ctx.finish(b);
  };
}

/** Index pairs clipped to each run; a pair past a domain is a warning, one outside every system an error. */
function slice(c: CompileContext): QueryFn<SeqSelection> {
  const q = c.query("0") as QueryFn<SeqSelection>;
  const raw = c.raw<number[]>("1") ?? [];
  const order = Array.from({ length: raw.length >> 1 }, (_, i) => i).sort((x, y) => raw[2 * x] - raw[2 * y]);
  const [starts, stops] = [Int32Array.from(order, (i) => raw[2 * i]), Int32Array.from(order, (i) => raw[2 * i + 1])];
  return (ctx) => {
    const s = q(ctx);
    const set = ctx.set;
    const b = ctx.builder();
    const inside = new Uint8Array(order.length);
    let backwards = -1;
    let clipK = -1;
    let clipP = -1;
    for (let g = 0; g < s.groupCount; g++) {
      b.begin();
      for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
        domainOf(set, s.kind[k], s.id[k], s.row[k], dom);
        const circle = circleOf(set, s.kind[k], s.id[k]);
        for (let p = 0; p < order.length; p++) {
          const p0 = starts[p];
          const p1 = stops[p];
          // A start after the stop: on a circle, the two runs across its origin (CM-10); elsewhere an error. Past a
          // circle's end it is outside it (LA-21).
          if (p0 > p1) {
            if (circle > 0 && (p0 >= circle || p1 > circle)) continue;
            if (circle === 0) {
              backwards = backwards < 0 ? p : backwards;
              continue;
            }
            inside[p] = 1;
            b.add(s.kind[k], s.id[k], s.row[k], Math.max(s.start[k], p0), Math.min(s.stop[k], circle));
            b.add(s.kind[k], s.id[k], s.row[k], Math.max(s.start[k], 0), Math.min(s.stop[k], p1));
            continue;
          }
          if (p0 < dom[1] && p1 > dom[0]) inside[p] = 1;
          if ((p0 < dom[0] && p1 > dom[0]) || (p0 < dom[1] && p1 > dom[1])) {
            clipK = k;
            clipP = p;
          }
          b.add(s.kind[k], s.id[k], s.row[k], Math.max(s.start[k], p0), Math.min(s.stop[k], p1));
        }
      }
      b.end();
    }
    const say = (k: number, p: number, what: string) =>
      `(slice ${nameOf(set, s.kind[k], s.id[k], s.row[k])} ${starts[p]} ${stops[p]}) ${what}: ${describe(set, s, k)}`;
    if (s.runCount > 0 && backwards >= 0) {
      const msg = say(0, backwards, "starts after its stop, which wraps only a cyclic sequence's origin");
      ctx.fail("seqql.slice.range", `${c.path}/args/1/${2 * order[backwards]}`, `${msg}; write the start first`);
    }
    const out = inside.indexOf(0);
    if (s.runCount > 0 && out >= 0)
      ctx.fail("seqql.slice.range", `${c.path}/args/1/${2 * order[out]}`, say(0, out, "is outside it"));
    if (clipK >= 0)
      ctx.warn(
        "seqql.slice.clipped",
        `${c.path}/args/1/${2 * order[clipP]}`,
        say(clipK, clipP, "reaches past it, and was clipped"),
      );
    return ctx.finish(b);
  };
}

export const BASE_RUNTIMES: readonly SymbolRuntime[] = [
  ...CORE,
  rt("seq.generator.all", () => all),
  rt("seq.generator.empty", () => (ctx) => SeqSelection.empty(ctx.set)),
  rt("seq.generator.systems", (c) => systems(c, c.raw<string>("kind") ?? "unit")),
  rt("seq.generator.rows", (c) => systems(c, "row")),
  rt("seq.generator.items", (c) => {
    const [of, where, path] = [c.raw<string>("of"), c.value("where"), `${c.path}/args/of`];
    const wants = "items takes :of an annotation";
    return (ctx) =>
      items(ctx, of, of === undefined ? -1 : annotationOf(ctx, of, path, Wants.annotation, wants), where, path, c.path);
  }),
  rt("seq.generator.values", values),
  rt("seq.generator.input", (c) => {
    const [name, path] = [c.raw<string>("0") as string, `${c.path}/args/0`];
    return (ctx) => {
      if (Object.hasOwn(ctx.inputs, name)) {
        if (ctx.inputs[name].groupCount === 0) ctx.because(c.path, `the app passed it empty`);
        return ctx.inputs[name];
      }
      const given = Object.keys(ctx.inputs);
      const near = suggest(name, given)[0];
      const msg = `no input "${name}" here; inputs given: ${given.join(", ") || "none"}${near ? `; did you mean "${near}"?` : ""}`;
      return ctx.fail("seqql.input.unknown", path, msg);
    };
  }),
  rt("seq.combinator.intersect", intersect),
  rt("seq.combinator.merge", merge),
  rt("seq.modifier.complement", complement),
  rt("seq.modifier.range", range),
  rt("seq.modifier.slice", slice),
];
