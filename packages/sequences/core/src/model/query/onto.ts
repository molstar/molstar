/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `onto`: each group's image in other coordinate systems, by overlap, along projection's routes
// (seqspec:model#projection). Within an alignment, a row and the columns (or another row, through them) through the
// row's map; rows and units through `RowBindings`, the columns through every bound row (a row bound to several
// sequences reaches each: `onto` names where it goes, so no binding is ambiguous); else the model's shortest path
// (`Mappings.path`, ambiguous paths skipped). Routes are built once per (source, target) system in a run; a walk reads
// the maps (`Coords`) with module visitors and allocates nothing per run.

import type { RowIndex } from "../common.js";
import { type CoordinateMap, Coords, CoordinateSystemKind as K } from "../coordinates.js";
import { Mappings, type PathStep, type RowBindings } from "../mapping.js";
import { KeyKind, type SeqSet } from "../set.js";
import type { CompileContext, QueryFn } from "./compile.js";
import type { QueryContext } from "./context.js";
import { addRuns, type SelectionBuilder, sameSystem } from "./kernels.js";
import { KIND_NAMES, resolve, rt, t, Wants, wrongKind } from "./runtime.js";
import type { SeqSelection } from "./selection.js";
import { domainOf, nameOf, unitsOf } from "./systems.js";

/** A leg: the input clipped to [lo, hi) of its system, then through `map` (read 1: from → to; -1: the inverse image). */
interface Step {
  readonly map: CoordinateMap;
  readonly read: number;
  readonly lo: number;
  readonly hi: number;
}
type Chain = readonly Step[];
interface Sys {
  readonly kind: number;
  readonly id: number;
  readonly row: number;
}
/** Why a source reached nothing (the first reason met), and the first ambiguous path. */
interface Why {
  say: string;
  ambiguous: string;
}

const dom = new Int32Array(2);
const steps: PathStep[] = [];
const sys = (kind: number, id: number, row = -1): Sys => ({ kind, id, row });
const same = (a: Sys, b: Sys) => a.kind === b.kind && a.id === b.id && a.row === b.row;
const alignmentOf = (s: Sys) => (s.kind === K.row || s.kind === K.columns ? s.id : -1);
const name = (set: SeqSet, s: Sys) => nameOf(set, s.kind, s.id, s.row);
const rowLength = (set: SeqSet, a: number, row: number) => set.alignments[a]?.rows.lengthOf(row as RowIndex) ?? 0;

/** Rows' maps while one source's routes are built (a row reached through many targets is mapped once). */
const rowMaps = new Map<number, CoordinateMap | null>();
/** A row's residues onto its alignment's columns (read 1), or back (-1); null when the row is not held. */
function rowStep(set: SeqSet, a: number, row: number, read: number, why: Why): Step | null {
  const al = set.alignments[a];
  const key = a * 0x100000000 + row;
  if (al && !rowMaps.has(key)) rowMaps.set(key, Coords.ofRow(al, row as RowIndex));
  const map = al ? (rowMaps.get(key) ?? null) : null;
  if (!al || !map) {
    why.say ||= `${nameOf(set, K.row, a, row)} is not loaded`;
    return null;
  }
  return { map, read, lo: 0, hi: read > 0 ? rowLength(set, a, row) : al.columnCount };
}
/** Binding entry e: its row onto its unit (`toUnit`), or the unit onto the row. */
function bindStep(set: SeqSet, b: RowBindings, a: number, e: number, toUnit: boolean): Step {
  const u = b.unit[e];
  const fromRow = b.rowIsFrom[e] === 1;
  const [lo, hi] = toUnit
    ? [0, rowLength(set, a, b.rows[b.entryRow[e]])]
    : [set.units.rangeStart[u], set.units.rangeStop[u]];
  return { map: b.map[e], read: fromRow === toUnit ? 1 : -1, lo, hi };
}
/** The shortest path's composed map; null when there is none or it is ambiguous (noted). */
function pathStep(set: SeqSet, S: Sys, T: Sys, why: Why): Step | null {
  const from = { kind: S.kind, id: S.id, row: S.row };
  const n = Mappings.path(set, from, { kind: T.kind, id: T.id, row: T.row }, steps);
  if (n === -2) why.ambiguous ||= `${name(set, S)} reaches ${name(set, T)} by paths that place it differently`;
  if (n <= 0) return null;
  const p = Mappings.compose(set, from, steps);
  domainOf(set, S.kind, S.id, S.row, dom);
  return { map: p.map, read: Coords.sameSystem(p.map.from, from) ? 1 : -1, lo: dom[0], hi: dom[1] };
}

/** Every chain from S onto T (its routes), appended to `out` after `prefix`. */
function routes(set: SeqSet, S: Sys, T: Sys, prefix: Step[], out: Step[][], why: Why): void {
  if (same(S, T)) {
    out.push(prefix);
    return;
  }
  const [sa, ta] = [alignmentOf(S), alignmentOf(T)];
  if (sa >= 0 && sa === ta) {
    const up = S.kind === K.row ? rowStep(set, sa, S.row, 1, why) : null;
    const down = T.kind === K.row ? rowStep(set, ta, T.row, -1, why) : null;
    if ((S.kind !== K.row || up) && (T.kind !== K.row || down)) out.push([...prefix, ...[up, down].filter((x) => !!x)]);
    return;
  }
  if (sa >= 0) {
    fromRows(set, S, T, prefix, out, why);
    return;
  }
  if (ta >= 0 && S.kind === K.unit && toRows(set, S, T, prefix, out, why)) return;
  const p = pathStep(set, S, T, why);
  if (p) out.push([...prefix, p]);
  else why.say ||= `${name(set, S)} has no route onto ${name(set, T)}`;
}

/** From a row (its own binding) or the columns (every bound row's) onto a system outside the alignment. */
function fromRows(set: SeqSet, S: Sys, T: Sys, prefix: Step[], out: Step[][], why: Why) {
  const a = S.id;
  const b = Mappings.bindingsOf(set, a);
  const rows = S.kind === K.row ? [S.row] : Array.from(b?.rows ?? []);
  const n0 = out.length;
  // The units bound to the rows first; when none is the target, the routes from them.
  for (let pass = 0; pass < 2 && out.length === n0; pass++) {
    for (const r of rows) {
      const k = b ? Mappings.boundRow(b, r) : -1;
      if (!b || k < 0) {
        why.say ||= `${nameOf(set, K.row, a, r)} is bound to no sequence`;
        continue;
      }
      const down = S.kind === K.columns ? rowStep(set, a, r, -1, why) : null;
      if (S.kind === K.columns && !down) continue;
      for (let e = b.start[k]; e < b.start[k + 1]; e++) {
        const direct = T.kind === K.unit && T.id === b.unit[e];
        if (pass === 0 && !direct) continue;
        const chain = [...prefix, ...(down ? [down] : []), bindStep(set, b, a, e, true)];
        if (pass === 0) out.push(chain);
        else routes(set, sys(K.unit, b.unit[e]), T, chain, out, why);
      }
    }
  }
  if (rows.length === 0) why.say ||= `no row of ${nameOf(set, K.columns, a, -1).split(" ")[0]} is bound to a sequence`;
}

/** From a unit bound to rows of T's alignment: its own row, else through its rows and the columns. False: unbound. */
function toRows(set: SeqSet, S: Sys, T: Sys, prefix: Step[], out: Step[][], why: Why): boolean {
  const b = Mappings.bindingsOf(set, T.id);
  const x = b ? Mappings.boundUnit(b, S.id) : -1;
  if (!b || x < 0) return false;
  const n0 = out.length;
  for (let pass = 0; pass < 2 && out.length === n0; pass++) {
    for (let i = b.unitStart[x]; i < b.unitStart[x + 1]; i++) {
      const e = b.byUnit[i];
      const r = b.rows[b.entryRow[e]];
      const direct = T.kind === K.row && T.row === r;
      if (pass === 0 && direct) out.push([...prefix, bindStep(set, b, T.id, e, false)]);
      if (pass === 1) routes(set, sys(K.row, T.id, r), T, [...prefix, bindStep(set, b, T.id, e, false)], out, why);
    }
  }
  return true;
}

// The walk: one chain onto one target, into the builder's open group; a level per step, made once.
const W = { b: null as SelectionBuilder | null, chain: [] as Chain, kind: 0, id: 0, row: 0, lo: 0, hi: 0 };
const LEVELS = Array.from({ length: 16 }, (_, i) => ({ i: i + 1 }));
function walk(i: number, s: number, e: number): void {
  const c = W.chain;
  if (i === c.length) {
    (W.b as SelectionBuilder).add(W.kind, W.id, W.row, s < W.lo ? W.lo : s, e > W.hi ? W.hi : e);
    return;
  }
  const st = c[i];
  const a = s < st.lo ? st.lo : s;
  const z = e > st.hi ? st.hi : e;
  if (a >= z) return;
  if (st.read > 0) Coords.forEachTo(st.map, a, z, next, LEVELS[i]);
  else Coords.forEachFrom(st.map, a, z, next, LEVELS[i]);
}
const next = (s: number, e: number, level: { i: number }) => walk(level.i, s, e);

const KINDS: Record<string, number> = { unit: K.unit, molecule: K.molecule, row: K.row, columns: K.columns };
/** The systems `onto` lands on: R's (an alignment its columns, or its rows with `:kind row`), or every one of `:kind`. */
function targetsOf(ctx: QueryContext, c: CompileContext, ref: string | undefined, kind: string | undefined): Sys[] {
  const [set, out] = [ctx.set, [] as Sys[]];
  const add = (k: number, id: number, row = -1) => {
    if (domainOf(set, k, id, row, dom)) out.push(sys(k, id, row));
  };
  const rows = (a: number) => {
    for (let r = 0; r < (set.alignments[a]?.rowCount ?? 0); r++) add(K.row, a, r);
  };
  if (ref === undefined) {
    const k = KINDS[kind as string];
    if (k === K.unit) for (let u = 0; u < set.units.count; u++) add(K.unit, u);
    if (k === K.molecule) for (let m = 0; m < set.molecules.length; m++) add(K.molecule, m);
    for (let a = 0; k >= K.row && a < set.alignments.length; a++) (k === K.row ? rows : add.bind(null, K.columns))(a);
    return out;
  }
  const [rPath, kPath] = [`${c.path}/args/1`, `${c.path}/args/kind`];
  const host = resolve(ctx, ref, rPath, Wants.system);
  // A selection's runs are in canonical order: a system's runs sit together.
  for (let k = 0; host && k < host.runCount; k++)
    if (k === 0 || !sameSystem(host, k - 1, k)) add(host.kind[k], host.id[k], host.row[k]);
  if (host) return out;
  const [kk, id, row] = [t.kind, t.id, t.row];
  if (kk === KeyKind.annotation || kk === KeyKind.mapping)
    wrongKind(ctx, ref, rPath, "onto takes a unit, molecule, row or alignment");
  if (kk === KeyKind.alignment) {
    if (kind === "row") rows(id);
    else if (kind === undefined || kind === "columns") add(K.columns, id);
    else {
      const say = `:kind ${kind} onto an alignment: it gives its columns, or its rows with :kind row; for the ${kind}s bound to its rows, leave R out: (onto Q :kind ${kind})`;
      ctx.fail("seqql.type", kPath, say);
    }
    return out;
  }
  if (kind !== undefined) {
    const alignment = nameOf(set, K.columns, id, -1).split(" ")[0];
    const say =
      kk === KeyKind.row
        ? `"${ref}" is one row already: :kind is for an alignment's rows or columns, (onto Q ${alignment} :kind row) for all its rows; it was ignored`
        : `"${ref}" is a ${KIND_NAMES[kk]}: :kind is for an alignment's rows or columns, or for every system of a kind without R; it was ignored`;
    ctx.warn("seqql.onto.kind", kPath, say);
  }
  const units = unitsOf(set);
  if (kk === KeyKind.unit) add(K.unit, id);
  else if (kk === KeyKind.row) add(K.row, id, row);
  else if (kk === KeyKind.complex)
    for (const u of [...(set.complexes[id]?.units ?? [])].sort((x, y) => x - y)) add(K.unit, u);
  else if (units.start[id + 1] === units.start[id]) add(K.molecule, id);
  else for (let k = units.start[id]; k < units.start[id + 1]; k++) add(K.unit, units.units[k]);
  return out;
}

/** A source's routes onto each target, and why it reached nothing, kept while the set and the targets stay. */
interface Routes {
  readonly chains: Chain[][];
  readonly say: string;
  readonly ambiguous: string;
}
const sameTargets = (a: readonly Sys[], b: readonly Sys[]) => a.length === b.length && a.every((x, i) => same(x, b[i]));

/** `(onto Q R :kind K :keep b)`: each group's image on R's systems; a group with none is dropped and counted. */
function onto(c: CompileContext): QueryFn<SeqSelection> {
  const q = c.query("0") as QueryFn<SeqSelection>;
  const [ref, kind, keep] = [c.raw<string>("1"), c.raw<string>("kind"), c.raw<boolean>("keep") === true];
  // Routes compose maps (a row's, a path's): built once per (set, targets), reused by every run until they change.
  const memo = { set: null as SeqSet | null, targets: [] as Sys[], routes: new Map<string, Routes>() };
  return (ctx) => {
    if (ref === undefined && kind === undefined) {
      const say =
        "onto needs where to go: a ref, (onto Q msa), or :kind for every system of a kind, (onto Q :kind unit)";
      ctx.fail("seqql.arg.missing", c.path, say);
    }
    const targets = targetsOf(ctx, c, ref === undefined ? undefined : String(ref), kind);
    const [set, s, b] = [ctx.set, q(ctx), ctx.builder()];
    if (memo.set !== set || !sameTargets(memo.targets, targets)) {
      memo.set = set;
      memo.targets = targets;
      memo.routes = new Map();
    }
    const why: Why = { say: "", ambiguous: "" };
    const [lo, hi] = [new Int32Array(targets.length), new Int32Array(targets.length)];
    for (let i = 0; i < targets.length; i++) {
      domainOf(set, targets[i].kind, targets[i].id, targets[i].row, dom);
      lo[i] = dom[0];
      hi[i] = dom[1];
    }
    let chains: Chain[][] = [];
    let unmapped = 0;
    W.b = b;
    for (let g = 0, last = -1; g < s.groupCount; g++) {
      b.begin();
      const before = b.runCount;
      for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
        if (last < 0 || !sameSystem(s, last, k)) {
          last = k;
          ctx.tick();
          const r = routesOf(set, memo.routes, sys(s.kind[k], s.id[k], s.row[k]), targets);
          why.say ||= r.say;
          why.ambiguous ||= r.ambiguous;
          chains = r.chains;
        }
        for (let i = 0; i < targets.length; i++) {
          const T = targets[i];
          W.kind = T.kind;
          W.id = T.id;
          W.row = T.row;
          W.lo = lo[i];
          W.hi = hi[i];
          const list = chains[i];
          for (let j = 0; j < list.length; j++) {
            W.chain = list[j];
            walk(0, s.start[k], s.stop[k]);
          }
        }
      }
      if (b.runCount === before) unmapped++;
      if (keep) addRuns(b, s, g);
      b.end();
    }
    W.b = null;
    const where = ref ?? `every ${kind}`;
    if (why.ambiguous) ctx.warn("seqql.onto.ambiguous", c.path, `${why.ambiguous}; it was skipped`);
    if (unmapped > 0) {
      const say = why.say || "their positions map to nothing there (gaps, inserts, or what the maps do not pair)";
      const n = `${unmapped} of ${s.groupCount} group${s.groupCount === 1 ? "" : "s"}`;
      ctx.warn("seqql.onto.unmapped", c.path, `${n} ${unmapped === 1 ? "has" : "have"} no image on ${where}: ${say}`);
    }
    return ctx.finish(b);
  };
}

/** A source system's routes onto every target, from the memo or built now. */
function routesOf(set: SeqSet, memo: Map<string, Routes>, S: Sys, targets: readonly Sys[]): Routes {
  const key = `${S.kind}:${S.id}:${S.row}`;
  let r = memo.get(key);
  if (!r) {
    const why: Why = { say: "", ambiguous: "" };
    let chains: Step[][][];
    try {
      chains = targets.map((T) => {
        const out: Step[][] = [];
        routes(set, S, T, [], out, why);
        return out;
      });
    } finally {
      // Cleared however the routes end: a throw left one set's row maps for the next to read (GI-14).
      rowMaps.clear();
    }
    r = { chains, say: why.say, ambiguous: why.ambiguous };
    memo.set(key, r);
  }
  return r;
}

export const ONTO_RUNTIMES = [rt("seq.modifier.onto", onto)];
