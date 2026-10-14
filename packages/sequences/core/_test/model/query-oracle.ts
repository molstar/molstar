/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A brute-force SeqQL evaluator for QR1's, QR2a's and QR2b's symbols (seqspec:seqql#labels-and-indices): every
// selection a list of groups of single positions, every symbol by its definition, position by position, with none of
// the runtime's kernels. Its labels follow seqspec:format#positions read from each position's own label
// (Numbering.labelOf), never through `rangeOf`. The runtime's results are compared against it (query.test.ts).

import { type Expression, isApply } from "@molstar/seq-core/query";
import {
  type Annotation,
  Annotations,
  ColumnRoles,
  type CoordinateMap,
  Coords,
  CoordinateSystemKind as K,
  KeyKind,
  Mappings,
  NO_LABEL,
  type Numbering,
  Numbering as Numberings,
  type PathStep,
  type PositionLabel,
  type RowBindings,
  type SeqSet,
  SeqSets,
  type Targets,
  type UnitId,
} from "../../src/model/index.js";
import type { SeqSelection } from "../../src/model/query/index.js";

/** One group: its item (-1, -1: none) and its positions, as sorted keys; a search hit's score. */
export interface OGroup {
  readonly ann: number;
  readonly item: number;
  readonly pos: readonly number[];
  readonly score?: number;
}
export class OracleError extends Error {
  readonly code: string;
  constructor(code: string) {
    super(code);
    this.code = code;
  }
}

const key = (kind: number, id: number, row: number, p: number) => ((kind * 1000 + id) * 1000 + row + 1) * 1e6 + p;
export const decode = (k: number) => {
  const p = k % 1e6;
  const s = (k - p) / 1e6;
  const row = (s % 1000) - 1;
  const id = Math.floor(s / 1000) % 1000;
  return { kind: Math.floor(s / 1e6), id, row, p, system: s };
};
const norm = (xs: Iterable<number>) => [...new Set(xs)].sort((a, b) => a - b);
const group = (pos: Iterable<number>, ann = -1, item = -1): OGroup => ({ ann, item, pos: norm(pos) });
const nonEmpty = (gs: OGroup[]) => gs.filter((g) => g.pos.length > 0);

/** A runtime selection as oracle groups, checking the runs' invariants on the way. */
export function groupsOf(s: SeqSelection): OGroup[] {
  const out: OGroup[] = [];
  for (let g = 0; g < s.groupCount; g++) {
    const pos: number[] = [];
    for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
      if (!(s.start[k] < s.stop[k])) throw new Error(`run ${k} is empty`);
      if (k > s.groupStart[g]) {
        const a = key(s.kind[k - 1], s.id[k - 1], s.row[k - 1], s.stop[k - 1]);
        if (!(a < key(s.kind[k], s.id[k], s.row[k], s.start[k]))) throw new Error(`runs ${k - 1}, ${k} out of order`);
      }
      for (let p = s.start[k]; p < s.stop[k]; p++) pos.push(key(s.kind[k], s.id[k], s.row[k], p));
    }
    if (pos.length === 0) throw new Error(`group ${g} is empty`);
    const score = s.groupScore[g];
    out.push({ ann: s.groupAnnotation[g], item: s.groupItem[g], pos, ...(Number.isNaN(score) ? {} : { score }) });
  }
  return out;
}

export class Oracle {
  readonly set: SeqSet;
  readonly inputs: Readonly<Record<string, OGroup[]>>;
  constructor(set: SeqSet, inputs: Readonly<Record<string, OGroup[]>> = {}) {
    this.set = set;
    this.inputs = inputs;
  }

  /** [lo, hi) of a system; hi -1 when none. */
  domain(kind: number, id: number, row: number): [number, number] {
    const { set } = this;
    const u = set.units;
    if (kind === K.unit) return u.alive[id] ? [u.rangeStart[id], u.rangeStop[id]] : [0, -1];
    if (kind === K.molecule) return [0, set.molecules[id]?.length ?? -1];
    const a = set.alignments[id];
    if (kind === K.row) return [0, a ? a.rows.lengthOf(row as never) : -1];
    return [0, a?.columnCount ?? -1];
  }
  whole(kind: number, id: number, row: number): OGroup {
    const [lo, hi] = this.domain(kind, id, row);
    return group(Array.from({ length: Math.max(0, hi - lo) }, (_, i) => key(kind, id, row, lo + i)));
  }
  unitsOver(m: number): number[] {
    const u = this.set.units;
    return Array.from({ length: u.count }, (_, i) => i).filter((i) => u.alive[i] && u.molecule[i] === m);
  }
  all(): OGroup[] {
    const { set } = this;
    const out = Array.from({ length: set.units.count }, (_, u) => this.whole(K.unit, u, -1));
    for (const a of set.alignments) for (let r = 0; a && r < a.rowCount; r++) out.push(this.whole(K.row, a.id, r));
    return nonEmpty(out);
  }

  /** Positions [s, e) of a target, a molecule's on every alive unit over it (or itself), clipped to domains. */
  target(tg: Targets, i: number, s: number, e: number): number[] {
    const out: number[] = [];
    const put = (kind: number, id: number, row: number) => {
      const [lo, hi] = this.domain(kind, id, row);
      for (let p = Math.max(s, lo); p < Math.min(e, hi); p++) out.push(key(kind, id, row, p));
    };
    if (tg.unit[i] >= 0) put(K.unit, tg.unit[i], -1);
    else if (tg.molecule[i] >= 0) {
      const units = this.unitsOver(tg.molecule[i]);
      for (const u of units) put(K.unit, u, -1);
      if (units.length === 0) put(K.molecule, tg.molecule[i], -1);
    } else if (tg.alignment[i] >= 0) put(tg.row[i] >= 0 ? K.row : K.columns, tg.alignment[i], tg.row[i]);
    return out;
  }
  items(a: Annotation): OGroup[] {
    if (a.dense || a.geometry === "matrix") throw new OracleError("seqql.ref.kind");
    const out: OGroup[] = [];
    for (let item = 0; item < a.itemCount; item++) {
      const pos: number[] = [];
      for (let r = 0; r < a.count; r++) {
        if (a.item[r] !== item) continue;
        if (a.geometry === "rows") pos.push(...this.target(a.target, r, 0, 1e5));
        else pos.push(...this.target(a.slots, a.slot[r], a.start[r], a.stop[r]));
        const p = a.pairs;
        if (p)
          pos.push(...this.target(p.bTarget, r, p.bStart[r] < 0 ? 0 : p.bStart[r], p.bStart[r] < 0 ? 1e5 : p.bStop[r]));
      }
      out.push(group(pos, a.id, item));
    }
    return nonEmpty(out);
  }

  ref(ref: string): OGroup[] {
    const { set } = this;
    const t = { kind: -1, id: -1, row: -1 };
    if (!set.byName(ref, t)) throw new OracleError("seqql.ref.unknown");
    if (t.kind === KeyKind.unit) return nonEmpty([this.whole(K.unit, t.id, -1)]);
    if (t.kind === KeyKind.molecule) {
      const units = this.unitsOver(t.id);
      return nonEmpty(units.length ? units.map((u) => this.whole(K.unit, u, -1)) : [this.whole(K.molecule, t.id, -1)]);
    }
    if (t.kind === KeyKind.complex) {
      const units = [...(set.complexes[t.id]?.units ?? [])].sort((a, b) => a - b);
      return nonEmpty(units.map((u) => this.whole(K.unit, u, -1)));
    }
    if (t.kind === KeyKind.row) return nonEmpty([this.whole(K.row, t.id, t.row)]);
    if (t.kind === KeyKind.alignment) return nonEmpty([this.whole(K.columns, t.id, -1)]);
    if (t.kind === KeyKind.annotation) return this.items(set.annotations[t.id] as Annotation);
    throw new OracleError("seqql.ref.kind");
  }

  /** A system's numbering by name, or its default; null: index + 1; undefined: it lacks that name. */
  numbering(kind: number, id: number, row: number, name: string | null): Numbering | null | undefined {
    const { set } = this;
    const n = name === null ? -1 : set.strings.find(name);
    const mol = kind === K.unit ? set.units.molecule[id] : kind === K.molecule ? id : -1;
    const own =
      kind === K.unit
        ? set.units.numberings[id]
        : kind === K.row
          ? set.alignments[id]?.numberingsOf(row as never)
          : null;
    const molecule = set.molecules[mol]?.numberings;
    if (name !== null) {
      if (kind === K.unit) return SeqSets.numberingOf(set, id as UnitId, n as never) ?? undefined;
      return (kind === K.row ? own : molecule)?.byName(n as never) ?? undefined;
    }
    const list = [...(own?.items ?? []), ...(molecule?.items ?? [])];
    return list.find((x) => x.isDefault) ?? list.find((x) => set.strings.get(x.name) === "ordinal") ?? null;
  }

  /** Whether position p's label lies in the range [a, z] by seqspec:seqql#overview: index order between present endpoints, an absent
   * endpoint placed by comparing integer part, then suffix (none, letters, `.n`). */
  inRange(n: Numbering | null, a: string, z: string): (p: number) => boolean {
    const { set } = this;
    const none = () => false;
    const lab: PositionLabel = { num: NO_LABEL, suffix: -1 as never, name: -1 as never };
    type L = [number, string];
    const at = (p: number): L | null => {
      if (!(n && Numberings.labelOf(n, p as never, lab))) return null;
      if (n.method === "names") return [0, set.strings.get(lab.name)];
      return [lab.num, lab.suffix >= 0 ? set.strings.get(lab.suffix) : ""];
    };
    const parse = (x: string): L | null => {
      if (n?.method === "names") return [0, x];
      const m = /^(-?[0-9]+)([A-Za-z]|\.[0-9]+)?$/.exec(x);
      return m ? [Number(m[1]), m[2] ?? ""] : null;
    };
    const [la, lz] = [parse(a), parse(z)];
    if (!la || !lz) return none;
    if (!n) return la[1] || lz[1] ? none : (p) => p + 1 >= la[0] && p + 1 <= lz[0];
    const rank = (x: string) =>
      !x ? 0 : x[0] === "." ? 1000 + Number(x.slice(1)) : x.toUpperCase().charCodeAt(0) - 64;
    const cmp = (x: L, y: L) => x[0] - y[0] || rank(x[1]) - rank(y[1]);
    const same = (x: L | null, y: L) => !!x && x[0] === y[0] && x[1] === y[1];
    let start = -1;
    let last = -1;
    for (let p = n.length - 1; p >= 0; p--) {
      if (same(at(p), la)) start = p;
      if (same(at(p), lz)) last = p;
    }
    if (n.method !== "names") {
      for (let p = 0; start < 0 && p < n.length; p++) if (at(p) && cmp(at(p) as L, la) >= 0) start = p;
      if (last < 0) for (let p = 0; p < n.length; p++) if (at(p) && cmp(at(p) as L, lz) <= 0) last = p;
    }
    if (start < 0 || last < start) return none;
    return (p) => p >= start && p <= last;
  }

  range(q: OGroup[], parts: string[], name: string | null): OGroup[] {
    const ends = parts.map((p) => {
      const m =
        /^(-?[0-9]+(?:[A-Za-z]|\.[0-9]+)?|[A-Za-z][A-Za-z0-9_]*)(?:-(-?[0-9]+(?:[A-Za-z]|\.[0-9]+)?|[A-Za-z][A-Za-z0-9_]*))?$/.exec(
          p,
        );
      return m ? [m[1], m[2] ?? m[1]] : [p, p];
    });
    let found = 0;
    let lacking = 0;
    const out = q.map((g) =>
      group(
        g.pos.filter((k) => {
          const { kind, id, row, p } = decode(k);
          const n = this.numbering(kind, id, row, name);
          if (n === undefined) {
            lacking++;
            return false;
          }
          found++;
          return ends.some(([a, z]) => this.inRange(n, a, z)(p));
        }),
      ),
    );
    if (found === 0 && lacking > 0) throw new OracleError("seqql.numbering.unknown");
    return nonEmpty(out);
  }

  slice(q: OGroup[], bounds: number[]): OGroup[] {
    const pairs = Array.from({ length: bounds.length / 2 }, (_, i) => [bounds[2 * i], bounds[2 * i + 1]]);
    const systems = new Set(q.flatMap((g) => g.pos.map((k) => decode(k).system)));
    for (const [p0, p1] of pairs) {
      const inside = [...systems].some((s) => {
        const { kind, id, row } = decode(s * 1e6);
        const [lo, hi] = this.domain(kind, id, row);
        return p0 < hi && p1 > lo;
      });
      if (systems.size > 0 && !inside) throw new OracleError("seqql.slice.range");
    }
    return nonEmpty(
      q.map((g) => group(g.pos.filter((k) => pairs.some(([a, b]) => decode(k).p >= a && decode(k).p < b)))),
    );
  }

  values(a: Annotation, lo: number, hi: number, above: number, below: number): OGroup[] {
    if (a.geometry !== "values") throw new OracleError("seqql.ref.kind");
    const ok = (v: number) => v >= lo && v <= hi && v > above && v < below;
    const pos: number[] = [];
    for (let s = 0; a.dense && s < a.dense.length; s++) {
      for (let i = 0; i < a.dense[s].length; i++)
        if (ok(Annotations.valueAt(a, s, i))) pos.push(...this.target(a.slots, s, i, i + 1));
    }
    for (let r = 0; !a.dense && a.sparse && r < a.count; r++)
      if (ok(a.sparse[r])) pos.push(...this.target(a.target, r, a.start[r], a.stop[r]));
    return runsOf(norm(pos));
  }

  systems(kind: string, of: Expression | undefined): OGroup[] {
    const { set } = this;
    const t = { kind: -1, id: -1, row: -1 };
    if (of !== undefined && !set.byName(String(of), t)) throw new OracleError("seqql.ref.unknown");
    const range = (n: number) => Array.from({ length: n }, (_, i) => i);
    const sorted = (c: number) => [...(set.complexes[c]?.units ?? [])].sort((a, b) => a - b);
    const rows = (a: number) => range(set.alignments[a]?.rowCount ?? 0).map((r) => [a, r]);
    const by: Record<string, Record<number, () => number[][]>> = {
      unit: {
        [-1]: () => range(set.units.count).map((u) => [u, -1]),
        [KeyKind.molecule]: () => this.unitsOver(t.id).map((u) => [u, -1]),
        [KeyKind.complex]: () => sorted(t.id).map((u) => [u, -1]),
        [KeyKind.unit]: () => [[t.id, -1]],
      },
      molecule: {
        [-1]: () => range(set.molecules.length).map((m) => [m, -1]),
        [KeyKind.molecule]: () => [[t.id, -1]],
        [KeyKind.unit]: () => [[set.units.molecule[t.id], -1]],
        [KeyKind.complex]: () => [...new Set(sorted(t.id).map((u) => set.units.molecule[u]))].map((m) => [m, -1]),
      },
      row: {
        [-1]: () => range(set.alignments.length).flatMap(rows),
        [KeyKind.alignment]: () => rows(t.id),
        [KeyKind.row]: () => [[t.id, t.row]],
      },
      columns: {
        [-1]: () => range(set.alignments.length).map((a) => [a, -1]),
        [KeyKind.alignment]: () => [[t.id, -1]],
        [KeyKind.row]: () => [[t.id, -1]],
      },
    };
    const k = { unit: K.unit, molecule: K.molecule, row: K.row, columns: K.columns }[kind] as number;
    const list = by[kind][of === undefined ? -1 : t.kind];
    if (!list) throw new OracleError("seqql.ref.kind");
    return nonEmpty(list().map(([id, row]) => this.whole(k, id, row)));
  }

  eval(e: Expression): OGroup[] {
    if (!isApply(e)) return this.ref(String(e));
    const arg = (k: string) => e.args?.[k] as Expression;
    const many = (k: string) => (e.args?.[k] ?? []) as Expression[];
    switch (e.head) {
      case "seq.generator.all":
        return this.all();
      case "seq.generator.empty":
        return [];
      case "seq.generator.input": {
        const x = this.inputs[String(arg("0"))];
        if (!x) throw new OracleError("seqql.input.unknown");
        return x;
      }
      case "seq.generator.items": {
        const of = arg("of");
        // QR2a: `:where` on the item slot, tested for every item in order (as the runtime does, empty ones included).
        const where = (a: Annotation) => {
          const keep = Array.from({ length: a.itemCount }, (_, item) => {
            this.slot = { ann: a.id, item };
            return arg("where") === undefined || this.test(arg("where")) === true;
          });
          return this.items(a).filter((g) => keep[g.item]);
        };
        if (of === undefined)
          return this.set.annotations.flatMap((a) => (a && !a.dense && a.geometry !== "matrix" ? where(a) : []));
        const t = { kind: -1, id: -1, row: -1 };
        if (!this.set.byName(String(of), t)) throw new OracleError("seqql.ref.unknown");
        if (t.kind !== KeyKind.annotation) throw new OracleError("seqql.ref.kind");
        const a = this.set.annotations[t.id] as Annotation;
        if (a.dense || a.geometry === "matrix") throw new OracleError("seqql.ref.kind");
        return where(a);
      }
      case "seq.generator.values": {
        const t = { kind: -1, id: -1, row: -1 };
        if (!this.set.byName(String(arg("of")), t)) throw new OracleError("seqql.ref.unknown");
        if (t.kind !== KeyKind.annotation) throw new OracleError("seqql.ref.kind");
        const num = (k: string, d: number) => (typeof arg(k) === "number" ? (arg(k) as number) : d);
        const a = this.set.annotations[t.id] as Annotation;
        return this.values(
          a,
          num("min", -Infinity),
          num("max", Infinity),
          num("above", -Infinity),
          num("below", Infinity),
        );
      }
      case "seq.generator.systems":
      case "seq.generator.rows":
        return this.systems(e.head.endsWith("rows") ? "row" : String(arg("kind") ?? "unit"), arg("of")).filter((g) => {
          const { kind, id, row } = decode(g.pos[0]);
          Object.assign(this.at, { kind, id, row, p: -1, g: null });
          return arg("where") === undefined || this.test(arg("where")) === true;
        });
      case "seq.combinator.intersect": {
        let a = this.eval(arg("0"));
        for (const x of many("1")) {
          const minus = isApply(x) && x.head === "seq.modifier.complement" && x.args?.in === undefined;
          // QR2b: a scan with no :in of its own reads the selection so far (seqspec:seqql#performance)
          const scan = isApply(x) && /^seq\.generator\.(letters|positions|search)$/.test(x.head) && !x.args?.in;
          if (scan) this.within = norm(a.flatMap((g) => g.pos));
          const u = new Set(this.eval(minus && isApply(x) ? (x.args?.["0"] as Expression) : x).flatMap((g) => g.pos));
          a = nonEmpty(a.map((g) => group(g.pos.filter((k) => u.has(k) !== minus))));
        }
        return a;
      }
      case "seq.combinator.merge": {
        const out: OGroup[] = [];
        for (const g of [arg("0"), ...many("1")].flatMap((x) => this.eval(x))) {
          const same = out.some((h) => h.ann === g.ann && h.item === g.item && h.pos.join() === g.pos.join());
          if (!same) out.push(g);
        }
        return out;
      }
      case "seq.modifier.complement": {
        const q = new Set(this.eval(arg("0")).flatMap((g) => g.pos));
        const u = arg("in") === undefined ? this.all() : this.eval(arg("in"));
        return nonEmpty([group(u.flatMap((g) => g.pos).filter((k) => !q.has(k)))]);
      }
      case "seq.modifier.range":
        return this.range(this.eval(arg("0")), many("1").map(String), (arg("numbering") as string) ?? null);
      case "seq.modifier.slice":
        return this.slice(this.eval(arg("0")), many("1") as number[]);
      default:
        return this.qr2a(e);
    }
  }

  // QR2a: onto, partners, expand, group, each with current, pick, first, ids, the item slot's tests.

  /** Inside `each`, the group being run on. */
  current: OGroup | null = null;
  /** The item a test reads: (annotation, item); -1 none. */
  slot = { ann: -1, item: -1 };

  qr2a(e: Expression & object): OGroup[] {
    if (!isApply(e)) throw new Error("not a call");
    const arg = (k: string) => e.args?.[k] as Expression;
    const many = (k: string) => (e.args?.[k] ?? []) as Expression[];
    switch (e.head) {
      case "seq.generator.current":
        if (!this.current) throw new OracleError("seqql.current");
        return [this.current];
      case "seq.modifier.expand": {
        const n = (arg("1") as number) ?? 0;
        const [l, r] = [(arg("left") as number) ?? n, (arg("right") as number) ?? n];
        return nonEmpty(
          this.eval(arg("0")).map((g) =>
            group(
              g.pos.flatMap((k) => {
                const { kind, id, row, p } = decode(k);
                const [lo, hi] = this.domain(kind, id, row);
                const out: number[] = [];
                for (let x = Math.max(lo, p - l); x < Math.min(hi, p + r + 1); x++) out.push(key(kind, id, row, x));
                return out;
              }),
            ),
          ),
        );
      }
      case "seq.modifier.group": {
        const all = norm(this.eval(arg("0")).flatMap((g) => g.pos));
        if (arg("by") === "run") return runsOf(all);
        if (arg("by") === "none") return nonEmpty([group(all)]);
        const by = new Map<number, number[]>();
        for (const k of all) by.set(decode(k).system, [...(by.get(decode(k).system) ?? []), k]);
        return [...by.values()].map((pos) => group(pos));
      }
      case "seq.modifier.each": {
        const out: OGroup[] = [];
        const prev = this.current;
        for (const g of this.eval(arg("0"))) {
          this.current = g;
          out.push(...this.eval(arg("1")));
        }
        this.current = prev;
        return out;
      }
      case "seq.filter.pick": {
        // QR2b: a group test's selection is the same for every group of one pick (computed once there,
        // seqspec:seqql#elements-and-groups)
        const memo = this.memo;
        this.memo = new Map();
        const out = this.eval(arg("0")).filter((g) => {
          this.slot = { ann: g.ann, item: g.item };
          const one = new Set(g.pos.map((k) => decode(k).system)).size === 1;
          const { kind, id, row } = one ? decode(g.pos[0]) : { kind: -1, id: -1, row: -1 };
          Object.assign(this.at, { kind, id, row, p: -1, g });
          return this.test(arg("where")) === true;
        });
        this.memo = memo;
        return out;
      }
      case "seq.filter.first":
        return this.eval(arg("0")).slice(0, (arg("1") as number) ?? 1);
      case "seq.filter.ids": {
        const ids = many("1").map(String);
        return this.eval(arg("0")).filter((g) => {
          const a = this.set.annotations[g.ann];
          const k = a && g.item >= 0 ? a.itemKey[g.item] : -1;
          return k >= 0 && ids.includes(this.set.strings.get(k));
        });
      }
      case "seq.modifier.onto":
        return this.onto(e);
      case "seq.modifier.partners":
        return this.partners(e);
      default:
        return this.qr2b(e);
    }
  }

  /** An item's positions, by its rows (the extent rule of `items`). */
  itemPositions(a: Annotation, item: number): readonly number[] {
    return this.items(a).find((g) => g.item === item)?.pos ?? [];
  }

  /** A test on the item slot: the core's logic, relations and `matches`, and the item's properties. */
  test(e: Expression): unknown {
    if (!isApply(e)) return e;
    const v = (k: string) => this.test(e.args?.[k] as Expression);
    const rest = () => ((e.args?.["1"] ?? []) as Expression[]).map((x) => this.test(x));
    const missing = (x: unknown) => x === undefined || x === null || (typeof x === "number" && Number.isNaN(x));
    const nums = (a: unknown, b: unknown) => typeof a === "number" && typeof b === "number";
    const a = this.set.annotations[this.slot.ann];
    const item = this.slot.item;
    const row = () => (a ? Array.from({ length: a.count }, (_, r) => r).find((r) => a.item[r] === item) : undefined);
    const col = (name: string) => {
      const r = row();
      return a && r !== undefined ? Annotations.rowValue(a, this.set, r, name) : undefined;
    };
    switch (e.head) {
      case "core.logic.and":
        return [e.args?.["0"] as Expression, ...((e.args?.["1"] ?? []) as Expression[])].every(
          (x) => this.test(x) === true,
        );
      case "core.logic.or":
        return [e.args?.["0"] as Expression, ...((e.args?.["1"] ?? []) as Expression[])].some(
          (x) => this.test(x) === true,
        );
      case "core.logic.not":
        return v("0") !== true;
      case "core.rel.eq":
        return ((x, y) => !missing(x) && x === y)(v("0"), v("1"));
      case "core.rel.neq":
        return ((x, y) => !missing(x) && !missing(y) && x !== y)(v("0"), v("1"));
      case "core.rel.lt":
        return ((x, y) => nums(x, y) && (x as number) < (y as number))(v("0"), v("1"));
      case "core.rel.gr":
        return ((x, y) => nums(x, y) && (x as number) > (y as number))(v("0"), v("1"));
      case "core.rel.lte":
        return ((x, y) => nums(x, y) && (x as number) <= (y as number))(v("0"), v("1"));
      case "core.rel.gre":
        return ((x, y) => nums(x, y) && (x as number) >= (y as number))(v("0"), v("1"));
      case "core.rel.in": {
        const x = v("0");
        return !missing(x) && rest().includes(x);
      }
      case "core.str.matches": {
        const [text, pattern] = [v("0"), v("1")];
        return typeof text === "string" && new RegExp(String(pattern), e.args?.case ? "" : "i").test(text);
      }
      case "seq.item.id":
        return a && item >= 0 && a.itemKey[item] >= 0 ? this.set.strings.get(a.itemKey[item]) : undefined;
      case "seq.item.label":
        return col("label");
      case "seq.item.category":
        return col("category");
      case "seq.item.value":
        return col("value");
      case "seq.item.type":
        return a && item >= 0 && a.semanticType >= 0 ? this.set.strings.get(a.semanticType) : undefined;
      case "seq.item.annotation":
        return a && item >= 0 ? this.set.strings.get(a.key) : undefined;
      case "seq.item.length":
        return a && item >= 0 ? this.itemPositions(a, item).length : undefined;
      case "seq.item.field": {
        const name = String(e.args?.["0"]);
        if (!a || item < 0) return undefined;
        const known =
          ColumnRoles[a.geometry][name] !== undefined || a.fields.some((f) => this.set.strings.get(f.name) === name);
        if (!known) throw new OracleError("seqql.field.unknown");
        return col(name);
      }
      default:
        return this.property(e);
    }
  }

  // QR2b: letters, positions, search; the position, system and group slots' properties.

  /** What a position, system or group test reads (the item slot is `slot`); -1, null: none. */
  at = { kind: -1, id: -1, row: -1, p: -1, g: null as OGroup | null };
  /** The group tests' selections within one evaluation of a pick. */
  memo = new Map<Expression, Set<number>>();
  /** An `and`'s selection so far, for the scan that is its next operand. */
  within: number[] | null = null;
  /** Letter-class tables for `search :classes`. */
  classes: Readonly<Record<string, Readonly<Record<string, string>>>> = {};

  /** A system's letters ("" not held: the columns, a lazy molecule). */
  letters(kind: number, id: number, row: number): string {
    const key = `${kind}:${id}:${row}`;
    if (!this.text.has(key)) this.text.set(key, this.read(kind, id, row));
    return this.text.get(key) as string;
  }
  private readonly text = new Map<string, string>();
  private read(kind: number, id: number, row: number): string {
    const { set } = this;
    if (kind === K.unit || kind === K.molecule) {
      const m = set.molecules[kind === K.unit ? set.units.molecule[id] : id];
      return m && !m.lazy ? new TextDecoder().decode(m.code) : "";
    }
    const rows = kind === K.row ? set.alignments[id]?.rows : undefined;
    const s = rows ? rows.slotOf(row as never) : -1;
    return rows && s >= 0
      ? new TextDecoder().decode(rows.residues.subarray(rows.residueStart[s], rows.residueStart[s + 1]))
      : "";
  }
  /** The positions a scan reads: its `:in`, else every unit (and, but for search, every row). */
  scope(e: Expression & object, rows: boolean): number[] {
    if (!isApply(e)) throw new Error("not a call");
    const within = this.within;
    this.within = null;
    if (e.args?.in !== undefined) return norm(this.eval(e.args.in as Expression).flatMap((g) => g.pos));
    if (within) return within;
    return norm(
      this.all()
        .filter((g) => rows || decode(g.pos[0]).kind === K.unit)
        .flatMap((g) => g.pos),
    );
  }
  /** Search hits in the runs of `pos`, each start tried: exact letters, a motif of classes, or a sticky regex. */
  hits(pos: number[], e: Expression & object): OGroup[] {
    if (!isApply(e)) throw new Error("not a call");
    const p = String(e.args?.["0"]);
    const table = e.args?.classes === undefined ? null : this.classes[String(e.args.classes)];
    if (e.args?.classes !== undefined && !table) throw new OracleError("seqql.classes.unknown");
    const re = e.args?.regex ? new RegExp(p, "iy") : null;
    const sets = [...p.toUpperCase()].map((c) => (table && Object.hasOwn(table, c) ? table[c] : c));
    const out: OGroup[] = [];
    for (const run of runsOf(pos)) {
      const { kind, id, row, p: first } = decode(run.pos[0]);
      const text = this.letters(kind, id, row).slice(first, first + run.pos.length);
      for (let i = 0; i < text.length; i++) {
        let n = -1;
        if (re) {
          re.lastIndex = i;
          n = re.exec(text)?.[0].length || -1;
        } else if (i + sets.length <= text.length && sets.every((s, j) => s.includes(text[i + j]))) n = sets.length;
        if (n > 0) out.push({ ...group(run.pos.slice(i, i + n)), score: 0 });
      }
    }
    return out;
  }

  qr2b(e: Expression & object): OGroup[] {
    if (!isApply(e)) throw new Error("not a call");
    const where = e.args?.where as Expression | undefined;
    switch (e.head) {
      case "seq.generator.letters": {
        const L = String(e.args?.["0"]).toUpperCase();
        const ok = (k: number) => {
          const { kind, id, row, p } = decode(k);
          const c = this.letters(kind, id, row)[p];
          return c !== undefined && L.includes(c);
        };
        return runsOf(this.scope(e, true).filter(ok));
      }
      case "seq.generator.positions":
        return runsOf(
          this.scope(e, true).filter((k) => {
            const { kind, id, row, p } = decode(k);
            Object.assign(this.at, { kind, id, row, p, g: null });
            return where === undefined || this.test(where) === true;
          }),
        );
      case "seq.generator.search":
        return this.hits(this.scope(e, false), e);
      default:
        throw new Error(`the oracle has no ${e.head}`);
    }
  }

  /** The position, system and group properties, read from `at`. */
  property(e: Expression & object): unknown {
    if (!isApply(e)) throw new Error("not a call");
    const { set } = this;
    const { kind, id, row, p, g } = this.at;
    const arg = (k: string) => e.args?.[k] as Expression;
    const str = (x: number) => (x >= 0 ? set.strings.get(x) : undefined);
    const mol =
      kind === K.unit ? set.molecules[set.units.molecule[id]] : kind === K.molecule ? set.molecules[id] : null;
    const aln = kind === K.row || kind === K.columns ? set.alignments[id] : null;
    /** Q's positions, evaluated without disturbing the test's slots. */
    const positionsOf = (x: Expression) => {
      const saved = [{ ...this.at }, this.slot] as const;
      const out = this.memo.get(x) ?? new Set(this.eval(x).flatMap((h) => h.pos));
      this.memo.set(x, out);
      [this.at, this.slot] = [saved[0], saved[1]];
      return out;
    };
    if (e.head.startsWith("seq.system.") && kind < 0) return undefined;
    if (e.head.startsWith("seq.group.") && !g) return undefined;
    switch (e.head) {
      case "seq.position.index":
        return p;
      case "seq.position.letter":
        return this.letters(kind, id, row)[p] ?? "";
      case "seq.position.label": {
        const n = this.numbering(kind, id, row, (arg("numbering") as string) ?? null);
        return n === undefined ? "" : n ? Numberings.format(n, p as never, set.strings) : String(p + 1);
      }
      case "seq.position.value_of": {
        const t = { kind: -1, id: -1, row: -1 };
        if (!set.byName(String(arg("0")), t)) throw new OracleError("seqql.ref.unknown");
        const a = set.annotations[t.id];
        if (t.kind !== KeyKind.annotation || !a || a.geometry !== "values") throw new OracleError("seqql.ref.kind");
        const sl = a.slots;
        const own = (s: number) =>
          kind === K.unit
            ? sl.unit[s] === id
            : kind === K.molecule
              ? sl.unit[s] < 0 && sl.molecule[s] === id
              : sl.alignment[s] === id && sl.row[s] === (kind === K.row ? row : -1);
        let s = sl.unit.findIndex((_, x) => own(x));
        if (s < 0 && kind === K.unit)
          s = sl.unit.findIndex((u, x) => u < 0 && sl.molecule[x] === set.units.molecule[id]);
        return s < 0 ? Number.NaN : Annotations.valueAt(a, s, p);
      }
      case "seq.system.ref":
        return kind === K.columns
          ? str(aln?.key ?? -1)
          : kind === K.row
            ? `${str(aln?.key ?? -1)}/${row}`
            : str(kind === K.unit ? set.units.key[id] : (mol?.key ?? -1));
      case "seq.system.kind":
        return ["molecule", "unit", "row", "columns"][kind];
      case "seq.system.label":
        return kind === K.unit
          ? str(set.units.label[id])
          : kind === K.molecule
            ? str(mol?.label ?? -1)
            : kind === K.row
              ? (aln?.labelText(row as never) ?? undefined)
              : undefined;
      case "seq.system.length": {
        const [lo, hi] = this.domain(kind, id, row);
        return Math.max(0, hi - lo);
      }
      case "seq.system.origin":
        return kind === K.unit ? str(set.units.origin[id]) : undefined;
      case "seq.system.score":
        return kind === K.unit ? set.units.score[id] : Number.NaN;
      case "seq.system.alphabet":
        return (mol?.alphabet ?? aln?.alphabet) || undefined;
      case "seq.system.field": {
        const name = String(arg("0"));
        if (arg("of") !== undefined) {
          const t = { kind: -1, id: -1, row: -1 };
          if (!set.byName(String(arg("of")), t)) throw new OracleError("seqql.ref.unknown");
          const a = set.annotations[t.id];
          if (t.kind !== KeyKind.annotation || !a || a.geometry !== "rows") throw new OracleError("seqql.ref.kind");
          if (!a.fields.some((f) => set.strings.get(f.name) === name)) throw new OracleError("seqql.field.unknown");
          const r = Array.from({ length: a.count }, (_, x) => x).find(
            (x) => kind === K.row && a.target.alignment[x] === id && a.target.row[x] === row,
          );
          return r === undefined ? undefined : Annotations.rowValue(a, set, r, name);
        }
        if (arg("ns") !== undefined) {
          const ns = String(arg("ns"));
          const own = kind === K.unit ? set.units.custom[id] : (mol?.custom ?? aln?.custom);
          const v = (own?.doc[ns] as Record<string, unknown> | undefined)?.[name];
          return v !== undefined || kind !== K.unit
            ? v
            : (mol?.custom.doc[ns] as Record<string, unknown> | undefined)?.[name];
        }
        if (name === "molecule") return mol ? str(mol.key) : undefined;
        if (name === "complex")
          return kind === K.unit ? str(set.complexes[set.units.complex[id]]?.key ?? -1) : undefined;
        if (name === "alignment") return aln ? str(aln.key) : undefined;
        return name === "row" && kind === K.row ? row : undefined;
      }
      case "seq.group.length":
        return g?.pos.length;
      case "seq.group.runs":
        return runsOf(g?.pos ?? []).length;
      case "seq.group.score":
        return g?.score ?? Number.NaN;
      case "seq.group.overlaps":
      case "seq.group.within": {
        const q = positionsOf(arg("0"));
        const n = e.head.endsWith("within") ? (arg("1") as number) : 0;
        const near = (k: number) => Array.from({ length: 2 * n + 1 }, (_, d) => k - n + d).some((x) => q.has(x));
        return (g?.pos ?? []).some(near);
      }
      case "seq.group.count": {
        const q = positionsOf(arg("0"));
        return (g?.pos ?? []).filter((k) => q.has(k)).length;
      }
      default:
        throw new Error(`the oracle's tests have no ${e.head}`);
    }
  }

  /** The systems `onto` lands on, as the runtime finds them (and its errors, in its order). */
  targets(ref: Expression | undefined, kind: string | undefined): Sys[] {
    const { set } = this;
    const out: Sys[] = [];
    const add = (k: number, id: number, row = -1) => {
      const [lo, hi] = this.domain(k, id, row);
      if (hi > lo && !out.some((x) => x.kind === k && x.id === id && x.row === row)) out.push({ kind: k, id, row });
    };
    const rows = (a: number) => {
      for (let r = 0; r < (set.alignments[a]?.rowCount ?? 0); r++) add(K.row, a, r);
    };
    if (ref === undefined && kind === undefined) throw new OracleError("seqql.arg.missing");
    if (ref === undefined) {
      if (kind === "unit") for (let u = 0; u < set.units.count; u++) add(K.unit, u);
      if (kind === "molecule") for (let m = 0; m < set.molecules.length; m++) add(K.molecule, m);
      for (let a = 0; a < set.alignments.length; a++)
        kind === "row" ? rows(a) : kind === "columns" ? add(K.columns, a) : 0;
      return out;
    }
    const t = { kind: -1, id: -1, row: -1 };
    if (!set.byName(String(ref), t)) throw new OracleError("seqql.ref.unknown");
    if (t.kind === KeyKind.annotation || t.kind === KeyKind.mapping) throw new OracleError("seqql.ref.kind");
    if (t.kind === KeyKind.alignment) {
      if (kind === "row") rows(t.id);
      else if (kind === undefined || kind === "columns") add(K.columns, t.id);
      else throw new OracleError("seqql.type");
    } else if (t.kind === KeyKind.unit) add(K.unit, t.id);
    else if (t.kind === KeyKind.row) add(K.row, t.id, t.row);
    else if (t.kind === KeyKind.complex)
      for (const u of [...(set.complexes[t.id]?.units ?? [])].sort((x, y) => x - y)) add(K.unit, u);
    else {
      const units = this.unitsOver(t.id);
      if (units.length) for (const u of units) add(K.unit, u);
      else add(K.molecule, t.id);
    }
    return out;
  }

  /** Every position of `T` position p of `S` lands on (the routes, read per position: `colOf`, `toAt`). */
  image(S: Sys, p: number, T: Sys): number[] {
    const { set } = this;
    const same = (a: Sys, b: Sys) => a.kind === b.kind && a.id === b.id && a.row === b.row;
    const alignmentOf = (x: Sys) => (x.kind === K.row || x.kind === K.columns ? x.id : -1);
    const [lo, hi] = this.domain(T.kind, T.id, T.row);
    const clip = (xs: number[]) => xs.filter((x) => x >= lo && x < hi);
    if (same(S, T)) return clip([p]);
    const [sa, ta] = [alignmentOf(S), alignmentOf(T)];
    const rowsOf = (a: number) => set.alignments[a]?.rows;
    /** Row r's residues at column c. */
    const residuesAt = (a: number, r: number, c: number) => {
      const out: number[] = [];
      for (let j = 0; j < (rowsOf(a)?.lengthOf(r as never) ?? 0); j++)
        if (rowsOf(a)?.colOf(r as never, j) === c) out.push(j);
      return out;
    };
    if (sa >= 0 && sa === ta) {
      const cols = S.kind === K.row ? [rowsOf(sa)?.colOf(S.row as never, p) ?? -1].filter((c) => c >= 0) : [p];
      return clip(T.kind === K.columns ? cols : cols.flatMap((c) => residuesAt(ta, T.row, c)));
    }
    /** Through a span-1 map from `from`'s side (1) or to's (-1): every element touching x. */
    const through = (m: CoordinateMap, x: number, read: number) => {
      if (read > 0) return Array.from({ length: m.span }, (_, ph) => Coords.toAt(m, x, ph)).filter((y) => y >= 0);
      const out: number[] = [];
      for (let y = m.fromStart; y < m.fromStop; y++) {
        for (let ph = 0; ph < m.span; ph++) if (Coords.toAt(m, y, ph) === x) out.push(y);
      }
      return [...new Set(out)];
    };
    if (sa >= 0) {
      const b = Mappings.bindingsOf(set, sa);
      const rows = S.kind === K.row ? [S.row] : Array.from(b?.rows ?? []);
      const entries: [number, number[]][] = [];
      for (const r of rows) {
        const k = b ? Mappings.boundRow(b, r) : -1;
        if (!b || k < 0) continue;
        const js = S.kind === K.row ? [p] : residuesAt(sa, r, p);
        for (let x = b.start[k]; x < b.start[k + 1]; x++) entries.push([x, js]);
      }
      const bb = b as RowBindings;
      const direct = entries.some(([x]) => T.kind === K.unit && bb.unit[x] === T.id);
      return norm(
        entries.flatMap(([x, js]) => {
          const u = bb.unit[x];
          const [ulo, uhi] = this.domain(K.unit, u, -1);
          const xs = js
            .flatMap((j) => through(bb.map[x], j, bb.rowIsFrom[x] ? 1 : -1))
            .filter((y) => y >= ulo && y < uhi);
          if (direct) return bb.unit[x] === T.id ? clip(xs) : [];
          return xs.flatMap((y) => this.image({ kind: K.unit, id: u, row: -1 }, y, T));
        }),
      );
    }
    const b = ta >= 0 && S.kind === K.unit ? Mappings.bindingsOf(set, ta) : null;
    const xu = b ? Mappings.boundUnit(b, S.id) : -1;
    if (b && xu >= 0) {
      const es = Array.from({ length: b.unitStart[xu + 1] - b.unitStart[xu] }, (_, i) => b.byUnit[b.unitStart[xu] + i]);
      const direct = es.some((x) => T.kind === K.row && b.rows[b.entryRow[x]] === T.row);
      return norm(
        es.flatMap((x) => {
          const r = b.rows[b.entryRow[x]];
          const len = rowsOf(ta)?.lengthOf(r as never) ?? 0;
          const js = through(b.map[x], p, b.rowIsFrom[x] ? -1 : 1).filter((j) => j >= 0 && j < len);
          if (direct) return r === T.row ? clip(js) : [];
          return js.flatMap((j) => this.image({ kind: K.row, id: ta, row: r }, j, T));
        }),
      );
    }
    const steps: PathStep[] = [];
    if (Mappings.path(set, S, T, steps) <= 0) return [];
    const m = Mappings.compose(set, S, steps).map;
    return clip(through(m, p, Coords.sameSystem(m.from, S) ? 1 : -1));
  }

  onto(e: Expression & object): OGroup[] {
    const args = (e as { args?: Record<string, Expression> }).args ?? {};
    const targets = this.targets(args["1"], args.kind as string | undefined);
    return nonEmpty(
      this.eval(args["0"]).map((g) => {
        const pos = g.pos.flatMap((k) => {
          const { kind, id, row, p } = decode(k);
          return targets.flatMap((T) => this.image({ kind, id, row }, p, T).map((x) => key(T.kind, T.id, T.row, x)));
        });
        return group(args.keep === true ? [...pos, ...g.pos] : pos);
      }),
    );
  }

  partners(e: Expression & object): OGroup[] {
    const args = (e as { args?: Record<string, Expression> }).args ?? {};
    const t = { kind: -1, id: -1, row: -1 };
    if (!this.set.byName(String(args.via), t)) throw new OracleError("seqql.ref.unknown");
    const a = this.set.annotations[t.id];
    if (t.kind !== KeyKind.annotation || !a?.pairs) throw new OracleError("seqql.ref.kind");
    const pc = a.pairs;
    const q = new Set(this.eval(args["0"]).flatMap((g) => g.pos));
    const side = (args.side as string) ?? "any";
    const end = (r: number, sd: number) => {
      const s = sd ? pc.bStart[r] : a.start[r];
      const z = sd ? pc.bStop[r] : a.stop[r];
      return this.target(
        sd ? pc.bTarget : a.slots,
        sd ? r : a.slot[r],
        s < 0 ? 0 : s,
        s < 0 ? 1e5 : Math.max(z, s + 1),
      );
    };
    const out: number[] = [];
    for (let r = 0; r < a.count; r++) {
      for (const sd of [0, 1]) {
        const other = 1 - sd;
        if (side !== "any" && side !== (other ? "b" : "a")) continue;
        if (!end(r, sd).some((k) => q.has(k))) continue;
        this.slot = { ann: a.id, item: a.item[r] };
        if (args.where !== undefined && this.test(args.where) !== true) continue;
        out.push(...end(r, other));
      }
    }
    return runsOf(norm(out));
  }
}

export interface Sys {
  readonly kind: number;
  readonly id: number;
  readonly row: number;
}

/** Sorted keys as one group per maximal run of consecutive positions in one system. */
export function runsOf(pos: readonly number[]): OGroup[] {
  const out: number[][] = [];
  for (const k of pos) {
    const last = out.at(-1);
    if (last && k === (last.at(-1) as number) + 1) last.push(k);
    else out.push([k]);
  }
  return out.map((p) => group(p));
}
