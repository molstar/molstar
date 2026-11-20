/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `(elements <items>)` (_spec/README.md "Queries", Language; locations-and-queries.md, The language): the schema
// selection, a Location. Its finest field decides how deep it goes: none below a unit, every matching unit's whole
// residues system; residue fields, those positions; an atom field, positions on the units' atoms systems. It is sugar
// over `filter` with its items' predicate (`itemsPredicate`), exactly; only speed differs: rows are grouped by their
// molecule and unit fields and each frame looks its rows up once, a row's positions are found from its `seq`, range
// or author number by index, and only the rest is tested per element (MolQL's `elementSchema`): rows plus matches,
// never rows × units.

import { type Atom, issueOf, type Node } from '@molstar/core/repr/tree';
import { checkLimit, type Ctx, type Resolve, tick } from './form-kit.js';
import { type Absent, coerce, compilePredicate, type Predicate } from './expressions.js';
import { QueryIssues } from './issues.js';
import { type Items, itemsLevel, RangeFields, readItems } from './items.js';
import { compareLabels, parseLabel } from './labels.js';
import { EmptyLocation, type Location, location, normalizeRanges } from './location.js';
import { findNumbering, labelIndex, snapLabel } from './numbering.js';
import { findProperty, levelRank, type PropertyDef } from './properties.js';
import { type Frame, FrameReads, framesOf, type Reader, stageFrame } from './structure.js';

const Unit = levelRank('unit');
const Residue = levelRank('residue');
const Atoms = levelRank('atom');

/** What a frame's elements are read with: its reads as `Readers`, or `'pending'` when one it needs has not arrived. */
export function readersOf(
  reads: FrameReads,
  uses: readonly PropertyDef[],
): ((d: PropertyDef) => Reader | undefined) | 'pending' {
  for (const d of uses) if (reads.reader(d) === 'pending') return 'pending';
  return (d) => reads.reader(d) as Reader | undefined;
}

/** One warning a property elements lacked (`seq.query.property.absent`), naming it and the count. */
export function reportAbsent(ctx: Ctx, absent: Absent) {
  for (const [name, n] of absent)
    ctx.issues.push(
      issueOf(
        QueryIssues,
        'seq.query.property.absent',
        `${n} element${n === 1 ? '' : 's'} lack ${name}: none of them matches`,
      ),
    );
}

/** Says once that a frame waits for a property's values. */
export function reportPending(ctx: Ctx, f: Frame, said: { said?: boolean }) {
  stageFrame(ctx, f);
  if (said.said) return;
  said.said = true;
  ctx.issues.push(
    issueOf(QueryIssues, 'seq.query.property.pending', 'A property’s values have not arrived', {
      address: f.system.address,
    }),
  );
}

/** A row's coarse fields (molecule and unit) and the rest, the rows indexed by their coarse values. */
interface RowIndex {
  /** By the set of coarse fields a row gives (its mask), the rows by those fields' values. */
  readonly groups: readonly { readonly fields: readonly PropertyDef[]; readonly rows: Map<string, number[]> }[];
}

const keyOf = (vs: readonly Atom[]) => JSON.stringify(vs.map((v) => [typeof v, v]));

function indexRows(items: Items, defs: ReadonlyMap<string, PropertyDef>): RowIndex {
  const coarse = items.fields.filter((f) => {
    const d = defs.get(f);
    return d && levelRank(d.level) <= Unit;
  });
  const groups = new Map<string, { fields: PropertyDef[]; rows: Map<string, number[]> }>();
  for (let r = 0; r < items.rows; r++) {
    const given = coarse.filter((f) => items.columns[f][r] !== null);
    const mask = given.join(' ');
    let g = groups.get(mask);
    if (!g) groups.set(mask, (g = { fields: given.map((f) => defs.get(f)!), rows: new Map() }));
    const key = keyOf(given.map((f) => items.columns[f][r]));
    let rows = g.rows.get(key);
    if (!rows) g.rows.set(key, (rows = []));
    rows.push(r);
  }
  return { groups: [...groups.values()] };
}

/** The rows a frame's coarse values select (a word set, `identifier`, by each of its words). */
function rowsFor(index: RowIndex, reads: FrameReads, absent: Absent): number[] {
  const out: number[] = [];
  for (const g of index.groups) {
    let keys: Atom[][] = [[]];
    let lacking = false;
    for (const d of g.fields) {
      const rd = reads.reader(d);
      const v = rd === 'pending' ? undefined : rd?.(0, -1);
      if (v === undefined) {
        absent.set(d.name, (absent.get(d.name) ?? 0) + 1);
        lacking = true;
        break;
      }
      const vs = (Array.isArray(v) ? v : [v]).map((x) => coerce(x as Atom, d.type) ?? null);
      keys = keys.flatMap((k) => vs.map((x) => [...k, x]));
    }
    if (lacking) continue;
    for (const k of keys) for (const r of g.rows.get(keyOf(k)) ?? []) out.push(r);
  }
  return [...new Set(out)].sort((a, b) => a - b);
}

/** A frame's author numbers: each integer's positions, built once a frame. */
function authIndex(f: Frame): Map<number, number[]> | undefined {
  const n = findNumbering(f.system, 'auth');
  if (!n) return undefined;
  const out = new Map<number, number[]>();
  for (let i = 0; i < f.system.length; i++) {
    const l = n.numbering.labelAt(i - n.shift);
    const int = l === undefined ? undefined : parseLabel(l).int;
    if (int === undefined) continue;
    let ps = out.get(int);
    if (!ps) out.set(int, (ps = []));
    ps.push(i);
  }
  return out;
}

/** The positions a row names by index (its `seq`, range, author number or author range), or undefined: every one. */
function rowPositions(
  ctx: Ctx,
  f: Frame,
  items: Items,
  r: number,
  auth: () => Map<number, number[]> | undefined,
  absent: Absent,
): number[] | undefined {
  const at = (k: string) => items.columns[k]?.[r] ?? null;
  const len = f.system.length;
  const start = f.system.parent?.start ?? 0;
  const seq = at('seq');
  const beg = at('beg');
  const end = at('end');
  let out: number[] | undefined;
  const meet = (ranges: number[]) => {
    out = out ? intersect(out, ranges) : ranges;
  };
  if (typeof seq === 'number') {
    const i = seq - 1 - start;
    meet(i >= 0 && i < len ? [i, i + 1] : []);
  }
  if (beg !== null || end !== null) {
    const b = Math.max(0, typeof beg === 'number' ? beg - 1 - start : 0);
    const e = Math.min(len, typeof end === 'number' ? end - start : len);
    meet(b < e ? [b, e] : []);
  }
  const authSeq = at('auth-seq');
  if (typeof authSeq === 'number') {
    const index = auth();
    if (!index) absent.set('auth-seq', (absent.get('auth-seq') ?? 0) + 1);
    const ps = index?.get(authSeq) ?? [];
    meet(ps.flatMap((p) => [p, p + 1]));
  }
  const ab = at('auth-beg');
  const ae = at('auth-end');
  if (ab !== null || ae !== null) {
    const n = findNumbering(f.system, 'auth');
    if (!n) meet([]);
    else {
      const end = (l: Atom, side: 'start' | 'end') => {
        const s = `${l}`;
        const i = labelIndex(n, f.system, s);
        return i >= 0 ? i : snapLabel(n, f.system, s, side);
      };
      const b = ab === null ? 0 : end(ab, 'start');
      const e = ae === null ? len - 1 : end(ae, 'end');
      // in the numbering's order, as a range's ends (EC-1 to EC-8): ends that snap past each other hold nothing
      const ordered = ab === null || ae === null || compareLabels(`${ab}`, `${ae}`) <= 0;
      meet(b >= 0 && e >= 0 && b <= e && ordered ? [b, e + 1] : []);
    }
  }
  if (out && f.residues) out = intersect(out, Array.from(f.residues));
  return out;
}

function intersect(a: readonly number[], b: readonly number[]): number[] {
  const out: number[] = [];
  let i = 0;
  let j = 0;
  while (i < a.length && j < b.length) {
    const s = Math.max(a[i], b[j]);
    const e = Math.min(a[i + 1], b[j + 1]);
    if (s < e) out.push(s, e);
    if (a[i + 1] < b[j + 1]) i += 2;
    else j += 2;
  }
  return out;
}

/** A row's fine fields (residue, atom; the range fields aside) as a predicate, compiled once a row. */
function finePredicate(ctx: Ctx, items: Items, r: number, defs: ReadonlyMap<string, PropertyDef>, absent: Absent) {
  const tests: Node[] = [];
  const ranks: number[] = [];
  for (const f of items.fields) {
    const d = defs.get(f);
    const v = items.columns[f][r];
    if (!d || v === null || levelRank(d.level) <= Unit || RangeFields[f]) continue;
    // a seq or author number named by index is not tested again
    if (f === 'seq' || f === 'auth-seq') continue;
    tests.push({ kind: '=', args: [{ kind: 'prop', args: [f] }, v] });
    ranks.push(levelRank(d.level));
  }
  // residue fields before atom fields: `and` stops at the first false, so an atom field lacks only where they matched
  const order = tests.map((_, i) => i).sort((x, y) => ranks[x] - ranks[y]);
  return tests.length
    ? compilePredicate(ctx.ext, { kind: 'and', args: order.map((i) => tests[i]) }, absent)
    : undefined;
}

/** The elements a frame gives: per row its positions by index, then its fine test per residue or atom. */
function walkFrame(
  ctx: Ctx,
  reads: FrameReads,
  items: Items,
  rows: readonly number[],
  level: number,
  preds: (r: number) => Predicate | undefined,
  absent: Absent,
  pending: { said?: boolean },
): { system: string; ranges: number[] } | undefined {
  const f = reads.frame;
  let authCache: Map<number, number[]> | undefined | null = null;
  const auth = () => (authCache === null ? (authCache = authIndex(f)) : authCache);
  // a frame without atoms (a polymer's unit) is still tested by residue: its atom fields then lack, counted only
  // where its residue fields matched, and it gives nothing
  const atoms = level >= Atoms ? reads.atoms() : undefined;
  const lacking = level >= Atoms && !atoms;
  const all = f.residues ? Array.from(f.residues) : [0, f.system.length];
  const out: number[] = [];
  for (const r of rows) {
    const positions = level >= Residue ? (rowPositions(ctx, f, items, r, auth, absent) ?? all) : all;
    const pred = preds(r);
    const rs = pred ? readersOf(reads, pred.uses) : undefined;
    if (rs === 'pending') {
      reportPending(ctx, f, pending);
      return undefined;
    }
    for (let k = 0; k < positions.length; k += 2)
      for (let p = positions[k]; p < positions[k + 1]; p++) {
        tick(ctx);
        if (!atoms) {
          if ((!pred || pred.test(rs!, p, -1)) && !lacking) out.push(p, p + 1);
          continue;
        }
        const [s, e] = atoms.atomsOf(p);
        for (let a = s; a < e; a++) {
          if (f.atoms && !holds(f.atoms, a)) continue;
          if (!pred || pred.test(rs!, p, a)) out.push(a, a + 1);
        }
      }
  }
  if (!out.length) return undefined;
  return { system: atoms ? atoms.address : f.system.address, ranges: Array.from(normalizeRanges(out)) };
}

function holds(ranges: Int32Array, i: number): boolean {
  for (let k = 0; k < ranges.length; k += 2) if (i >= ranges[k] && i < ranges[k + 1]) return true;
  return false;
}

/** The property of each field of items (none for a range field). */
function fieldDefs(ctx: Ctx, items: Items): Map<string, PropertyDef> {
  const out = new Map<string, PropertyDef>();
  for (const f of items.fields) {
    const d = RangeFields[f] ? undefined : findProperty(ctx.ext.properties, f);
    if (d && 'def' in d) out.set(f, d.def);
  }
  return out;
}

/** `elements` over frames: the Location its items name. */
export function elementsOver(ctx: Ctx, frames: readonly Frame[], items: Items): Location {
  checkLimit(ctx, items.rows);
  const defs = fieldDefs(ctx, items);
  const level = itemsLevel(ctx.ext.properties, items);
  const index = indexRows(items, defs);
  const absent: Absent = new Map();
  const pending: { said?: boolean } = {};
  const preds = new Map<number, Predicate | undefined>();
  const predOf = (r: number) => {
    if (!preds.has(r)) preds.set(r, finePredicate(ctx, items, r, defs, absent));
    return preds.get(r);
  };
  const parts: { system: string; ranges: number[] }[] = [];
  for (const f of frames) {
    tick(ctx);
    const reads = new FrameReads(ctx, f);
    const rows = rowsFor(index, reads, absent);
    if (!rows.length) continue;
    const part = walkFrame(ctx, reads, items, rows, level, predOf, absent, pending);
    if (part) parts.push(part);
  }
  reportAbsent(ctx, absent);
  return parts.length ? location({ parts }) : EmptyLocation;
}

export const elements: Resolve = (ctx, form, prefix) => {
  const items = readItems(ctx.ext.properties, form);
  return { type: 'location', location: elementsOver(ctx, framesOf(ctx, prefix, 'elements'), items) };
};
