/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The structure forms beside `elements` (_spec/README.md "Queries", Language): the generator `groups` (Keys: one a
// group, by atom, residue, unit or molecule, each level's test reading that level's properties), the generic `filter
// :test` (the prefix's type, a subset: per position of a Location, per key of Keys, per unit of nodes) and the
// modifiers `expand :to` (atoms to their residues, residues to their units or their molecule's units) and `linked :via
// :type` (the prefix and what a `pairs` annotation's rows join to it, never a projection through them: links are never
// maps). Every form is monadic over its prefix's frames.

import { type Atom, printAddress, type Selected, type Value } from '@molstar/core/repr/tree';
import { readersOf, reportAbsent, reportPending } from './elements-form.js';
import { evaluate } from './evaluate.js';
import { type Absent, compilePredicate, type Predicate } from './expressions.js';
import {
  checkLimit,
  checkParams,
  type Ctx,
  fail,
  param,
  prefixNodes,
  type Resolve,
  required,
  stage,
  tick,
} from './form-kit.js';
import { keysExtent } from './key-extent.js';
import { type Key, keySet } from './keys.js';
import { EmptyLocation, type Location, location, overlapsRanges, unionLocations } from './location.js';
import { levelRank, type PropertyLevel, PropertyLevels } from './properties.js';
import { fromNamed, type QueryResult } from './result.js';
import { type Frame, FrameReads, framesOf } from './structure.js';

const Residue = levelRank('residue');
const Atoms = levelRank('atom');

/** The positions of a frame's elements at a level: its residues, or its atoms (each with its residue). */
function* elementsOf(reads: FrameReads, atomLevel: boolean): Generator<[number, number]> {
  const f = reads.frame;
  const residues = f.residues ?? Int32Array.of(0, f.system.length);
  const atoms = atomLevel ? reads.atoms() : undefined;
  for (let k = 0; k < residues.length; k += 2)
    for (let r = residues[k]; r < residues[k + 1]; r++) {
      if (!atoms) {
        yield [r, -1];
        continue;
      }
      const [s, e] = atoms.atomsOf(r);
      for (let a = s; a < e; a++) if (!f.atoms || overlapsRanges(f.atoms, a, a + 1)) yield [r, a];
    }
}

/** A Location's frames tested element by element: those whose test is true. */
export function filterFrames(ctx: Ctx, frames: readonly Frame[], pred: Predicate, absent: Absent): Location {
  const pending: { said?: boolean } = {};
  const parts: { system: string; ranges: number[] }[] = [];
  for (const f of frames) {
    const reads = new FrameReads(ctx, f);
    // a frame is read at its own level: an atoms part's atoms, else its residues
    const atomLevel = f.atoms !== undefined;
    const rs = readersOf(reads, pred.uses);
    if (rs === 'pending') {
      reportPending(ctx, f, pending);
      continue;
    }
    const out: number[] = [];
    for (const [r, a] of elementsOf(reads, atomLevel)) {
      tick(ctx);
      if (pred.test(rs, r, a)) out.push(atomLevel ? a : r, (atomLevel ? a : r) + 1);
    }
    if (out.length) parts.push({ system: atomLevel ? reads.atoms()!.address : f.system.address, ranges: out });
  }
  return parts.length ? location({ parts }) : EmptyLocation;
}

/** The finest level a prefix's elements are: atoms when any part is on an atoms system. */
const frameLevel = (frames: readonly Frame[]) => (frames.some((f) => f.atoms) ? Atoms : Residue);

const filter: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['test']);
  const test = required(form, 'test');
  const absent: Absent = new Map();
  const r = prefix.kind === 'nodes' ? undefined : fromNamed(prefix);
  let out: QueryResult;
  if (r?.type === 'location' || r?.type === 'keys') {
    const frames = framesOf(ctx, prefix, 'filter');
    // an empty prefix has no level to read past: nothing to test, nothing refused
    const pred = compilePredicate(
      ctx.ext,
      test,
      absent,
      frames.length ? frameLevel(frames) : undefined,
      'A filter’s test',
    );
    const loc = filterFrames(ctx, frames, pred, absent);
    out =
      r.type === 'location'
        ? { type: 'location', location: loc }
        : { type: 'keys', keys: keysWithin(ctx, r.keys.keys, loc) };
  } else {
    // nodes: a unit (or a molecule) passes when its test does, read at the unit level
    const pred = compilePredicate(ctx.ext, test, absent, levelRank('unit'), 'A filter’s test over nodes');
    const pending: { said?: boolean } = {};
    const kept: Selected[] = [];
    for (const m of prefixNodes(prefix, 'filter')) {
      const frames = framesOf(ctx, { kind: 'nodes', matches: [m] }, 'filter');
      if (frames.some((f) => passes(ctx, f, pred, pending))) kept.push(m);
    }
    checkLimit(ctx, kept.length);
    out = { type: 'nodes', nodes: kept };
  }
  reportAbsent(ctx, absent);
  return out;
};

/**
 * Whether a frame's unit-level test passes; a property it reads pending stages the frame and says so (EC-169), the
 * node then not held until its values arrive, never silently dropped.
 */
function passes(ctx: Ctx, f: Frame, pred: Predicate, pending: { said?: boolean }): boolean {
  tick(ctx);
  const rs = readersOf(new FrameReads(ctx, f), pred.uses);
  if (rs === 'pending') {
    reportPending(ctx, f, pending);
    return false;
  }
  return pred.test(rs, 0, -1);
}

/** The keys whose extent the location holds. */
function keysWithin(ctx: Ctx, keys: readonly Key[], loc: Location): ReturnType<typeof keySet> {
  return keySet(
    keys.filter((k) => {
      const e = keysExtent(ctx, keySet([k]));
      return e.parts.some((p) => {
        const at = loc.parts.find((q) => q.system === p.system);
        return at && p.ranges.length > 0 && overlapsRanges(at.ranges, p.ranges[0], p.ranges[p.ranges.length - 1]);
      });
    }),
  );
}

// ---- groups ----

const GroupLevels: readonly PropertyLevel[] = ['atom', 'residue', 'unit', 'molecule'];

const groups: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['molecule-test', 'unit-test', 'residue-test', 'atom-test', 'by']);
  const by = param(form, 'by') ?? 'residue';
  if (typeof by !== 'string' || !GroupLevels.includes(by as PropertyLevel))
    fail('seq.query.form.param', 'A form of kind groups takes :by atom, residue, unit or molecule');
  const absent: Absent = new Map();
  const tests = PropertyLevels.map((level) => {
    const t = param(form, `${level}-test`);
    return t === undefined
      ? undefined
      : compilePredicate(ctx.ext, t, absent, levelRank(level), `A form of kind groups’ :${level}-test`);
  });
  const byRank = levelRank(by as PropertyLevel);
  const atomLevel = byRank === Atoms || tests[Atoms] !== undefined;
  const fineLevel = Math.max(byRank, tests[Atoms] ? Atoms : tests[Residue] ? Residue : -1);
  const pending: { said?: boolean } = {};
  const out: Key[] = [];
  const seen = new Set<string>();
  const push = (owner: string, key: string | number) => {
    const id = `${owner} ${typeof key} ${key}`;
    if (seen.has(id)) return;
    seen.add(id);
    out.push({ owner, key });
  };
  for (const f of framesOf(ctx, prefix, 'groups')) {
    tick(ctx);
    const reads = new FrameReads(ctx, f);
    const uses = tests.flatMap((t) => t?.uses ?? []);
    const rs = readersOf(reads, uses);
    if (rs === 'pending') {
      reportPending(ctx, f, pending);
      continue;
    }
    if (tests[0] && !tests[0].test(rs, 0, -1)) continue;
    if (tests[1] && !tests[1].test(rs, 0, -1)) continue;
    const coarse = () => {
      if (byRank === levelRank('unit')) {
        if (f.complex && f.unit) push(f.complex, f.unit.id);
      } else if (f.complex && f.unit?.entity !== undefined) push(f.complex, f.unit.entity);
    };
    if (fineLevel < Residue) {
      coarse();
      continue;
    }
    const atoms = atomLevel ? reads.atoms() : undefined;
    if (atomLevel && !atoms) {
      absent.set('atoms', (absent.get('atoms') ?? 0) + 1);
      continue;
    }
    let any = false;
    for (const [r, a] of elementsOf(reads, atomLevel)) {
      tick(ctx);
      if (tests[Residue] && !tests[Residue].test(rs, r, a)) continue;
      if (tests[Atoms] && !tests[Atoms].test(rs, r, a)) continue;
      any = true;
      if (byRank === Atoms) push(atoms!.address, a);
      else if (byRank === Residue) push(f.system.address, r);
      else break;
    }
    if (any && byRank < Residue) coarse();
  }
  reportAbsent(ctx, absent);
  const keys = keySet(out);
  checkLimit(ctx, keys.keys.length);
  return { type: 'keys', keys };
};

// ---- modifiers ----

const ExpandLevels = ['residue', 'unit', 'molecule'];

const expand: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['to']);
  const to = required(form, 'to');
  if (typeof to !== 'string' || !ExpandLevels.includes(to))
    fail('seq.query.form.param', 'A form of kind expand takes :to residue, unit or molecule');
  if (prefix.kind === 'nodes') fail('seq.query.form.prefix', 'A form of kind expand applies to a location or keys');
  const parts: { system: string; ranges: number[] }[] = [];
  for (const f of framesOf(ctx, prefix, 'expand')) {
    tick(ctx);
    if (to === 'residue') {
      parts.push({ system: f.system.address, ranges: Array.from(f.residues ?? [0, f.system.length]) });
      continue;
    }
    const units =
      to === 'molecule' && f.complex ? (ctx.model.units?.(f.complex) ?? []).filter((u) => u.of === f.of) : undefined;
    if (!units?.length) {
      parts.push({ system: f.system.address, ranges: [0, f.system.length] });
      continue;
    }
    for (const u of units) {
      const s = ctx.model.system(u.address);
      if (s) parts.push({ system: u.address, ranges: [0, s.length] });
    }
  }
  return { type: 'location', location: location({ parts }) };
};

/** A param of words: one, or a list's. */
function wordsOf(v: Value | undefined, kind: string, key: string): ReadonlySet<Atom> | undefined {
  if (v === undefined) return undefined;
  const list = (Array.isArray(v) ? v : [v]) as Atom[];
  if (!list.every((x) => typeof x === 'string' || typeof x === 'number'))
    fail('seq.query.form.param', `A form of kind ${kind} takes words for :${key}`);
  return new Set(list.map(String));
}

const linked: Resolve = (ctx, form, prefix) => {
  checkParams(form, ['via', 'type']);
  const types = wordsOf(param(form, 'type'), 'linked', 'type');
  const via = evaluate(ctx, required(form, 'via'), prefix);
  if (via.type !== 'nodes')
    return fail('seq.query.form.type', 'A form of kind linked follows the pairs nodes :via names');
  if (prefix.kind === 'nodes') fail('seq.query.form.prefix', 'A form of kind linked applies to a location or keys');
  // the prefix's residues: an atoms part read as the residues holding its atoms
  const frames = framesOf(ctx, prefix, 'linked');
  let held: Location = location({
    parts: frames.map((f) => ({ system: f.system.address, ranges: Array.from(f.residues ?? [0, f.system.length]) })),
  });
  const joined: { system: string; ranges: number[] }[] = [];
  for (const m of via.nodes) {
    const address = printAddress(m.address);
    const a = ctx.model.annotation(address);
    if (!a) return fail('seq.query.form.prefix', 'A form of kind linked follows a pairs annotation', { address });
    const read = a.items();
    stage(ctx, read.pending);
    const byKey = new Map<string | number, { system: string; start: number; end: number }[]>();
    for (const it of read.items) {
      tick(ctx);
      const t = it.fields?.type;
      if (types && !(typeof t === 'string' && types.has(t))) continue;
      let rows = byKey.get(it.key);
      if (!rows) byKey.set(it.key, (rows = []));
      rows.push({ system: it.system ?? a.target, start: it.start, end: it.end });
    }
    for (const rows of byKey.values()) {
      const touches = rows.some((row) => {
        const part = held.parts.find((p) => p.system === row.system);
        return part && overlapsRanges(part.ranges, row.start, row.end);
      });
      if (touches) for (const row of rows) joined.push({ system: row.system, ranges: [row.start, row.end] });
    }
  }
  held = unionLocations(held, location({ parts: joined }));
  return { type: 'location', location: held };
};

/** The structure forms beside `elements`, by kind. */
export const StructureForms = { groups, filter, expand, linked } as const;
