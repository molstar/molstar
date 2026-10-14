/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Queries (seqspec:libraries#interactions) over a kernel pairs annotation, built on `Pairs.partnersOf` and
// `Pairs.touching`.
// Each row's kind is read once per annotation (a runtime property); after that the queries allocate nothing per
// call: module-level visitors and scratch (grown when too small, kept), a fresh scratch only when re-entered.
// They take the set (an addition to the spec): `category` is a NameId, read through the set's strings.

import {
  type Annotation,
  CustomData,
  type CustomProperty,
  MoleculeKind,
  Pairs,
  type SeqIndex,
  type SeqSet,
  type UnitId,
} from "@molstar/seq-core/model";
import { Kind, type KindMask, kindByName } from "./vocabulary.js";

/** Per row of a pairs annotation, its kind, read from `category` once. */
const RowKinds: CustomProperty<Annotation, Uint8Array> = {
  name: "interactions:kinds",
  compute(a, set) {
    const out = new Uint8Array(a.count).fill(Kind.other);
    for (let r = 0; a.category && r < a.count; r++) out[r] = kindByName(set.strings.get(a.category[r]));
    return out;
  },
};
const rowKinds = (set: SeqSet, a: Annotation) => CustomData.get(a.custom, RowKinds, a, set) as Uint8Array;

function kindOf(set: SeqSet, pairs: Annotation, row: number): number {
  return rowKinds(set, pairs)[row];
}

// contactsOf: (unit, index + 1, kind) packed into one number per contact, sorted, runs merged into masks.
const UNIT = 2 ** 24;
interface Contacts {
  keys: Float64Array;
  n: number;
  kinds: Uint8Array;
  mask: KindMask;
}
const contacts: Contacts = { keys: new Float64Array(256), n: 0, kinds: new Uint8Array(0), mask: 0 };
let contactsDepth = 0;

function collect(unit: UnitId, index: SeqIndex, row: number, s: Contacts) {
  const kind = s.kinds[row];
  if (!(s.mask & (1 << kind))) return;
  if (s.n === s.keys.length) {
    const grown = new Float64Array(s.n * 2);
    grown.set(s.keys);
    s.keys = grown;
  }
  s.keys[s.n++] = (unit * UNIT + index + 1) * 16 + kind;
}

/** Shell sort of `a[0, n)` in place (a typed array's own sort would need a subarray view). */
function sortPrefix(a: Float64Array, n: number) {
  for (let gap = n >> 1; gap > 0; gap >>= 1) {
    for (let i = gap; i < n; i++) {
      const x = a[i];
      let j = i;
      for (; j >= gap && a[j - gap] > x; j -= gap) a[j] = a[j - gap];
      a[j] = x;
    }
  }
}

/**
 * Residues that contact a ligand unit, with the kinds per residue (only kinds in `kinds`), by unit then position;
 * index -1: a whole-unit endpoint.
 */
function contactsOf<C>(
  set: SeqSet,
  pairs: Annotation,
  ligand: UnitId,
  kinds: KindMask,
  visit: (unit: UnitId, index: SeqIndex, kinds: KindMask, ctx: C) => void,
  ctx: C,
): void {
  const s = contactsDepth++ === 0 ? contacts : { keys: new Float64Array(64), n: 0, kinds: contacts.kinds, mask: 0 };
  try {
    s.n = 0;
    s.kinds = rowKinds(set, pairs);
    s.mask = kinds;
    Pairs.partnersOf(pairs, ligand, collect, s);
    sortPrefix(s.keys, s.n);
    for (let i = 0; i < s.n; ) {
      const at = Math.floor(s.keys[i] / 16);
      let mask = 0;
      for (; i < s.n && Math.floor(s.keys[i] / 16) === at; i++) mask |= 1 << (s.keys[i] % 16);
      visit(Math.floor(at / UNIT) as UnitId, ((at % UNIT) - 1) as SeqIndex, mask, ctx);
    }
  } finally {
    contactsDepth--;
  }
}

// atomsContacting: the other endpoint of each pair touching the residue, on ligand units.
interface Atoms {
  set: SeqSet;
  pairs: Annotation;
  kinds: Uint8Array;
  visit: (ligand: UnitId, atom: SeqIndex, kind: number, ctx: unknown) => void;
  ctx: unknown;
  /** The pair row being visited, and the other endpoint's atoms: [start, stop). */
  row: number;
  start: number;
  stop: number;
}
const atoms = {} as Atoms;
let atomsDepth = 0;

/** One unit of the other endpoint (`Pairs.endUnits`): its ligand atoms, or the whole unit (atom -1). */
function emitUnit(unit: UnitId, s: Atoms) {
  const { units, molecules } = s.set;
  if (molecules[units.molecule[unit]]?.kind !== MoleculeKind.ligand) return;
  const kind = s.kinds[s.row];
  if (s.start < 0) s.visit(unit, -1 as SeqIndex, kind, s.ctx);
  else for (let i = s.start; i < s.stop; i++) s.visit(unit, i as SeqIndex, kind, s.ctx);
}

function touched(row: number, side: number, s: Atoms) {
  const a = s.pairs;
  const pc = a.pairs;
  if (!pc) return;
  const other = side ? 0 : 1;
  s.row = row;
  s.start = other ? pc.bStart[row] : a.start[row];
  s.stop = Math.max(other ? pc.bStop[row] : a.stop[row], s.start + 1);
  Pairs.endUnits(a, row, other, emitUnit, s);
}

/** Ligand atoms that contact one residue, with the kind of each pair (an atom -1: a whole-unit endpoint). */
function atomsContacting<C>(
  set: SeqSet,
  pairs: Annotation,
  unit: UnitId,
  index: SeqIndex,
  visit: (ligand: UnitId, atom: SeqIndex, kind: number, ctx: C) => void,
  ctx: C,
): void {
  const s = atomsDepth++ === 0 ? atoms : ({} as Atoms);
  try {
    s.set = set;
    s.pairs = pairs;
    s.kinds = rowKinds(set, pairs);
    s.visit = visit as Atoms["visit"];
    s.ctx = ctx;
    Pairs.touching(pairs, unit, index, touched, s);
  } finally {
    atomsDepth--;
    s.ctx = s.set = s.pairs = null as never; // the shared state keeps no set alive
  }
}

export const Interactions = { kindOf, contactsOf, atomsContacting };
