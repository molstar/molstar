/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Molecules, complexes and units resolved (_spec/README.md "Resolve", Systems; Molecules and units): a molecule's
// residues and atoms as two systems joined by a group map that applies at a path's ends only, a unit as a system on
// its molecule's positions over its range (its position i the molecule's `start + i`, so its steps up and down are the
// paths' free unit steps), copies derived from equal `of` and range, a sequence's modified index. Pure functions of a
// node and what it names, so the static model and the overlay build the same entries.

import { CategoryColumn, isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import { type SequenceAlphabet, standardComponent } from '@molstar/core/util/sequence-constants';
import { inlineSeq, SeqRuntime, seqToValue, type Tables } from '@molstar/seqx-core/data';
import { unitRange } from '@molstar/seqx-core/lang';
import { ImplicitLinear } from './linear.js';
import type { AccessPlan, SystemEntry } from './model.js';
import type { Edge } from './paths.js';

/** A molecule's atoms system's address: its residues system's, then `(atoms)`. */
export function atomsSystemOf(system: string): string {
  return `${system}/(atoms)`;
}

/** A molecule (`molecule` node): its class, name, residues and identifiers; its atoms when it has them. */
export interface MoleculeEntry {
  readonly address: string;
  /** Its class: `ligand` when the node names none. */
  readonly class: string;
  readonly name?: string;
  /** A component code a residue (undefined: the molecule is a source's, its residues not held). */
  readonly residues?: readonly string[];
  readonly identifiers: readonly { readonly namespace: string; readonly value: string }[];
  readonly atoms?: AtomsEntry;
}

/**
 * An atoms system and the group map onto its residues system (`of`): residue `r` owns the atoms
 * `[offsets[r], offsets[r + 1])`, contiguous by construction (`atoms.residue` non-decreasing). Never a path step: it
 * applies at a path's ends (`collapseAtoms`, `expandAtoms`).
 */
export interface AtomsEntry {
  readonly address: string;
  /** The residues system it groups onto: a molecule, or a unit of one. */
  readonly of: string;
  /** `n + 1` offsets for `n` residues, from 0. */
  readonly offsets: Int32Array;
  /** Each atom's name and element (the molecule's atoms from `first` on). */
  readonly names: { get(i: number): string | null };
  readonly elements: { get(i: number): string | null };
  /** Its first atom among its molecule's (0 for a molecule's own; a unit's from its range's first residue). */
  readonly first: number;
}

/** A sequence's modified index (_spec/README.md "Resolve", Molecules and units): its sparse rows, sorted. */
export interface ModifiedIndex {
  /** Each row's position, in order; rows at one position are alternatives, the first giving the letter. */
  readonly position: Int32Array;
  readonly component: readonly string[];
  /** The positions each row's component covers (1 but for one standing for several parents). */
  readonly length: Int32Array;
}

/** A unit of a complex: a system on its molecule's positions `[start, end)`, counted from 0. */
export interface UnitEntry {
  readonly address: string;
  readonly complex: string;
  /** Its id in its complex: what links' `a_unit` and `b_unit` name. */
  readonly id: string;
  /** The sequence or molecule it instantiates. */
  readonly of: string;
  readonly start: number;
  /** Undefined while its molecule's length is unknown (a source's) and it names no range. */
  readonly end?: number;
  /** `polymer` for a unit of a sequence, else its molecule's class. */
  readonly class: string;
  readonly chain?: string;
  readonly authChain?: string;
  readonly operator?: string;
  readonly label?: string;
  readonly group?: string;
  /** Its observed runs on its own positions, `[start, end, …]`; undefined when all were observed. */
  readonly observed?: Int32Array;
  /** Its copy group's key: units of one complex with equal `of` and range share it (EC-65), whatever their order. */
  readonly copy: string;
}

/** A complex: its units in order, its copy groups and the addresses of its links (`pairs` children). */
export interface ComplexEntry {
  readonly address: string;
  readonly name?: string;
  readonly units: readonly UnitEntry[];
  /** Copy groups by key, each its units' addresses in the complex's order; a group of one is a unit with no copy. */
  readonly copies: ReadonlyMap<string, readonly string[]>;
  readonly links: readonly string[];
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

const str = (v: Value | undefined): string | undefined => (typeof v === 'string' ? v : undefined);

type Strings = { get(i: number): string | null; readonly length: number };

/** A column of strings as normalized (a category column), or as a plain list (a tree read without normalizing). */
function strings(v: Value | undefined): Strings | undefined {
  if (v instanceof CategoryColumn) return v;
  if (Array.isArray(v)) return { get: (i) => (typeof v[i] === 'string' ? (v[i] as string) : null), length: v.length };
  return undefined;
}

const ints = (v: Value | undefined): ArrayLike<number> | undefined =>
  isColumn(v) && !(v instanceof CategoryColumn)
    ? (v as ArrayLike<number>)
    : Array.isArray(v)
      ? (v as number[])
      : undefined;

/** Inline tables as a `seq` source of one tile over the whole system, named by its address. */
function inlinePlan(address: string, system: string, length: number | undefined, tables: Tables): AccessPlan {
  const value = seqToValue(inlineSeq({ name: system, length, cyclic: false }, tables));
  return { kind: 'inline', source: { address, format: 'seq', value }, runtime: SeqRuntime, system, shape: 'letters' };
}

/** A molecule's entry and its systems: its residues (its address) and, with atoms, `<molecule>/(atoms)`. */
export function moleculeEntry(
  address: string,
  node: Node,
  source: string | undefined,
  known?: number,
): { molecule: MoleculeEntry; systems: SystemEntry[] } {
  const res = strings(p(node, 'residues'));
  const residues = res ? Array.from({ length: res.length }, (_, i) => res.get(i) ?? '') : undefined;
  const name = str(p(node, 'name'));
  const ids = p(node, 'identifiers') as Record<string, Value> | undefined;
  const ns = strings(ids?.namespace);
  const vs = ids?.value;
  const identifiers = ns
    ? Array.from({ length: ns.length }, (_, i) => ({
        namespace: ns.get(i) ?? '',
        value: Array.isArray(vs) ? String(vs[i] ?? '') : '',
      }))
    : [];
  const length = residues ? residues.length : known;
  const label = name ?? node.id ?? address;
  const system: SystemEntry = {
    address,
    kind: 'molecule',
    ...(length !== undefined ? { length } : {}),
    cyclic: false,
    name: label,
    numberings: ['index', 'ordinal'],
    ...(source ? { source } : {}),
    access: source
      ? { kind: 'source', source, system: label, shape: 'letters' }
      : inlinePlan(address, address, length, {}),
  };
  const systems = [system];
  let atoms: AtomsEntry | undefined;
  const a = p(node, 'atoms') as Record<string, Value> | undefined;
  const residue = ints(a?.residue);
  const names = strings(a?.name);
  const elements = strings(a?.element);
  if (residues && residue && names && elements) {
    const offsets = new Int32Array(residues.length + 1);
    // residue r's atoms start at the first atom of a residue ≥ r (non-decreasing, checked)
    let k = 0;
    for (let r = 0; r <= residues.length; r++) {
      while (k < residue.length && residue[k] < r) k++;
      offsets[r] = k;
    }
    atoms = { address: atomsSystemOf(address), of: address, offsets, names, elements, first: 0 };
    systems.push(atomsSystem(atoms, label));
  }
  const molecule: MoleculeEntry = {
    address,
    class: str(p(node, 'class')) ?? 'ligand',
    ...(name !== undefined ? { name } : {}),
    ...(residues ? { residues } : {}),
    identifiers,
    ...(atoms ? { atoms } : {}),
  };
  return { molecule, systems };
}

function atomsSystem(atoms: AtomsEntry, name: string): SystemEntry {
  const length = atoms.offsets[atoms.offsets.length - 1];
  return {
    address: atoms.address,
    kind: 'atoms',
    length,
    cyclic: false,
    name: `${name} atoms`,
    numberings: ['index', 'ordinal'],
    access: inlinePlan(atoms.address, atoms.address, length, {}),
  };
}

/** A sequence's modified index from its `modified` table, or undefined without one. */
export function modifiedIndex(sequence: Node): ModifiedIndex | undefined {
  const m = p(sequence, 'modified') as Record<string, Value> | undefined;
  const pos = ints(m?.position);
  const comp = strings(m?.component);
  if (!pos || !comp) return undefined;
  const len = ints(m?.length);
  const n = Math.min(pos.length, comp.length);
  return {
    position: Int32Array.from({ length: n }, (_, i) => pos[i]),
    component: Array.from({ length: n }, (_, i) => comp.get(i) ?? ''),
    length: Int32Array.from({ length: n }, (_, i) => Math.max(1, len ? len[i] : 1)),
  };
}

/** The last row whose position is at most `p`, or −1. */
function lastAtOrBefore(position: Int32Array, p: number): number {
  let lo = 0;
  let hi = position.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (position[mid] <= p) lo = mid + 1;
    else hi = mid;
  }
  return lo - 1;
}

/**
 * The non-standard components at a position: every row covering it (alternatives in the source's order), or none
 * (the position has its letter's standard component).
 */
export function modifiedAt(index: ModifiedIndex, p: number): readonly string[] {
  let i = lastAtOrBefore(index.position, p);
  if (i < 0) return [];
  // the rows at that position, from the first (rows are in order; a length stops before the next row's position)
  const at = index.position[i];
  while (i > 0 && index.position[i - 1] === at) i--;
  const out: string[] = [];
  for (; i < index.position.length && index.position[i] === at; i++)
    if (at + index.length[i] > p) out.push(index.component[i]);
  return out;
}

/**
 * A position's components: its modified rows, else its letter's standard component (undefined where the letter has
 * none, or the letters are not held).
 */
export function componentsAt(
  index: ModifiedIndex | undefined,
  alphabet: SequenceAlphabet | undefined,
  letters: string | undefined,
  p: number,
): readonly string[] {
  const own = index ? modifiedAt(index, p) : [];
  if (own.length) return own;
  const std = alphabet && letters && p < letters.length ? standardComponent(alphabet, letters[p]) : undefined;
  return std ? [std] : [];
}

/** A modified index on a unit's positions: its molecule's rows shifted by the unit's start and clipped to it. */
export function shiftModified(index: ModifiedIndex, start: number, end: number): ModifiedIndex {
  const position: number[] = [];
  const component: string[] = [];
  const length: number[] = [];
  for (let i = 0; i < index.position.length; i++) {
    const s = Math.max(index.position[i], start);
    const e = Math.min(index.position[i] + index.length[i], end);
    if (e <= s) continue;
    position.push(s - start);
    component.push(index.component[i]);
    length.push(e - s);
  }
  return { position: Int32Array.from(position), component, length: Int32Array.from(length) };
}

/** What a unit's entry reads of what it instantiates: its system and, for a molecule, its entry. */
export interface UnitOf {
  readonly system: SystemEntry;
  readonly molecule?: MoleculeEntry;
}

/**
 * A unit's entry and its systems: the unit (its molecule's positions over its range, from 0; its declared numberings'
 * ids) and, for a molecule with atoms, `<unit>/(atoms)`: the molecule's atoms of the unit's residues.
 */
export function unitEntry(
  address: string,
  complex: string,
  node: Node,
  of: UnitOf,
): { unit: UnitEntry; systems: SystemEntry[]; atoms?: AtomsEntry } {
  const range = unitRange(node);
  const start = range ? range[0] : 0;
  const end = range ? range[1] : of.system.length;
  const length = end === undefined ? undefined : Math.max(0, end - start);
  const id = node.id ?? '';
  const obs = p(node, 'observed') as Record<string, Value> | undefined;
  const os = ints(obs?.start);
  const oe = ints(obs?.end);
  let observed: Int32Array | undefined;
  if (os && oe) {
    observed = new Int32Array(2 * Math.min(os.length, oe.length));
    for (let i = 0; 2 * i < observed.length; i++) {
      observed[2 * i] = os[i];
      observed[2 * i + 1] = oe[i];
    }
  }
  const declared: string[] = [];
  for (const c of node.children ?? [])
    if (c.kind === 'numbering' && c.id !== undefined && c.id !== 'index' && c.id !== 'ordinal') declared.push(c.id);
  const label = str(p(node, 'label'));
  const s = of.system;
  // a unit's letters are its sequence's over its range: inline, sliced; a source's, its sequence's plan, read on the
  // sequence's positions through the unit's step (a layer on the unit reads its parent's, `parent`), pending while
  // unread, never drawn empty (EC-169)
  const letters = s.letters?.slice(start, end);
  const access: AccessPlan =
    letters !== undefined
      ? inlinePlan(address, address, length, { letters: { form: 'letters', start: 0, letters } })
      : s.access.kind === 'source'
        ? s.access
        : inlinePlan(address, address, length, {});
  const system: SystemEntry = {
    address,
    kind: 'unit',
    ...(s.alphabet ? { alphabet: s.alphabet } : {}),
    ...(length !== undefined ? { length } : {}),
    cyclic: false,
    ...(letters !== undefined ? { letters } : {}),
    ...(s.source ? { source: s.source } : {}),
    name: label ?? id,
    numberings: ['index', 'ordinal', ...declared],
    access,
    parent: { system: s.address, start },
  };
  const unit: UnitEntry = {
    address,
    complex,
    id,
    of: s.address,
    start,
    ...(end !== undefined ? { end } : {}),
    class: s.kind === 'molecule' ? (of.molecule?.class ?? 'ligand') : 'polymer',
    ...pick(node, { chain: 'chain', authChain: 'auth_chain', operator: 'operator', label: 'label', group: 'group' }),
    ...(observed ? { observed } : {}),
    copy: `${s.address}|${start}|${end ?? ''}`,
  };
  const systems = [system];
  const ma = of.molecule?.atoms;
  if (!ma || end === undefined) return { unit, systems };
  const first = ma.offsets[start];
  const offsets = Int32Array.from({ length: end - start + 1 }, (_, r) => ma.offsets[start + r] - first);
  const atoms: AtomsEntry = {
    address: atomsSystemOf(address),
    of: address,
    offsets,
    names: shifted(ma.names, first),
    elements: shifted(ma.elements, first),
    first,
  };
  systems.push(atomsSystem(atoms, system.name));
  return { unit, systems, atoms };
}

function shifted(c: { get(i: number): string | null }, by: number): { get(i: number): string | null } {
  return by === 0 ? c : { get: (i) => c.get(i + by) };
}

function pick<K extends string>(node: Node, keys: Readonly<Record<K, string>>): Partial<Record<K, string>> {
  const out: Partial<Record<K, string>> = {};
  for (const k of Object.keys(keys) as K[]) {
    const v = str(p(node, keys[k]));
    if (v !== undefined) out[k] = v;
  }
  return out;
}

/**
 * A complex's entry from its units in order: copy groups derived by key (EC-65), the same groups for any order of the
 * units; its links' addresses.
 */
export function complexEntry(address: string, node: Node, units: readonly UnitEntry[], links: readonly string[]) {
  const copies = new Map<string, string[]>();
  for (const u of units) {
    let g = copies.get(u.copy);
    if (!g) copies.set(u.copy, (g = []));
    g.push(u.address);
  }
  const name = str(p(node, 'name'));
  const entry: ComplexEntry = { address, ...(name !== undefined ? { name } : {}), units, copies, links };
  return entry;
}

/**
 * A unit's two steps (_spec/README.md "Resolve", Paths): up to its molecule (position i to `start + i`) and down from
 * it, free in the search and under the copy rule (no path climbs to a molecule and down to another unit of it).
 */
export function unitEdges(unit: UnitEntry): readonly [Edge, Edge] {
  const up = new ImplicitLinear(unit.start, 1, 0, unit.end === undefined ? Infinity : unit.end - unit.start);
  return [
    { map: unit.address, from: unit.address, to: unit.of, backward: false, step: 0, linear: up, unit: 'up' },
    {
      map: unit.address,
      from: unit.of,
      to: unit.address,
      backward: true,
      step: 0,
      linear: up.inverse(),
      unit: 'down',
    },
  ];
}

// ---- the atoms group map, at a path's ends ----

/** The residue that owns an atom: the last `r` with `offsets[r] ≤ atom`. */
function residueOf(offsets: Int32Array, atom: number): number {
  let lo = 0;
  let hi = offsets.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >>> 1;
    if (offsets[mid] <= atom) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/**
 * Atoms' ranges collapsed onto the residues they touch (as nucleotides onto residues, EC-18): `[start, end, …]` on the
 * atoms system to `[start, end, …]` on its residues system (unsorted, possibly overlapping).
 */
export function collapseAtoms(atoms: AtomsEntry, ranges: ArrayLike<number>, out: number[] = []): number[] {
  out.length = 0;
  const n = atoms.offsets[atoms.offsets.length - 1];
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const s = Math.max(0, ranges[i]);
    const e = Math.min(n, ranges[i + 1]);
    if (e <= s) continue;
    out.push(residueOf(atoms.offsets, s), residueOf(atoms.offsets, e - 1) + 1);
  }
  return out;
}

/** Residues' ranges expanded onto their atoms' runs: `[start, end, …]` on the residues system to its atoms system. */
export function expandAtoms(atoms: AtomsEntry, ranges: ArrayLike<number>, out: number[] = []): number[] {
  out.length = 0;
  const last = atoms.offsets.length - 1;
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const s = Math.max(0, Math.min(last, ranges[i]));
    const e = Math.max(0, Math.min(last, ranges[i + 1]));
    if (atoms.offsets[e] > atoms.offsets[s]) out.push(atoms.offsets[s], atoms.offsets[e]);
  }
  return out;
}
