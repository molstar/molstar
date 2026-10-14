/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Units, complexes, copy groups and the set (seqspec:model#identity-and-lineage).

import { type Alignment, type AlignmentInput, Alignments } from './alignment.js';
import { type Annotation, type AnnotationInput, Annotations } from './annotation.js';
import { type ComplexId, CustomData, Issue, type NameId, StringTable, type UnitId } from './common.js';
import { type Mapping, type MappingInput, Mappings, type RowBindings } from './mapping.js';
import { type Molecule, type MoleculeInput, Molecules } from './molecule.js';
import { Numbering, type NumberingInput, type NumberingSet } from './numbering.js';
import type { SortedArray, SortedRanges } from './ordered-set.js';

/** Struct-of-arrays, one row per unit id; removed ids are tombstones (`alive` 0). */
export interface UnitTable {
  readonly count: number;
  readonly alive: Uint8Array;
  /** NameId: its ref ("1brs/A", "pose3", "p53/unit"). */
  readonly key: Int32Array;
  /** MoleculeId. */
  readonly molecule: Int32Array;
  /** ComplexId, -1: none. */
  readonly complex: Int32Array;
  /** The shown subset, [rangeStart, rangeStop). */
  readonly rangeStart: Int32Array;
  readonly rangeStop: Int32Array;
  /** null: all of the range. */
  readonly observed: readonly (SortedRanges | null)[];
  /** GroupId. */
  readonly group: Int32Array;
  /** NameId, open ("sequence", "experimental", "docked", "<ns>:<name>"); -1: none. */
  readonly origin: Int32Array;
  /** NaN: none. */
  readonly score: Float64Array;
  /** NameId. */
  readonly label: Int32Array;
  /** Unit-level (seqspec:model#numberings): they win over the molecule's of the same name. */
  readonly numberings: readonly NumberingSet[];
  readonly custom: readonly CustomData[];
}

/** Units over one molecule with an equal range. */
export interface CopyGroups {
  readonly count: number;
  /** UnitId of the first member in data order. */
  readonly first: Int32Array;
  readonly members: readonly SortedArray<UnitId>[];
}

export interface Complex {
  readonly id: ComplexId;
  /** Its ref ("1brs"): what a provider's refs are under. */
  readonly key: NameId;
  readonly label: NameId;
  /** UnitId, in the complex's order. */
  readonly units: Int32Array;
  readonly custom: CustomData;
}

export const KeyKind = { molecule: 0, unit: 1, alignment: 2, row: 3, mapping: 4, annotation: 5, complex: 6 } as const;
/** Reused by the caller. */
export interface KeyTarget {
  kind: number;
  id: number;
  row: number;
}

/**
 * What validation resolves refs against: a set, or a host's resolver over its build inputs' keys, so data is validated
 * before it is assembled. `molecules` by id: their kinds, so a ligand's key is no sequence's.
 */
export interface KeyResolver {
  byName(key: string, out: KeyTarget): boolean;
  readonly molecules: readonly ({ readonly kind: number } | null)[];
  /**
   * A ref that may name a dynamic node of a streamed source, not held yet (seqspec:sources#dynamic-nodes): valid, kind
   * unknown.
   */
  dynamic?(key: string): boolean;
}

export interface SeqSet {
  /** One per build; kept by withChanges. */
  readonly lineage: number;
  readonly strings: StringTable;
  /** By MoleculeId; null: removed. */
  readonly molecules: readonly (Molecule | null)[];
  readonly units: UnitTable;
  readonly complexes: readonly (Complex | null)[];
  readonly groups: CopyGroups;
  readonly alignments: readonly (Alignment | null)[];
  readonly mappings: readonly (Mapping | null)[];
  readonly annotations: readonly (Annotation | null)[];
  /** What a ref names; `"R/3"` names row 3 of alignment R. */
  byKey(key: NameId, out: KeyTarget): boolean;
  /** As `byKey`, by the ref's text: a row ref need not be interned. Allocation-free. */
  byName(key: string, out: KeyTarget): boolean;
  readonly custom: CustomData;
}

export interface UnitInput {
  key: string;
  molecule: string;
  range?: [number, number];
  observed?: [number, number][];
  label?: string;
  origin?: string;
  score?: number;
  numberings?: NumberingInput[];
  custom?: object;
}
/** Units by key, in order. */
export interface ComplexInput {
  key: string;
  units: string[];
  label?: string;
  custom?: object;
}
export interface SetInput {
  molecules: MoleculeInput[];
  units?: UnitInput[];
  complexes?: ComplexInput[];
  alignments?: AlignmentInput[];
  mappings?: MappingInput[];
  annotations?: AnnotationInput[];
}

interface SetImpl extends SeqSet {
  /** NameId → id * 8 + KeyKind. */
  readonly keys: ReadonlyMap<number, number>;
  /** A removed or evicted key's id * 8 + KeyKind: given back when it is upserted again (the stable-id rule). */
  readonly gone: ReadonlyMap<number, number>;
  /** The row → bound spaces tables, by AlignmentId (seqspec:model#locations-and-loci). */
  readonly bindings: ReadonlyMap<number, RowBindings>;
}

function byKey(this: SetImpl, key: NameId, out: KeyTarget): boolean {
  const v = this.keys.get(key);
  out.row = -1;
  if (v === undefined) return key >= 0 && rowRef(this, this.strings.get(key), out);
  out.kind = v & 7;
  out.id = v >> 3;
  return true;
}

function byName(this: SetImpl, key: string, out: KeyTarget): boolean {
  const id = this.strings.find(key);
  return id >= 0 && this.keys.has(id) ? this.byKey(id, out) : rowRef(this, key, out);
}

/** The row a row ref names (`msa/3`: 3): the digits after its last `/`; -1: none. Data mode reads it too (`Systems`). */
export function rowOfRef(s: string): number {
  const slash = s.lastIndexOf('/');
  if (slash <= 0 || slash === s.length - 1) return -1;
  let row = 0;
  for (let i = slash + 1; i < s.length; i++) {
    const d = s.charCodeAt(i) - 48;
    if (d < 0 || d > 9) return -1;
    row = row * 10 + d;
  }
  return row;
}

/** `"<alignment key>/<k>"`: row k, without slicing or interning (alignments are tens). */
function rowRef(set: SetImpl, s: string, out: KeyTarget): boolean {
  const row = rowOfRef(s);
  if (row < 0) return false;
  const slash = s.lastIndexOf('/');
  for (let a = 0; a < set.alignments.length; a++) {
    const aln = set.alignments[a];
    if (!aln || row >= aln.rowCount) continue;
    const k = set.strings.get(aln.key);
    if (k.length !== slash || !s.startsWith(k)) continue;
    out.kind = KeyKind.row;
    out.id = a;
    out.row = row;
    return true;
  }
  return false;
}

let nextLineage = 1;

/** The unit table's typed columns but the derived `group`; `observed`, `numberings` and `custom` are plain lists. */
const TYPED = {
  alive: Uint8Array,
  key: Int32Array,
  molecule: Int32Array,
  complex: Int32Array,
  rangeStart: Int32Array,
  rangeStop: Int32Array,
  origin: Int32Array,
  score: Float64Array,
  label: Int32Array,
} as const;
const LISTS = ['observed', 'numberings', 'custom'] as const;
type Draft = { [K in keyof typeof TYPED]: number[] } & {
  observed: (SortedRanges | null)[];
  numberings: NumberingSet[];
  custom: CustomData[];
};

function draftOf(t: UnitTable): Draft {
  const d: Record<string, unknown[]> = {};
  for (const k of Object.keys(TYPED) as (keyof typeof TYPED)[]) d[k] = Array.from(t[k]);
  for (const k of LISTS) d[k] = [...t[k]];
  return d as unknown as Draft;
}

/** The table, with copy groups: equal (molecule, rangeStart, rangeStop); an own unit is a group of its own. */
function freeze(d: Draft, sequence: number): { units: UnitTable; groups: CopyGroups } {
  const t: Record<string, unknown> = { count: d.key.length, group: new Int32Array(d.key.length).fill(-1) };
  for (const k of Object.keys(TYPED) as (keyof typeof TYPED)[]) t[k] = (TYPED[k] as typeof Int32Array).from(d[k]);
  for (const k of LISTS) t[k] = d[k];
  const units = t as unknown as UnitTable & { group: Int32Array };
  const groupOf = new Map<string, number>();
  const members: number[][] = [];
  for (let u = 0; u < units.count; u++) {
    if (!units.alive[u]) continue;
    const k =
      units.origin[u] === sequence ? `#${u}` : `${units.molecule[u]}:${units.rangeStart[u]}:${units.rangeStop[u]}`;
    if (!groupOf.has(k)) {
      groupOf.set(k, members.length);
      members.push([]);
    }
    units.group[u] = groupOf.get(k) as number;
    members[units.group[u]].push(u);
  }
  const first = Int32Array.from(members, (m) => m[0]);
  return { units, groups: { count: members.length, first, members: members.map((m) => Int32Array.from(m)) } };
}

function emptySet(): SetImpl {
  const draft = Object.fromEntries([...Object.keys(TYPED), ...LISTS].map((k) => [k, []])) as unknown as Draft;
  return {
    ...freeze(draft, -1),
    molecules: [],
    complexes: [],
    alignments: [],
    mappings: [],
    annotations: [],
    lineage: nextLineage++,
    strings: new StringTable(),
    keys: new Map(),
    gone: new Map(),
    bindings: new Map(),
    byKey,
    byName,
    custom: CustomData.create(),
  };
}

function apply(
  base: SetImpl,
  upsert: Partial<SetInput>,
  remove: readonly string[],
  at: string,
  evict: readonly string[] = [],
) {
  const issues: Issue[] = [];
  const { strings } = base;
  const keys = new Map(base.keys);
  // Copied on the first change: most versions remove and revive nothing.
  let gone = base.gone;
  const tomb = (key: number, v: number) => {
    if (gone === base.gone) gone = new Map(gone);
    (gone as Map<number, number>).set(key, v);
  };
  /** Keys evicted, and their dependents evicted with them: quietly. */
  const evicted = new Set<number>();
  /** Unit ids given back to a key removed earlier: in no complex. */
  const revived = new Set<number>();
  const molecules = base.molecules.slice();
  const complexes = base.complexes.slice();
  const lists = {
    alignments: base.alignments.slice(),
    mappings: base.mappings.slice(),
    annotations: base.annotations.slice(),
  };
  const seen = new Set<number>();
  const lengthChanged = new Set<number>();
  const upserted = new Set<number>();
  const ownUnits: [UnitInput, string][] = [];
  /** Per key (NameId) of a molecule, unit or alignment: the path of the upsert or removal that changed it. */
  const changedBy = new Map<number, string>();
  const sequence = strings.intern('sequence');

  /** Keys given back their old id by `claim`: no longer gone once their build succeeds (the sweep below). */
  const claimed: number[] = [];
  /**
   * The existing id of `key` with this kind, or the id it had before it was removed; -1 when new; undefined (with an
   * error) on a clash.
   */
  const claim = (key: string, kind: number, path: string): number | undefined => {
    const id = strings.intern(key);
    const v = keys.get(id);
    if (seen.has(id) || (v !== undefined && (v & 7) !== kind)) {
      issues.push(Issue.error('model.key.duplicate', path, `Key "${key}" names something else already.`));
      return undefined;
    }
    seen.add(id);
    if (v !== undefined) return v >> 3;
    const g = gone.get(id);
    if (g === undefined || (g & 7) !== kind) return -1;
    // Left in `gone` until the build succeeds: a failed one keeps the old id for the next try (CM-07).
    claimed.push(id);
    if (kind === KeyKind.unit) revived.add(g >> 3);
    return g >> 3;
  };
  const find = (key: string, kind: number) => {
    const v = keys.get(strings.find(key));
    return v !== undefined && (v & 7) === kind ? v >> 3 : -1;
  };

  (upsert.molecules ?? []).forEach((input, k) => {
    const path = `${at}/molecules/${k}`;
    const prev = claim(input.key, KeyKind.molecule, path);
    if (prev === undefined) return;
    const id = prev < 0 ? molecules.length : prev;
    const m = Molecules.build(input, id, strings, issues, path);
    if (!m) return;
    if (prev >= 0 && molecules[prev]?.length !== m.length) lengthChanged.add(id);
    changedBy.set(m.key, path);
    molecules[id] = m;
    keys.set(m.key, id * 8 + KeyKind.molecule);
    if (input.ownUnit)
      ownUnits.push([{ key: `${input.key}/unit`, molecule: input.key, label: input.label ?? input.key }, path]);
  });

  // The unit table is re-frozen only when units change: set again after removals.
  let touched = !!upsert.units?.length || !!upsert.complexes?.length || ownUnits.length + lengthChanged.size > 0;
  const d = draftOf(base.units);
  const kill = (u: number, issue: Issue | null, path: string) => {
    d.alive[u] = 0;
    keys.delete(d.key[u]);
    tomb(d.key[u], u * 8 + KeyKind.unit);
    changedBy.set(d.key[u], path);
    if (issue) issues.push(issue);
  };

  const putUnit = (input: UnitInput, path: string, own: boolean) => {
    const err = (code: string, m: string) => issues.push(Issue.error(code, path, `Unit "${input.key}": ${m}`));
    const prev = claim(input.key, KeyKind.unit, path);
    if (prev === undefined) return;
    const mol = molecules[find(input.molecule, KeyKind.molecule)];
    if (!mol) return err('model.unit.molecule', `no molecule "${input.molecule}".`);
    const [start, stop] = input.range ?? [0, mol.length];
    if (mol.cyclic && start > stop) {
      const m = `range [${start}, ${stop}) wraps the origin; a unit's range does not yet: show it with a block's range.`;
      return err('model.unit.range', m);
    }
    if (!(Number.isInteger(start) && Number.isInteger(stop) && start >= 0 && start < stop && stop <= mol.length)) {
      return err('model.unit.range', `range [${start}, ${stop}) is not within 0..${mol.length}.`);
    }
    let observed: Int32Array | null = Int32Array.from((input.observed ?? []).flat());
    for (let i = 0; i < observed.length; i += 2) {
      const lo = i ? observed[i - 1] : start;
      if (!(observed[i] >= lo && observed[i] < observed[i + 1] && observed[i + 1] <= stop)) {
        return err('model.unit.observed', `observed ranges must be sorted, disjoint and within [${start}, ${stop}).`);
      }
    }
    if (!input.observed || (observed.length === 2 && observed[0] === start && observed[1] === stop)) observed = null;
    if (!own && input.origin === 'sequence')
      return err('model.unit.origin', `origin "sequence" is reserved for own units.`);
    const u = prev < 0 ? d.key.length : prev;
    const p = `${path}/numberings`;
    d.alive[u] = 1;
    d.key[u] = strings.intern(input.key);
    d.molecule[u] = mol.id;
    d.complex[u] = prev < 0 || revived.has(u) ? -1 : d.complex[u];
    d.rangeStart[u] = start;
    d.rangeStop[u] = stop;
    d.observed[u] = observed;
    d.origin[u] = own ? sequence : input.origin == null ? -1 : strings.intern(input.origin);
    d.score[u] = input.score ?? Number.NaN;
    d.label[u] = strings.intern(input.label ?? input.key);
    d.numberings[u] = Numbering.buildSet(input.numberings, mol.length, strings, issues, p, false);
    d.custom[u] = CustomData.create(input.custom);
    keys.set(d.key[u], u * 8 + KeyKind.unit);
    changedBy.set(d.key[u], path);
    upserted.add(u);
  };
  for (const [input, path] of ownUnits) putUnit(input, path, true);
  (upsert.units ?? []).forEach((input, k) => {
    putUnit(input, `${at}/units/${k}`, false);
  });

  (upsert.complexes ?? []).forEach((input, k) => {
    const path = `${at}/complexes/${k}`;
    const prev = claim(input.key, KeyKind.complex, path);
    if (prev === undefined) return;
    const id = prev < 0 ? complexes.length : prev;
    for (let u = 0; u < d.complex.length; u++) if (d.complex[u] === id) d.complex[u] = -1;
    const units: number[] = [];
    input.units.forEach((key, j) => {
      const u = find(key, KeyKind.unit);
      if (u < 0 || d.complex[u] >= 0) {
        const m = u < 0 ? `no unit "${key}".` : `unit "${key}" belongs to a complex already.`;
        issues.push(Issue.error('model.complex.unit', `${path}/units/${j}`, `Complex "${input.key}": ${m}`));
        return;
      }
      d.complex[u] = id;
      units.push(u);
    });
    const key = strings.intern(input.key);
    complexes[id] = {
      id: id as ComplexId,
      key,
      label: strings.intern(input.label ?? input.key),
      units: Int32Array.from(units),
      custom: CustomData.create(input.custom),
    };
    keys.set(key, id * 8 + KeyKind.complex);
  });

  (upsert.alignments ?? []).forEach((input, k) => {
    const path = `${at}/alignments/${k}`;
    const prev = claim(input.key, KeyKind.alignment, path);
    if (prev === undefined) return;
    const id = prev < 0 ? lists.alignments.length : prev;
    const a = Alignments.build(input, id, strings, issues, path);
    if (!a) return;
    lists.alignments[id] = a;
    keys.set(a.key, id * 8 + KeyKind.alignment);
    changedBy.set(a.key, path);
  });

  // Mappings are built once the endpoints are final (assemble, below).
  const pending: Pending<MappingInput>[] = [];
  (upsert.mappings ?? []).forEach((input, k) => {
    const path = `${at}/mappings/${k}`;
    const prev = claim(input.key, KeyKind.mapping, path);
    if (prev === undefined) return;
    const id = prev < 0 ? Math.max(lists.mappings.length, ...pending.map((p) => p.id + 1)) : prev;
    pending.push({ id, input, path });
    keys.set(strings.intern(input.key), id * 8 + KeyKind.mapping);
  });
  // Annotations likewise, once mappings are (their targets are final then).
  const annotations: Pending<AnnotationInput>[] = [];
  (upsert.annotations ?? []).forEach((input, k) => {
    const path = `${at}/annotations/${k}`;
    const prev = claim(input.key, KeyKind.annotation, path);
    if (prev === undefined) return;
    const id = prev < 0 ? Math.max(lists.annotations.length, ...annotations.map((p) => p.id + 1)) : prev;
    annotations.push({ id, input, path });
    keys.set(strings.intern(input.key), id * 8 + KeyKind.annotation);
  });

  // Evicted keys go as removed ones do, an unknown one quietly.
  const drop = (key: string, k: number, quiet: boolean) => {
    const name = strings.find(key);
    const v = keys.get(name);
    if (v === undefined) {
      if (!quiet) issues.push(Issue.warning('model.key.unknown', `/remove/${k}`, `No key "${key}" to remove.`));
      return;
    }
    if (quiet) evicted.add(name);
    keys.delete(name);
    tomb(name, v);
    changedBy.set(name, quiet ? `/evict/${k}` : `/remove/${k}`);
    const kind = v & 7;
    const id = v >> 3;
    touched ||= kind === KeyKind.molecule || kind === KeyKind.unit || kind === KeyKind.complex;
    if (kind === KeyKind.molecule) molecules[id] = null;
    else if (kind === KeyKind.unit) d.alive[id] = 0;
    else if (kind === KeyKind.complex) complexes[id] = null;
    else
      lists[kind === KeyKind.alignment ? 'alignments' : kind === KeyKind.mapping ? 'mappings' : 'annotations'][id] =
        null;
  };
  for (let k = 0; k < remove.length; k++) drop(remove[k], k, false);
  for (let k = 0; k < evict.length; k++) drop(evict[k], k, true);

  // Units over removed molecules go; units over a molecule whose length changed must come with it,
  // unless a full range with no numberings, which follows the new length. Issues point at that change.
  for (let u = 0; u < d.key.length; u++) {
    if (!d.alive[u]) continue;
    const mol = molecules[d.molecule[u]];
    const path = changedBy.get((mol ?? base.molecules[d.molecule[u]])?.key ?? -1) ?? '';
    const unit = `Unit "${strings.get(d.key[u])}"`;
    if (!mol) {
      const ev = evicted.size > 0 && evicted.has(base.molecules[d.molecule[u]]?.key ?? -1);
      if (ev) evicted.add(d.key[u]);
      const quiet = ev || d.origin[u] === sequence;
      kill(u, quiet ? null : Issue.warning('model.unit.orphaned', path, `${unit} lost its molecule.`), path);
      continue;
    }
    if (lengthChanged.has(mol.id) && !upserted.has(u)) {
      const full = d.rangeStart[u] === 0 && d.rangeStop[u] === base.molecules[mol.id]?.length && !d.observed[u];
      if (full && d.numberings[u].count === 0) {
        d.rangeStop[u] = mol.length;
        changedBy.set(d.key[u], path);
      } else kill(u, Issue.warning('model.unit.stale', path, `${unit} no longer fits its molecule.`), path);
    }
    if (d.complex[u] >= 0 && !complexes[d.complex[u]]) d.complex[u] = -1;
  }
  for (let c = 0; c < complexes.length; c++) {
    const cx = complexes[c];
    if (cx?.units.some((u) => !d.alive[u])) complexes[c] = { ...cx, units: cx.units.filter((u) => d.alive[u]) };
  }

  // A key that came back under its old id is no longer gone; one whose build failed still is.
  for (const id of claimed) {
    if (!keys.has(id) || !gone.has(id)) continue;
    if (gone === base.gone) gone = new Map(gone);
    (gone as Map<number, number>).delete(id);
  }

  const table = touched ? freeze(d, sequence) : base;
  const set = {
    lineage: base.lineage,
    strings,
    molecules: same(molecules, base.molecules),
    units: table.units,
    complexes: same(complexes, base.complexes),
    groups: table.groups,
    alignments: same(lists.alignments, base.alignments),
    mappings: lists.mappings as readonly (Mapping | null)[],
    annotations: lists.annotations as readonly (Annotation | null)[],
    keys,
    gone,
    bindings: base.bindings,
    byKey,
    byName,
    custom: base.custom,
  };
  /** Whether a dependent is over an evicted key (a row ref by its alignment): evicted with it, quietly. */
  const isEvicted = (over: readonly string[]) => {
    for (let k = 0; evicted.size > 0 && k < over.length; k++)
      if (evicted.has(strings.find(over[k])) || evicted.has(strings.find(over[k].replace(/\/[0-9]+$/, ''))))
        return true;
    return false;
  };
  /** The change a dependent's issue points at: the first of its keys (a row ref by its alignment) that changed. */
  const causeOf = (over: readonly string[]) => {
    for (const k of over) {
      const at = changedBy.get(strings.find(k)) ?? changedBy.get(strings.find(k.replace(/\/[0-9]+$/, '')));
      if (at !== undefined) return at;
    }
    return changedBy.values().next().value ?? '';
  };
  const mappings: Dependents<Mapping, MappingInput> = {
    name: 'Mapping',
    kind: KeyKind.mapping,
    ...Mappings.dependents,
  };
  const tombs = { tomb, evicted, isEvicted };
  // What mappings and annotations are over unchanged (an annotation-only version: a stream's dynamic nodes), none is
  // checked again: the version costs its changes, not the set's size (ST1b-b).
  const stable = set.molecules === base.molecules && set.units === base.units && set.alignments === base.alignments;
  assemble(set, lists.mappings, pending, keys, issues, mappings, causeOf, tombs, stable);
  set.mappings = same(lists.mappings, base.mappings);
  const anns: Dependents<Annotation, AnnotationInput> = {
    name: 'Annotation',
    kind: KeyKind.annotation,
    ...Annotations.dependents,
  };
  assemble(set, lists.annotations, annotations, keys, issues, anns, causeOf, tombs, stable);
  set.gone = gone;
  set.annotations = same(lists.annotations, base.annotations);
  const unchanged = set.mappings === base.mappings && set.units === base.units && set.alignments === base.alignments;
  if (!unchanged) set.bindings = Mappings.bind(set);
  return { set: set as SeqSet, issues };
}

type Pending<I> = { id: number; input: I; path: string };
/** Tombstones and evictions, shared by the dependents' passes. */
interface Gone {
  tomb(key: number, v: number): void;
  readonly evicted: Set<number>;
  isEvicted(over: readonly string[]): boolean;
}
/** A kind of dependent: over endpoints or targets that a change may rebuild or remove. */
interface Dependents<T, I> {
  readonly name: string;
  readonly kind: number;
  build(set: SeqSet, input: I, id: number, issues: Issue[], path: string): T | null;
  /** Against the set's current version: 0 unchanged, 1 to rebuild, 2 orphaned. */
  check(set: SeqSet, entry: T): number;
  inputOf(entry: T): I;
  /** The keys it is over. */
  keysOf(entry: T): readonly string[];
}

/**
 * Dependents after a change: `pending` built, then every other entry checked against this version (`revalidate`):
 * kept, rebuilt, or removed as orphaned; none when what they are over is `stable` (the same molecules, units and
 * alignments). A failed upsert keeps the previous version, which is revalidated too.
 * Issues of a rebuild or an orphan point at the change that caused them (`causeOf`). `list` is updated in
 * place; `keys` loses removed keys.
 */
function assemble<T extends { readonly key: NameId }, I extends { readonly key: string }>(
  set: SeqSet,
  list: (T | null)[],
  pending: readonly Pending<I>[],
  keys: Map<number, number>,
  issues: Issue[],
  d: Dependents<T, I>,
  causeOf: (keys: readonly string[]) => string,
  gone: Gone,
  stable: boolean,
) {
  const revalidate = (id: number) => {
    const e = list[id];
    const status = e ? d.check(set, e) : 0;
    if (!e || status === 0) return;
    const input = d.inputOf(e);
    const local: Issue[] = [];
    const next = status === 1 ? d.build(set, input, id, local, '') : null;
    const lost = d.kind === KeyKind.mapping ? 'an endpoint' : 'its targets';
    const code = `model.${d.name.toLowerCase()}.orphaned`;
    // Over an evicted node: evicted with it, nothing reported.
    const quiet = status === 2 && gone.isEvicted(d.keysOf(e));
    if (quiet) gone.evicted.add(e.key);
    else if (status === 2) local.push(Issue.warning(code, '', `${d.name} "${input.key}" lost ${lost}.`));
    const cause = causeOf(d.keysOf(e));
    for (const issue of local) issues.push({ ...issue, path: cause });
    list[id] = next;
    if (!next) {
      keys.delete(e.key);
      gone.tomb(e.key, id * 8 + d.kind);
    }
  };
  const fresh = new Set(pending.map((p) => p.id));
  for (let id = 0; !stable && id < list.length; id++) if (!fresh.has(id)) revalidate(id);
  for (const p of pending) {
    const key = set.strings.intern(p.input.key);
    if (keys.get(key) !== p.id * 8 + d.kind) continue;
    const built = d.build(set, p.input, p.id, issues, p.path);
    if (built) list[p.id] = built;
    else if (list[p.id]) revalidate(p.id);
    else keys.delete(key);
  }
  for (let id = 0; id < list.length; id++) list[id] ??= null;
}

function same<T>(next: readonly T[], prev: readonly T[]): readonly T[] {
  return next.length === prev.length && next.every((x, i) => x === prev[i]) ? prev : next;
}

/** A unit's numbering by name: its own, else its molecule's. */
function numberingOf(set: SeqSet, unit: UnitId, name: NameId): Numbering | null {
  return (
    set.units.numberings[unit]?.byName(name) ?? set.molecules[set.units.molecule[unit]]?.numberings.byName(name) ?? null
  );
}

export const SeqSets = {
  build: (input: SetInput) => apply(emptySet(), input, [], ''),
  /**
   * Upserts by key, then removes and evicts keys; keeps the lineage and shares every table the changes do not touch.
   * A removed key upserted again gets its old id back; an evicted one's dependents go with it, quietly.
   */
  withChanges: (set: SeqSet, changes: { upsert?: Partial<SetInput>; remove?: string[]; evict?: string[] }) =>
    apply(set as SetImpl, changes.upsert ?? {}, changes.remove ?? [], '/upsert', changes.evict),
  numberingOf,
};
