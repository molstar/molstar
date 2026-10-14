/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coordinate systems as a query reads them (seqspec:seqql#elements-and-groups, seqspec:seqql#labels-and-indices): a
// system's domain, name, numbering and labels, and the alive units over a molecule. Edges: called per system, never per
// position.

import type { NameId, UnitId } from '../common.js';
import { CoordinateSystemKind as K } from '../coordinates.js';
import { LABEL, NO_LABEL, Numbering, type NumberingSet, type PositionLabel, UNKNOWN_SUFFIX } from '../numbering.js';
import { type SeqSet, SeqSets } from '../set.js';

/** Writes the system's indices [lo, hi) to `out`; false when it has none held (a removed unit, a row not held). */
export function domainOf(set: SeqSet, kind: number, id: number, row: number, out: Int32Array): boolean {
  const u = set.units;
  const a = set.alignments[id];
  let lo = 0;
  let hi = -1;
  if (kind === K.unit && u.alive[id]) {
    lo = u.rangeStart[id];
    hi = u.rangeStop[id];
  } else if (kind === K.molecule) hi = set.molecules[id]?.length ?? -1;
  else if (kind === K.row) hi = a ? a.rows.lengthOf(row as never) : -1;
  else if (kind === K.columns) hi = a?.columnCount ?? -1;
  out[0] = lo;
  out[1] = hi;
  return hi > lo;
}

/**
 * The circle a system runs on: a cyclic molecule's length when the system is the whole molecule (it, or a unit over
 * all of it), else 0. On a circle a range, a slice and an expand wrap the origin (seqspec:seqql#labels-and-indices).
 */
export function circleOf(set: SeqSet, kind: number, id: number): number {
  const mol = kind === K.unit ? set.units.molecule[id] : kind === K.molecule ? id : -1;
  const m = mol >= 0 ? set.molecules[mol] : undefined;
  if (!m?.cyclic) return 0;
  if (kind === K.unit && (set.units.rangeStart[id] !== 0 || set.units.rangeStop[id] !== m.length)) return 0;
  return m.length;
}

/** Its ref as a person writes it: `1brs/A`, `p53`, `msa/3`; columns: `msa columns`. */
export function nameOf(set: SeqSet, kind: number, id: number, row: number): string {
  const s = set.strings;
  if (kind === K.unit) return s.get(set.units.key[id]);
  if (kind === K.molecule) return s.get(set.molecules[id]?.key ?? -1);
  if (kind === K.complex) return s.get(set.complexes[id]?.key ?? -1);
  const a = s.get(set.alignments[id]?.key ?? -1);
  return kind === K.row ? `${a}/${row}` : `${a} columns`;
}

const ordinal = (set: SeqSet) => set.strings.find('ordinal') as NameId;

/**
 * The numbering `name` of a system, or its default (null `name`): its own default, its molecule's, else `ordinal`.
 * null: none (labels are index + 1); undefined: the system has no numbering by that name.
 */
export function numberingOf(set: SeqSet, kind: number, id: number, row: number, name: string | null) {
  const lists =
    kind === K.unit
      ? [set.units.numberings[id], set.molecules[set.units.molecule[id]]?.numberings]
      : kind === K.molecule
        ? [set.molecules[id]?.numberings]
        : kind === K.row
          ? [set.alignments[id]?.numberingsOf(row as never)]
          : [];
  if (name !== null) {
    const n = set.strings.find(name) as NameId;
    const exact = kind === K.unit ? SeqSets.numberingOf(set, id as UnitId, n) : lists[0]?.byName(n);
    if (exact) return exact;
    const under = numberingsUnder(set, lists, name);
    return under.length === 1 ? under[0] : undefined;
  }
  for (const l of lists) for (const n of l?.items ?? []) if (n.isDefault) return n;
  return lists.at(-1)?.byName(ordinal(set)) ?? null;
}

/**
 * A namespace alone (`uniprot`) names the numberings under it (`uniprot:P04637`, as loaders write them; CM-02):
 * distinct by name, the unit's before its molecule's. A name with a `:` names nothing more than itself.
 */
function numberingsUnder(set: SeqSet, lists: readonly (NumberingSet | null | undefined)[], ns: string): Numbering[] {
  if (ns.includes(':')) return [];
  const out = new Map<string, Numbering>();
  for (const l of lists)
    for (const x of l?.items ?? []) {
      const name = set.strings.get(x.name);
      if (name.startsWith(`${ns}:`) && !out.has(name)) out.set(name, x);
    }
  return [...out.values()];
}

/** The numberings a namespace alone names on a system (more than one: the name is ambiguous there), by name. */
export function numberingNamesUnder(set: SeqSet, kind: number, id: number, row: number, ns: string): string[] {
  return numberingNames(set, kind, id, row).filter((n) => !ns.includes(':') && n.startsWith(`${ns}:`));
}

/** The names of a system's numberings, for messages. */
export function numberingNames(set: SeqSet, kind: number, id: number, row: number): string[] {
  const own = kind === K.unit ? set.units.numberings[id]?.items : [];
  const mol = kind === K.unit ? set.units.molecule[id] : kind === K.molecule ? id : -1;
  const rows = kind === K.row ? set.alignments[id]?.numberingsOf(row as never)?.items : [];
  const all = [...(own ?? []), ...(set.molecules[mol]?.numberings.items ?? []), ...(rows ?? [])];
  return [...new Set(all.map((n) => set.strings.get(n.name)))];
}

/** A position's label in a numbering (null: index + 1); "" when unlabeled. */
export const labelAt = (set: SeqSet, n: Numbering | null, index: number) =>
  n ? Numbering.format(n, index as never, set.strings) : String(index + 1);

/**
 * A written label in a numbering (seqspec:format#positions): a name in a `names` numbering, else an integer and a
 * suffix.
 */
export function toLabel(set: SeqSet, n: Numbering | null, text: string, out: PositionLabel): boolean {
  out.num = NO_LABEL;
  out.suffix = out.name = -1 as NameId;
  out.rank = undefined;
  if (n?.method === 'names') {
    out.name = set.strings.find(text) as NameId;
    return out.name >= 0;
  }
  const m = LABEL.exec(text);
  if (!m) return false;
  out.num = Number(m[1]);
  if (!m[2]) return true;
  // Found, never interned: a query reads its set (CM-11). A suffix no string holds is no label's, placed by its order.
  const id = set.strings.find(m[2]);
  out.suffix = id >= 0 ? (id as NameId) : UNKNOWN_SUFFIX;
  if (id < 0) out.rank = Numbering.suffixOrder(m[2]);
  return true;
}

const byMolecule = new WeakMap<SeqSet, { start: Int32Array; units: Int32Array }>();
/** Alive units over each molecule, in id order: built once per set. */
export function unitsOf(set: SeqSet): { readonly start: Int32Array; readonly units: Int32Array } {
  let x = byMolecule.get(set);
  if (!x) {
    const u = set.units;
    const start = new Int32Array(set.molecules.length + 1);
    for (let k = 0; k < u.count; k++) if (u.alive[k]) start[u.molecule[k] + 1]++;
    for (let m = 0; m < set.molecules.length; m++) start[m + 1] += start[m];
    const units = new Int32Array(start[set.molecules.length]);
    const next = start.slice(0, -1);
    for (let k = 0; k < u.count; k++) if (u.alive[k]) units[next[u.molecule[k]]++] = k;
    x = { start, units };
    byMolecule.set(set, x);
  }
  return x;
}
