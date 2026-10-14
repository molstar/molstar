/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A set's coordinate systems as their readers name, bound and number them (seqspec:model#coordinate-systems), one copy
// for every reader: SeqQL's labels and messages, and a viewer's own system keys, descriptors and rulers. A system is a
// (kind, id, row) of `CoordinateSystemKind`. Edge rate: called per system, never per position.

import { CoordinateSystemKind as K } from "./coordinates.js";
import type { Numbering } from "./numbering.js";
import { circleOf, domainOf, nameOf, numberingNames, numberingOf } from "./query/systems.js";
import { rowOfRef, type SeqSet } from "./set.js";

/**
 * The key naming a system in its set, as the set's `byName` reads it back: `1brs/A`, `p53`, an alignment's columns by
 * the alignment's own key (`msa`), a row `msa/3`. null: a removed unit, or a kind no key names (an item, an external
 * system). Not `nameOf`, which names columns for a person (`msa columns`).
 */
function refOf(set: SeqSet, kind: number, id: number, row: number): string | null {
  const s = (k: number | undefined) => (k === undefined ? null : set.strings.get(k));
  if (kind === K.unit) return set.units.alive[id] ? s(set.units.key[id]) : null;
  if (kind === K.molecule) return s(set.molecules[id]?.key);
  if (kind === K.complex) return s(set.complexes[id]?.key);
  const a = s(set.alignments[id]?.key);
  return kind === K.columns ? a : kind === K.row && a !== null ? `${a}/${row}` : null;
}

/** The row a row ref names (`msa/3`: 3), as the set's `byName` reads it; -1: none. */
const rowOf = rowOfRef;

/**
 * A system's numberings by name, in the order asked (a namespace alone names the one under it, as `numberingOf`), the
 * names it has none by left out; none asked: its default, else none (its labels are its indices + 1).
 */
function numberingsOf(set: SeqSet, kind: number, id: number, row: number, names: readonly string[]): Numbering[] {
  if (names.length === 0) {
    const n = numberingOf(set, kind, id, row, null);
    return n?.isDefault ? [n] : [];
  }
  const out: Numbering[] = [];
  for (const name of names) {
    const n = numberingOf(set, kind, id, row, name);
    if (n) out.push(n);
  }
  return out;
}

/** Coordinate systems of a set: their keys, names, domains and numberings. */
export const Systems = {
  refOf,
  rowOf,
  /** Its name for a person, in messages: its key; an alignment's columns `msa columns`. */
  nameOf,
  /** Writes its indices [lo, hi) to `out`; false when it holds none (a removed unit, a row not held). */
  domainOf,
  /** The circle it runs on: a cyclic molecule's length when it is the whole molecule, else 0. */
  circleOf,
  /** One numbering by name, else its default (null name); null: none (its indices + 1); undefined: no such name. */
  numberingOf,
  numberingsOf,
  /** The distinct names of its numberings (a unit's own first, then its molecule's), for messages and descriptors. */
  numberingNames,
};
