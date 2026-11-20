/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the molecule forms read of a model (_spec/README.md "Queries", Forms; "Resolve", Molecules and units): a
// system's components, a complex's units, a molecule, an atoms system. One implementation over the entries of
// `molecules.ts`, so the static model and the overlay answer alike.

import { type Node, parseAddress } from '@molstar/core/repr/tree';
import type { Numbering, QueryAtoms, QueryComponents, QueryMolecule, QueryUnit } from '@molstar/seqx-core/query';
import type { SystemEntry } from './model.js';
import {
  type AtomsEntry,
  atomsSystemOf,
  type ComplexEntry,
  componentsAt,
  type ModifiedIndex,
  type MoleculeEntry,
  shiftModified,
  type UnitEntry,
} from './molecules.js';

/** The molecule entries of a model or an overlay, by address. */
export interface Molecular {
  readonly molecules: ReadonlyMap<string, MoleculeEntry>;
  readonly complexes: ReadonlyMap<string, ComplexEntry>;
  readonly units: ReadonlyMap<string, UnitEntry>;
  /** Atoms systems' group maps, by the atoms system's address. */
  readonly atoms: ReadonlyMap<string, AtomsEntry>;
  /** Sequences' modified indices, by the sequence's address. */
  readonly modified: ReadonlyMap<string, ModifiedIndex>;
  /** Molecule and unit nodes, by address (what a node form selects). */
  readonly nodes: ReadonlyMap<string, Node>;
  /** Dynamic units' declared numberings (inline labels), by the unit's address, then id: the overlay's. */
  readonly numberings?: ReadonlyMap<string, ReadonlyMap<string, Numbering>>;
  /** Dynamic compounds (a source's `compounds` node), by code: the overlay's. */
  readonly compounds?: ReadonlyMap<string, CompoundEntry>;
}

export const EmptyMolecular: Molecular = {
  molecules: new Map(),
  complexes: new Map(),
  units: new Map(),
  atoms: new Map(),
  modified: new Map(),
  nodes: new Map(),
};

/** A compound described by a `compound` node: its name, class and parents. */
export interface CompoundEntry {
  readonly id: string;
  readonly name?: string;
  readonly class?: string;
  readonly parents?: readonly string[];
}

const Shifted = new WeakMap<UnitEntry, ModifiedIndex>();

/** A sequence's modified index, or a unit's of one (its sequence's shifted and clipped, made once). */
export function modifiedOf(m: Molecular, system: string): ModifiedIndex | undefined {
  const own = m.modified.get(system);
  if (own) return own;
  const u = m.units.get(system);
  const of = u && m.modified.get(u.of);
  if (!u || !of || u.end === undefined) return undefined;
  let s = Shifted.get(u);
  if (!s) Shifted.set(u, (s = shiftModified(of, u.start, u.end)));
  return s;
}

/** A system's components: a sequence's (or a unit's of one) beside its letters, a molecule's residues (or a unit's). */
export function componentsOf(
  m: Molecular,
  system: SystemEntry | undefined,
  systemOf: (address: string) => SystemEntry | undefined,
): QueryComponents | undefined {
  if (!system || system.length === undefined) return undefined;
  const length = system.length;
  const unit = system.kind === 'unit' ? m.units.get(system.address) : undefined;
  const molecule = m.molecules.get(unit ? unit.of : system.address);
  if (molecule) {
    const start = unit?.start ?? 0;
    const residues = molecule.residues;
    return {
      length,
      polymer: false,
      at: (p) => {
        const code = residues?.[start + p];
        return code ? [code] : [];
      },
    };
  }
  if (system.kind !== 'sequence' && system.kind !== 'unit') return undefined;
  const index = modifiedOf(m, system.address);
  const letters = system.letters ?? (unit ? systemOf(unit.of)?.letters?.slice(unit.start, unit.end) : undefined);
  return {
    length,
    polymer: true,
    // a source's letters are not held here: a form reading standard components waits
    ...(letters === undefined && system.source ? { pending: true } : {}),
    ...(index ? { modified: index } : {}),
    at: (p) => componentsAt(index, system.alphabet, letters, p),
  };
}

/** A complex's units as the node forms read them, or undefined when the address is no complex. */
export function unitsOf(m: Molecular, complex: string): readonly QueryUnit[] | undefined {
  const c = m.complexes.get(complex);
  if (!c) return undefined;
  return c.units.map((u) => ({
    address: u.address,
    node: m.nodes.get(u.address) ?? { kind: 'unit', id: u.id },
    id: u.id,
    of: u.of,
    ...(u.chain !== undefined ? { chain: u.chain } : {}),
    ...(u.authChain !== undefined ? { authChain: u.authChain } : {}),
    ...(u.operator !== undefined ? { operator: u.operator } : {}),
    ...entityOf(u.of),
  }));
}

/** A source's molecule's entity id: its dynamic node's key (`<source>/(node :key 3)`); none for a spec's molecule. */
function entityOf(of: string): { entity?: string } {
  const s = parseAddress(of).pattern?.segments;
  const last = s?.[s.length - 1];
  return last && typeof last !== 'string' && last.kind === 'node' ? { entity: String(last.params?.key) } : {};
}

/** The complexes at or under an address, in the order they were read. */
export function complexesUnder(m: Molecular, owner: string): string[] {
  const out: string[] = [];
  for (const c of m.complexes.keys()) if (c === owner || c.startsWith(owner === '/' ? '/' : `${owner}/`)) out.push(c);
  return out;
}

/** A molecule as `(molecule …)` reads it. */
export function moleculeOf(m: Molecular, address: string): QueryMolecule | undefined {
  const e = m.molecules.get(address);
  if (!e) return undefined;
  return {
    address,
    node: m.nodes.get(address) ?? { kind: 'molecule' },
    class: e.class,
    identifiers: e.identifiers.map((i) => `${i.namespace}:${i.value}`),
  };
}

/** The atoms system of a molecule or a unit (or an atoms system itself). */
export function atomsOf(m: Molecular, system: string): QueryAtoms | undefined {
  const a = m.atoms.get(system) ?? m.atoms.get(atomsSystemOf(system));
  if (!a) return undefined;
  const { offsets } = a;
  const n = offsets.length - 1;
  return {
    address: a.address,
    length: offsets[n],
    names: a.names,
    element: a.elements,
    of: a.of,
    residue: (atom) => {
      // the last offset at or before the atom: residues own contiguous runs (`AtomsEntry.offsets`)
      let lo = 0;
      let hi = n - 1;
      while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if (offsets[mid] <= atom) lo = mid;
        else hi = mid - 1;
      }
      return lo;
    },
    atomsOf: (r) => [offsets[r], offsets[r + 1]],
  };
}
