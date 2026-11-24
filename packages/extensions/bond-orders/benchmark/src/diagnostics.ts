/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { type Structure, Unit } from '@molstar/model/model/structure';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { BondType } from '@molstar/model/model/structure/model/types';

/** Atoms written but not judged: tier C's caps (plan §4.2). */
export type Skip = (element: ElementIndex) => boolean;
const none: Skip = () => false;

/** Highest bond-order sum, for the elements whose valence does not expand (a stated or assigned charge included). */
const MaxValence: Record<string, number> = { H: 1, B: 4, C: 4, N: 4, O: 2, F: 1 };

export interface BondRow {
  key: string;
  order: number;
  aromatic: boolean;
  /** the atoms' positions in the mol block (0-based): `toMolBlock` writes the atomic units' elements in order */
  atoms: [number, number];
}

/** The covalent bonds of a structure by atom-name pair (`a-b`, `a < b`): the L0 comparison (plan §5). */
export function bondRows(structure: Structure, skip = none): Map<string, BondRow> {
  const rows = new Map<string, BondRow>();
  let offset = 0;
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const { a, b, edgeProps } = unit.bonds;
    for (let t = 0; t < a.length; t++) {
      if (a[t] > b[t] || !BondType.isCovalent(edgeProps.flags[t])) continue;
      if (skip(unit.elements[a[t]]) || skip(unit.elements[b[t]])) continue;
      const na = label_atom_id.value(unit.elements[a[t]]);
      const nb = label_atom_id.value(unit.elements[b[t]]);
      const key = na < nb ? `${na}-${nb}` : `${nb}-${na}`;
      rows.set(key, {
        key,
        order: edgeProps.order[t],
        aromatic: BondType.is(edgeProps.flags[t], BondType.Names.aromatic),
        atoms: [offset + a[t], offset + b[t]],
      });
    }
    offset += unit.elements.length;
  }
  return rows;
}

export interface OrderDiff {
  /** `a-b expected→perceived` */
  bonds: string[];
  /** the same bonds by their atoms' positions in the mol block, for the report's highlights (plan §4.3) */
  atoms: [number, number][];
}

/** The bonds whose order differs. */
export function orderDiff(expected: Map<string, BondRow>, perceived: Map<string, BondRow>): OrderDiff {
  const diff: OrderDiff = { bonds: [], atoms: [] };
  for (const [key, e] of expected) {
    const p = perceived.get(key);
    if (!p || p.order === e.order) continue;
    diff.bonds.push(`${key} ${e.order}→${p.order}`);
    diff.atoms.push(e.atoms);
  }
  return diff;
}

/**
 * Checks of a perceived structure independent of the reference (plan §5.2): an atom over its valence, an atom of an
 * aromatic bond with more than one multiple bond, an aromatic C with none (its ring left without its double).
 */
export function diagnostics(structure: Structure, skip = none): string[] {
  const found: string[] = [];
  // bonds between units (a covalent ligand's link) count in their atoms' valence
  const inter = new Map<string, { sum: number; multiple: number }>();
  for (const e of structure.interUnitBonds.edges) {
    if (!BondType.isCovalent(e.props.flag)) continue;
    const key = `${e.unitA}:${e.indexA}`;
    const v = inter.get(key) ?? { sum: 0, multiple: 0 };
    v.sum += e.props.order;
    if (e.props.order > 1) v.multiple++;
    inter.set(key, v);
  }
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const { label_atom_id, type_symbol } = unit.model.atomicHierarchy.atoms;
    const { offset, edgeProps } = unit.bonds;
    for (let i = 0; i < unit.elements.length; i++) {
      if (skip(unit.elements[i])) continue;
      const across = inter.get(`${unit.id}:${i}`);
      let sum = across?.sum ?? 0;
      let multiple = across?.multiple ?? 0;
      let aromatic = false;
      for (let t = offset[i]; t < offset[i + 1]; t++) {
        if (!BondType.isCovalent(edgeProps.flags[t])) continue;
        sum += edgeProps.order[t];
        if (edgeProps.order[t] > 1) multiple++;
        if (BondType.is(edgeProps.flags[t], BondType.Names.aromatic)) aromatic = true;
      }
      const name = label_atom_id.value(unit.elements[i]);
      const element = type_symbol.value(unit.elements[i]);
      const max = MaxValence[element];
      if (max !== void 0 && sum > max) found.push(`${name}: valence ${sum}`);
      if (aromatic && multiple > 1) found.push(`${name}: aromatic, ${multiple} multiple bonds`);
      if (aromatic && multiple === 0 && element === 'C') found.push(`${name}: aromatic C without a double`);
    }
  }
  return found;
}
