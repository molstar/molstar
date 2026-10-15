/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '@molstar/core/data/db';
import type { Model } from '@molstar/model/model/structure';
import type { SchemaItem } from '@molstar/model/model/structure/structure/element/schema';
import type { RowOperator } from '@molstar/model/model/structure/structure/element/schema-selection/selections';

/** What rows can refer to, of one atom (a missing label_seq_id or auth_seq_id is `undefined`, a missing insertion code is `''`) */
interface AtomRecord {
  /** Position of the chain and of the residue in the file */
  iChain: number;
  iRes: number;
  label_asym_id: string;
  auth_asym_id: string;
  label_entity_id: string;
  label_seq_id: number | undefined;
  auth_seq_id: number | undefined;
  pdbx_PDB_ins_code: string;
  residue_index: number;
  label_comp_id: string;
  auth_comp_id: string;
  label_atom_id: string;
  auth_atom_id: string;
  type_symbol: string;
  atom_id: number;
  atom_index: number;
}

const defined = (x: unknown) => x !== undefined && x !== null;
const equals = (required: unknown, value: unknown) => !defined(required) || required === value;
/** `value` is within the bounds which are given (a value which is missing is not within any bound) */
function within(beg: number | undefined | null, end: number | undefined | null, value: number | undefined) {
  if (defined(beg) && (value === undefined || value < beg!)) return false;
  if (defined(end) && (value === undefined || value > end!)) return false;
  return true;
}

/** An independent and fast definition of "an atom satisfies a row" for tests: what every atom has is read from the model once,
 * so that many rows can be tested against all atoms (the production `atomQualifies` reads the columns for every pair).
 * The rule: every criterion of the row which is given has to hold, and a value which is missing (no `label_seq_id`) satisfies no criterion about it. */
export function createRowOracle(model: Model) {
  const h = model.atomicHierarchy;
  const n = h.atoms._rowCount;
  const { Present } = Column.ValueKind;
  const atoms: AtomRecord[] = new Array(n);
  for (let iAtom = 0; iAtom < n; iAtom++) {
    const iChain = h.chainAtomSegments.index[iAtom];
    const iRes = h.residueAtomSegments.index[iAtom];
    atoms[iAtom] = {
      iChain,
      iRes,
      label_asym_id: h.chains.label_asym_id.value(iChain),
      auth_asym_id: h.chains.auth_asym_id.value(iChain),
      label_entity_id: h.chains.label_entity_id.value(iChain),
      label_seq_id:
        h.residues.label_seq_id.valueKind(iRes) === Present ? h.residues.label_seq_id.value(iRes) : undefined,
      auth_seq_id: h.residues.auth_seq_id.valueKind(iRes) === Present ? h.residues.auth_seq_id.value(iRes) : undefined,
      pdbx_PDB_ins_code: h.residues.pdbx_PDB_ins_code.value(iRes),
      residue_index: h.residueSourceIndex.value(iRes),
      label_comp_id: h.atoms.label_comp_id.value(iAtom),
      auth_comp_id: h.atoms.auth_comp_id.value(iAtom),
      label_atom_id: h.atoms.label_atom_id.value(iAtom),
      auth_atom_id: h.atoms.auth_atom_id.value(iAtom),
      type_symbol: h.atoms.type_symbol.value(iAtom),
      atom_id: model.atomicConformation.atomId.value(iAtom),
      atom_index: h.atomSourceIndex.value(iAtom),
    };
  }

  // The residues of every auth chain (all the chains which have an `auth_asym_id`) in the order of the file (all atoms of a residue are the same here)
  const residuesOfAuthChain = new Map<
    string,
    { iRes: number; auth_seq_id: number | undefined; ins: string; water: boolean }[]
  >();
  for (const a of atoms) {
    const list = residuesOfAuthChain.get(a.auth_asym_id) ?? [];
    residuesOfAuthChain.set(a.auth_asym_id, list);
    if (list.length === 0 || list[list.length - 1].iRes !== a.iRes)
      list.push({
        iRes: a.iRes,
        auth_seq_id: a.auth_seq_id,
        ins: a.pdbx_PDB_ins_code,
        water: a.label_comp_id === 'HOH',
      });
  }

  /** For a row with a range with insertion codes: the residues of each auth chain from the first residue of the beginning to the last of the end
   * (a bound without a code is the number: the first residue with it, the last), by the definition on the residues of the auth chain.
   * An auth chain without a bound, or with the bounds in reverse order, has none. `undefined` for a row without such a range. */
  function insCodeRangeOf(row: SchemaItem): Map<string, [number, number]> | 'nothing' | undefined {
    if (!defined(row.beg_auth_ins_code) && !defined(row.end_auth_ins_code)) return undefined;
    if (
      (defined(row.beg_auth_ins_code) && !defined(row.beg_auth_seq_id)) ||
      (defined(row.end_auth_ins_code) && !defined(row.end_auth_seq_id))
    )
      return 'nothing';
    const result = new Map<string, [number, number]>();
    for (const [authAsymId, residues] of residuesOfAuthChain) {
      // (a bound is a residue which is not a water, and a water only if no other residue has it)
      const find = (seq: number, ins: string | undefined | null) => {
        const matches = residues.filter((r) => r.auth_seq_id === seq && (!defined(ins) || r.ins === ins));
        const notWaters = matches.filter((r) => !r.water);
        return (notWaters.length > 0 ? notWaters : matches).map((r) => r.iRes);
      };
      let lo = residues[0].iRes;
      let hi = residues[residues.length - 1].iRes;
      if (defined(row.beg_auth_seq_id)) {
        const found = find(row.beg_auth_seq_id, row.beg_auth_ins_code);
        if (found.length === 0) continue;
        lo = Math.min(...found);
      }
      if (defined(row.end_auth_seq_id)) {
        const found = find(row.end_auth_seq_id, row.end_auth_ins_code);
        if (found.length === 0) continue;
        hi = Math.max(...found);
      }
      if (lo <= hi) result.set(authAsymId, [lo, hi]);
    }
    return result;
  }

  /** Does the atom satisfy the row, for units with this operator? (the most selective criteria first) */
  function satisfies(
    row: SchemaItem,
    a: AtomRecord,
    operator: RowOperator,
    insRange: ReturnType<typeof insCodeRangeOf>,
  ): boolean {
    if (insRange === 'nothing') return false;
    if (insRange) {
      const range = insRange.get(a.auth_asym_id);
      if (!range || a.iRes < range[0] || a.iRes > range[1]) return false;
    }
    if (defined(row.instance_id) && row.instance_id !== operator.instanceId) return false;
    if (defined(row.operator_name) && row.operator_name !== operator.name) return false;
    return (
      equals(row.label_asym_id, a.label_asym_id) &&
      equals(row.label_seq_id, a.label_seq_id) &&
      equals(row.label_atom_id, a.label_atom_id) &&
      equals(row.auth_asym_id, a.auth_asym_id) &&
      equals(row.auth_seq_id, a.auth_seq_id) &&
      equals(row.label_entity_id, a.label_entity_id) &&
      equals(row.pdbx_PDB_ins_code, a.pdbx_PDB_ins_code) &&
      within(row.beg_label_seq_id, row.end_label_seq_id, a.label_seq_id) &&
      (insRange !== undefined || within(row.beg_auth_seq_id, row.end_auth_seq_id, a.auth_seq_id)) &&
      equals(row.residue_index, a.residue_index) &&
      equals(row.label_comp_id, a.label_comp_id) &&
      equals(row.auth_comp_id, a.auth_comp_id) &&
      equals(row.auth_atom_id, a.auth_atom_id) &&
      equals(defined(row.type_symbol) ? row.type_symbol!.toUpperCase() : undefined, a.type_symbol) &&
      equals(row.atom_id, a.atom_id) &&
      equals(row.atom_index, a.atom_index)
    );
  }

  return {
    nAtoms: n,
    /** The index of the last of `rows` which the atom satisfies, -1 for none (later rows win), for every atom */
    rowOfEachAtom(rows: readonly SchemaItem[], operator: RowOperator): number[] {
      const out = new Array<number>(n).fill(-1);
      for (let iRow = 0; iRow < rows.length; iRow++) {
        const row = rows[iRow];
        const insRange = insCodeRangeOf(row);
        for (let iAtom = 0; iAtom < n; iAtom++) {
          if (satisfies(row, atoms[iAtom], operator, insRange)) out[iAtom] = iRow;
        }
      }
      return out;
    },
    /** The atoms (ascending) which satisfy at least one of `rows` */
    atomsOf(rows: readonly SchemaItem[], operator: RowOperator): number[] {
      return this.rowOfEachAtom(rows, operator).flatMap((v, i) => (v >= 0 ? [i] : []));
    },
  };
}
