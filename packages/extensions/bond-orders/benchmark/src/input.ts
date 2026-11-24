/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { createAndFillArray } from '@molstar/core/data/db/column-helpers';
import { CifBlock, CifCategory, type CifField } from '@molstar/io/reader/cif/data-model';

/** What perception must not see (plan §3.3), removed on the parsed block, never on the file. */
export interface Removals {
  hydrogens: boolean;
  charges: boolean;
}

/** A CCD component without its H (atoms and their bonds) and without its charges, as asked. */
export function ccdInput(block: CifBlock, removals: Removals): CifBlock {
  const atoms = block.categories.chem_comp_atom;
  const bonds = block.categories.chem_comp_bond;
  const heavy = new Set<string>();
  const atomRows: number[] = [];
  const atomId = atoms.getField('atom_id')!;
  const typeSymbol = atoms.getField('type_symbol')!;
  for (let i = 0; i < atoms.rowCount; i++) {
    if (removals.hydrogens && isHydrogen(typeSymbol.str(i))) continue;
    atomRows.push(i);
    heavy.add(atomId.str(i));
  }

  const categories: { [name: string]: CifCategory } = { ...block.categories };
  categories.chem_comp_atom = viewCategory(atoms, atomRows, removals.charges ? ['charge'] : []);
  if (bonds) {
    const id1 = bonds.getField('atom_id_1')!;
    const id2 = bonds.getField('atom_id_2')!;
    const bondRows: number[] = [];
    for (let i = 0; i < bonds.rowCount; i++) {
      if (heavy.has(id1.str(i)) && heavy.has(id2.str(i))) bondRows.push(i);
    }
    // no bond left between heavy atoms (a lone CH3 of the CCD): no category, as Mol* reads a bond table's first row
    if (bondRows.length) categories.chem_comp_bond = viewCategory(bonds, bondRows, []);
    else delete categories.chem_comp_bond;
  }
  const names = block.categoryNames.filter((n) => categories[n]);
  return CifBlock(names, categories, block.header, block.saveFrames);
}

/**
 * An entry for tier C (plan §3.2, §3.3): model 1, the first altloc of each residue, without H (as asked) and without
 * charges. Its `chem_comp_bond` stays: the connectivity, whose orders `forceCompute` resets.
 */
export function entryInput(block: CifBlock, removals: Removals): CifBlock {
  const atoms = block.categories.atom_site;
  const typeSymbol = atoms.getField('type_symbol')!;
  const model = atoms.getField('pdbx_PDB_model_num');
  const alt = atoms.getField('label_alt_id');
  const asym = atoms.getField('label_asym_id')!;
  const seq = atoms.getField('auth_seq_id');
  const ins = atoms.getField('pdbx_PDB_ins_code');
  const firstModel = model?.str(0);
  const altOf = new Map<string, string>();
  const rows: number[] = [];
  for (let i = 0; i < atoms.rowCount; i++) {
    if (model && model.str(i) !== firstModel) continue;
    if (removals.hydrogens && isHydrogen(typeSymbol.str(i))) continue;
    const a = alt?.str(i) ?? '';
    if (a && a !== '.' && a !== '?') {
      const residue = `${asym.str(i)} ${seq?.str(i)} ${ins?.str(i)}`;
      const first = altOf.get(residue);
      if (first === void 0) altOf.set(residue, a);
      else if (first !== a) continue;
    }
    rows.push(i);
  }
  const categories: { [name: string]: CifCategory } = { ...block.categories };
  categories.atom_site = viewCategory(atoms, rows, removals.charges ? ['pdbx_formal_charge'] : []);
  return CifBlock([...block.categoryNames], categories, block.header, block.saveFrames);
}

export function isHydrogen(typeSymbol: string) {
  const s = typeSymbol.toUpperCase();
  return s === 'H' || s === 'D';
}

/** `category` restricted to `rows`, without the fields `drop` (read as absent). */
function viewCategory(category: CifCategory, rows: number[], drop: string[]): CifCategory {
  const fieldNames = category.fieldNames.filter((f) => !drop.includes(f));
  const fields: { [name: string]: CifField } = {};
  for (const name of fieldNames) fields[name] = viewField(category.getField(name)!, rows);
  return CifCategory(category.name, rows.length, fieldNames, fields);
}

function viewField(field: CifField, rows: number[]): CifField {
  const rowCount = rows.length;
  const str = (row: number) => field.str(rows[row]);
  const int = (row: number) => field.int(rows[row]);
  const float = (row: number) => field.float(rows[row]);
  return {
    binaryEncoding: void 0,
    isDefined: field.isDefined,
    rowCount,
    str,
    int,
    float,
    valueKind: (row) => field.valueKind(rows[row]),
    areValuesEqual: (a, b) => field.areValuesEqual(rows[a], rows[b]),
    toStringArray: (params) => createAndFillArray(rowCount, str, params),
    toIntArray: (params) => createAndFillArray(rowCount, int, params),
    toFloatArray: (params) => createAndFillArray(rowCount, float, params),
  };
}
