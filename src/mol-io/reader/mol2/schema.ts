/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../../mol-data/db';

// Full format https://zhanggroup.org/DockRMSD/mol2.pdf
// there are many records but for now ignore (pass over) all but the following
// @<TRIPOS>MOLECULE
// @<TRIPOS>ATOM
// @<TRIPOS>BOND
// @<TRIPOS>SUBSTRUCTURE
// @<TRIPOS>CRYSIN
//
// note that the format is not a fixed column format but white space separated

export interface Mol2Molecule {
    mol_name: string
    num_atoms: number
    num_bonds: number
    num_subst: number
    num_feat: number
    num_sets: number
    mol_type: string
    charge_type: string
    status_bits: string
    mol_comment: string
}

export interface Mol2Atoms {
    count: number,

    atom_id: Column<number>,
    atom_name: Column<string>,
    x: Column<number>,
    y: Column<number>,
    z: Column<number>,
    atom_type: Column<string>,

    // optional in the format, assign UndefinedColumn if not available
    subst_id: Column<number>,
    subst_name: Column<string>,
    charge: Column<number>,
    status_bits: Column<string>
}

export interface Mol2Bonds {
    count: number,

    bond_id: Column<number>,
    origin_atom_id: Column<number>,
    target_atom_id: Column<number>,
    bond_type: Column<string>,

    // optional in the format, assign UndefinedColumn if not available
    status_bits: Column<string>
}

export interface Mol2Substructure {
    count: number,

    subst_id: Column<number>,
    subst_name: Column<string>,
    root_atom: Column<number>,

    // optional in the format, assign UndefinedColumn if not available
    subst_type: Column<string>,
    dict_type: Column<string>,
    chain: Column<string>,
    sub_type: Column<string>,
    inter_bonds: Column<number>,
    status_bits: Column<string>
}

export interface Mol2Crysin {
    a: number
    b: number
    c: number
    alpha: number
    beta: number
    gamma: number
    spaceGroup: number
    setting: number
}

export interface Mol2Structure {
    molecule: Readonly<Mol2Molecule>,
    atoms: Readonly<Mol2Atoms>,
    bonds: Readonly<Mol2Bonds>,
    substructures?: Readonly<Mol2Substructure>,
    crysin?: Readonly<Mol2Crysin>
}

export interface Mol2File {
    name: string
    structures: Mol2Structure[]
}