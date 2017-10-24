/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../common/column'

// Full format http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
// there are many records but for now ignore (pass over) all but the following
// @<TRIPOS>MOLECULE
// @<TRIPOS>ATOM
// @<TRIPOS>BOND
//
// note that the format is not a fixed column format but white space separated

export interface Molecule {
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

export interface Atoms {
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
    status_bit: Column<string>
}

export interface Bonds {
    count: number,

    bond_id: Column<number>,
    origin_atom_id: Column<number>,
    target_atom_id: Column<number>,
    bond_type: Column<string>,

    // optional in the format, assign UndefinedColumn if not available
    status_bits: Column<string>
}

export interface Structure {
    molecule: Readonly<Molecule>,
    atoms: Readonly<Atoms>,
    bonds: Readonly<Bonds>
}

export interface File {
    structures: Structure[]
}