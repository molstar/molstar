/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Augment Agent
 */

/**
 * ZML (Zipped Molecular Lot) file format. See the `ama.s3t.internal.zml`
 * reference Python implementation. Supported versions:
 *
 *  - V1: `molsys.json` + `pos.npy`.
 *  - V3: `header.json` + columnar topology buffers + `pos.f64`.
 *
 * Force field entries (`ff*.pkl`, `ff_*`) are ignored. V0 (DMS + pickle)
 * and V2 (full pickle) are not supported.
 */

export const ZmlSupportedVersions: ReadonlyArray<number> = [1, 3];

export const ZmlV1FileNames = {
    version: 'version',
    molsys: 'molsys.json',
    positions: 'pos.npy',
} as const;

export const ZmlV3FileNames = {
    version: 'version',
    header: 'header.json',
    Z: 'Z.i16',
    formalCharges: 'fc.i16',
    atomIds: 'atom_ids.i32',
    resIds: 'res_ids.i32',
    resAtomOffsets: 'res_atom_offsets.i32',
    molResOffsets: 'mol_res_offsets.i32',
    bonds: 'bonds.i32',
    bondOrders: 'bond_orders.i8',
    positions: 'pos.f64',
} as const;

/** JSON header of a V3 archive. */
export interface ZmlV3Header {
    readonly name: string;
    readonly box: ReadonlyArray<ReadonlyArray<number>> | null;
    readonly props: { readonly [k: string]: number | string };
    readonly atom_names: ReadonlyArray<string>;
    readonly atom_props: ReadonlyArray<{ readonly [k: string]: number | string }>;
    readonly res_names: ReadonlyArray<string>;
    readonly res_icodes: ReadonlyArray<string>;
    readonly res_props: ReadonlyArray<{ readonly [k: string]: number | string }>;
    readonly mol_names: ReadonlyArray<string>;
    readonly mol_chain_ids: ReadonlyArray<string>;
    readonly mol_props: ReadonlyArray<{ readonly [k: string]: number | string }>;
}

export interface ZmlAtom {
    readonly id: number;
    readonly name: string;
    readonly atomic_number: number;
    readonly formal_charge: number;
    readonly props: { readonly [k: string]: number | string };
}

export interface ZmlResidue {
    readonly id: number;
    readonly name: string;
    readonly insertion_code: string;
    readonly atoms: ReadonlyArray<ZmlAtom>;
    readonly props: { readonly [k: string]: number | string };
}

export interface ZmlMol {
    readonly name: string;
    readonly chain_id: string;
    readonly comtype: string;
    readonly residues: ReadonlyArray<ZmlResidue>;
    readonly props: { readonly [k: string]: number | string };
}

/** `[atomIndexA, atomIndexB, order]` with indices into the flat atom list. */
export type ZmlBond = readonly [number, number, number];

export interface ZmlMolSys {
    readonly name: string;
    readonly molecules: ReadonlyArray<ZmlMol>;
    readonly bonds: ReadonlyArray<ZmlBond>;
    /** 3x3 box vectors in Angstrom, or null for no periodic box. */
    readonly box: ReadonlyArray<ReadonlyArray<number>> | null;
    readonly props: { readonly [k: string]: number | string };
}

export interface ZmlFile {
    readonly name: string;
    readonly molsys: ZmlMolSys;
    /** x/y/z interleaved, Angstrom, length == 3 * atomCount. */
    readonly positions: Float64Array;
    readonly atomCount: number;
}
