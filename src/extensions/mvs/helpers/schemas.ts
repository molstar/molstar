/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from '../../../mol-data/db';
import { pickObjectKeys } from '../../../mol-util/object';
import { Choice } from '../../../mol-util/param-choice';

const { str, int } = Column.Schema;


/** Names of allowed MVS annotation schemas (values for the annotation schema parameter) */
export type MVSAnnotationSchema = Choice.Values<typeof MVSAnnotationSchema>
export const MVSAnnotationSchema = new Choice(
    {
        whole_structure: 'Whole Structure',
        entity: 'Entity',
        chain: 'Chain (label*)',
        auth_chain: 'Chain (auth*)',
        residue: 'Residue (label*)',
        auth_residue: 'Residue (auth*)',
        residue_range: 'Residue range (label*)',
        auth_residue_range: 'Residue range (auth*)',
        atom: 'Atom (label*)',
        auth_atom: 'Atom (auth*)',
        all_atomic: 'All atomic selectors',
    },
    'all_atomic',
);

/** Represents a set of criteria for selection of atoms in a model (in `all_atomic` schema).
 * Missing/undefined values mean that we do not care about that specific atom property. */
export type MVSAnnotationRow = Partial<Table.Row<typeof AllAtomicCifAnnotationSchema>>


/** Get CIF schema definition for given annotation schema name */
export function getCifAnnotationSchema<K extends MVSAnnotationSchema>(schemaName: K): Pick<typeof AllAtomicCifAnnotationSchema, (typeof FieldsForSchemas)[K][number]> {
    return pickObjectKeys(AllAtomicCifAnnotationSchema, FieldsForSchemas[schemaName]);
}


/** Definition of `all_atomic` schema for CIF (other atomic schemas are subschemas of this one) */
const AllAtomicCifAnnotationSchema = {
    /** Tag for grouping multiple annotation rows with the same `group_id` (e.g. to show one label for two chains);
     * if the `group_id` is not given, each row is processed separately */
    group_id: str,

    label_entity_id: str,
    label_asym_id: str,
    auth_asym_id: str,

    label_seq_id: int,
    auth_seq_id: int,
    pdbx_PDB_ins_code: str,
    /** Minimum label_seq_id (inclusive) */
    beg_label_seq_id: int,
    /** Maximum label_seq_id (inclusive) */
    end_label_seq_id: int,
    /** Minimum auth_seq_id (inclusive) */
    beg_auth_seq_id: int,
    /** Maximum auth_seq_id (inclusive) */
    end_auth_seq_id: int,
    label_comp_id: str,
    auth_comp_id: str,
    // residue_index: int, // 0-based residue index in the source file // TODO this is defined in Python builder but not supported by Molstar yet

    /** Atom name like 'CA', 'N', 'O'... */
    label_atom_id: str,
    /** Atom name like 'CA', 'N', 'O'... */
    auth_atom_id: str,
    /** Element symbol like 'H', 'He', 'Li', 'Be' (case-insensitive)... */
    type_symbol: str,
    /** Unique atom identifier across conformations (_atom_site.id) */
    atom_id: int,
    /** 0-based index of the atom in the source data */
    atom_index: int,
    /** Instance identifier to distinguish instances of the same chain created by applying different symmetry operators,
     * like 'ASM-X0-1' for assemblies or '1_555' for crystals */
    instance_id: str,
} satisfies Table.Schema;

/** Allowed fields (i.e. CIF columns or JSON keys) for each annotation schema
 * (other fields will just be ignored) */
const FieldsForSchemas = {
    whole_structure: ['group_id'],
    entity: ['group_id', 'label_entity_id'],
    chain: ['group_id', 'label_entity_id', 'label_asym_id'],
    auth_chain: ['group_id', 'auth_asym_id'],
    residue: ['group_id', 'label_entity_id', 'label_asym_id', 'label_seq_id'],
    auth_residue: ['group_id', 'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code'],
    residue_range: ['group_id', 'label_entity_id', 'label_asym_id', 'beg_label_seq_id', 'end_label_seq_id'],
    auth_residue_range: ['group_id', 'auth_asym_id', 'beg_auth_seq_id', 'end_auth_seq_id'],
    atom: ['group_id', 'label_entity_id', 'label_asym_id', 'label_seq_id', 'label_atom_id', 'type_symbol', 'atom_id', 'atom_index'],
    auth_atom: ['group_id', 'auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_atom_id', 'type_symbol', 'atom_id', 'atom_index'],
    all_atomic: Object.keys(AllAtomicCifAnnotationSchema) as (keyof typeof AllAtomicCifAnnotationSchema)[],
} satisfies { [schema in MVSAnnotationSchema]: (keyof typeof AllAtomicCifAnnotationSchema)[] };
