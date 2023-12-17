/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import * as iots from 'io-ts';
import { HexColor, isHexColorString } from '../../helpers/utils';
import { ValueFor, float, int, list, literal, str, tuple, union } from '../generic/params-schema';


/** `format` parameter values for `parse` node in MVS tree */
export const ParseFormatT = literal('mmcif', 'bcif', 'pdb');
export type ParseFormatT = ValueFor<typeof ParseFormatT>

/** `format` parameter values for `parse` node in Molstar tree */
export const MolstarParseFormatT = literal('cif', 'pdb');
export type MolstarParseFormatT = ValueFor<typeof MolstarParseFormatT>

/** `kind` parameter values for `structure` node in MVS tree */
export const StructureTypeT = literal('model', 'assembly', 'symmetry', 'symmetry_mates');

/** `selector` parameter values for `component` node in MVS tree */
export const ComponentSelectorT = literal('all', 'polymer', 'protein', 'nucleic', 'branched', 'ligand', 'ion', 'water');

/** `selector` parameter values for `component` node in MVS tree */
export const ComponentExpressionT = iots.partial({
    label_entity_id: str,
    label_asym_id: str,
    auth_asym_id: str,
    label_seq_id: int,
    auth_seq_id: int,
    pdbx_PDB_ins_code: str,
    beg_label_seq_id: int,
    end_label_seq_id: int,
    beg_auth_seq_id: int,
    end_auth_seq_id: int,
    label_atom_id: str,
    auth_atom_id: str,
    type_symbol: str,
    atom_id: int,
    atom_index: int,
});

/** `type` parameter values for `representation` node in MVS tree */
export const RepresentationTypeT = literal('ball_and_stick', 'cartoon', 'surface');

/** `schema` parameter values for `*_from_uri` and `*_from_source` nodes in MVS tree */
export const SchemaT = literal('whole_structure', 'entity', 'chain', 'auth_chain', 'residue', 'auth_residue', 'residue_range', 'auth_residue_range', 'atom', 'auth_atom', 'all_atomic');

/** `format` parameter values for `*_from_uri` nodes in MVS tree */
export const SchemaFormatT = literal('cif', 'bcif', 'json');

/** Parameter values for vector params, e.g. `position` */
export const Vector3 = tuple([float, float, float]);

/** Parameter values for matrix params, e.g. `rotation` */
export const Matrix = list(float);

/** `color` parameter values for `color` node in MVS tree */
export const HexColorT = new iots.Type<HexColor>(
    'HexColor',
    ((value: any) => typeof value === 'string') as any,
    (value, ctx) => isHexColorString(value) ? { _tag: 'Right', right: value } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid hex color string` }] },
    value => value
);

/** `color` parameter values for `color` node in MVS tree */
export const ColorNamesT = literal('white', 'gray', 'black', 'red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'magenta');

/** `color` parameter values for `color` node in MVS tree */
export const ColorT = union([HexColorT, ColorNamesT]);
