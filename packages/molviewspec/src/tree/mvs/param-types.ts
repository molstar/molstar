/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { HexColor, ColorName } from '../../util/helpers';
import { ValueFor, float, int, list, literal, str, tuple, union } from '../generic/field-schema';
import { ColorNames } from '../../util/color';


/** `format` parameter values for `parse` node in MVS tree */
export const ParseFormatT = literal('mmcif', 'bcif', 'pdb', 'map');
export type ParseFormatT = ValueFor<typeof ParseFormatT>

/** `format` parameter values for `parse` node in Molstar tree */
export const MolstarParseFormatT = literal('cif', 'pdb', 'map');
export type MolstarParseFormatT = ValueFor<typeof MolstarParseFormatT>

/** `kind` parameter values for `structure` node in MVS tree */
export const StructureTypeT = literal('model', 'assembly', 'symmetry', 'symmetry_mates');

/** `selector` parameter values for `component` node in MVS tree */
export const ComponentSelectorT = literal('all', 'polymer', 'protein', 'nucleic', 'branched', 'ligand', 'ion', 'water', 'coarse');

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
    label_comp_id: str,
    auth_comp_id: str,
    // residue_index: int, // 0-based residue index in the source file // TODO this is defined in Python builder but not supported by Molstar yet
    label_atom_id: str,
    auth_atom_id: str,
    type_symbol: str,
    atom_id: int,
    atom_index: int,
});
export type ComponentExpressionT = ValueFor<typeof ComponentExpressionT>

/** `schema` parameter values for `*_from_uri` and `*_from_source` nodes in MVS tree */
export const SchemaT = literal('whole_structure', 'entity', 'chain', 'auth_chain', 'residue', 'auth_residue', 'residue_range', 'auth_residue_range', 'atom', 'auth_atom', 'all_atomic');

/** `format` parameter values for `*_from_uri` nodes in MVS tree */
export const SchemaFormatT = literal('cif', 'bcif', 'json');

/** Parameter values for vector params, e.g. `position` */
export const Vector3 = tuple([float, float, float]);
export type Vector3 = ValueFor<typeof Vector3>

/** Parameter values for matrix params, e.g. `rotation` */
export const Matrix = list(float);

/** Primitives-related types */
export const PrimitiveComponentExpressionT = iots.partial({ structure_ref: str, expression_schema: SchemaT, expressions: list(ComponentExpressionT) });
export type PrimitiveComponentExpressionT = ValueFor<typeof PrimitiveComponentExpressionT>
export const PrimitivePositionT = iots.union([Vector3, ComponentExpressionT, PrimitiveComponentExpressionT]);
export type PrimitivePositionT = ValueFor<typeof PrimitivePositionT>

export const FloatList = list(float);
export const IntList = list(int);
export const StrList = list(str);


/** `color` parameter values for `color` node in MVS tree */
export const HexColorT = new iots.Type<HexColor>(
    'HexColor',
    ((value: any) => typeof value === 'string') as any,
    (value, ctx) => HexColor.is(value) ? { _tag: 'Right', right: value } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid hex color string` }] },
    value => value
);

/** `color` parameter values for `color` node in MVS tree */
export const ColorNameT = new iots.Type<ColorName>(
    'ColorName',
    ((value: any) => typeof value === 'string') as any,
    (value, ctx) => ColorName.is(value) ? { _tag: 'Right', right: value } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid hex color string` }] },
    value => value
);

/** `color` parameter values for `color` node in MVS tree */
export const ColorNamesT = literal(...Object.keys(ColorNames) as (keyof typeof ColorNames)[]);

/** `color` parameter values for `color` node in MVS tree */
export const ColorT = union([ColorNameT, HexColorT]);
export type ColorT = ValueFor<typeof ColorT>

/** Type helpers */
export function isVector3(x: any): x is Vector3 {
    return !!x && Array.isArray(x) && x.length === 3 && typeof x[0] === 'number';
}

export function isPrimitiveComponentExpressions(x: any): x is PrimitiveComponentExpressionT {
    return !!x && Array.isArray(x.expressions);
}

export function isComponentExpression(x: any): x is ComponentExpressionT {
    return !!x && typeof x === 'object' && !x.expressions;
}
