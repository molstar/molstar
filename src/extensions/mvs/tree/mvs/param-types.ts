/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { ColorNames } from '../../../../mol-util/color/names';
import { ColorName, HexColor } from '../../helpers/utils';
import { ValueFor, bool, dict, float, int, list, literal, nullable, obj, partial, str, tuple, union } from '../generic/field-schema';


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
export const ComponentExpressionT = partial({
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
export const PrimitiveComponentExpressionT = partial({ structure_ref: str, expression_schema: SchemaT, expressions: list(ComponentExpressionT) });
export type PrimitiveComponentExpressionT = ValueFor<typeof PrimitiveComponentExpressionT>
export const PrimitivePositionT = union([Vector3, ComponentExpressionT, PrimitiveComponentExpressionT]);
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
    (value, ctx) => ColorName.is(value) ? { _tag: 'Right', right: value } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid color name` }] },
    value => value
);

/** `color` parameter values for `color` node in MVS tree */
export const ColorNamesT = literal(...Object.keys(ColorNames) as (keyof ColorNames)[]);

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


export const ColorListNameT = literal(
    // Color lists from https://observablehq.com/@d3/color-schemes (definitions: https://colorbrewer2.org/export/colorbrewer.js)
    // Sequential single-hue
    'Blues', 'Greens', 'Greys', 'Oranges', 'Purples', 'Reds',
    // Sequential multi-hue
    'BuGn', 'BuPu', 'GnBu', 'OrRd', 'PuBuGn', 'PuBu', 'PuRd', 'RdPu', 'YlGnBu', 'YlGn', 'YlOrBr', 'YlOrRd',
    'Cividis', 'Viridis', 'Inferno', 'Magma', 'Plasma', 'Warm', 'Cool', 'CubehelixDefault', 'Turbo',
    // Diverging
    'BrBG', 'PRGn', 'PiYG', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
    // Cyclical
    'Rainbow', 'Sinebow',
    // Categorical
    'Category10', 'Observable10', 'Tableau10',
    'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Dark2', 'Paired', 'Accent',

    // Additional list, not standard for visualization in general, but commonly used for structures
    'Chainbow',
);
export type ColorListNameT = ValueFor<typeof ColorListNameT>;

export const ColorMappingNameT = literal('ElementSymbol', 'ResidueName', 'ResidueProperties');
// TODO add meaningful options
// TODO decide on naming (ResidueName vs JmolResidueName, ResidueProperties vs ClustalResidueProperties)
// TODO would it make sense to have a switch for case-insensitive
export type ColorMappingNameT = ValueFor<typeof ColorMappingNameT>;

export const CategoricalPaletteNameT = union([ColorListNameT, ColorMappingNameT]);
export type CategoricalPaletteNameT = ValueFor<typeof CategoricalPaletteNameT>;

export const CategoricalPalette = iots.intersection([
    obj({ kind: literal('categorical') }),
    partial({
        colors: union([
            CategoricalPaletteNameT,
            list(ColorT),
            dict(str, ColorT),
        ]),
        /** Color to use when a) `color` is a dictionary and given key is not present, or b) `color` is a list or a named palette and there are more real values than listed values and `repeat_color_list` is not true. */
        missing_color: ColorT,
        /** Repeat color list once all colors are depleted (only applies if `color` is a list or a named palette). */
        repeat_color_list: bool,
        /** Sort real values before assigning colors from a list or named palette. */
        sort: nullable(literal('lexical', 'numeric')),
        /** Sort direction (only applies if `sort` is provided). */
        sort_direction: literal('ascending', 'descending'),
    }),
]);
export type CategoricalPalette = ValueFor<typeof CategoricalPalette>;

// TODO consider spreading the palette param directly into color_from_uri/color_from_source params (though this will be tricky) or achieve smart error messages

export const Palette = CategoricalPalette;
// export const Palette = union([CategoricalPalette, DiscretePalette, ContinuousPalette]);

// Draft from https://docs.google.com/document/d/1p9yePdtvO8RzYQ90jEdCHM5sMqxFpy4f8DbleXagRDE/edit?tab=t.0

// class GradientPalette:
//     kind: Literal["gradient"] = "gradient"
//     # either uniformly distributed or explicitly scaled between 0, 1
//     # [('red', 0), ('green', 0.2), ('blue', 1)]
//     stops: list[tuple[ColorT, float]] | list[tuple[ColorT, float, float]] | list[ColorT] | None
//     stop_value_kind: Literal["normalized", "explicit"]
//     name: PalleteNameT | None

//     value_domain: tuple[float, float] | None = None  # min, max | none <=> auto


// class DiscretePalette:
//     kind: Literal["discrete"] = "discrete"
//     # either uniformly distributed or explicitly scaled between 0, 1
//     # [('red', 0), ('green', 0.2), ('blue', 1)]
//     stops: list[tuple[ColorT, float]] | list[tuple[ColorT, float, float]] | list[ColorT] | None
//     stop_value_kind: Literal["normalized", "explicit"]
//     name: PalleteNameT | None

//     value_domain: tuple[float, float] | None = None  # min, max | none <=> auto


// class CategoricalPalette:
//     kind: Literal["categorical"] = "categorical"

//     colors: dict[Any, ColorT] | list[ColorT] | None
//     name: PalleteNameT | None

//     missing_color: ColorT | None = None # applied when color is missing dict[Any, ColorT]
//     sort_values: Literal["ascending", "descending"] | None = None
//     sort_kind: Literal["lexical", "numeric"] | None = None
