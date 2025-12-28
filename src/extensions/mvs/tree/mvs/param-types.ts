/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as iots from 'io-ts';
import { ColorNames } from '../../../../mol-util/color/names';
import { ValueFor, bool, dict, float, int, list, literal, nullable, object, partial, str, tuple, union } from '../generic/field-schema';

/** `format` parameter values for `parse` node in MVS tree */
export type ParseFormatT =
    // trajectory
    | 'mmcif' | 'bcif' | 'pdb' | 'pdbqt' | 'gro' | 'xyz' | 'mol' | 'sdf' | 'mol2' | 'lammpstrj'
    // coordinates
    | 'xtc' | 'nctraj' | 'dcd' | 'trr'
    // topology
    | 'psf' | 'prmtop' | 'top'
    // volumes
    | 'map' | 'dx' | 'dxbin'
export const ParseFormatT = literal<ParseFormatT>(
    // trajectory
    'mmcif',
    'bcif', // +volumes
    'pdb',
    'pdbqt',
    'gro',
    'xyz',
    'mol',
    'sdf',
    'mol2',
    'lammpstrj', // + coordinates
    // coordinates
    'xtc',
    'nctraj',
    'dcd',
    'trr',
    // topology
    'psf',
    'prmtop',
    'top',
    // volumes
    'map',
    'dx',
    'dxbin',
);

/** `format` parameter values for `parse` node in Molstar tree */
export type MolstarParseFormatT =
    // trajectory
    | 'cif' | 'pdb' | 'pdbqt' | 'gro' | 'xyz' | 'mol' | 'sdf' | 'mol2' | 'lammpstrj'
    // coordinates
    | 'xtc' | 'nctraj' | 'dcd' | 'trr'
    // topology
    | 'psf' | 'prmtop' | 'top'
    // volumes
    | 'map' | 'dx' | 'dxbin'
export const MolstarParseFormatT = literal<MolstarParseFormatT>(
    // trajectory
    'cif', // +volumes
    'pdb',
    'pdbqt',
    'gro',
    'xyz',
    'mol',
    'sdf',
    'mol2',
    'lammpstrj',
    // coordinates
    'xtc',
    'nctraj',
    'dcd',
    'trr',
    // topology
    'psf',
    'prmtop',
    'top',
    // volumes
    'map',
    'dx',
    'dxbin',
);

/** `kind` parameter values for `structure` node in MVS tree */
export type StructureTypeT = 'model' | 'assembly' | 'symmetry' | 'symmetry_mates';
export const StructureTypeT = literal<StructureTypeT>('model', 'assembly', 'symmetry', 'symmetry_mates');

/** `selector` parameter values for `component` node in MVS tree */
export type ComponentSelectorT = 'all' | 'polymer' | 'protein' | 'nucleic' | 'branched' | 'ligand' | 'ion' | 'water' | 'coarse';
export const ComponentSelectorT = literal<ComponentSelectorT>('all', 'polymer', 'protein', 'nucleic', 'branched', 'ligand', 'ion', 'water', 'coarse');

/** `selector` parameter values for `component` node in MVS tree */
const _ComponentExpressionT = partial({
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
    /** 0-based residue index in the source file */
    residue_index: int, // TODO this is defined in Python builder but not supported by Molstar yet
    label_atom_id: str,
    auth_atom_id: str,
    type_symbol: str,
    atom_id: int,
    atom_index: int,
    /** Instance identifier to distinguish instances of the same chain created by applying different symmetry operators,
     * like 'ASM-X0-1' for assemblies or '1_555' for crystals */
    instance_id: str,
});
/** `selector` parameter values for `component` node in MVS tree */
export interface ComponentExpressionT extends ValueFor<typeof _ComponentExpressionT> { }
export const ComponentExpressionT: iots.Type<ComponentExpressionT> = _ComponentExpressionT;


/** `schema` parameter values for `*_from_uri` and `*_from_source` nodes in MVS tree */
export type SchemaT = 'whole_structure' | 'entity' | 'chain' | 'auth_chain' | 'residue' | 'auth_residue' | 'residue_range' | 'auth_residue_range' | 'atom' | 'auth_atom' | 'all_atomic';
export const SchemaT = literal<SchemaT>('whole_structure', 'entity', 'chain', 'auth_chain', 'residue', 'auth_residue', 'residue_range', 'auth_residue_range', 'atom', 'auth_atom', 'all_atomic');

/** `format` parameter values for `*_from_uri` nodes in MVS tree */
export type SchemaFormatT = 'cif' | 'bcif' | 'json';
export const SchemaFormatT = literal<SchemaFormatT>('cif', 'bcif', 'json');

/** Parameter values for vector params, e.g. `position` */
export type Vector3 = [number, number, number];
export const Vector3: iots.Type<Vector3> = tuple([float, float, float]);

/** Parameter values for matrix params, e.g. `rotation` */
export const Matrix = list(float); // TODO impl custom types Matrix3x3 and Matrix4x4

export type LabelAttachments = 'bottom-left' | 'bottom-center' | 'bottom-right' | 'middle-left' | 'middle-center' | 'middle-right' | 'top-left' | 'top-center' | 'top-right';
export const LabelAttachments = literal<LabelAttachments>('bottom-left', 'bottom-center', 'bottom-right', 'middle-left', 'middle-center', 'middle-right', 'top-left', 'top-center', 'top-right');

/** Primitives-related types */
const _PrimitiveComponentExpressionT = partial({
    structure_ref: str,
    expression_schema: SchemaT,
    expressions: list(ComponentExpressionT),
});
/** Primitives-related types */
export interface PrimitiveComponentExpressionT extends ValueFor<typeof _PrimitiveComponentExpressionT> { }
export const PrimitiveComponentExpressionT: iots.Type<PrimitiveComponentExpressionT> = _PrimitiveComponentExpressionT;

export const PrimitivePositionT = union(Vector3, ComponentExpressionT, PrimitiveComponentExpressionT);
export type PrimitivePositionT = ValueFor<typeof PrimitivePositionT>

export const FloatList = list(float);
export const IntList = list(int);
export const StrList = list(str);


/** Hexadecimal color string, e.g. '#FF1100' (the type matches more than just valid HexColor strings) */
export type HexColorT = `#${string}`;
export const HexColorT = new iots.Type<HexColorT>(
    'HexColor',
    ((value: any) => typeof value === 'string') as any,
    (value, ctx) => isHexColorT(value) ? { _tag: 'Right', right: value } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid hex color string` }] },
    value => value
);
/** Regular expression matching a hexadecimal color string, e.g. '#FF1100' or '#f10' */
const hexColorRegex = /^#([0-9A-F]{3}){1,2}$/i;
/** Decide if a string is a valid hexadecimal color string (6-digit or 3-digit, e.g. '#FF1100' or '#f10') */
function isHexColorT(str: any): str is HexColorT {
    return typeof str === 'string' && hexColorRegex.test(str);
}

/** Named color string (e.g. 'red') for `color` parameter values for `color` node in MVS tree */
export const ColorNameT = new iots.Type<ColorNameT>(
    'ColorName',
    ((value: any) => typeof value === 'string') as any,
    (value, ctx) => isColorNameT(value) ? { _tag: 'Right', right: value } : { _tag: 'Left', left: [{ value: value, context: ctx, message: `"${value}" is not a valid color name` }] },
    value => value
);
export type ColorNameT =
    | 'aliceblue' | 'antiquewhite' | 'aqua' | 'aquamarine' | 'azure' | 'beige' | 'bisque' | 'black'
    | 'blanchedalmond' | 'blue' | 'blueviolet' | 'brown' | 'burlywood' | 'cadetblue' | 'chartreuse'
    | 'chocolate' | 'coral' | 'cornflower' | 'cornflowerblue' | 'cornsilk' | 'crimson' | 'cyan' | 'darkblue'
    | 'darkcyan' | 'darkgoldenrod' | 'darkgray' | 'darkgreen' | 'darkgrey' | 'darkkhaki' | 'darkmagenta'
    | 'darkolivegreen' | 'darkorange' | 'darkorchid' | 'darkred' | 'darksalmon' | 'darkseagreen' | 'darkslateblue'
    | 'darkslategray' | 'darkslategrey' | 'darkturquoise' | 'darkviolet' | 'deeppink' | 'deepskyblue' | 'dimgray'
    | 'dimgrey' | 'dodgerblue' | 'firebrick' | 'floralwhite' | 'forestgreen' | 'fuchsia' | 'gainsboro'
    | 'ghostwhite' | 'gold' | 'goldenrod' | 'gray' | 'green' | 'greenyellow' | 'grey' | 'honeydew' | 'hotpink'
    | 'indianred' | 'indigo' | 'ivory' | 'khaki' | 'laserlemon' | 'lavender' | 'lavenderblush' | 'lawngreen'
    | 'lemonchiffon' | 'lightblue' | 'lightcoral' | 'lightcyan' | 'lightgoldenrod' | 'lightgoldenrodyellow'
    | 'lightgray' | 'lightgreen' | 'lightgrey' | 'lightpink' | 'lightsalmon' | 'lightseagreen' | 'lightskyblue'
    | 'lightslategray' | 'lightslategrey' | 'lightsteelblue' | 'lightyellow' | 'lime' | 'limegreen' | 'linen'
    | 'magenta' | 'maroon' | 'maroon2' | 'maroon3' | 'mediumaquamarine' | 'mediumblue' | 'mediumorchid' | 'mediumpurple'
    | 'mediumseagreen' | 'mediumslateblue' | 'mediumspringgreen' | 'mediumturquoise' | 'mediumvioletred' | 'midnightblue'
    | 'mintcream' | 'mistyrose' | 'moccasin' | 'navajowhite' | 'navy' | 'oldlace' | 'olive' | 'olivedrab' | 'orange'
    | 'orangered' | 'orchid' | 'palegoldenrod' | 'palegreen' | 'paleturquoise' | 'palevioletred' | 'papayawhip'
    | 'peachpuff' | 'peru' | 'pink' | 'plum' | 'powderblue' | 'purple' | 'purple2' | 'purple3' | 'rebeccapurple'
    | 'red' | 'rosybrown' | 'royalblue' | 'saddlebrown' | 'salmon' | 'sandybrown' | 'seagreen' | 'seashell' | 'sienna'
    | 'silver' | 'skyblue' | 'slateblue' | 'slategray' | 'slategrey' | 'snow' | 'springgreen' | 'steelblue' | 'tan'
    | 'teal' | 'thistle' | 'tomato' | 'turquoise' | 'violet' | 'wheat' | 'white' | 'whitesmoke' | 'yellow' | 'yellowgreen'
/** Decide if a string is a valid named color string */
function isColorNameT(str: any): str is ColorNameT {
    return str in ColorNames;
}

/** `color` parameter values for `color` node in MVS tree */
export type ColorT = ColorNameT | HexColorT;
export const ColorT: iots.Type<ColorT> = union(ColorNameT, HexColorT);


// Type helpers

export function isVector3(x: any): x is Vector3 {
    return !!x && Array.isArray(x) && x.length === 3 && typeof x[0] === 'number';
}

export function isPrimitiveComponentExpressions(x: any): x is PrimitiveComponentExpressionT {
    return !!x && Array.isArray(x.expressions);
}

export function isComponentExpression(x: any): x is ComponentExpressionT {
    return !!x && typeof x === 'object' && !x.expressions;
}


export type ColorListNameT =
    // Color lists from https://observablehq.com/@d3/color-schemes (definitions: https://colorbrewer2.org/export/colorbrewer.js)
    // Sequential single-hue
    | 'Reds' | 'Oranges' | 'Greens' | 'Blues' | 'Purples' | 'Greys'
    // Sequential multi-hue
    | 'OrRd' | 'BuGn' | 'PuBuGn' | 'GnBu' | 'PuBu' | 'BuPu' | 'RdPu' | 'PuRd' | 'YlOrRd' | 'YlOrBr' | 'YlGn' | 'YlGnBu'
    | 'Magma' | 'Inferno' | 'Plasma' | 'Viridis' | 'Cividis' | 'Turbo' | 'Warm' | 'Cool' | 'CubehelixDefault'
    // Cyclical
    | 'Rainbow' | 'Sinebow'
    // Diverging
    | 'RdBu' | 'RdGy' | 'PiYG' | 'BrBG' | 'PRGn' | 'PuOr' | 'RdYlGn' | 'RdYlBu' | 'Spectral'
    // Categorical
    | 'Category10' | 'Observable10' | 'Tableau10'
    | 'Set1' | 'Set2' | 'Set3' | 'Pastel1' | 'Pastel2' | 'Dark2' | 'Paired' | 'Accent'
    // Additional lists, not standard for visualization in general, but commonly used for structures
    | 'Chainbow'
export const ColorListNameT = literal<ColorListNameT>(
    // Color lists from https://observablehq.com/@d3/color-schemes (definitions: https://colorbrewer2.org/export/colorbrewer.js)
    // Sequential single-hue
    'Reds', 'Oranges', 'Greens', 'Blues', 'Purples', 'Greys',
    // Sequential multi-hue
    'OrRd', 'BuGn', 'PuBuGn', 'GnBu', 'PuBu', 'BuPu', 'RdPu', 'PuRd', 'YlOrRd', 'YlOrBr', 'YlGn', 'YlGnBu',
    'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Turbo', 'Warm', 'Cool', 'CubehelixDefault',
    // Cyclical
    'Rainbow', 'Sinebow',
    // Diverging
    'RdBu', 'RdGy', 'PiYG', 'BrBG', 'PRGn', 'PuOr', 'RdYlGn', 'RdYlBu', 'Spectral',
    // Categorical
    'Category10', 'Observable10', 'Tableau10',
    'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Dark2', 'Paired', 'Accent',

    // Additional lists, not standard for visualization in general, but commonly used for structures
    'Chainbow',
);

export type ColorDictNameT = 'ElementSymbol' | 'ResidueName' | 'ResidueProperties' | 'SecondaryStructure';
export const ColorDictNameT = literal<ColorDictNameT>('ElementSymbol', 'ResidueName', 'ResidueProperties', 'SecondaryStructure');


const _CategoricalPalette = object(
    {
        kind: literal('categorical'),
    },
    // Optionals:
    {
        colors: union(
            ColorListNameT,
            ColorDictNameT,
            list(ColorT),
            dict(str, ColorT),
        ),
        /** Repeat color list once all colors are depleted (only applies if `colors` is a list or a color list name). */
        repeat_color_list: bool,
        /** Sort actual annotation values before assigning colors from a list (none = take values in order of their first occurrence). */
        sort: literal('none', 'lexical', 'numeric'),
        /** Sort direction. */
        sort_direction: literal('ascending', 'descending'),
        /** Treat annotation values as case-insensitive strings. */
        case_insensitive: bool,
        /** Color to use when a) `colors` is a dictionary (or a color dictionary name) and given key is not present, or b) `colors` is a list (or a color list name) and there are more actual annotation values than listed colors and `repeat_color_list` is not true. */
        missing_color: nullable(ColorT),
    }
);
export interface CategoricalPalette extends ValueFor<typeof _CategoricalPalette> { }
export const CategoricalPalette: iots.Type<CategoricalPalette> = _CategoricalPalette;

export const CategoricalPaletteDefaults: Required<CategoricalPalette> = {
    kind: 'categorical',
    colors: 'Category10', // this is also default for categorical in Matplotlib
    repeat_color_list: false,
    sort: 'none',
    sort_direction: 'ascending',
    case_insensitive: false,
    missing_color: null,
};


const _DiscretePalette = object(
    {
        kind: literal('discrete'),
    },
    // Optionals:
    {
        /** Define colors for the discrete color palette and optionally corresponding checkpoints.
         * Checkpoints refer to the values normalized to interval [0, 1] if `mode` is `"normalized"` (default), or to the values directly if `mode` is `"absolute"`.
         * If checkpoints are not provided, they will created automatically (uniformly distributed over interval [0, 1]).
         * If 1 checkpoint is provided for each color, then the color applies to values from this checkpoint (inclusive) until the next listed checkpoint (exclusive); the last color applies until Infinity.
         * If 2 checkpoints are provided for each color, then the color applies to values from the first until the second checkpoint (inclusive); null means +/-Infinity; if ranges overlap, the later listed takes precedence.
         */
        colors: union(
            ColorListNameT,
            list(ColorT),
            list(tuple([ColorT, float])),
            list(tuple([nullable(ColorT), nullable(float), nullable(float)])),
        ),
        /** Reverse order of `colors` list. Only has effect when `colors` is a color list name or a color list without explicit checkpoints. */
        reverse: bool,
        /** Defines whether the annotation values should be normalized before assigning color based on checkpoints in `colors` (`x_normalized = (x - x_min) / (x_max - x_min)`, where `[x_min, x_max]` are either `value_domain` if provided, or the lowest and the highest value encountered in the annotation). Default is `"normalized"`. */
        mode: literal('normalized', 'absolute'),
        /** Defines `x_min` and `x_max` for normalization of annotation values. Either can be `null`, meaning that minimum/maximum of the actual values will be used. Only used when `mode` is `"normalized"`. */
        value_domain: tuple([nullable(float), nullable(float)]),
    }
);
export interface DiscretePalette extends ValueFor<typeof _DiscretePalette> { }
export const DiscretePalette: iots.Type<DiscretePalette> = _DiscretePalette;

export const DiscretePaletteDefaults: Required<DiscretePalette> = {
    kind: 'discrete',
    colors: 'YlGn', // YlGn was selected as default because (a) Matplotlib's default Viridis looks ugly in 3D and (b) YlGn does not contain white, so it's easier to see that it's doing something even when values are in wrong range
    reverse: false,
    mode: 'normalized',
    value_domain: [null, null],
};


const _ContinuousPalette = object(
    {
        kind: literal('continuous'),
    },
    // Optionals:
    {
        /** Define colors for the continuous color palette and optionally corresponding checkpoints (i.e. annotation values that are mapped to each color).
         * Checkpoints refer to the values normalized to interval [0, 1] if `mode` is `"normalized"` (default), or to the values directly if `mode` is `"absolute"`.
         * If checkpoints are not provided, they will created automatically (uniformly distributed over interval [0, 1]). */
        colors: union(
            ColorListNameT,
            list(ColorT),
            list(tuple([ColorT, float])),
        ),
        /** Reverse order of `colors` list. Only has effect when `colors` is a color list name or a color list without explicit checkpoints. */
        reverse: bool,
        /** Defines whether the annotation values should be normalized before assigning color based on checkpoints in `colors` (`x_normalized = (x - x_min) / (x_max - x_min)`, where `[x_min, x_max]` are either `value_domain` if provided, or the lowest and the highest value encountered in the annotation). Default is `"normalized"`. */
        mode: literal('normalized', 'absolute'),
        /** Defines `x_min` and `x_max` for normalization of annotation values. Either can be `null`, meaning that minimum/maximum of the actual values will be used. Only used when `mode` is `"normalized"`. */
        value_domain: tuple([nullable(float), nullable(float)]),
        /** Color to use for values below the lowest checkpoint. 'auto' means color of the lowest checkpoint. */
        underflow_color: nullable(union(literal('auto'), ColorT)),
        /** Color to use for values above the highest checkpoint. 'auto' means color of the highest checkpoint. */
        overflow_color: nullable(union(literal('auto'), ColorT)),
    }
);
export interface ContinuousPalette extends ValueFor<typeof _ContinuousPalette> { }
export const ContinuousPalette: iots.Type<ContinuousPalette> = _ContinuousPalette;

export const ContinuousPaletteDefaults: Required<ContinuousPalette> = {
    kind: 'continuous',
    colors: 'YlGn', // YlGn was selected as default because (a) Matplotlib's default Viridis looks ugly in 3D and (b) YlGn does not contain white, so it's easier to see that it's doing something even when values are in wrong range
    reverse: false,
    mode: 'normalized',
    value_domain: [null, null],
    underflow_color: null,
    overflow_color: null,
};

// TODO consider spreading the palette param directly into color_from_uri/color_from_source params (though this will be tricky)
// TODO consider implementing some kind of recursion for object-typed params to achieve smart error messages and default value handling

export const Palette = union(CategoricalPalette, DiscretePalette, ContinuousPalette);
