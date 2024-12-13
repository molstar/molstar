/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { bool, float, int, literal, mapping, nullable, OptionalField, RequiredField, str, union } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema, ValuesFor } from '../generic/params-schema';
import type { MVSNode } from './mvs-tree';
import { ColorT, FloatList, IntList, PrimitivePositionT, StrList } from './param-types';


// TODO ensure every null color/label_color falls back to the parent 'primitives' node color/label_color + mention in docs

const _LineBase = {
    start: RequiredField(PrimitivePositionT, 'Start of this line.'),
    end: RequiredField(PrimitivePositionT, 'End of this line.'),
    thickness: OptionalField(float, 0.05, 'Thickness of this line.'),
    dash_length: OptionalField(nullable(float), null, 'Length of each dash.'),
    color: OptionalField(nullable(ColorT), null, 'Color of the line. If not specified, the primitives group color is used.'),
    // TODO naming inconsistency thickness vs radius in LinesParams
    // TODO thickness seems to be in 3D coordinates, while LinesParams.line_radius is in screen coordinates (does not get thicker with zoom)
};

const MeshParams = {
    vertices: RequiredField(FloatList, '3N length array of floats with vertex position (x1, y1, z1, ...).'),
    indices: RequiredField(IntList, '3N length array of indices into vertices that form triangles (t1_1, t1_2, t1_3, ...).'),
    triangle_groups: OptionalField(nullable(IntList), null, 'Assign a number to each triangle to group them.'),
    group_colors: OptionalField(mapping(int, ColorT), {}, 'Assign a color to each group. If not assigned, default primitives group color is used.'),
    group_tooltips: OptionalField(mapping(int, str), {}, 'Assign an optional tooltip to each group.'),
    color: OptionalField(nullable(ColorT), null, 'Default color of the triangles and wireframe.'),
    tooltip: OptionalField(nullable(str), null, 'Tooltip shown when hovering over the mesh. Assigned `group_tooltips` take precedence.'),
    show_triangles: OptionalField(bool, true, 'Determine whether to render triangles of the mesh.'),
    show_wireframe: OptionalField(bool, false, 'Determine whether to render wireframe of the mesh.'),
    wireframe_radius: OptionalField(float, 1, 'Wireframe line radius'),
    wireframe_color: OptionalField(nullable(ColorT), null, 'Wireframe color, uses triangle/group colors when not set'),
    // TODO docstrings
};

const LinesParams = {
    vertices: RequiredField(FloatList, '3N length array of floats with vertex position (x1, y1, z1, ...)'),
    indices: RequiredField(IntList, '2N length array of indices into vertices that form lines (l1_1, ll1_2, ...)'),
    line_groups: OptionalField(nullable(IntList), null, 'Assign a number to each triangle to group them.'),
    group_colors: OptionalField(mapping(int, ColorT), {}, 'Assign a color to each group. If not assigned, default primitives group color is used. Takes precedence over line_colors.'),
    group_tooltips: OptionalField(mapping(int, str), {}, 'Assign an optional tooltip to each group.'),
    group_radius: OptionalField(mapping(int, float), {}, 'Assign an optional radius to each group. Take precedence over line_radius.'),
    line_colors: OptionalField(nullable(StrList), null, 'Assign a color to each line.'),
    tooltip: OptionalField(nullable(str), null, 'Tooltip shown when hovering over the lines. Assigned group_tooltips take precedence.'),
    color: OptionalField(nullable(ColorT), null, 'Default color of the lines.'),
    line_radius: OptionalField(float, 1, 'Line radius'),
    // TODO line_radius - would make more sense to call it just `radius`?
    // TODO remove line_colors, instead have implicit groups?
};

const LineParams = {
    ..._LineBase,
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering on the line.'),
};

const DistanceMeasurementParams = {
    ..._LineBase,
    label_template: OptionalField(nullable(str), null, 'Template used to construct the label. Use {{distance}} as placeholder for the distance.'),
    label_size: OptionalField(union([float, literal('auto')]), 'auto', 'Size of the label. Auto scales it by the distance.'),
    label_auto_size_scale: OptionalField(float, 0.1, 'Scaling factor for auto size.'),
    label_auto_size_min: OptionalField(float, 0, 'Minimum size for auto size.'),
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label.'),
    // TODO think about merging 'line' and 'distance_measurement', clearly 'line' is just a 'distance_measurement' with label_template=null
    // TODO dash_start, gap_length are in Python builder function but do not exist in MVS nodes!
    // TODO {{distance}} being in Angstroms might be quite restrictive for future usecases (cell-scale objects) (could have param distance_unit: number|'A'|'nm'|'um'... ?)
    // TODO label_size could have just null instead of 'auto'
    // TODO label_size, label_auto_size_min -> what size in what units? text height in Angstroms or smth else?
};

const PrimitiveLabelParams = {
    position: RequiredField(PrimitivePositionT, 'Position of this label.'),
    text: RequiredField(str, 'The label.'),
    label_size: OptionalField(float, 1, 'Size of the label.'),
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label.'),
    label_offset: OptionalField(float, 0, 'Camera-facing offset to prevent overlap with geometry.'),
};

export const MVSPrimitiveParams = UnionParamsSchema('kind', {
    'mesh': SimpleParamsSchema(MeshParams),
    'lines': SimpleParamsSchema(LinesParams),
    'line': SimpleParamsSchema(LineParams),
    'distance_measurement': SimpleParamsSchema(DistanceMeasurementParams),
    'label': SimpleParamsSchema(PrimitiveLabelParams),
});

export type MVSPrimitive = ValuesFor<typeof MVSPrimitiveParams>
export type MVSPrimitiveKind = MVSPrimitive['kind']
export type MVSPrimitiveOptions = MVSNode<'primitives'>['params']
export type MVSPrimitiveParams<T extends MVSPrimitiveKind = MVSPrimitiveKind> = Extract<MVSPrimitive, { kind: T }>
