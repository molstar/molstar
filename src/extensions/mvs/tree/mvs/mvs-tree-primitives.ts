/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { bool, float, int, mapping, nullable, OptionalField, RequiredField, str } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema, ValuesFor } from '../generic/params-schema';
import type { MVSNode } from './mvs-tree';
import { ColorT, FloatList, IntList, PrimitivePositionT } from './param-types';


const _TubeBase = {
    start: RequiredField(PrimitivePositionT, 'Start of this tube.'),
    end: RequiredField(PrimitivePositionT, 'End of this tube.'),
    radius: OptionalField(float, 0.05, 'Tube radius (in Angstroms).'),
    dash_length: OptionalField(nullable(float), null, 'Length of each dash.'),
    color: OptionalField(nullable(ColorT), null, 'Color of the tube. If not specified, the primitives group color is used.'),
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
    wireframe_width: OptionalField(float, 1, 'Wireframe line width'),
    wireframe_color: OptionalField(nullable(ColorT), null, 'Wireframe color, uses triangle/group colors when not set'),
    // TODO docstrings
    // TODO avoid pseudo-defaults in primitives.ts
};

const LinesParams = {
    vertices: RequiredField(FloatList, '3N length array of floats with vertex position (x1, y1, z1, ...)'),
    indices: RequiredField(IntList, '2N length array of indices into vertices that form lines (l1_1, ll1_2, ...)'),
    line_groups: OptionalField(nullable(IntList), null, 'Assign a number to each triangle to group them.'),
    group_colors: OptionalField(mapping(int, ColorT), {}, 'Assign a color to each group. If not assigned, default primitives group color is used. Takes precedence over line_colors.'),
    group_tooltips: OptionalField(mapping(int, str), {}, 'Assign an optional tooltip to each group.'),
    group_width: OptionalField(mapping(int, float), {}, 'Assign an optional line width to each group. Take precedence over `width`.'),
    tooltip: OptionalField(nullable(str), null, 'Tooltip shown when hovering over the lines. Assigned group_tooltips take precedence.'),
    color: OptionalField(nullable(ColorT), null, 'Default color of the lines.'),
    width: OptionalField(float, 1, 'Line width (in screen-space units)'),
};

const TubeParams = {
    ..._TubeBase,
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering on the tube.'),
};

const DistanceMeasurementParams = {
    ..._TubeBase,
    label_template: OptionalField(nullable(str), null, 'Template used to construct the label. Use {{distance}} as placeholder for the distance.'),
    label_size: OptionalField(nullable(float), null, 'Size of the label (text height in Angstroms). If not provided (null), size will be computed relative to the distance (see label_auto_size_scale, label_auto_size_min).'),
    label_auto_size_scale: OptionalField(float, 0.1, 'Scaling factor for auto size.'),
    label_auto_size_min: OptionalField(float, 0, 'Minimum size for auto size.'),
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label.'),
};

const PrimitiveLabelParams = {
    position: RequiredField(PrimitivePositionT, 'Position of this label.'),
    text: RequiredField(str, 'The label.'),
    label_size: OptionalField(float, 1, 'Size of the label.'),
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label.'),
    label_offset: OptionalField(float, 0, 'Camera-facing offset to prevent overlap with geometry.'),
};

export const MVSPrimitiveParams = UnionParamsSchema(
    'kind',
    'Kind of geometrical primitive',
    {
        'mesh': SimpleParamsSchema(MeshParams),
        'lines': SimpleParamsSchema(LinesParams),
        'tube': SimpleParamsSchema(TubeParams),
        'distance_measurement': SimpleParamsSchema(DistanceMeasurementParams),
        'label': SimpleParamsSchema(PrimitiveLabelParams),
    },
);

export type MVSPrimitive = ValuesFor<typeof MVSPrimitiveParams>
export type MVSPrimitiveKind = MVSPrimitive['kind']
export type MVSPrimitiveOptions = MVSNode<'primitives'>['params']
export type MVSPrimitiveParams<T extends MVSPrimitiveKind = MVSPrimitiveKind> = Extract<MVSPrimitive, { kind: T }>
