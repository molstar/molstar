/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { bool, float, int, mapping, nullable, OptionalField, RequiredField, str } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { ColorT, FloatList, IntList, PrimitivePositionT } from './param-types';


const _TubeBase = {
    /** Start point of the tube. */
    start: RequiredField(PrimitivePositionT, 'Start point of the tube.'),
    /** End point of the tube. */
    end: RequiredField(PrimitivePositionT, 'End point of the tube.'),
    /** Tube radius (in Angstroms). */
    radius: OptionalField(float, 0.05, 'Tube radius (in Angstroms).'),
    /** Length of each dash and gap between dashes. If not specified (null), draw full line. */
    dash_length: OptionalField(nullable(float), null, 'Length of each dash and gap between dashes. If not specified (null), draw full line.'),
    /** Color of the tube. If not specified, uses the parent primitives group `color`. */
    color: OptionalField(nullable(ColorT), null, 'Color of the tube. If not specified, uses the parent primitives group `color`.'),
};

const MeshParams = {
    /** 3*n_vertices length array of floats with vertex position (x1, y1, z1, ...). */
    vertices: RequiredField(FloatList, '3*n_vertices length array of floats with vertex position (x1, y1, z1, ...).'),
    /** 3*n_triangles length array of indices into vertices that form triangles (t1_1, t1_2, t1_3, ...). */
    indices: RequiredField(IntList, '3*n_triangles length array of indices into vertices that form triangles (t1_1, t1_2, t1_3, ...).'),
    /** Assign a number to each triangle to group them. If not specified, each triangle is considered a separate group (triangle i = group i). */
    triangle_groups: OptionalField(nullable(IntList), null, 'Assign a number to each triangle to group them. If not specified, each triangle is considered a separate group (triangle i = group i).'),
    /** Assign a color to each group. Where not assigned, uses `color`. */
    group_colors: OptionalField(mapping(int, ColorT), {}, 'Assign a color to each group. Where not assigned, uses `color`.'),
    /** Assign a tooltip to each group. Where not assigned, uses `tooltip`. */
    group_tooltips: OptionalField(mapping(int, str), {}, 'Assign a tooltip to each group. Where not assigned, uses `tooltip`.'),
    /** Color of the triangles and wireframe. Can be overwritten by `group_colors`. If not specified, uses the parent primitives group `color`. */
    color: OptionalField(nullable(ColorT), null, 'Color of the triangles and wireframe. Can be overwritten by `group_colors`. If not specified, uses the parent primitives group `color`.'),
    /** Tooltip shown when hovering over the mesh. Can be overwritten by `group_tooltips`. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip shown when hovering over the mesh. Can be overwritten by `group_tooltips`. If not specified, uses the parent primitives group `tooltip`.'),
    /** Determine whether to render triangles of the mesh. */
    show_triangles: OptionalField(bool, true, 'Determine whether to render triangles of the mesh.'),
    /** Determine whether to render wireframe of the mesh. */
    show_wireframe: OptionalField(bool, false, 'Determine whether to render wireframe of the mesh.'),
    /** Wireframe line width (in screen-space units). */
    wireframe_width: OptionalField(float, 1, 'Wireframe line width (in screen-space units).'),
    /** Wireframe color. If not specified, uses `group_colors`. */
    wireframe_color: OptionalField(nullable(ColorT), null, 'Wireframe color. If not specified, uses `group_colors`.'),
};

const LinesParams = {
    /** 3*n_vertices length array of floats with vertex position (x1, y1, z1, ...). */
    vertices: RequiredField(FloatList, '3*n_vertices length array of floats with vertex position (x1, y1, z1, ...).'),
    /** 2*n_lines length array of indices into vertices that form lines (l1_1, l1_2, ...). */
    indices: RequiredField(IntList, '2*n_lines length array of indices into vertices that form lines (l1_1, l1_2, ...).'),
    /** Assign a number to each triangle to group them. If not specified, each line is considered a separate group (line i = group i). */
    line_groups: OptionalField(nullable(IntList), null, 'Assign a number to each triangle to group them. If not specified, each line is considered a separate group (line i = group i).'),
    /** Assign a color to each group. Where not assigned, uses `color`. */
    group_colors: OptionalField(mapping(int, ColorT), {}, 'Assign a color to each group. Where not assigned, uses `color`.'),
    /** Assign a tooltip to each group. Where not assigned, uses `tooltip`. */
    group_tooltips: OptionalField(mapping(int, str), {}, 'Assign a tooltip to each group. Where not assigned, uses `tooltip`.'),
    /** Assign a line width to each group. Where not assigned, uses `width`. */
    group_widths: OptionalField(mapping(int, float), {}, 'Assign a line width to each group. Where not assigned, uses `width`.'),
    /** Color of the lines. Can be overwritten by `group_colors`. If not specified, uses the parent primitives group `color`. */
    color: OptionalField(nullable(ColorT), null, 'Color of the lines. Can be overwritten by `group_colors`. If not specified, uses the parent primitives group `color`.'),
    /** Tooltip shown when hovering over the lines. Can be overwritten by `group_tooltips`. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip shown when hovering over the lines. Can be overwritten by `group_tooltips`. If not specified, uses the parent primitives group `tooltip`.'),
    /** Line width (in screen-space units). Can be overwritten by `group_widths`. */
    width: OptionalField(float, 1, 'Line width (in screen-space units). Can be overwritten by `group_widths`.'),
};

const TubeParams = {
    ..._TubeBase,
    /** Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`.'),
};

const DistanceMeasurementParams = {
    ..._TubeBase,
    /** Template used to construct the label. Use {{distance}} as placeholder for the distance. */
    label_template: OptionalField(str, '{{distance}}', 'Template used to construct the label. Use {{distance}} as placeholder for the distance.'),
    /** Size of the label (text height in Angstroms). If not specified, size will be relative to the distance (see label_auto_size_scale, label_auto_size_min). */
    label_size: OptionalField(nullable(float), null, 'Size of the label (text height in Angstroms). If not specified, size will be relative to the distance (see label_auto_size_scale, label_auto_size_min).'),
    /** Scaling factor for relative size. */
    label_auto_size_scale: OptionalField(float, 0.1, 'Scaling factor for relative size.'),
    /** Minimum size for relative size. */
    label_auto_size_min: OptionalField(float, 0, 'Minimum size for relative size.'),
    /** Color of the label. If not specified, uses the parent primitives group `label_color`. */
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label. If not specified, uses the parent primitives group `label_color`.'),
};

const PrimitiveLabelParams = {
    /** Position of this label. */
    position: RequiredField(PrimitivePositionT, 'Position of this label.'),
    /** The label. */
    text: RequiredField(str, 'The label.'),
    /** Size of the label (text height in Angstroms). */
    label_size: OptionalField(float, 1, 'Size of the label (text height in Angstroms).'),
    /** Color of the label. If not specified, uses the parent primitives group `label_color`. */
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label. If not specified, uses the parent primitives group `label_color`.'),
    /** Camera-facing offset to prevent overlap with geometry. */
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
