/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { bool, dict, float, int, nullable, OptionalField, RequiredField, str, union } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { ColorT, FloatList, IntList, PrimitivePositionT, Vector3 } from './param-types';


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
    group_colors: OptionalField(dict(int, ColorT), {}, 'Assign a color to each group. Where not assigned, uses `color`.'),
    /** Assign a tooltip to each group. Where not assigned, uses `tooltip`. */
    group_tooltips: OptionalField(dict(int, str), {}, 'Assign a tooltip to each group. Where not assigned, uses `tooltip`.'),
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
    group_colors: OptionalField(dict(int, ColorT), {}, 'Assign a color to each group. Where not assigned, uses `color`.'),
    /** Assign a tooltip to each group. Where not assigned, uses `tooltip`. */
    group_tooltips: OptionalField(dict(int, str), {}, 'Assign a tooltip to each group. Where not assigned, uses `tooltip`.'),
    /** Assign a line width to each group. Where not assigned, uses `width`. */
    group_widths: OptionalField(dict(int, float), {}, 'Assign a line width to each group. Where not assigned, uses `width`.'),
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

const ArrowParams = {
    /** Start point of the tube. */
    start: RequiredField(PrimitivePositionT, 'Start point of the arrow.'),
    /** End point of the tube. */
    end: OptionalField(nullable(PrimitivePositionT), null, 'End point of the arrow.'),
    /** If specified, the endpoint is computed as start + direction. */
    direction: OptionalField(nullable(Vector3), null, 'If specified, the endpoint is computed as start + direction.'),
    /** Length of the arrow. If unset, the distance between start and end is used. */
    length: OptionalField(nullable(float), null, 'Length of the arrow. If unset, the distance between start and end is used.'),
    /** Draw a cap at the start of the arrow. */
    show_start_cap: OptionalField(bool, false, 'Draw a cap at the start of the arrow.'),
    /** Length of the start cap. */
    start_cap_length: OptionalField(nullable(float), null, 'Length of the start cap. If not provided, will be 2 * start_cap_radius.'),
    /** Radius of the start cap. */
    start_cap_radius: OptionalField(nullable(float), null, 'Radius of the start cap. If not provided, will be 2 * tube_radius.'),
    /** Draw an arrow at the end of the arrow. */
    show_end_cap: OptionalField(bool, false, 'Draw a cap at the end of the arrow.'),
    /** Height of the arrow at the end. */
    end_cap_length: OptionalField(nullable(float), null, 'Length of the end cap. If not provided, will be 2 * end_cap_radius.'),
    /** Radius of the arrow at the end. */
    end_cap_radius: OptionalField(nullable(float), null, 'Radius of the end cap. If not provided, will be 2 * tube_radius.'),
    /** Draw a tube connecting the start and end points. */
    show_tube: OptionalField(bool, true, 'Draw a tube connecting the start and end points.'),
    /** Tube radius (in Angstroms). */
    tube_radius: OptionalField(float, 0.05, 'Tube radius (in Angstroms).'),
    /** Length of each dash and gap between dashes. If not specified (null), draw full line. */
    tube_dash_length: OptionalField(nullable(float), null, 'Length of each dash and gap between dashes. If not specified (null), draw full line.'),
    /** Color of the tube. If not specified, uses the parent primitives group `color`. */
    color: OptionalField(nullable(ColorT), null, 'Color of the tube. If not specified, uses the parent primitives group `color`.'),
    /** Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering over the arrow. If not specified, uses the parent primitives group `tooltip`.'),
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

const AngleMeasurementParams = {
    /** Point A. */
    a: RequiredField(PrimitivePositionT, 'Point A.'),
    /** Point B. */
    b: RequiredField(PrimitivePositionT, 'Point B.'),
    /** Point C. */
    c: RequiredField(PrimitivePositionT, 'Point C.'),
    /** Template used to construct the label. Use {{angle}} as placeholder for the angle in radians. */
    label_template: OptionalField(str, '{{angle}}', 'Template used to construct the label. Use {{angle}} as placeholder for the angle in radians.'),
    /** Size of the label (text height in Angstroms). If not specified, size will be relative to the distance (see label_auto_size_scale, label_auto_size_min). */
    label_size: OptionalField(nullable(float), null, 'Size of the label (text height in Angstroms). If not specified, size will be relative to the distance (see label_auto_size_scale, label_auto_size_min).'),
    /** Scaling factor for relative size. */
    label_auto_size_scale: OptionalField(float, 0.33, 'Scaling factor for relative size.'),
    /** Minimum size for relative size. */
    label_auto_size_min: OptionalField(float, 0, 'Minimum size for relative size.'),
    /** Color of the label. If not specified, uses the parent primitives group `label_color`. */
    label_color: OptionalField(nullable(ColorT), null, 'Color of the label. If not specified, uses the parent primitives group `label_color`.'),
    /** Draw vectors between (a, b) and (b, c). */
    show_vector: OptionalField(bool, true, 'Draw vectors between (a, b) and (b, c).'),
    /** Color of the vectors. */
    vector_color: OptionalField(nullable(ColorT), null, 'Color of the vectors.'),
    /** Radius of the vectors. */
    vector_radius: OptionalField(float, 0.05, 'Radius of the vectors.'),
    /** Draw a filled circle section representing the angle. */
    show_section: OptionalField(bool, true, 'Draw a filled circle section representing the angle.'),
    /** Color of the angle section. If not specified, the primitives group color is used. */
    section_color: OptionalField(nullable(ColorT), null, 'Color of the angle section. If not specified, the primitives group color is used.'),
    /** Radius of the angle section. In angstroms. */
    section_radius: OptionalField(nullable(float), null, 'Radius of the angle section. In angstroms.'),
    /** Factor to scale the radius of the angle section. Ignored if section_radius is set. */
    section_radius_scale: OptionalField(float, 0.33, 'Factor to scale the radius of the angle section. Ignored if section_radius is set.'),
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

const EllipseParams = {
    /** Color of the primitive. If not specified, uses the parent primitives group `color`. */
    color: OptionalField(nullable(ColorT), null, 'Color of the ellipse. If not specified, uses the parent primitives group `color`.'),
    /** If true, ignores radius_minor/magnitude of the minor axis */
    as_circle: OptionalField(bool, false, 'If true, ignores radius_minor/magnitude of the minor axis.'),
    /** ellipse center. */
    center: RequiredField(PrimitivePositionT, 'The center of the ellipse.'),
    /** Major axis of this ellipse. */
    major_axis: OptionalField(nullable(Vector3), null, 'Major axis of this ellipse.'),
    /** Minor axis of this ellipse. */
    minor_axis: OptionalField(nullable(Vector3), null, 'Minor axis of this ellipse.'),
    /** Major axis endpoint. If specified, overrides major axis to be major_axis_endpoint - center. */
    major_axis_endpoint: OptionalField(nullable(PrimitivePositionT), null, 'Major axis endpoint. If specified, overrides major axis to be major_axis_endpoint - center.'),
    /** Minor axis endpoint. If specified, overrides minor axis to be minor_axis_endpoint - center. */
    minor_axis_endpoint: OptionalField(nullable(PrimitivePositionT), null, 'Minor axis endpoint. If specified, overrides minor axis to be minor_axis_endpoint - center.'),
    /** Radius of the major axis. If unset, the length of the major axis is used. */
    radius_major: OptionalField(nullable(float), null, 'Radius of the major axis. If unset, the length of the major axis is used.'),
    /** Radius of the minor axis. If unset, the length of the minor axis is used. */
    radius_minor: OptionalField(nullable(float), null, 'Radius of the minor axis. If unset, the length of the minor axis is used.'),
    /** Start of the arc. In radians */
    theta_start: OptionalField(float, 0, 'Start of the arc. In radians'),
    /** End of the arc. In radians */
    theta_end: OptionalField(float, 2 * Math.PI, 'End of the arc. In radians'),
    /** Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`.'),
};

const EllipsoidParams = {
    /** Color of the primitive. If not specified, uses the parent primitives group `color`. */
    color: OptionalField(nullable(ColorT), null, 'Color of the ellipsoid. If not specified, uses the parent primitives group `color`.'),
    /** Ellipsoid center. */
    center: RequiredField(PrimitivePositionT, 'The center of the ellipsoid.'),
    /** Major axis of this ellipsoid. */
    major_axis: OptionalField(nullable(Vector3), null, 'Major axis of this ellipsoid.'),
    /** Minor axis of this ellipsoid. */
    minor_axis: OptionalField(nullable(Vector3), null, 'Minor axis of this ellipsoid.'),
    /** Major axis endpoint. If specified, overrides major axis to be major_axis_endpoint - center. */
    major_axis_endpoint: OptionalField(nullable(PrimitivePositionT), null, 'Major axis endpoint. If specified, overrides major axis to be major_axis_endpoint - center.'),
    /** Minor axis endpoint. If specified, overrides minor axis to be minor_axis_endpoint - center. */
    minor_axis_endpoint: OptionalField(nullable(PrimitivePositionT), null, 'Minor axis endpoint. If specified, overrides minor axis to be minor_axis_endpoint - center.'),
    /** Radii of the ellipsoid along each axis. */
    radius: OptionalField(nullable(union(Vector3, float)), null, 'Radii of the ellipsoid along each axis.'),
    /** Added to the radii of the ellipsoid along each axis. */
    radius_extent: OptionalField(nullable(union(Vector3, float)), null, 'Added to the radii of the ellipsoid along each axis.'),
    /** Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`.'),
};

const BoxParams = {
    /** The center of the box. */
    center: RequiredField(PrimitivePositionT, 'The center of the box.'),
    /** The width, the height, and the depth of the box. Added to the bounding box determined by the center. */
    extent: OptionalField(nullable(Vector3), null, 'The width, the height, and the depth of the box. Added to the bounding box determined by the center.'),
    /** Determine whether to render the faces of the box. */
    show_faces: OptionalField(bool, true, 'Determine whether to render the faces of the box.'),
    /** Color of the box faces. */
    face_color: OptionalField(nullable(ColorT), null, 'Color of the box faces.'),
    /** Determine whether to render the edges of the box. */
    show_edges: OptionalField(bool, false, 'Determine whether to render the edges of the box.'),
    /** Radius of the box edges. In angstroms. */
    edge_radius: OptionalField(float, 0.1, 'Radius of the box edges. In angstroms.'),
    /** Color of the box edges. */
    edge_color: OptionalField(nullable(ColorT), null, 'Color of the edges.'),
    /** Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`. */
    tooltip: OptionalField(nullable(str), null, 'Tooltip to show when hovering over the tube. If not specified, uses the parent primitives group `tooltip`.'),
};

export const MVSPrimitiveParams = UnionParamsSchema(
    'kind',
    'Kind of geometrical primitive',
    {
        'mesh': SimpleParamsSchema(MeshParams),
        'lines': SimpleParamsSchema(LinesParams),
        'tube': SimpleParamsSchema(TubeParams),
        'arrow': SimpleParamsSchema(ArrowParams),
        'distance_measurement': SimpleParamsSchema(DistanceMeasurementParams),
        'angle_measurement': SimpleParamsSchema(AngleMeasurementParams),
        'label': SimpleParamsSchema(PrimitiveLabelParams),
        'ellipse': SimpleParamsSchema(EllipseParams),
        'ellipsoid': SimpleParamsSchema(EllipsoidParams),
        'box': SimpleParamsSchema(BoxParams),
    },
);
