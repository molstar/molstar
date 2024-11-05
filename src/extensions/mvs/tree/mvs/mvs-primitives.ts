/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, int, literal, mapping, nullable, obj, str, union, ValueFor } from '../generic/params-schema';
import type { MVSNode } from './mvs-tree';
import { ColorT, FloatList, IntList, PrimitivePositionT, StrList } from './param-types';

// TODO: Figure out validation and default values for these

const _LineBase = {
    start: PrimitivePositionT,
    end: PrimitivePositionT,
    thickness: nullable(float),
    color: nullable(ColorT),
    dash_length: nullable(float),
};

const _TransformBase = {
    rotation: FloatList,
    translation: nullable(FloatList)
};

const CylinderParams = obj({
    kind: literal('cylinder'),
    rotation: FloatList,
    // color
    // color: Optional[ColorT] = Field(description="Default color for the box.")
    // radius_top: float = Field(description="The radius of the top of the cylinder top. Radius equal to zero will yield a cone.")
    // radius_bottom: float = Field(description="The radius of the bottom of the cylinder. Radius equal to zero will yield a reversed cone.")
    // height: float = Field(description="The height of the cone.")
    // # TODO: meaning of the following two? Check Sebastian's answers in some of the PRs.
    // theta_start: float = Field(description="TODO")
    // theta_length: float = Field(description="TODO")
    // # TODO: type for rotation as field
    // rotation: Optional[Mat3[float]] = Field(
    //     description="9d vector describing the rotation, in a column major (j * 3 + i indexing) format, this is equivalent to Fortran-order in numpy, to be multiplied from the left",
    // )
    // bottom_cap: bool = Field(description="Determine whether to cap the top of the cylinder.")
    // top_cap: bool = Field(description="Determine whether to cap the bottom of the cylinder.")
    
});

const BoxParams = obj({
    kind: literal('box'),
    ..._TransformBase,
    color: nullable(ColorT),
    center: PrimitivePositionT,
    extent: FloatList,
    scaling: nullable(FloatList),
    as_edges: nullable(bool),
    edge_radius: nullable(float),
    box_groups: nullable(IntList),
});

const MeshParams = obj({
    kind: literal('mesh'),
    vertices: FloatList,
    indices: IntList,
    triangle_colors: nullable(StrList),
    triangle_groups: nullable(IntList),
    group_colors: nullable(mapping(int, ColorT)),
    group_tooltips: nullable(mapping(int, str)),
    tooltip: nullable(str),
    show_triangles: nullable(bool),
    show_wireframe: nullable(bool),
    color: nullable(ColorT),
    wireframe_radius: nullable(float),
    wireframe_color: nullable(ColorT),
});

const LinesParams = obj({
    kind: literal('lines'),
    vertices: FloatList,
    indices: IntList,
    line_colors: nullable(StrList),
    line_groups: nullable(IntList),
    group_colors: nullable(mapping(int, ColorT)),
    group_tooltips: nullable(mapping(int, str)),
    group_radius: nullable(mapping(int, float)),
    tooltip: nullable(str),
    color: nullable(ColorT),
    line_radius: nullable(float),
});

const LineParams = obj({
    kind: literal('line'),
    ..._LineBase,
    tooltip: nullable(str),
});

const DistanceMeasurementParams = obj({
    kind: literal('distance_measurement'),
    ..._LineBase,
    label_template: nullable(str),
    label_size: nullable(union([float, literal('auto')])),
    label_auto_size_scale: nullable(float),
    label_auto_size_min: nullable(float),
    label_color: nullable(ColorT),
});

const PrimitiveLabelParams = obj({
    kind: literal('label'),
    position: PrimitivePositionT,
    text: str,
    label_size: nullable(float),
    label_color: nullable(ColorT),
    label_offset: nullable(float),
});

export const MVSPrimitiveParams = union([MeshParams, BoxParams, CylinderParams, LinesParams, LineParams, DistanceMeasurementParams, PrimitiveLabelParams]);

export type MVSPrimitive = ValueFor<typeof MVSPrimitiveParams>
export type MVSPrimitiveKind = MVSPrimitive['kind']
export type MVSPrimitiveOptions = MVSNode<'primitives'>['params']
export type MVSPrimitiveParams<T extends MVSPrimitiveKind = MVSPrimitiveKind> = Extract<MVSPrimitive, { kind: T }>