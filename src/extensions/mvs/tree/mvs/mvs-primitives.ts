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
    dash_start: nullable(float),
    dash_length: nullable(float),
    gap_length: nullable(float),
};

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
    wireframe_radius: nullable(float),
    wireframe_color: nullable(ColorT),
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

export const MVSPrimitiveParams = union([MeshParams, LineParams, DistanceMeasurementParams, PrimitiveLabelParams]);

export type MVSPrimitive = ValueFor<typeof MVSPrimitiveParams>
export type MVSPrimitiveKind = MVSPrimitive['kind']
export type MVSPrimitiveOptions = MVSNode<'primitives'>['params']
export type MVSPrimitiveParams<T extends MVSPrimitiveKind = MVSPrimitiveKind> = Extract<MVSPrimitive, { kind: T }>