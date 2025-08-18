/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, int, list, OptionalField, RequiredField, str, union, nullable, literal, ValueFor } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { NodeFor, ParamsOfKind, SubtreeOfKind, TreeFor, TreeSchema } from '../generic/tree-schema';
import { ColorT, ContinuousPalette, DiscretePalette, Matrix, Vector3 } from '../mvs/param-types';

const Easing = literal(
    'linear',
    'bounce-in', 'bounce-out', 'bounce-in-out',
    'circle-in', 'circle-out', 'circle-in-out',
    'cubic-in', 'cubic-out', 'cubic-in-out',
    'exp-in', 'exp-out', 'exp-in-out',
    'quad-in', 'quad-out', 'quad-in-out',
    'sin-in', 'sin-out', 'sin-in-out',
);

export type MVSAnimationEasing = ValueFor<typeof Easing>;

const _Noise = {
    noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the interpolated value.')
    // support cummulative noise?
};

const _Common = {
    target_ref: RequiredField(str, 'Reference to the node.'),
    property: RequiredField(union(str, list(union(str, int))), 'Value accessor.'),
    start_ms: OptionalField(float, 0, 'Start time of the transition in milliseconds.'),
    duration_ms: RequiredField(float, 'Duration of the transition in milliseconds.'),
};

const _Frequency = {
    frequency: OptionalField(int, 1, 'Determines how many times the interpolation loops. Current T = frequency * t mod 1.'),
    alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
};

const _Easing = {
    easing: OptionalField(Easing, 'linear', 'Easing function to use for the transition.'),
};

const ScalarInterpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    start: OptionalField(nullable(union(float, list(float))), null, 'Start value. If a list of values is provided, each element will be interpolated separately. If unset, parent state value is used.'),
    end: OptionalField(nullable(union(float, list(float))), null, 'End value. If a list of values is provided, each element will be interpolated separately. If unset, only noise is applied.'),
    ..._Noise,
};

const Vec3Interpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    start: OptionalField(nullable(list(float)), null, 'Start value. If unset, parent state value is used. Must be array of length 3N (x1, y1, z1, x2, y2, z2, ...).'),
    end: OptionalField(nullable(list(float)), null, 'End value. Must be array of length 3N (x1, y1, z1, x2, y2, z2, ...). If unset, only noise is applied.'),
    spherical: OptionalField(bool, false, 'Whether to use spherical interpolation.'),
    ..._Noise,
};

const RotationMatrixInterpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    start: OptionalField(nullable(Matrix), null, 'Start value. If unset, parent state value is used.'),
    end: OptionalField(nullable(Matrix), null, 'End value. If unset, only noise is applied.'),
    ..._Noise,
};

const ColorInterpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    start: OptionalField(nullable(ColorT), null, 'Start value. If unset, parent state value is used.'),
    end: OptionalField(nullable(ColorT), null, 'End value.'),
    palette: OptionalField(nullable(union(DiscretePalette, ContinuousPalette)), null, 'Palette to sample colors from. Overrides start and end values.'),
};

const TransformationMatrixInterpolation = {
    ..._Common,
    pivot: OptionalField(nullable(Vector3), null, 'Pivot point for rotation and scale.'),
    rotation_start: OptionalField(nullable(Matrix), null, 'Start rotation value. If unset, parent state value is used.'),
    rotation_end: OptionalField(nullable(Matrix), null, 'End rotation value. If unset, only noise is applied.'),
    rotation_noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the rotation.'),
    rotation_easing: OptionalField(Easing, 'linear', 'Easing function to use for the rotation.'),
    rotation_frequency: OptionalField(int, 1, 'Determines how many times the rotation interpolation loops. Current T = frequency * t mod 1.'),
    rotation_alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
    translation_start: OptionalField(nullable(Vector3), null, 'Start translation value. If unset, parent state value is used.'),
    translation_end: OptionalField(nullable(Vector3), null, 'End translation value. If unset, only noise is applied.'),
    translation_noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the translation.'),
    translation_easing: OptionalField(Easing, 'linear', 'Easing function to use for the translation.'),
    translation_frequency: OptionalField(int, 1, 'Determines how many times the translation interpolation loops. Current T = frequency * t mod 1.'),
    translation_alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
    scale_start: OptionalField(nullable(Vector3), null, 'Start scale value. If unset, parent state value is used.'),
    scale_end: OptionalField(nullable(Vector3), null, 'End scale value. If unset, only noise is applied.'),
    scale_noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the scale.'),
    scale_easing: OptionalField(Easing, 'linear', 'Easing function to use for the scale.'),
    scale_frequency: OptionalField(int, 1, 'Determines how many times the scale interpolation loops. Current T = frequency * t mod 1.'),
    scale_alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
};

export const MVSAnimationSchema = TreeSchema({
    rootKind: 'animation',
    nodes: {
        animation: {
            description: 'Animation root node',
            parent: [],
            params: SimpleParamsSchema({
                frame_time_ms: OptionalField(float, 1000 / 60, 'Frame time in milliseconds'),
                duration_ms: OptionalField(nullable(float), null, 'Total duration of the animation. If not specified, computed as maximum of all transitions.'),
                autoplay: OptionalField(bool, true, 'Determines whether the animation should autoplay when a snapshot is loaded'),
                loop: OptionalField(bool, false, 'Determines whether the animation should loop when it reaches the end'),
                include_camera: OptionalField(bool, false, 'Determines whether the camera state should be included in the animation'),
                include_canvas: OptionalField(bool, false, 'Determines whether the canvas state should be included in the animation'),
            }),
        },
        interpolate: {
            description: 'This node enables interpolating between values',
            parent: ['animation'],
            params: UnionParamsSchema(
                'kind',
                'Interpolation kind',
                {
                    scalar: SimpleParamsSchema(ScalarInterpolation),
                    vec3: SimpleParamsSchema(Vec3Interpolation),
                    rotation_matrix: SimpleParamsSchema(RotationMatrixInterpolation),
                    transform_matrix: SimpleParamsSchema(TransformationMatrixInterpolation),
                    color: SimpleParamsSchema(ColorInterpolation),
                },
            )
        }
    }
});

export type MVSAnimationKind = keyof typeof MVSAnimationSchema.nodes
export type MVSAnimationNode<TKind extends MVSAnimationKind = MVSAnimationKind> = NodeFor<typeof MVSAnimationSchema, TKind>
export type MVSAnimationTree = TreeFor<typeof MVSAnimationSchema>
export type MVSAnimationNodeParams<TKind extends MVSAnimationKind> = ParamsOfKind<MVSAnimationTree, TKind>
export type MVSAnimationSubtree<TKind extends MVSAnimationKind = MVSAnimationKind> = SubtreeOfKind<MVSAnimationTree, TKind>