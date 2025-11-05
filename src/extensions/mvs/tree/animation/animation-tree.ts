/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, int, list, OptionalField, RequiredField, str, union, nullable, literal, ValueFor, dict } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { NodeFor, ParamsOfKind, SubtreeOfKind, TreeFor, TreeSchema } from '../generic/tree-schema';
import { ColorT, ContinuousPalette, DiscretePalette, Matrix, Vector3 } from '../mvs/param-types';

type Easing =
    | 'linear'
    | 'bounce-in' | 'bounce-out' | 'bounce-in-out'
    | 'circle-in' | 'circle-out' | 'circle-in-out'
    | 'cubic-in' | 'cubic-out' | 'cubic-in-out'
    | 'exp-in' | 'exp-out' | 'exp-in-out'
    | 'quad-in' | 'quad-out' | 'quad-in-out'
    | 'sin-in' | 'sin-out' | 'sin-in-out'
const Easing = literal<Easing>(
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
    /** Magnitude of the noise to apply to the interpolated value. */
    noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the interpolated value.')
    // support cummulative noise?
};

const _Common = {
    /** Reference to the node. */
    target_ref: RequiredField(str, 'Reference to the node.'),
    /** Value accessor. */
    property: RequiredField(union(str, list(union(str, int))), 'Value accessor.'),
    /** Start time of the transition in milliseconds. */
    start_ms: OptionalField(float, 0, 'Start time of the transition in milliseconds.'),
    /** Duration of the transition in milliseconds. */
    duration_ms: RequiredField(float, 'Duration of the transition in milliseconds.'),
};

const _Frequency = {
    /** Determines how many times the interpolation loops. Current T = frequency * t mod 1. */
    frequency: OptionalField(int, 1, 'Determines how many times the interpolation loops. Current T = frequency * t mod 1.'),
    /** Whether to alternate the direction of the interpolation for frequency > 1. */
    alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
};

const _Easing = {
    /** Easing function to use for the transition. */
    easing: OptionalField(Easing, 'linear', 'Easing function to use for the transition.'),
};

const ScalarInterpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    /** Start value. If a list of values is provided, each element will be interpolated separately. If unset, parent state value is used. */
    start: OptionalField(nullable(union(float, list(float))), null, 'Start value. If a list of values is provided, each element will be interpolated separately. If unset, parent state value is used.'),
    /** End value. If a list of values is provided, each element will be interpolated separately. If unset, only noise is applied. */
    end: OptionalField(nullable(union(float, list(float))), null, 'End value. If a list of values is provided, each element will be interpolated separately. If unset, only noise is applied.'),
    /** Whether to round the values to the closest integer. Useful for example for trajectory animation. */
    discrete: OptionalField(bool, false, 'Whether to round the values to the closest integer. Useful for example for trajectory animation.'),
    ..._Noise,
};

const Vec3Interpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    /** Start value. If unset, parent state value is used. Must be array of length 3N (x1, y1, z1, x2, y2, z2, ...). */
    start: OptionalField(nullable(list(float)), null, 'Start value. If unset, parent state value is used. Must be array of length 3N (x1, y1, z1, x2, y2, z2, ...).'),
    /** End value. Must be array of length 3N (x1, y1, z1, x2, y2, z2, ...). If unset, only noise is applied. */
    end: OptionalField(nullable(list(float)), null, 'End value. Must be array of length 3N (x1, y1, z1, x2, y2, z2, ...). If unset, only noise is applied.'),
    /** Whether to use spherical interpolation. */
    spherical: OptionalField(bool, false, 'Whether to use spherical interpolation.'),
    ..._Noise,
};

const RotationMatrixInterpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    /** Start value. If unset, parent state value is used. */
    start: OptionalField(nullable(Matrix), null, 'Start value. If unset, parent state value is used.'),
    /** End value. If unset, only noise is applied. */
    end: OptionalField(nullable(Matrix), null, 'End value. If unset, only noise is applied.'),
    ..._Noise,
};

const ColorInterpolation = {
    ..._Common,
    ..._Frequency,
    ..._Easing,
    /** Start value. If unset, parent state value is used. */
    start: OptionalField(union(nullable(ColorT), dict(union(int, str), ColorT)), null, 'Start value. If unset, parent state value is used.'),
    /** End value. */
    end: OptionalField(union(nullable(ColorT), dict(union(int, str), ColorT)), null, 'End value.'),
    /** Palette to sample colors from. Overrides start and end values. */
    palette: OptionalField(nullable(union(DiscretePalette, ContinuousPalette)), null, 'Palette to sample colors from. Overrides start and end values.'),
};

const TransformationMatrixInterpolation = {
    ..._Common,
    /** Pivot point for rotation and scale. */
    pivot: OptionalField(nullable(Vector3), null, 'Pivot point for rotation and scale.'),
    /** Start rotation value. If unset, parent state value is used. */
    rotation_start: OptionalField(nullable(Matrix), null, 'Start rotation value. If unset, parent state value is used.'),
    /** End rotation value. If unset, only noise is applied */
    rotation_end: OptionalField(nullable(Matrix), null, 'End rotation value. If unset, only noise is applied.'),
    /** Magnitude of the noise to apply to the rotation. */
    rotation_noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the rotation.'),
    /** Easing function to use for the rotation. */
    rotation_easing: OptionalField(Easing, 'linear', 'Easing function to use for the rotation.'),
    /** Determines how many times the rotation interpolation loops. Current T = frequency * t mod 1. */
    rotation_frequency: OptionalField(int, 1, 'Determines how many times the rotation interpolation loops. Current T = frequency * t mod 1.'),
    /** Whether to alternate the direction of the interpolation for frequency > 1. */
    rotation_alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
    /** Start translation value. If unset, parent state value is used. */
    translation_start: OptionalField(nullable(Vector3), null, 'Start translation value. If unset, parent state value is used.'),
    /** End translation value. If unset, only noise is applied. */
    translation_end: OptionalField(nullable(Vector3), null, 'End translation value. If unset, only noise is applied.'),
    /** Magnitude of the noise to apply to the translation. */
    translation_noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the translation.'),
    /** Easing function to use for the translation. */
    translation_easing: OptionalField(Easing, 'linear', 'Easing function to use for the translation.'),
    /** Determines how many times the translation interpolation loops. Current T = frequency * t mod 1. */
    translation_frequency: OptionalField(int, 1, 'Determines how many times the translation interpolation loops. Current T = frequency * t mod 1.'),
    /** Whether to alternate the direction of the interpolation for frequency > 1. */
    translation_alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
    /** Start scale value. If unset, parent state value is used. */
    scale_start: OptionalField(nullable(Vector3), null, 'Start scale value. If unset, parent state value is used.'),
    /** End scale value. If unset, only noise is applied. */
    scale_end: OptionalField(nullable(Vector3), null, 'End scale value. If unset, only noise is applied.'),
    /** Magnitude of the noise to apply to the scale. */
    scale_noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the scale.'),
    /** Easing function to use for the scale. */
    scale_easing: OptionalField(Easing, 'linear', 'Easing function to use for the scale.'),
    /** Determines how many times the scale interpolation loops. Current T = frequency * t mod 1. */
    scale_frequency: OptionalField(int, 1, 'Determines how many times the scale interpolation loops. Current T = frequency * t mod 1.'),
    /** Whether to alternate the direction of the interpolation for frequency > 1. */
    scale_alternate_direction: OptionalField(bool, false, 'Whether to alternate the direction of the interpolation for frequency > 1.'),
};

export const MVSAnimationSchema = TreeSchema({
    rootKind: 'animation',
    nodes: {
        animation: {
            description: 'Animation root node',
            parent: [],
            params: SimpleParamsSchema({
                /** Frame time in milliseconds. */
                frame_time_ms: OptionalField(float, 1000 / 60, 'Frame time in milliseconds.'),
                /** Total duration of the animation. If not specified, computed as maximum of all transitions. */
                duration_ms: OptionalField(nullable(float), null, 'Total duration of the animation. If not specified, computed as maximum of all transitions.'),
                /** Determines whether the animation should autoplay when a snapshot is loaded */
                autoplay: OptionalField(bool, true, 'Determines whether the animation should autoplay when a snapshot is loaded.'),
                /** Determines whether the animation should loop when it reaches the end. */
                loop: OptionalField(bool, false, 'Determines whether the animation should loop when it reaches the end.'),
                /** Determines whether the camera state should be included in the animation. */
                include_camera: OptionalField(bool, false, 'Determines whether the camera state should be included in the animation.'),
                /** Determines whether the canvas state should be included in the animation. */
                include_canvas: OptionalField(bool, false, 'Determines whether the canvas state should be included in the animation.'),
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