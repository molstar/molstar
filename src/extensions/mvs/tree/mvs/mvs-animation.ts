/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, int, list, OptionalField, RequiredField, str, union, nullable, literal, ValueFor } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { NodeFor, ParamsOfKind, SubtreeOfKind, TreeFor, TreeSchema } from '../generic/tree-schema';
import { Matrix, Vector3 } from './param-types';

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

const _Common = {
    target_ref: RequiredField(str, 'Reference to the node.'),
    property: RequiredField(union(str, list(union(str, int))), 'Value accessor.'),
    start_ms: OptionalField(float, 0, 'Start time of the transition in milliseconds.'),
    duration_ms: RequiredField(float, 'End time of the transition in milliseconds.'),
    easing: OptionalField(Easing, 'linear', 'Easing function to use for the transition.'),
    // TODO: support cummulative noise?
    noise_magnitude: OptionalField(float, 0, 'Magnitude of the noise to apply to the transition.'),
};

const ScalarTransition = {
    ..._Common,
    from: OptionalField(nullable(float), null, 'Start value. If unset, source value is used.'),
    to: RequiredField(float, 'End value for the transition.'),
};

const Vec3Transition = {
    ..._Common,
    from: OptionalField(nullable(Vector3), null, 'Start value. If unset, source value is used.'),
    to: RequiredField(Vector3, 'End value for the transition.'),
    spherical: OptionalField(bool, false, 'Whether to use spherical interpolation.'),
};

const RotationMatrixTransition = {
    ..._Common,
    from: OptionalField(nullable(Matrix), null, 'Start value. If unset, source value is used.'),
    to: RequiredField(Matrix, 'End value for the transition.'),
};

export const MVSAnimationSchema = TreeSchema({
    rootKind: 'animation',
    nodes: {
        animation: {
            description: 'Animation root node',
            parent: [],
            params: SimpleParamsSchema({
                frame_time_ms: OptionalField(float, 1000 / 60, 'Frame time in milliseconds'),
                autoplay: OptionalField(bool, true, 'Determines whether the animation should autoplay when a snapshot is loaded'),
                loop: OptionalField(bool, false, 'Determines whether the animation should loop when it reaches the end'),
                include_camera: OptionalField(bool, false, 'Determines whether the camera state should be included in the animation'),
                include_canvas: OptionalField(bool, false, 'Determines whether the canvas state should be included in the animation'),
            }),
        },
        interpolate: {
            description: 'This node enables transitions between different states',
            parent: ['animation'],
            params: UnionParamsSchema(
                'type',
                'Transition type',
                {
                    scalar: SimpleParamsSchema(ScalarTransition),
                    vec3: SimpleParamsSchema(Vec3Transition),
                    rotation_matrix: SimpleParamsSchema(RotationMatrixTransition),
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