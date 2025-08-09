/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, int, list, OptionalField, RequiredField, str, union, nullable, literal, ValueFor } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { NodeFor, ParamsOfKind, SubtreeOfKind, TreeFor, TreeSchema } from '../generic/tree-schema';

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

const ScalarTransition = {
    target_ref: RequiredField(str, 'Reference to the node.'),
    property: RequiredField(union(str, list(union(str, int))), 'Scalar value accessor.'),
    start_value: OptionalField(nullable(float), null, 'Start value. If unset, source value is used.'),
    end_value: RequiredField(float, 'End value for the transition.'),
    start_ms: OptionalField(float, 0, 'Start time of the transition in milliseconds.'),
    end_ms: RequiredField(float, 'End time of the transition in milliseconds.'),
    easing: OptionalField(Easing, 'linear', 'Easing function to use for the transition.'),
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
        transition: {
            description: 'This node enables transitions between different states',
            parent: ['animation'],
            params: UnionParamsSchema(
                'type',
                'Transition type',
                {
                    scalar: SimpleParamsSchema(ScalarTransition),
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