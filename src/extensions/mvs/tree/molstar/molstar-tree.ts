/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { omitObjectKeys, pickObjectKeys } from '../../../../mol-util/object';
import { RequiredField, bool, str } from '../generic/field-schema';
import { SimpleParamsSchema } from '../generic/params-schema';
import { NodeFor, ParamsOfKind, SubtreeOfKind, TreeFor, TreeSchema } from '../generic/tree-schema';
import { FullMVSTreeSchema } from '../mvs/mvs-tree';
import { MolstarParseFormatT } from '../mvs/param-types';


/** Schema for `MolstarTree` (intermediate tree representation between `MVSTree` and a real Molstar state) */
export const MolstarTreeSchema = TreeSchema({
    rootKind: 'root',
    nodes: {
        ...FullMVSTreeSchema.nodes,
        download: {
            ...FullMVSTreeSchema.nodes.download,
            params: SimpleParamsSchema({
                ...FullMVSTreeSchema.nodes.download.params.fields,
                is_binary: RequiredField(bool, 'Specifies whether file is downloaded as bytes array or string'),
            }),
        },
        parse: {
            ...FullMVSTreeSchema.nodes.parse,
            params: SimpleParamsSchema({
                format: RequiredField(MolstarParseFormatT, 'File format'),
            }),
        },
        /** Auxiliary node corresponding to Molstar's CoordinatesFrom*. */
        coordinates: {
            description: "Auxiliary node corresponding to Molstar's CoordinatesFrom*.",
            parent: ['parse'],
            params: SimpleParamsSchema({
                format: RequiredField(MolstarParseFormatT, 'File format'),
            }),
        },
        /** Auxiliary node corresponding to Molstar's TrajectoryFrom*. */
        trajectory: {
            description: "Auxiliary node corresponding to Molstar's TrajectoryFrom*.",
            parent: ['parse'],
            params: SimpleParamsSchema({
                format: RequiredField(MolstarParseFormatT, 'File format'),
                ...pickObjectKeys(FullMVSTreeSchema.nodes.structure.params.fields, ['block_header', 'block_index'] as const),
            }),
        },
        /** Auxiliary node corresponding to Molstar's TrajectoryFrom*. */
        trajectory_with_coordinates: {
            description: "Auxiliary node corresponding to assigning a separate coordinates to a trajectory.",
            parent: ['model'],
            params: SimpleParamsSchema({
                coordinates_ref: RequiredField(str, 'Coordinates reference'),
            }),
        },
        /** Auxiliary node corresponding to Molstar's ModelFromTrajectory. */
        model: {
            description: "Auxiliary node corresponding to Molstar's ModelFromTrajectory.",
            parent: ['trajectory', 'trajectory_with_coordinates'],
            params: SimpleParamsSchema(
                pickObjectKeys(FullMVSTreeSchema.nodes.structure.params.fields, ['model_index'] as const)
            ),
        },
        /** Auxiliary node corresponding to Molstar's StructureFromModel. */
        structure: {
            ...FullMVSTreeSchema.nodes.structure,
            parent: ['model'],
            params: SimpleParamsSchema(
                omitObjectKeys(FullMVSTreeSchema.nodes.structure.params.fields, ['block_header', 'block_index', 'model_index', 'coordinates_ref'] as const)
            ),
        },
    }
});


/** Node kind in a `MolstarTree` */
export type MolstarKind = keyof typeof MolstarTreeSchema.nodes;

/** Node in a `MolstarTree` */
export type MolstarNode<TKind extends MolstarKind = MolstarKind> = NodeFor<typeof MolstarTreeSchema, TKind>

/** Params for a specific node kind in a `MolstarTree` */
export type MolstarNodeParams<TKind extends MolstarKind> = ParamsOfKind<MolstarTree, TKind>

/** Intermediate tree representation between `MVSTree` and a real Molstar state */
export type MolstarTree = TreeFor<typeof MolstarTreeSchema>

/** Any subtree in a `MolstarTree` (e.g. its root doesn't need to be 'root') */
export type MolstarSubtree<TKind extends MolstarKind = MolstarKind> = SubtreeOfKind<MolstarTree, TKind>
