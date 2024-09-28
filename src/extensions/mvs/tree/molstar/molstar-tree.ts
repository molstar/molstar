/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { omitObjectKeys, pickObjectKeys } from '../../../../mol-util/object';
import { OptionalField, RequiredField, bool } from '../generic/params-schema';
import { NodeFor, TreeFor, TreeSchema } from '../generic/tree-schema';
import { FullMVSTreeSchema } from '../mvs/mvs-tree';
import { MolstarParseFormatT, RawVolumeOptionsT, RawVolumeSourceT } from '../mvs/param-types';


/** Schema for `MolstarTree` (intermediate tree representation between `MVSTree` and a real Molstar state) */
export const MolstarTreeSchema = TreeSchema({
    rootKind: 'root',
    nodes: {
        ...FullMVSTreeSchema.nodes,
        download: {
            ...FullMVSTreeSchema.nodes.download,
            params: {
                ...FullMVSTreeSchema.nodes.download.params,
                is_binary: RequiredField(bool),
            },
        },
        parse: {
            ...FullMVSTreeSchema.nodes.parse,
            params: {
                format: RequiredField(MolstarParseFormatT),
            },
        },
        raw_volume: {
            ...FullMVSTreeSchema.nodes.raw_volume,
            params: {
                source: RequiredField(RawVolumeSourceT),
                options: OptionalField(RawVolumeOptionsT)
            }
        },
        /** Auxiliary node corresponding to Molstar's TrajectoryFrom*. */
        trajectory: {
            description: "Auxiliary node corresponding to Molstar's TrajectoryFrom*.",
            parent: ['parse'],
            params: {
                format: RequiredField(MolstarParseFormatT),
                ...pickObjectKeys(FullMVSTreeSchema.nodes.structure.params, ['block_header', 'block_index'] as const),
            },
        },
        /** Auxiliary node corresponding to Molstar's ModelFromTrajectory. */
        model: {
            description: "Auxiliary node corresponding to Molstar's ModelFromTrajectory.",
            parent: ['trajectory'],
            params: pickObjectKeys(FullMVSTreeSchema.nodes.structure.params, ['model_index'] as const),
        },
        /** Auxiliary node corresponding to Molstar's StructureFromModel. */
        structure: {
            ...FullMVSTreeSchema.nodes.structure,
            parent: ['model'],
            params: omitObjectKeys(FullMVSTreeSchema.nodes.structure.params, ['block_header', 'block_index', 'model_index'] as const),
        },
    }
});


/** Node kind in a `MolstarTree` */
export type MolstarKind = keyof typeof MolstarTreeSchema.nodes;

/** Node in a `MolstarTree` */
export type MolstarNode<TKind extends MolstarKind = MolstarKind> = NodeFor<typeof MolstarTreeSchema, TKind>

/** Intermediate tree representation between `MVSTree` and a real Molstar state */
export type MolstarTree = TreeFor<typeof MolstarTreeSchema>
