/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
    createTreeSchema,
    TreeSchema,
} from 'molviewspec/tree/generic/tree-schema';
import { MVSTree, MVSTreeSchema } from 'molviewspec/tree/mvs/mvs-tree';

/** MolstarTree extends MVSTree with additional Molstar-specific nodes */
export type MolstarTree = MVSTree;

/** Schema for MolstarTree (currently same as MVS) */
export const MolstarTreeSchema: TreeSchema<MolstarTree> = MVSTreeSchema;

/** Union of all possible node kinds in MolstarTree */
export type MolstarKind = MVSTree['kind'];

/** Union of all possible nodes in MolstarTree */
export type MolstarNode<K extends MolstarKind = MolstarKind> = Extract<
MVSTree,
{ kind: K }
>;

/** Union of all possible subtrees in MolstarTree */
export type MolstarSubtree<K extends MolstarKind = MolstarKind> =
  MolstarNode<K>;

/** Node parameters for a specific node kind */
export type MolstarNodeParams<K extends MolstarKind> =
  MolstarNode<K> extends { params: infer P } ? P : never;
