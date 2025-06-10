/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  MVSTree,
  MVSTreeSchema,
  MVSKind,
  MVSNode,
  MVSNodeParams,
} from "molviewspec/tree/mvs/mvs-tree";

/** MolstarTree extends MVSTree with additional Molstar-specific nodes */
export type MolstarTree = MVSTree;

/** Schema for MolstarTree (currently same as MVS) */
export const MolstarTreeSchema = MVSTreeSchema;

/** Union of all possible node kinds in MolstarTree */
export type MolstarKind = MVSKind;

/** Union of all possible nodes in MolstarTree */
export type MolstarNode<K extends MolstarKind = MolstarKind> = MVSNode<K>;

/** Union of all possible subtrees in MolstarTree */
export type MolstarSubtree<K extends MolstarKind = MolstarKind> = MVSNode<K>;

/** Node parameters for a specific node kind */
export type MolstarNodeParams<K extends MolstarKind> = MVSNodeParams<K>;
