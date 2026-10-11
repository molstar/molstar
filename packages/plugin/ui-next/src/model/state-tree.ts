/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import type { HomeAction } from './home.js';

/** A cell of the state tree (`plugin.state.data`): what the transform made, and which transform. */
export interface StateTreeNode {
  ref: string;
  label: string;
  /** The object's description ("3,624 atoms", "1 model"). */
  description?: string;
  icon: IconName;
  /** The transformer's display name ("Download", "Parse CIF", "Structure from Model"). */
  transform: string;
  /** Decorators folded into the cell (A → A transforms: "Custom Structure Properties"), by display name. */
  decorators?: readonly string[];
  children?: readonly StateTreeNode[];
  /** Undefined when the object has no visual. */
  visible?: boolean;
  status?: 'ok' | 'pending' | 'error';
  error?: string;
}

/**
 * The selected cell: its updates (Mol*'s "Update …": the cell's own transform first, then the decorators folded into it,
 * the A → A transforms like custom structure properties that Mol* shows as part of the cell, not as cells of their
 * own), and the state actions applicable to it (Mol*'s "Apply Action"), each a form with an apply button, like Home's.
 */
export interface StateTreeInspection {
  node: StateTreeNode;
  updates: readonly HomeAction[];
  actions: readonly HomeAction[];
}

/** The state tree for advanced use (implemented over `plugin.state.data`). */
export interface StateTreeModel {
  roots: ReadonlyState<readonly StateTreeNode[]>;
  selected: ReadonlyState<string | undefined>;
  select(ref: string | undefined): void;
  inspect(ref: string): StateTreeInspection | undefined;
  toggleVisible(ref: string): void;
  remove(ref: string): void;
}
