/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import type { ParamsFormModel } from '../forms/params-form.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export type SceneNodeKind =
  | 'structure'
  | 'component'
  | 'representation'
  | 'volume'
  | 'volume-representation'
  | 'measurements'
  | 'measurement';

/**
 * One row of the scene tree: a structure, a component, a representation, a volume, a measurement. Implemented over the
 * structure/volume hierarchy managers, keyed by state-tree ref, so a node keeps its identity across hierarchy rebuilds;
 * each field is its own state, so a row updates only where it changed.
 */
export interface SceneNode {
  /** The state-tree ref. */
  id: string;
  kind: SceneNodeKind;
  icon: IconName;
  label: ReadonlyState<string>;
  /** A short second line ("3,624 atoms", "Cartoon · Chain ID"). */
  description: ReadonlyState<string | undefined>;
  isHidden: ReadonlyState<boolean>;
  isExpanded: ReadonlyState<boolean>;
  children: ReadonlyState<readonly SceneNode[]>;
  toggleVisibility(): void;
  toggleExpanded(): void;
  /** Highlights it in 3D while hovered. */
  highlight(on: boolean): void;
  focus(): void;
  remove(): void;
}

/** What the inspector shows for a node: its parameters as a form (changes apply to the state tree). */
export interface SceneInspection {
  node: SceneNode;
  /** The node's own parameters (a representation's type and theme, a structure's model and assembly, …). */
  form: ParamsFormModel<any> | undefined;
  /** Section slot for contributed sections (`inspector:structure`, `inspector:component`, …). */
  slot: string;
}

export interface SceneModel {
  /** Top-level nodes: structures, volumes, measurements. */
  roots: ReadonlyState<readonly SceneNode[]>;
  /** The node picked for the inspector. */
  picked: ReadonlyState<SceneNode | undefined>;
  pick(node: SceneNode | undefined): void;
  inspect(node: SceneNode): SceneInspection;
}
