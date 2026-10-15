/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

export type SceneNodeKind =
  | 'structure'
  | 'component'
  | 'representation'
  | 'volume'
  | 'volume-representation'
  | 'unit-cell'
  | 'contents'
  | 'entity'
  | 'chain'
  | 'measurements'
  | 'measurement'
  /** A row an extension adds under a structure (volume streaming, assembly symmetry). */
  | 'extension';

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
  /**
   * Its children are computed when first expanded (a structure's contents: entities, then chains), so a mesoscale
   * structure with thousands of chains costs nothing until someone looks. Expandable before it has any.
   */
  lazy?: boolean;
  /** Computes the children (once; `children` updates when they arrive). */
  load?(): void;
  toggleVisibility(): void;
  toggleExpanded(): void;
  /** Highlights it in 3D while hovered. */
  highlight(on: boolean): void;
  focus(): void;
  remove(): void;
  /** Can it be removed (a structure's contents, entities and chains, cannot)? True when absent. */
  removable?: boolean;
  /**
   * The row's menu (the "…"): actions on this node, as Mol*'s component menu has them: add representation, set
   * coloring, modify by selection (include, subtract, intersect), select this, edit label. Entries with `items` open
   * inline.
   */
  menu?: readonly SceneMenuItem[];
}

export interface SceneMenuItem {
  id: string;
  label: string;
  icon?: IconName;
  /** Destructive (Remove): set apart at the end, in the danger color. */
  danger?: boolean;
  run?(): void;
  items?: readonly SceneMenuItem[];
}

/** What the inspector shows for a node: its parameters as a form (changes apply to the state tree). */
export interface SceneInspection {
  node: SceneNode;
  /** The node's own parameters (a representation's type and theme, a structure's model and assembly, …). */
  form: ParamsFormModel<any> | undefined;
  /** Section slot for contributed sections (`inspector:structure`, `inspector:component`, …). */
  slot: string;
  /** Read-only facts shown compactly above the form (atoms, chains; a map's grid, voxel size, σ). */
  facts?: readonly (readonly [label: string, value: string])[];
  /** Things to do from here, as plain buttons under the facts (Measurements: measure the selection, in its tab). */
  actions?: readonly SceneMenuItem[];
}

export interface SceneModel {
  /** Top-level nodes: structures, volumes, measurements. */
  roots: ReadonlyState<readonly SceneNode[]>;
  /** The node picked for the inspector. */
  picked: ReadonlyState<SceneNode | undefined>;
  pick(node: SceneNode | undefined): void;
  inspect(node: SceneNode): SceneInspection;
}

/** A representation preset (Mol*'s structure representation preset registry): its name, group, description. */
export interface StructurePreset {
  id: string;
  label: string;
  description?: string;
}

/**
 * The representation presets for the Preset menu, grouped as the registry groups them (ungrouped first: Empty,
 * Automatic; then Basic, Miscellaneous, and what extensions register), and the one in effect.
 */
export interface PresetsModel {
  groups: readonly { label?: string; presets: readonly StructurePreset[] }[];
  current: ReadonlyState<string | undefined>;
  apply(id: string): void;
}

/**
 * Quick styles (Mol*'s Quick Styles panel): one click restyles everything, by representation preset or by visual
 * style. Commands underneath (`structure.preset`, `canvas.style`), so the command menu offers them too.
 */
export interface QuickStylesModel {
  representations: readonly { id: string; label: string }[];
  styles: readonly { id: string; label: string }[];
  applyRepresentation(id: string): void;
  applyStyle(id: string): void;
}
