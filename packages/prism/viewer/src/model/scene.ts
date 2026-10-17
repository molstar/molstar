/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { Observable } from 'rxjs';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import type { RepresentationQuickModel } from './representation-quick.js';
import type { VolumeLevelsModel } from './volume-levels.js';
import type { VolumeTransferModel } from './volume-transfer.js';

export type SceneNodeKind =
  | 'structure'
  | 'component'
  | 'representation'
  | 'volume'
  | 'volume-representation'
  | 'particles'
  | 'particles-representation'
  | 'measurements'
  | 'measurement'
  | 'unit-cell'
  | 'contents'
  | 'entity'
  | 'chain'
  /** A structure's focus group: the focus behavior's components (viewer-spec.md §2.1, always there). */
  | 'focus'
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
  /** A color shown before the description (a volume's isosurface: which one is which); undefined: a theme's colors. */
  swatch?: ReadonlyState<number | undefined>;
  isHidden: ReadonlyState<boolean>;
  isExpanded: ReadonlyState<boolean>;
  children: ReadonlyState<readonly SceneNode[]>;
  /**
   * Its children are computed when first expanded (a structure's contents: entities, then chains), so a mesoscale
   * structure with thousands of chains costs nothing until someone looks. Expandable before it has any.
   */
  lazy?: boolean;
  /**
   * Computes the children (`children` updates when they arrive); from then on they follow the scene. With `peek` (a
   * search looking in), they are computed while it stays lazy: they do not follow the scene, and are computed again at
   * the next look after it changed.
   */
  load?(peek?: boolean): void;
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
  /**
   * Whether it has a menu, without building it (a menu read when it opens may be built then, from the state: a row
   * asks this when it is made, and a virtualized tree makes rows as it scrolls). Absent: `menu` is read.
   */
  readonly hasMenu?: boolean;
}

export interface SceneMenuItem {
  id: string;
  label: string;
  icon?: IconName;
  /** Destructive (Remove): set apart at the end, in the danger color. */
  danger?: boolean;
  run?(): void;
  items?: readonly SceneMenuItem[];
  /** Its items are icons on its own row, always shown (Modify by selection: include, subtract, intersect), not a fold. */
  inline?: true;
  /** Marked as the current choice (the representation's type in its Type menu). */
  checked?: boolean;
  /** Unavailable while true (an action on the selection with nothing selected). */
  disabled?: Observable<boolean>;
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
  /** Gives it another name (a component: `component.rename`); the inspector's title becomes editable. */
  rename?(label: string): void;
  /** A map's levels, above the form (a volume, an isosurface, streamed density: one per map). */
  levels?: readonly VolumeLevelsModel[];
  /** A direct volume's transfer function (opacity and color over the map's values), above the form. */
  transfer?: readonly VolumeTransferModel[];
  /** A structure representation's quick settings (type, color, opacity, size), above the form. */
  quick?: RepresentationQuickModel;
  /** A structure's contents: its counts, and its entities to open in the tree. */
  contents?: ContentsSummary;
}

/** What a structure is made of, at a glance. */
export interface ContentsSummary {
  /** Entities, chains, polymer residues, atoms (models when several), each a figure. */
  readonly counts: readonly { readonly label: string; readonly value: string }[];
  readonly entities: readonly ContentsEntity[];
  /** Opens an entity's row in the tree (Contents opened) and picks it. */
  reveal(entity: string): void;
}

export interface ContentsEntity {
  readonly id: string;
  readonly name: string;
  /** mmCIF's: polymer, non-polymer, branched, water, … */
  readonly type: string;
  /** Its chains (author ids). */
  readonly chains: readonly string[];
  readonly atoms: number;
}

export interface SceneModel {
  /** Top-level nodes: structures, volumes, measurements. */
  roots: ReadonlyState<readonly SceneNode[]>;
  /** The node picked for the inspector. */
  picked: ReadonlyState<SceneNode | undefined>;
  pick(node: SceneNode | undefined): void;
  inspect(node: SceneNode): SceneInspection;
  /**
   * The picked node's inspection, built anew when what it shows is no longer what it was built from (a structure's
   * facts and contents after an assembly change, a volume's levels and transfer functions as its representations
   * change); its form stays the same. Absent: `inspect(picked)` holds.
   */
  inspection?: ReadonlyState<SceneInspection | undefined>;
  /** Picks a node by id, its ancestors expanded so the tree shows where it is; false while it is not in the tree. */
  reveal?(id: string): boolean;
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
