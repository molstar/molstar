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
import type { ObservableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import type { Structure } from '@molstar/model/model/structure';
import type { SequenceChain } from '../sequence/chain-filter.js';
import type { ChainFilter } from '../sequence/tracks.js';
import type { RepresentationQuickModel } from './representation-quick.js';
import type { SceneChanges } from './scene-changes.js';
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
  | 'extension'
  /** A structure's sequence (sequence.md §8.1): its eye is the band, its children its tracks. */
  | 'sequence'
  | 'sequence-track'
  /** A structure's components, grouped (what is drawn of it): open at first. */
  | 'components'
  /**
   * A structure's other rows, grouped (mostly off: the unit cell, volume streaming, the extensions' rows that ask for
   * it, `structureOfGroup`; the PAE, assembly symmetry, the membrane): open at first.
   */
  | 'other';

/** The id of a structure's Sequence node (`reveal` takes it). */
export const sequenceNodeId = (structure: string) => `${structure}#sequence`;
/** The ids of a structure's Components and Other groups. */
export const componentsNodeId = (structure: string) => `${structure}#components`;
export const otherNodeId = (structure: string) => `${structure}#other`;
/**
 * The structure (its state ref) a group (Components, Other) is of; undefined for any other node. An extension's row in
 * the Other group (a membrane, symmetry axes, the PAE) asks for `node.kind === 'other'` and reads its structure here.
 */
export function structureOfGroup(node: SceneNode): string | undefined {
  if (node.kind !== 'components' && node.kind !== 'other') return undefined;
  return node.id.slice(0, node.id.lastIndexOf('#'));
}

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
  /** False: it shows or hides nothing as one (Entities, Other), so its row has no eye. Absent: an eye. */
  readonly hasEye?: boolean;
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
  /** A structure built from a trajectory with several frames: its frames, shown and played (a row of its own). */
  readonly trajectory?: SceneTrajectory;
}

/**
 * A structure's frames (an NMR ensemble's models, a simulation's frames) as the Structure panel shows them: a row of
 * its own under the structure (the frame, a scrubber, ‹ ›, play) and its controls in the inspector. It acts through the
 * API (`trajectory.frame`, `trajectory.play` of this structure, `trajectory.stop`).
 */
export interface SceneTrajectory {
  /** What the frames are called: `models` (an ensemble) or `frames` (a simulation's coordinates). */
  readonly unit: 'models' | 'frames';
  readonly count: ReadonlyState<number>;
  /** The frame shown, from 0. */
  readonly index: ReadonlyState<number>;
  /** This structure's frames are playing. */
  readonly playing: ReadonlyState<boolean>;
  /** How it plays (`trajectory.play`'s mode and seconds per pass), kept while the structure is. */
  readonly mode: ObservableState<'loop' | 'palindrome' | 'once'>;
  readonly seconds: ObservableState<number>;
  /** Shows a frame (from 0): at most one in flight, the latest wins (a scrubber's drag). */
  seek(index: number): void;
  /** Moves by `by` frames, wrapping around. */
  step(by: number): void;
  togglePlay(): void;
  /** Its row is on screen (`true`) or no longer: the canvas' bar steps back while it is (ui/shell.tsx). */
  markShown(on: boolean): void;
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

/**
 * A row of the scene tree's node menus that an extension adds (`ui.nodeActions`: the clipping feature's Clip this…),
 * on the nodes of its kinds that pass `when`; after the node's own rows, before Remove.
 */
export interface SceneNodeAction {
  id: string;
  label: string;
  icon?: IconName;
  order?: number;
  kinds: readonly SceneNodeKind[];
  when?(node: SceneNode): boolean;
  run(node: SceneNode): void;
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
  /** A Sequence node's controls: the band's style and letters, Add track, the chains (sequence.md §8.1). */
  sequence?: SequenceInspectorModel;
}

/**
 * A Sequence node's inspector (sequence.md §8.1), made once per node and following the edits: what the band shows, its
 * style and letter colors (the band's: the same on every node), Add track, and the chains the structure's sequence
 * shows (the presets, each chain by its entity, every copy apart).
 */
export interface SequenceInspectorModel {
  /** What the band shows of it, in a few words ("Tracks · auto · 4 tracks"; "Not shown"). */
  readonly status: ReadonlyState<string>;
  readonly style: ReadonlyState<'tracks' | 'text'>;
  setStyle(style: 'tracks' | 'text'): void;
  readonly colors: ReadonlyState<string | null>;
  /** The letter color schemes: a name for people and a swatch of its colors (none: not painted). */
  readonly schemes: ReadonlyState<
    readonly { readonly id: string | null; readonly label: string; readonly swatch?: string }[]
  >;
  setColors(scheme: string | null): void;
  /** A track can be added here (the structure the band shows, a provider out): Add track opens the command menu. */
  readonly canAdd: ReadonlyState<boolean>;
  addTrack(): void;
  /** The structure's entities and their chains (new when the structure under the node is another). */
  readonly entities: ReadonlyState<readonly SequenceChainsEntity[]>;
  /** The chains shown, by name. */
  readonly shown: ReadonlyState<ReadonlySet<string>>;
  /** The preset in effect; undefined: chains chosen one by one. */
  readonly preset: ReadonlyState<SequenceChainsPreset | undefined>;
  apply(preset: SequenceChainsPreset): void;
  /** Shows or leaves out one chain. */
  toggleChain(name: string): void;
  /** Shows every chain of an entity, or none when every one is shown. */
  toggleEntity(id: string): void;
  /** Every copy a block of its own (an assembly's). */
  readonly each: ReadonlyState<boolean>;
  setEach(on: boolean): void;
}

export type SequenceChainsPreset = 'all' | 'none' | 'polymer';

export interface SequenceChainsEntity {
  readonly id: string;
  /** Its description, one line ("Cellular retinoic acid binding protein type II"). */
  readonly name: string;
  readonly kind: 'polymer' | 'branched' | 'ligand';
  /** Its chain instances: `name` as the sequence names the unit, `label` with the author's id where it differs. */
  readonly chains: readonly { readonly name: string; readonly label: string }[];
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
  /**
   * The scene's structure as one change stream (scene-changes.ts): the nodes whose children or expansion changed, once
   * per Mol* transaction or toggle. A view of the whole tree listens to this, not to each node. Absent (a fixture): the
   * tree follows each shown node's `children` and `isExpanded`.
   */
  changes?: ReadonlyState<SceneChanges>;
  /**
   * The tree filtered by the selection (the search's non-text part, a chip in its field): only the rows that hold it.
   * Absent (a fixture without it): no such filter.
   */
  selectionFilter?: SceneSelectionFilter;
  /** The rows the extensions add to the node menus (`ui.nodeActions`); absent: none can be added. */
  nodeActions?: ObservableState<readonly SceneNodeAction[]>;
}

/** Which rows hold the current selection: structures, their components and representations, Contents' entities and chains. */
export interface SceneSelectionFilter {
  readonly on: ReadonlyState<boolean>;
  /** On; off when the selection is cleared (nothing would show), and it cannot go on with nothing selected. */
  set(on: boolean): void;
  /** Told when the selection changes while on: the tree searches again. */
  readonly version: ReadonlyState<number>;
  holds(node: SceneNode): boolean;
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

/**
 * Each structure's sequence as the structure tree shows it (sequence.md §8.1): a Sequence node under every structure
 * with polymers, its eye the band, its children the tracks. Made over the sequence model (`sources/sequence-scene.ts`):
 * no Seq* here, and the tree builds no Mol* cell for it.
 */
export interface SceneSequence {
  /** Told when any structure's row may read otherwise (what the band follows, a track, a filter). */
  readonly version: ReadonlyState<number>;
  /** A structure's row; undefined: it has no polymer (no Sequence node). */
  row(structure: string, data: Structure): SceneSequenceRow | undefined;
  /** The eye: the band shows this structure, or no sequence when it does. */
  toggle(structure: string): void;
  /** Follows the rule again (Auto): `sequence.follow` with none. */
  followAutomatically(): void;
  setTrackHidden(structure: string, data: Structure, track: string, hidden: boolean): void;
  addTrack(structure: string, data: Structure, provider: string): void;
  removeTrack(structure: string, data: Structure, track: string): void;
  setChains(structure: string, data: Structure, filter: ChainFilter | undefined): void;
  /** The colors' legend (the band's). */
  legend?(): void;
  /** The band's style and letter colors (one for the band, sequence.md §7.1): `sequence.style`, `sequence.colors`. */
  setStyle(style: 'tracks' | 'text'): void;
  setColors(scheme: string | null): void;
  /** The letter color schemes, by name (null: no colors), their names for people and swatches (once the sequence is in). */
  schemes(): readonly { readonly id: string | null; readonly label: string; readonly swatch?: string }[];
  /** Opens the command menu at Add track, for the structure the band shows (its commands act on it). */
  openAddTrack?(): void;
}

export interface SceneSequenceRow {
  /** The band shows this structure. */
  readonly shown: boolean;
  /** The person chose it (`sequence.follow`); false: the rule's pick (Auto) or not shown. */
  readonly chosen: boolean;
  /** "3 tracks · auto · chains A, B", "Not shown · 3 tracks". */
  readonly description: string;
  readonly tracks: readonly SceneSequenceTrack[];
  /** What Add track offers. */
  readonly addable: readonly { readonly id: string; readonly label: string }[];
  /** The chain instances it can show, and the filter on them. */
  readonly chains: readonly SequenceChain[];
  readonly filter: ChainFilter | undefined;
  /** The band's style and letter colors (the same on every row: they are the band's). */
  readonly style: 'tracks' | 'text';
  readonly colors: string | null;
}

export interface SceneSequenceTrack {
  readonly ref: string;
  readonly label: string;
  readonly hidden: boolean;
  readonly removable: boolean;
  /** What it draws, in a few words ("Letters · Zappo"). */
  readonly description?: string;
  /** Its colorings of the 3D view (the structure shown only). */
  readonly actions: readonly { readonly id: string; readonly label: string; readonly icon?: IconName; run(): void }[];
}
