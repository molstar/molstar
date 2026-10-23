/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type EmptiedComponent, emptiedComponentsOf, hiddenContents } from '../api/contents.js';
import { StructureFocusRepresentationTags } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { InitVolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { SpacegroupCell } from '@molstar/core/math/geometry';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { type State, type StateObjectCell, StateSelection, StateTree } from '@molstar/core/state';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import type { Model } from '@molstar/model/model/structure';
import { Volume } from '@molstar/model/model/volume';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import {
  type Structure,
  StructureElement,
  StructureProperties,
  StructureSelection,
  type Unit,
} from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import type {
  ModelUnitcellRef,
  StructureComponentRef,
  StructureRef,
  StructureRepresentationRef,
} from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type {
  ParticleListRef,
  ParticleRepresentationRef,
} from '@molstar/plugin/state/manager/particles/hierarchy-state';
import type { VolumeRef, VolumeRepresentationRef } from '@molstar/plugin/state/manager/volume/hierarchy-state';
import { type IconName, type ObservableState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerApi } from '../api/area.js';
import type {
  ContentsEntity,
  ContentsSummary,
  SceneInspection,
  SceneMenuItem,
  SceneModel,
  SceneNode,
  SceneNodeAction,
  SceneNodeKind,
  SceneSelectionFilter,
  SceneTrajectory,
} from '../model/scene.js';
import { SceneChangeLog, sameChildren } from '../model/scene-changes.js';
import type { VolumeLevelUnit } from '../model/volume-levels.js';
import { namesOf } from '../api/scope.js';
import { currentTransform } from '../api/structure.js';
import { VolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/behavior';
import {
  CreateVolumeStreamingBehavior,
  CreateVolumeStreamingInfo,
} from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import {
  gridText,
  type StreamingState,
  shortNumber,
  streamingOf,
  streamingSummary,
  voxelText,
  volumeStats,
} from '../api/volume.js';
import { type MeasuredCell, MeasureKinds, measurementLocis, readMeasurements } from './measurements.js';
import { type TrajectorySource, trajectorySourceOf } from './trajectory.js';
import { derive } from '@molstar/prism-ui/reactive/state';
import { structureNames } from './naming.js';
import { fromChanged, onceATask, type SourceScope } from './state.js';
import type { Observable } from 'rxjs';
import { createVolumeLevels, mapLevels, streamingLevels, type VolumeLevelsSource } from './volume-levels.js';
import { createVolumeTransfer, type VolumeTransferSource } from './volume-transfer.js';
import {
  createRepresentationQuick,
  representationChoices,
  type RepresentationQuickSource,
} from './representation-quick.js';
import type { RepresentationQuickChange } from '../model/representation-quick.js';
import { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';
import type { BondProviderRegistry } from '@molstar/plugin/state/bonds/registry';

// The scene as the Structure tab shows it (state.md §1, row 2; view-models.md §3): one record per node, keyed by state
// ref, so a node keeps its identity across hierarchy rebuilds, and each field its own state, so a row updates only where
// it changed. A structure's contents (entities, then chains) are computed when first expanded. Every action goes through
// the API, bound to the person (api.md §11): the tree changes state only by operations.

export const ReprIcons: Readonly<Record<string, IconName>> = {
  cartoon: 'cartoon',
  'ball-and-stick': 'ball-and-stick',
  spacefill: 'spacefill',
  'molecular-surface': 'surface',
  'gaussian-surface': 'surface',
  'gaussian-volume': 'surface',
  putty: 'cartoon',
  backbone: 'cartoon',
  line: 'ball-and-stick',
  label: 'tag',
  carbohydrate: 'shapes',
  ellipsoid: 'spacefill',
};

/** The operation that applies an inspector form, by the kind of representation it edits. */
const UpdateOperations: Partial<Record<SceneNodeKind, string>> = {
  representation: 'representation.update',
  'volume-representation': 'volume.update',
  'particles-representation': 'particles.update',
};

/** What a node id is kept by: the part before its `#` (a structure's contents, entities and chains), else itself. */
const ownerOf = (id: string): string => {
  const i = id.indexOf('#');
  return i < 0 ? id : id.slice(0, i);
};

/**
 * The source's nodes by ref: a node put in is given the source's change log. It knows the ids put in since the last
 * sync's clean-up (`fresh`) and the ids each owner keeps (`under`), so the clean-up looks only at what may have gone.
 */
class NodeMap extends Map<string, Node> {
  /** Ids put in since `fresh` was last cleared. */
  readonly fresh = new Set<string>();
  readonly #byOwner = new Map<string, Set<string>>();
  constructor(
    readonly log: SceneChangeLog,
    /** The rows the extensions add to the node menus, each node reads them when its menu opens */
    readonly actions?: ReadonlyState<readonly SceneNodeAction[]>,
  ) {
    super();
  }
  override set(key: string, node: Node): this {
    node.log = this.log;
    node.added = this.actions;
    if (!super.has(key)) {
      this.fresh.add(key);
      const owner = ownerOf(key);
      if (owner !== key) {
        let ids = this.#byOwner.get(owner);
        if (!ids) this.#byOwner.set(owner, (ids = new Set()));
        ids.add(key);
      }
    }
    return super.set(key, node);
  }
  override delete(key: string): boolean {
    if (!super.delete(key)) return false;
    this.fresh.delete(key);
    const owner = ownerOf(key);
    if (owner !== key) {
      const ids = this.#byOwner.get(owner);
      ids?.delete(key);
      if (ids?.size === 0) this.#byOwner.delete(owner);
    }
    return true;
  }
  /** The ids kept by `owner` (`owner#…`). */
  under(owner: string): ReadonlySet<string> | undefined {
    return this.#byOwner.get(owner);
  }
}

/**
 * What a structure's rows were built from (its hierarchy refs' cells, transforms and objects, its name, its hidden
 * contents), compared in place at each sync: the same values in the same order are "as before", and the rows are kept.
 * No allocation while it is the same.
 */
class Signature {
  readonly #values: unknown[] = [];
  #at = 0;
  #same = true;
  begin(): void {
    this.#at = 0;
    this.#same = true;
  }
  add(value: unknown): void {
    const at = this.#at++;
    if (this.#same && at < this.#values.length && Object.is(this.#values[at], value)) return;
    this.#same = false;
    this.#values[at] = value;
  }
  /** Whether every value was as before, and as many. */
  end(): boolean {
    if (this.#at !== this.#values.length) {
      this.#same = false;
      this.#values.length = this.#at;
    }
    return this.#same;
  }
}

/** A structure's rows as last built (`syncRows`): its node, the ids they hold, its focus parts. */
interface StructureRows {
  readonly signature: Signature;
  node?: Node;
  /** The ids of its rows (what the clean-up keeps). */
  ids: readonly string[];
  /** Its focus components' rows (under the focus row) and the focus parts with atoms. */
  focused: readonly SceneNode[];
  parts: readonly StructureComponentRef[];
  /** Its volume streaming's root and behavior cells, while it streams. */
  streaming?: { readonly root: string; readonly behavior: string };
}

const NoEmptied: readonly EmptiedComponent[] = [];

class Node implements SceneNode {
  /** Where its structural changes are told (the source's log; set when it is put in the source's nodes). */
  log?: SceneChangeLog;
  readonly #label: ObservableState<string>;
  readonly #description = observableState<string | undefined>(undefined);
  readonly #hidden = observableState(false);
  readonly #expanded: ObservableState<boolean>;
  readonly #children = observableState<readonly SceneNode[]>([]);
  #swatch?: ObservableState<number | undefined>;
  /** A structure's frames (several): `structureNode` sets it, and takes it away when it has one model again. */
  trajectory?: SceneTrajectory;
  icon: IconName;
  #menu?: readonly SceneMenuItem[];
  /** Builds the menu when read, for one that follows the state (a representation's current type marked). */
  menuSource?: () => readonly SceneMenuItem[] | undefined;
  /** The extensions' rows for every node (`ui.nodeActions`): this node's are those of its kind that take it */
  added?: ReadonlyState<readonly SceneNodeAction[]>;
  removable?: boolean;
  lazy?: boolean;
  load?: (peek?: boolean) => void;

  constructor(
    readonly id: string,
    readonly kind: SceneNodeKind,
    icon: IconName,
    label: string,
    expanded: boolean,
    readonly actions: {
      visibility?(hidden: boolean): void;
      highlight?(on: boolean): void;
      focus?(): void;
      remove?(): void;
    },
  ) {
    this.icon = icon;
    this.#label = observableState(label);
    this.#expanded = observableState(expanded);
  }

  get menu(): readonly SceneMenuItem[] | undefined {
    const own = this.menuSource?.() ?? this.#menu;
    const added = this.#addedItems();
    return added.length ? [...(own ?? []), ...added] : own;
  }
  get hasMenu(): boolean {
    return this.menuSource !== undefined || (this.#menu?.length ?? 0) > 0 || this.#addedItems().length > 0;
  }
  #addedItems(): SceneMenuItem[] {
    const actions = this.added?.value;
    if (!actions?.length) return [];
    return actions
      .filter((a) => a.kinds.includes(this.kind) && (a.when?.(this) ?? true))
      .sort((a, b) => (a.order ?? 0) - (b.order ?? 0))
      .map((a) => ({ id: a.id, label: a.label, ...(a.icon ? { icon: a.icon } : {}), run: () => a.run(this) }));
  }
  set menu(items: readonly SceneMenuItem[] | undefined) {
    this.#menu = items;
  }
  get label(): ReadonlyState<string> {
    return this.#label.asReadonly;
  }
  get description(): ReadonlyState<string | undefined> {
    return this.#description.asReadonly;
  }
  get isHidden(): ReadonlyState<boolean> {
    return this.#hidden.asReadonly;
  }
  get isExpanded(): ReadonlyState<boolean> {
    return this.#expanded.asReadonly;
  }
  get children(): ReadonlyState<readonly SceneNode[]> {
    return this.#children.asReadonly;
  }
  get swatch(): ReadonlyState<number | undefined> | undefined {
    return this.#swatch?.asReadonly;
  }

  /** The row's color (a volume's isosurface); only rows that have one get it. */
  setSwatch(color: number | undefined): void {
    if (!this.#swatch) this.#swatch = observableState(color);
    else this.#swatch.set(color);
  }

  set(label: string, description: string | undefined, hidden: boolean): void {
    this.#label.set(label);
    this.#description.set(description);
    this.#hidden.set(hidden);
  }

  setChildren(children: readonly SceneNode[]): void {
    // A sync makes the arrays again: the same nodes in order are no change (told to nothing, the tree not reflowed).
    if (sameChildren(this.#children.value, children)) return;
    this.#children.set(children);
    this.log?.touch(this.id);
  }

  toggleVisibility(): void {
    this.actions.visibility?.(!this.#hidden.value);
  }
  toggleExpanded(): void {
    const log = this.log;
    const run = () => {
      this.#expanded.set(!this.#expanded.value);
      log?.touch(this.id);
      if (this.#expanded.value && this.lazy) this.load?.();
    };
    // Opened with its children loaded: one change.
    if (log) log.batch(run);
    else run();
  }
  highlight(on: boolean): void {
    this.actions.highlight?.(on);
  }
  focus(): void {
    this.actions.focus?.();
  }
  remove(): void {
    this.actions.remove?.();
  }
}

export interface SceneSource extends SceneModel {
  /** Rebuilds from the hierarchy (called on its changes; tests may call it). */
  sync(): void;
}

/**
 * The scene over the plugin. `api` is bound to the person (its calls are the person's steps); `scope` holds the
 * subscriptions.
 */
/** Opens a command in the command menu with values (a component's own row: the component filled in). */
export type OpenCommand = (id: string, values: Readonly<Record<string, unknown>>) => void;

export function createSceneSource(
  scope: SourceScope,
  plugin: PluginContext,
  api: ViewerApi,
  openCommand?: OpenCommand,
): SceneSource {
  // The scene's structure as one stream, told once a sync (a Mol* transaction) or a toggle (model/scene-changes.ts).
  const log = new SceneChangeLog();
  const nodeActions = observableState<readonly SceneNodeAction[]>([]);
  const nodes = new NodeMap(log, nodeActions.asReadonly);
  // Nothing selected: what acts on the selection is unavailable (a component's Modify by selection).
  const selection = plugin.managers.structure.selection;
  // Read once a task (a select over many structures tells `changed` once a structure, and each read sums the whole
  // selection).
  const noSelection = fromChanged(scope, () => selection.stats.elementCount === 0, onceATask(selection.events.changed));
  const roots = observableState<readonly SceneNode[]>([]);
  const picked = observableState<SceneNode | undefined>(undefined);
  // The tree filtered by the selection (the selection bar's Show in structure tree): which rows hold it, asked by the
  // tree's search per row and kept until the selection changes; followed only while the filter is on (idle is idle).
  const filterOn = observableState(false);
  const filterVersion = observableState(0);
  const holding = new Map<string, boolean>();
  scope.add(
    selection.events.changed.subscribe(() => {
      holding.clear();
      if (!filterOn.value) return;
      // Cleared: the filter goes with it (nothing would show)
      if (selection.stats.elementCount === 0) filterOn.set(false);
      else filterVersion.set(filterVersion.value + 1);
    }),
  );
  const selectedIn = (structureRef: string | undefined) => {
    const data = structureRef ? dataOf(plugin, structureRef) : undefined;
    const l = data && selection.getLoci(data);
    return l && StructureElement.Loci.is(l) && !StructureElement.Loci.isEmpty(l) ? l : undefined;
  };
  const holdsNow = (node: SceneNode): boolean => {
    const refs = plugin.managers.structure.hierarchy.current.refs;
    switch (node.kind) {
      case 'structure':
        return !!selectedIn(node.id);
      case 'contents':
        return !!selectedIn(node.id.slice(0, node.id.indexOf('#')));
      case 'component':
      case 'representation': {
        const ref = refs.get(node.id);
        const component = ref?.kind === 'structure-representation' ? ref.component : ref;
        if (component?.kind !== 'structure-component') return false;
        const sel = selectedIn(component.structure.cell.transform.ref);
        const data = component.cell.obj?.data;
        return !!sel && !!data && !StructureElement.Loci.isEmpty(StructureElement.Loci.remap(sel, data));
      }
      case 'entity':
      case 'chain': {
        const at = node.id.indexOf('#');
        const structureRef = node.id.slice(0, at);
        const sel = selectedIn(structureRef);
        const data = dataOf(plugin, structureRef);
        if (!sel || !data) return false;
        const [, entity, chain] = node.id.slice(at + 1).split(':');
        return StructureElement.Loci.areIntersecting(
          sel,
          contentsLoci(data, entity, node.kind === 'chain' ? chain : undefined),
        );
      }
    }
    return false;
  };
  const selectionFilter: SceneSelectionFilter = {
    on: filterOn.asReadonly,
    version: filterVersion.asReadonly,
    set(on) {
      if (on && selection.stats.elementCount === 0) on = false;
      holding.clear();
      if (filterOn.value !== on) filterOn.set(on);
    },
    holds(node) {
      let h = holding.get(node.id);
      if (h === undefined) holding.set(node.id, (h = holdsNow(node)));
      return h;
    },
  };
  // The frames of the structures that have several (made when the first such structure is shown; its states listen to
  // the plugin only while a row or the inspector reads them).
  const trajectories = () => trajectorySourceOf(plugin);
  /** The frames of particle lists read from a trajectory (Simularium), by list ref */
  const particleTrajectories = new Map<string, ParticleSceneTrajectory>();
  // Inspections by node id, each with what it was built from (`inspectionKey`): a sync refreshes what they hold, and
  // builds one anew when that changed. Only the picked node's is kept once the pick changes (every sync refreshes the
  // ones kept); `inspection` is the picked node's, the inspector's (view-models.md §4.3).
  const inspections = new Map<string, SceneInspection>();
  const inspectionKeys = new Map<string, unknown>();
  const inspection = observableState<SceneInspection | undefined>(undefined);
  /** The structure a unit cell row is under (its operations name the structure). */
  const unitCellOwners = new Map<string, string>();
  /** The structure a volume streaming row is under, and the streaming's behavior cell (its params change). */
  const streamingOwners = new Map<string, string>();
  // The view streaming last had that was not off (and its radius), which the eye brings back.
  const streamingShown = new Map<string, { view: string; radius?: number }>();
  const streamingBehaviors = new Map<string, string>();
  // The structures whose streaming the person turned on: the streaming's node is picked once made, so its controls are in
  // the inspector at once (as a unit cell's when its eye turns it on).
  const pickStreaming = new Set<string>();
  /** How map levels are shown and typed (σ or map units), one choice for every map. */
  const levelUnit = observableState<VolumeLevelUnit>('sigma');
  // The ids to refresh once a sync is done (a volume and each of its representations ask for it).
  let refreshing: Set<string> | undefined;
  const refreshLevels = (id: string) => {
    if (refreshing) {
      refreshing.add(id);
      return;
    }
    const inspection = inspections.get(id);
    for (const l of (inspection?.levels ?? []) as VolumeLevelsSource[]) l.refresh();
    for (const t of (inspection?.transfer ?? []) as VolumeTransferSource[]) t.refresh();
  };

  const structureNode = (s: StructureRef, name: string): Node => {
    const ref = s.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'structure', 'structure', name, true, {
        visibility: (hide) => void api.run('structure.visibility', { structure: ref, visible: !hide }),
        highlight: (on) =>
          on ? void api.run('highlight.set', { target: 'all', structure: ref }) : void api.run('highlight.clear'),
        focus: () => void api.run('focus.set', { target: 'all', structure: ref }),
        remove: () => void api.run('structure.remove', { structure: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const data = s.cell.obj?.data;
    // Several frames (an ensemble, a simulation): its row says how many, and a row under it shows and plays them.
    const frames = s.model?.trajectory?.cell.obj?.data.frameCount ?? 1;
    const unit = FromCoordinates.has(s.model?.trajectory?.cell.transform.transformer.definition.name ?? '')
      ? 'frames'
      : 'models';
    // Gained, lost or changed (its unit): told as a change of the node, so the tree flattens it again (its frames' row).
    const before = n.trajectory;
    if (frames < 2) n.trajectory = undefined;
    else if (!n.trajectory || n.trajectory.unit !== unit)
      n.trajectory = sceneTrajectory(trajectories(), api, ref, unit);
    if (n.trajectory !== before) n.log?.touch(n.id);
    n.set(
      name,
      data ? `${describeStructure(data)}${frames > 1 ? ` · ${frames.toLocaleString('en')} ${unit}` : ''}` : undefined,
      !!s.cell.state.isHidden,
    );
    // Its menu follows what the structure offers (its assemblies, a crystal's symmetry, streaming on or off): built when
    // read (the row's ⋯ opened), from the structure as it is then. Built every sync, it was every operation's offer and
    // every preset's check, for every structure, at every change (structures × definitions a sync).
    if (!n.menuSource) {
      n.menuSource = () => {
        const now = plugin.managers.structure.hierarchy.current.refs.get(ref);
        return now?.kind === 'structure' ? structureMenu(plugin, api, now) : undefined;
      };
    }
    // An open inspector follows the source and the transform (Undo, superposition), without reporting it as an edit.
    if (data) inspections.get(ref)?.form?.load(structureValues(s, data));
    return n;
  };

  /** The unit cell of a structure's model (v5's model-level Unit Cell row; parity.md §2). */
  const unitCellNode = (u: ModelUnitcellRef, structureRef: string): Node => {
    const ref = u.cell.transform.ref;
    unitCellOwners.set(ref, structureRef);
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'unit-cell', 'box', 'Unit cell', false, {
        visibility: (hide) => {
          // Turned on, it is picked: its controls in the inspector.
          if (!hide) picked.set(nodes.get(ref));
          void api.run('structure.unit-cell', { structure: unitCellOwners.get(ref), visible: !hide });
        },
      });
      n.removable = false;
      nodes.set(ref, n);
    }
    n.set('Unit cell', u.cell.obj?.description, !!u.cell.state.isHidden);
    inspections.get(ref)?.form?.load(unitCellValues(u.cell));
    return n;
  };

  /**
   * Volume streaming on a structure (parity.md §2, v5's Volume Streaming node): the view and channels in its inspector,
   * the eye turns the view off and back on, remove turns streaming off.
   */
  const streamingNode = (on: StreamingState, structureRef: string): Node => {
    const ref = on.root;
    streamingOwners.set(ref, structureRef);
    streamingBehaviors.set(on.behavior, ref);
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'extension', 'volume', 'Volume streaming', false, {
        visibility: (hide) => {
          if (!hide) picked.set(nodes.get(ref));
          const owner = streamingOwners.get(ref);
          const last = streamingShown.get(ref);
          const view = hide ? 'off' : (last?.view ?? 'auto');
          void api.run('volume.streaming-view', {
            structure: owner,
            view,
            ...(!hide && last?.radius !== undefined ? { radius: last.radius } : {}),
          });
        },
        remove: () => void api.run('volume.streaming', { structure: streamingOwners.get(ref), on: false }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const view = on.params.entry.params.view;
    if (view.name !== 'off') {
      const radius = (view.params as { radius?: unknown } | undefined)?.radius;
      streamingShown.set(ref, { view: view.name, ...(typeof radius === 'number' ? { radius } : {}) });
    }
    n.set('Volume streaming', streamingSummary(on), view.name === 'off');
    inspections.get(ref)?.form?.load(streamingValues(on));
    refreshLevels(ref);
    if (pickStreaming.delete(structureRef)) picked.set(n);
    return n;
  };

  /**
   * Volume streaming not on yet, for a structure that has density to stream (v5's Volume Streaming row, parity.md §2):
   * a hidden row whose eye turns streaming on.
   */
  const streamingOffNode = (structureRef: string): Node => {
    const id = `${structureRef}#streaming-off`;
    let n = nodes.get(id);
    if (!n) {
      n = new Node(id, 'extension', 'volume', 'Volume streaming', false, {
        visibility: (hide) => {
          if (hide) return;
          // The streaming's own node takes its place, picked (`streamingNode`): this row's inspector would stay on a
          // row that is gone.
          pickStreaming.add(structureRef);
          void api.run('volume.streaming', { structure: structureRef, on: true }).then((r) => {
            if (!r.ok) pickStreaming.delete(structureRef);
          });
        },
      });
      n.removable = false;
      nodes.set(id, n);
    }
    n.set('Volume streaming', 'off', true);
    return n;
  };

  const componentNode = (c: StructureComponentRef, structureRef: string): Node => {
    const ref = c.cell.transform.ref;
    let n = nodes.get(ref);
    // A bare row (made by `emptiedNode` with no component row before it, as after a restore) has no actions: make it whole.
    if (!n || !n.actions.visibility) {
      const loci = () => {
        const data = c.cell.obj?.data;
        return data ? Structure_toLoci(data) : undefined;
      };
      n = new Node(ref, 'component', 'layers', c.cell.obj?.label ?? 'Component', false, {
        visibility: (hide) =>
          void api.run('component.visibility', { component: ref, structure: structureRef, visible: !hide }),
        highlight: (on) => {
          const l = loci();
          if (on && l) void api.run('highlight.set', { target: { loci: l } });
          else void api.run('highlight.clear');
        },
        focus: () => {
          const l = loci();
          if (l) void api.run('focus.set', { target: { loci: l } });
        },
        remove: () => void api.run('component.remove', { component: ref, structure: structureRef }),
      });
      n.removable = true;
      n.menu = componentMenu(plugin, api, ref, structureRef, loci, openCommand, c.cell.obj?.label ?? '', noSelection);
      nodes.set(ref, n);
    }
    const reprs = c.representations.map((r) => reprLabel(plugin, r)).join(', ');
    const atoms = c.cell.obj?.data.elementCount;
    n.set(
      (c.cell.obj?.label ?? 'Component').replace(FocusPrefix, ''),
      [reprs, atoms !== undefined ? `${atoms.toLocaleString('en')} atoms` : undefined].filter(Boolean).join(' · ') ||
        undefined,
      !!c.cell.state.isHidden,
    );
    return n;
  };

  const reprNode = (r: StructureRepresentationRef): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    const type = (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name ?? '';
    if (!n) {
      n = new Node(ref, 'representation', ReprIcons[type] ?? 'cartoon', reprLabel(plugin, r), false, {
        visibility: (hide) => void api.run('representation.visibility', { node: ref, visible: !hide }),
        // This one by its ref: its type may have changed since, and another of that type may sit beside it.
        remove: () => void api.run('representation.remove', { node: ref }),
      });
      n.removable = true;
      const node = n;
      n.menuSource = () => representationMenu(plugin, api, ref, () => picked.set(node));
      nodes.set(ref, n);
    }
    n.icon = ReprIcons[type] ?? 'cartoon';
    n.set(reprLabel(plugin, r), colorLabel(plugin, r), !!r.cell.state.isHidden);
    // An open inspector follows the cell (Undo, another surface), without reporting it as an edit.
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    (inspections.get(ref)?.quick as RepresentationQuickSource | undefined)?.refresh();
    return n;
  };

  /** A structure's contents: entities, then their chains; computed when first expanded. */
  const contentsNode = (s: StructureRef, structureRef: string): Node => {
    const id = `${structureRef}#contents`;
    let n = nodes.get(id);
    if (!n) {
      const node = new Node(id, 'contents', 'list-tree', 'Contents', false, {});
      node.lazy = true;
      node.removable = false;
      n = node;
      nodes.set(id, n);
    }
    // The structure as it is now: a load computes the rows, and once they are shown each sync refreshes them (what is
    // hidden, a changed structure). A search's look (`peek`) leaves it lazy: computed once until the next sync, so a
    // search does not make every later sync compute every structure's rows.
    const node = n;
    let peeked = false;
    node.load = (peek) => {
      if (peek && peeked) return;
      peeked = !!peek;
      if (!peek) node.lazy = false;
      const data = s.cell.obj?.data;
      // Its entities and their chains: one change.
      if (data) log.batch(() => node.setChildren(entities(plugin, data, structureRef, nodes, api)));
    };
    if (!node.lazy) node.load();
    const data = s.cell.obj?.data;
    n.set('Contents', data ? `${data.model.entities.data._rowCount} entities` : undefined, false);
    return n;
  };

  // The focus: one row at the top level, after the structures (the focus is one, whichever structure it is on), always
  // there so the tree does not shift as it comes and goes; its rows are the focus behavior's components (the target, its
  // surroundings) while something is focused.
  const FocusId = 'scene#focus';
  const focusNode = (
    parts: readonly { component: StructureComponentRef; structure: string }[],
    name?: string,
  ): Node => {
    let n = nodes.get(FocusId);
    if (!n) {
      const target = () => {
        const data = plugin.managers.structure.focus.current?.loci;
        return data && !StructureElement.Loci.isEmpty(data) ? data : undefined;
      };
      n = new Node(FocusId, 'focus', 'focus', 'Focus', false, {
        // The parts are one step (one undo brings back all of them).
        visibility: (hide) => {
          const parts = focusPartsNow;
          if (parts.length === 0) return;
          void api.group('Focus visibility', async (g) => {
            for (const { component, structure } of parts) {
              await g.run('component.visibility', {
                component: component.cell.transform.ref,
                structure,
                visible: !hide,
              });
            }
          });
        },
        highlight: (on) => {
          const l = target();
          if (on && l) void api.run('highlight.set', { target: { loci: l } });
          else void api.run('highlight.clear');
        },
        focus: () => {
          const l = target();
          if (l) void api.run('focus.set', { target: { loci: l } });
        },
      });
      n.removable = false;
      n.menu = [{ id: 'clear', label: 'Clear focus', icon: 'x', run: () => void api.run('focus.clear') }];
      nodes.set(FocusId, n);
    }
    focusPartsNow = parts;
    const atoms = parts.find((p) =>
      p.component.cell.transform.tags?.includes(StructureFocusRepresentationTags.TargetSel),
    )?.component.cell.obj?.data.elementCount;
    n.set(
      'Focus',
      atoms ? `${atoms.toLocaleString('en')} ${atoms === 1 ? 'atom' : 'atoms'}${name ? ` · ${name}` : ''}` : 'none',
      parts.length > 0 && parts.every((p) => !!p.component.cell.state.isHidden),
    );
    return n;
  };
  let focusPartsNow: readonly { component: StructureComponentRef; structure: string }[] = [];

  // A component left empty by hidden contents (the hierarchy leaves it out): its row stays, hidden, saying why; showing
  // the contents again brings it back as it was.
  const emptiedNode = (ref: string, label: string): Node => {
    let n = nodes.get(ref);
    if (!n || n.kind !== 'component') {
      n = new Node(ref, 'component', 'layers', label, false, {});
      n.removable = false;
      nodes.set(ref, n);
    }
    n.set(label, 'hidden in Contents', true);
    n.setChildren([]);
    return n;
  };

  const measurementsNode = (): Node => {
    const id = '#measurements';
    let n = nodes.get(id);
    if (!n) {
      n = new Node(id, 'measurements', 'ruler', 'Measurements', true, {});
      n.removable = false;
      nodes.set(id, n);
    }
    return n;
  };

  const measurementNode = (m: MeasuredCell): Node => {
    let n = nodes.get(m.ref);
    if (!n) {
      const ref = m.ref;
      const locis = () => measurementLocis(plugin, ref);
      n = new Node(ref, 'measurement', MeasureKinds[m.kind].icon, m.label, false, {
        visibility: (hide) => void api.run('measure.visibility', { measurement: ref, visible: !hide }),
        highlight: (on) => {
          const l = locis()[0];
          if (on && l) void api.run('highlight.set', { target: { loci: l } });
          else void api.run('highlight.clear');
        },
        focus: () => {
          const l = locis()[0];
          if (l) void api.run('focus.set', { target: { loci: l } });
        },
        remove: () => void api.run('measure.remove', { measurement: ref }),
      });
      // Its eye and remove need the measurements feature's operations.
      n.removable = !!api.describe('measure.remove');
      nodes.set(ref, n);
    }
    n.set(m.value ? `${m.value}` : m.label, m.value ? m.label : MeasureKinds[m.kind].title, m.hidden);
    return n;
  };

  // Volumes and particles: top-level nodes like structures, their representations as children. Eyes and
  // removes go through `volume.*` and `particles.*`.
  const volumeNode = (v: VolumeRef, name: string): Node => {
    const ref = v.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'volume', 'volume', name, true, {
        visibility: (hide) => void api.run('volume.visibility', { volume: ref, visible: !hide }),
        remove: () => void api.run('volume.remove', { volume: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const data = v.cell.obj?.data;
    const stats = data && volumeStats(data);
    n.set(name, stats ? `${gridText(stats)} · ${voxelText(stats)}` : undefined, !!v.cell.state.isHidden);
    refreshLevels(ref);
    return n;
  };

  const volumeReprNode = (r: VolumeRepresentationRef, volumeRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'volume-representation', 'volume', '', false, {
        visibility: (hide) => void api.run('volume.visibility', { node: ref, visible: !hide }),
        remove: () => void api.run('volume.remove', { node: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const p = r.cell.transform.params as
      | {
          type?: { params?: { isoValue?: Volume.IsoValue } };
          colorTheme?: { name: string; params?: { value?: number } };
        }
      | undefined;
    const iso = p?.type?.params?.isoValue;
    n.set(
      r.cell.obj?.label ?? 'Representation',
      iso ? Volume.IsoValue.toString(iso) : undefined,
      !!r.cell.state.isHidden,
    );
    n.setSwatch(p?.colorTheme?.name === 'uniform' ? p.colorTheme.params?.value : undefined);
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    refreshLevels(ref);
    refreshLevels(volumeRef);
    return n;
  };

  const particlesNode = (l: ParticleListRef, name: string): Node => {
    const ref = l.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'particles', 'shapes', name, true, {
        visibility: (hide) => void api.run('particles.visibility', { particles: ref, visible: !hide }),
        remove: () => void api.run('particles.remove', { particles: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const data = particleListData(plugin.state.data, ref);
    n.set(name, data ? particlesText(data) : undefined, !!l.cell.state.isHidden);
    // Built when the row's ⋯ opens, from the list as it is then (its representations, its unit cell).
    if (!n.menuSource) n.menuSource = () => particlesMenu(plugin, api, ref, openCommand);
    inspections.get(ref)?.form?.load(l.cell.transform.params as never);
    // A trajectory's frames (Simularium): a row under it shows and plays them, as a structure's (`particles.frame`)
    const frames = particleFrames(plugin, l.cell);
    const before = n.trajectory;
    if (!frames) {
      n.trajectory = undefined;
      particleTrajectories.delete(ref);
    } else {
      let t = particleTrajectories.get(ref);
      if (!t) particleTrajectories.set(ref, (t = particleTrajectory(plugin, api, ref)));
      t.refresh(frames);
      n.trajectory = t;
    }
    if (n.trajectory !== before) n.log?.touch(n.id);
    return n;
  };

  /** A particle list's unit cell (its periodic box): its eye and inspector through `particles.unit-cell`. */
  const particlesUnitCellNode = (cell: StateObjectCell, listRef: string): Node => {
    const ref = cell.transform.ref;
    unitCellOwners.set(ref, listRef);
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'unit-cell', 'box', 'Unit cell', false, {
        visibility: (hide) => {
          if (!hide) picked.set(nodes.get(ref));
          void api.run('particles.unit-cell', { particles: listRef, visible: !hide });
        },
      });
      n.removable = false;
      nodes.set(ref, n);
    }
    n.set('Unit cell', undefined, !!cell.state.isHidden);
    inspections.get(ref)?.form?.load(unitCellValues(cell));
    return n;
  };

  const particlesReprNode = (r: ParticleRepresentationRef, listRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'particles-representation', 'shapes', '', false, {
        visibility: (hide) => void api.run('particles.visibility', { node: ref, visible: !hide }),
        remove: () => void api.run('particles.remove', { node: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const p = r.cell.transform.params as
      | { type?: { name: string }; colorTheme?: { name: string; params?: { value?: number } } }
      | undefined;
    const { themes } = plugin.representation.particles;
    const theme = p?.colorTheme?.name;
    const themeLabel =
      theme && themes.colorThemeRegistry.has(theme) ? themes.colorThemeRegistry.get(theme).label : theme;
    n.set(r.cell.obj?.label ?? 'Representation', themeLabel, !!r.cell.state.isHidden);
    n.setSwatch(theme === 'uniform' ? p?.colorTheme?.params?.value : undefined);
    if (!n.menuSource)
      n.menuSource = () => particlesRepresentationMenu(plugin, api, ref, () => picked.set(nodes.get(ref)));
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    return n;
  };

  // The rows of each structure as last built, kept while what they were built from is the same (`Signature`): a sync
  // builds only the structures that changed, so loading N structures one by one is N builds, not N² (state.md §4;
  // AGENTS.md, bulk operations). The clean-up follows: how many rows hold each id, and only the ids that lost theirs
  // (or are new) are looked at.
  const structureRows = new Map<string, StructureRows>();
  const seenCount = new Map<string, number>();
  /** The ids held by the rows after the structures (the focus, volumes, particles, measurements) at the last sync. */
  let restIds: readonly string[] = [];
  /** Structures built anew at the next sync whatever their signature (a streaming cell came or went under one). */
  const stale = new Set<string>();
  const markStale = (ref: string | undefined) => {
    const tree = plugin.state.data.tree;
    for (let at = ref; at; ) {
      if (structureRows.has(at)) {
        stale.add(at);
        return;
      }
      const parent: string | undefined = tree.transforms.get(at)?.parent;
      if (!parent || parent === at) return;
      at = parent;
    }
  };

  /** Writes what a structure's rows are built from into its signature; whether it is as before. */
  const signed = (
    rows: StructureRows,
    s: StructureRef,
    name: string,
    unitCell: ModelUnitcellRef | undefined,
    emptied: readonly EmptiedComponent[],
  ): boolean => {
    const state = plugin.state.data;
    const tree = state.tree;
    const sig = rows.signature;
    const ref = s.cell.transform.ref;
    sig.begin();
    sig.add(name);
    sig.add(s.cell);
    sig.add(s.cell.transform);
    sig.add(s.cell.obj);
    sig.add(!!s.cell.state.isHidden);
    // Its children come and go (volume streaming's info, a component).
    sig.add(tree.children.get(ref));
    sig.add(s.transform?.cell.transform);
    sig.add(s.model?.cell.obj);
    sig.add(s.model?.trajectory?.cell.transform);
    sig.add(s.model?.trajectory?.cell.obj);
    sig.add(s.volumeStreaming?.cell.transform);
    sig.add(s.volumeStreaming?.cell.obj);
    if (rows.streaming) {
      sig.add(tree.children.get(rows.streaming.root));
      sig.add(state.cells.get(rows.streaming.root)?.obj);
      sig.add(tree.transforms.get(rows.streaming.behavior));
    }
    sig.add(unitCell?.cell.transform);
    sig.add(unitCell?.cell.obj);
    sig.add(!!unitCell?.cell.state.isHidden);
    sig.add(hiddenContents(plugin, ref));
    sig.add(s.components.length);
    for (const c of s.components) {
      sig.add(c.cell);
      sig.add(c.cell.transform);
      sig.add(c.cell.obj);
      sig.add(!!c.cell.state.isHidden);
      sig.add(c.representations.length);
      for (const r of c.representations) {
        sig.add(r.cell);
        sig.add(r.cell.transform);
        sig.add(r.cell.obj);
        sig.add(!!r.cell.state.isHidden);
      }
    }
    sig.add(emptied.length);
    for (const e of emptied) {
      sig.add(e.ref);
      sig.add(e.label);
      sig.add(tree.children.get(e.parent));
    }
    return sig.end();
  };

  /** Builds a structure's rows (its node, its children); the ids they hold. */
  const buildStructure = (
    rows: StructureRows,
    s: StructureRef,
    name: string,
    unitCell: ModelUnitcellRef | undefined,
    emptied: readonly EmptiedComponent[],
  ): string[] => {
    const ref = s.cell.transform.ref;
    const ids: string[] = [];
    const sn = structureNode(s, name);
    ids.push(ref);
    const contents = contentsNode(s, ref);
    ids.push(contents.id);
    const kids: SceneNode[] = [contents];
    const focused: SceneNode[] = [];
    for (const c of s.components) {
      const isFocus = isFocusPart(c);
      // An emptied focus component (the focus cleared) has no row.
      if (isFocus && !c.cell.obj?.data.elementCount) continue;
      const cn = componentNode(c, ref);
      if (isFocus) cn.removable = false;
      ids.push(cn.id);
      const reprs = c.representations.map((r) => {
        const rn = reprNode(r);
        ids.push(rn.id);
        return rn;
      });
      cn.setChildren(reprs);
      (isFocus ? focused : kids).push(cn);
    }
    // Components emptied by hidden contents, where they stand among the others (the state tree's order).
    if (emptied.length) {
      const order = new Map<string, number>();
      let i = 0;
      for (const child of plugin.state.data.tree.children.get(emptied[0].parent).values()) order.set(child, i++);
      for (const e of emptied) {
        const en = emptiedNode(e.ref, e.label);
        ids.push(en.id);
        kids.push(en);
      }
      // The contents row first; the components by their place in the tree (one not there first, as `indexOf`'s −1).
      const [first, ...rest] = kids;
      rest.sort((a, b) => (order.get(a.id) ?? -1) - (order.get(b.id) ?? -1));
      kids.splice(0, kids.length, first, ...rest);
    }
    rows.parts = focusParts(plugin, ref).filter((c) => !!c.cell.obj?.data.elementCount);
    rows.focused = focused;
    // The model's unit cell, under the first of its structures (another structure of the same model shares it).
    if (unitCell) {
      const un = unitCellNode(unitCell, ref);
      ids.push(un.id);
      kids.push(un);
    }
    const streaming = streamingOf(plugin, s);
    rows.streaming = streaming ? { root: streaming.root, behavior: streaming.behavior } : undefined;
    if (streaming) {
      const vn = streamingNode(streaming, ref);
      ids.push(vn.id);
      kids.push(vn);
    } else if (canStream(plugin, s)) {
      const vn = streamingOffNode(ref);
      ids.push(vn.id);
      kids.push(vn);
    }
    sn.setChildren(kids);
    rows.node = sn;
    return ids;
  };

  const syncRows = () => {
    const list = plugin.managers.structure.hierarchy.current.structures;
    const names = structureNames(list);
    const top: SceneNode[] = [];
    const allFocused: SceneNode[] = [];
    const allFocusParts: { component: StructureComponentRef; structure: string }[] = [];
    let focusedOn: string | undefined;
    // The ids whose last holder let go during this sync: the clean-up's candidates.
    const dropped: string[] = [];
    const hold = (ids: readonly string[]) => {
      for (const id of ids) seenCount.set(id, (seenCount.get(id) ?? 0) + 1);
    };
    const release = (ids: readonly string[]) => {
      for (const id of ids) {
        const count = (seenCount.get(id) ?? 0) - 1;
        if (count > 0) seenCount.set(id, count);
        else {
          seenCount.delete(id);
          dropped.push(id);
        }
      }
    };
    const present = new Set<string>();
    for (const s of list) present.add(s.cell.transform.ref);
    // Read once for every structure.
    const emptiedBy = emptiedComponentsOf(plugin, present);
    // The unit cells placed so far (each under the first structure of its model).
    const placed = new Set<string>();
    for (const s of list) {
      const ref = s.cell.transform.ref;
      const name = names.get(ref) ?? ref;
      const emptied = emptiedBy.get(ref) ?? NoEmptied;
      const u = s.model?.unitcell;
      const unitCell = u?.cell.obj?.data?.repr && !placed.has(u.cell.transform.ref) ? u : undefined;
      if (unitCell) placed.add(unitCell.cell.transform.ref);
      let rows = structureRows.get(ref);
      if (!rows) {
        rows = { signature: new Signature(), ids: [], focused: [], parts: [] };
        structureRows.set(ref, rows);
      }
      const same =
        signed(rows, s, name, unitCell, emptied) && !stale.has(ref) && !!rows.node && nodes.get(ref) === rows.node;
      if (same) {
        // Its streaming's levels follow the density as it arrives.
        if (rows.streaming) refreshLevels(rows.streaming.root);
      } else {
        const ids = buildStructure(rows, s, name, unitCell, emptied);
        hold(ids);
        release(rows.ids);
        rows.ids = ids;
        // What it was built from now (its streaming read anew).
        signed(rows, s, name, unitCell, emptied);
      }
      if (rows.parts.length) focusedOn = name;
      for (const component of rows.parts) allFocusParts.push({ component, structure: ref });
      for (const f of rows.focused) allFocused.push(f);
      top.push(rows.node!);
    }
    for (const [ref, rows] of structureRows) {
      if (present.has(ref)) continue;
      release(rows.ids);
      structureRows.delete(ref);
    }
    stale.clear();
    const rest: string[] = [];
    if (list.length) {
      // Named by its structure when there are several.
      const fn = focusNode(allFocusParts, list.length > 1 ? focusedOn : undefined);
      rest.push(fn.id);
      fn.setChildren(allFocused);
      top.push(fn);
    }
    const volumes = plugin.managers.volume.hierarchy.current.volumes;
    const volumeLabels = namesOf(volumes);
    for (const v of volumes) {
      const ref = v.cell.transform.ref;
      const vn = volumeNode(v, volumeLabels.get(ref) ?? ref);
      rest.push(ref);
      vn.setChildren(
        v.representations.map((r) => {
          const rn = volumeReprNode(r, ref);
          rest.push(rn.id);
          return rn;
        }),
      );
      top.push(vn);
    }
    const lists = plugin.managers.particles.hierarchy.current.lists;
    const listLabels = namesOf(lists);
    for (const l of lists) {
      const ref = l.cell.transform.ref;
      const ln = particlesNode(l, listLabels.get(ref) ?? ref);
      rest.push(ref);
      const unitCell = particleUnitCellOf(plugin.state.data, ref);
      ln.setChildren([
        ...l.representations.map((r) => {
          const rn = particlesReprNode(r, ref);
          rest.push(rn.id);
          return rn;
        }),
        ...(unitCell ? [particlesUnitCellNode(unitCell, ref)] : []).map((un) => {
          rest.push(un.id);
          return un;
        }),
      ]);
      top.push(ln);
    }
    // Measurements (they may span structures): one group after the data, a row each; their eyes and removes go
    // through `measure.*` (the measurements feature), their hover and focus through the measured atoms.
    // The group is there once something is loaded, empty or not: its inspector says how to measure.
    const measured = readMeasurements(plugin);
    if (measured.length || top.length) {
      const gn = measurementsNode();
      rest.push(gn.id);
      gn.set('Measurements', measurementCounts(measured), false);
      gn.setChildren(
        measured.map((m) => {
          const mn = measurementNode(m);
          rest.push(mn.id);
          return mn;
        }),
      );
      top.push(gn);
    }
    hold(rest);
    release(restIds);
    restIds = rest;
    // The clean-up: a node stays while a row holds its id or its owner's (lazy contents' entity and chain nodes stay
    // with their structure). Only what may have lost that is looked at: the ids let go, the ids they own, new ones.
    const gone = (id: string) => !seenCount.has(id) && !seenCount.has(ownerOf(id));
    const drop = (id: string) => {
      const n = nodes.get(id);
      if (!n || !gone(id)) return;
      nodes.delete(id);
      inspections.delete(id);
      inspectionKeys.delete(id);
      unitCellOwners.delete(id);
      particleTrajectories.delete(id);
      streamingOwners.delete(id);
      streamingShown.delete(id);
      if (picked.value === n) picked.set(undefined);
    };
    for (const id of dropped) {
      drop(id);
      const owned = nodes.under(id);
      if (owned) for (const o of [...owned]) drop(o);
    }
    for (const id of [...nodes.fresh]) drop(id);
    nodes.fresh.clear();
    // An inspection built from what has changed since (a structure's assembly or model, a volume's representations) is
    // built anew: its facts, contents, levels and transfer functions are fixed when built. Its form stays (it follows
    // the cell as it is, and the inspector keeps it, a gesture on it going on).
    for (const [id, key] of inspectionKeys) {
      const n = nodes.get(id);
      if (!n || Object.is(inspectionKey(plugin, n), key)) continue;
      const form = inspections.get(id)?.form;
      const fresh = build(n);
      inspections.set(id, form && fresh.form ? { ...fresh, form } : fresh);
      inspectionKeys.set(id, inspectionKey(plugin, n));
      if (picked.value === n) inspection.set(inspections.get(id));
    }
    // The same roots (a sync makes the array again): no change, the tree not flattened again.
    if (!sameChildren(roots.value, top)) roots.set(top);
  };

  const sync = () => {
    const outer = refreshing;
    const ids = outer ?? new Set<string>();
    refreshing = ids;
    try {
      log.batch(syncRows);
    } finally {
      if (!outer) {
        refreshing = undefined;
        for (const id of ids) refreshLevels(id);
      }
    }
  };

  const build = (n: SceneNode) =>
    inspect(plugin, api, n, unitCellOwners.get(n.id) ?? streamingOwners.get(n.id), levelUnit, (row) => picked.set(row));
  const inspectNode = (n: SceneNode): SceneInspection => {
    let built = inspections.get(n.id);
    if (!built) {
      built = build(n);
      inspections.set(n.id, built);
      inspectionKeys.set(n.id, inspectionKey(plugin, n));
    }
    return built;
  };

  scope.add(
    picked.asReadonly.subscribe((n) => {
      for (const id of inspections.keys()) {
        if (id === n?.id) continue;
        inspections.delete(id);
        inspectionKeys.delete(id);
      }
      inspection.set(n && inspectNode(n));
    }),
  );

  // Whether particle frames play follows the animation (played, stopped), not the frame shown
  scope.add(
    plugin.managers.animation.events.updated.subscribe(() => {
      for (const t of particleTrajectories.values()) t.refreshPlaying();
    }),
  );
  scope.add(plugin.managers.structure.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.volume.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.particles.hierarchy.behaviors.selection.subscribe(() => sync()));
  // The measurements (one added, removed, or given new params: its label, its unit): the manager says its state
  // changed on every state change, so the rows sync only when a list holds other cells, or a cell another transform or
  // object (a cell keeps its identity through an update).
  const MeasurementKinds = ['labels', 'distances', 'angles', 'dihedrals', 'orientations', 'planes'] as const;
  let measured: readonly unknown[] = [];
  scope.add(
    plugin.managers.structure.measurement.behaviors.state.subscribe((state) => {
      const now: unknown[] = [];
      for (const k of MeasurementKinds) {
        const list = state[k] as readonly { transform: unknown; obj?: unknown }[];
        now.push(list.length);
        for (const c of list) now.push(c, c.transform, c.obj);
      }
      const same = now.length === measured.length && now.every((x, i) => x === measured[i]);
      measured = now;
      if (!same) sync();
    }),
  );
  // Streaming's cells are not in the hierarchies: their coming and going rebuilds, a params change refreshes the row.
  const isStreaming = (t: { transformer: unknown } | undefined) =>
    t?.transformer === CreateVolumeStreamingInfo || t?.transformer === CreateVolumeStreamingBehavior;
  scope.add(
    plugin.state.data.events.cell.created.subscribe((e) => {
      if (!isStreaming(e.cell.transform)) return;
      // Its structure's rows are built anew (the hierarchy may not have it yet).
      markStale(e.ref);
      sync();
    }),
  );
  scope.add(
    plugin.state.data.events.cell.removed.subscribe((e) => {
      if (!streamingOwners.has(e.ref) && !streamingBehaviors.has(e.ref)) return;
      streamingBehaviors.delete(e.ref);
      markStale(e.parent);
      sync();
    }),
  );
  scope.add(
    plugin.state.data.events.object.updated.subscribe((e) => {
      // A channel's visual under the behavior: new density arrived (the levels' histogram).
      const parent = plugin.state.data.tree.transforms.get(e.ref)?.parent;
      const above = parent !== undefined ? streamingBehaviors.get(parent) : undefined;
      if (above) refreshLevels(above);
      const root = streamingBehaviors.get(e.ref);
      const owner = root ? streamingOwners.get(root) : undefined;
      const s = owner ? plugin.managers.structure.hierarchy.current.refs.get(owner) : undefined;
      const on = s?.kind === 'structure' ? streamingOf(plugin, s) : undefined;
      if (on) streamingNode(on, owner!);
    }),
  );
  // A streaming behavior given new params (`volume.streaming-view`) takes them in place: no object event, only a new
  // transform. Its row follows its transform, read on each state change (one lookup per streaming row), not the
  // measurement manager's state, which said "changed" on every state change only by accident.
  const streamingSeen = new Map<string, unknown>();
  scope.add(
    plugin.state.data.events.changed.subscribe(() => {
      if (!streamingBehaviors.size) return;
      for (const [behavior, root] of streamingBehaviors) {
        const t = plugin.state.data.tree.transforms.get(behavior);
        if (!t || streamingSeen.get(behavior) === t) continue;
        streamingSeen.set(behavior, t);
        const owner = streamingOwners.get(root);
        const s = owner ? plugin.managers.structure.hierarchy.current.refs.get(owner) : undefined;
        const on = s?.kind === 'structure' ? streamingOf(plugin, s) : undefined;
        if (on) streamingNode(on, owner!);
      }
      for (const behavior of streamingSeen.keys())
        if (!streamingBehaviors.has(behavior)) streamingSeen.delete(behavior);
    }),
  );
  scope.add(
    plugin.state.data.events.cell.stateUpdated.subscribe((e) => {
      const n = nodes.get(e.ref);
      if (n) n.set(n.label.value, n.description.value, !!e.cell.state.isHidden);
    }),
  );

  return {
    roots: roots.asReadonly,
    picked: picked.asReadonly,
    pick: (n) => picked.set(n === picked.value ? undefined : n),
    selectionFilter,
    nodeActions,
    reveal(id) {
      const path = (list: readonly SceneNode[]): SceneNode[] | undefined => {
        for (const n of list) {
          if (n.id === id) return [n];
          const below = path(n.children.value);
          if (below) return [n, ...below];
        }
        return undefined;
      };
      const found = path(roots.value);
      if (!found) return false;
      // The path opened in one change.
      log.batch(() => {
        for (const n of found.slice(0, -1)) if (!n.isExpanded.value) n.toggleExpanded();
      });
      picked.set(found[found.length - 1]);
      return true;
    },
    changes: log.changes,
    inspection: inspection.asReadonly,
    inspect: inspectNode,
    sync,
  };
}

/**
 * What a node's inspection is built from, where the inspection fixes it when built: a structure's facts and its
 * contents' summary come from the structure (an assembly or a model changes it), a volume's levels and transfer functions
 * from its representations and their types, an isosurface's or a direct volume's from its type. Others: nothing.
 */
function inspectionKey(plugin: PluginContext, node: SceneNode): unknown {
  switch (node.kind) {
    case 'structure':
      return dataOf(plugin, node.id);
    case 'contents':
      return dataOf(plugin, node.id.slice(0, node.id.indexOf('#')));
    case 'volume': {
      let key = '';
      for (const r of volumeRefOf(plugin, node.id)?.representations ?? []) {
        key += `${r.cell.transform.ref}:${typeName(r.cell.transform.params)};`;
      }
      return key;
    }
    case 'volume-representation':
      return typeName(plugin.state.data.cells.get(node.id)?.transform.params);
    case 'particles':
      // its facts follow the list as its representations see it (a target attached)
      return particleListData(plugin.state.data, node.id);
  }
  return undefined;
}

function typeName(params: unknown): string {
  return (params as { type?: { name: string } } | undefined)?.type?.name ?? '';
}

/** A volume of the hierarchy by its ref. */
function volumeRefOf(plugin: PluginContext, ref: string): VolumeRef | undefined {
  return plugin.managers.volume.hierarchy.current.volumes.find((v) => v.cell.transform.ref === ref);
}

function Structure_toLoci(s: Structure): StructureElement.Loci {
  return StructureElement.Loci.all(s);
}

/** The transformer that adds coordinates (a simulation's frames) to a model: its frames are `frames`, not `models`. */
// (MolViewSpec's `trajectory-with-coordinates` too: the same trajectory, built from an MVS coordinates node)
const FromCoordinates = new Set(['trajectory-from-model-and-coordinates', 'trajectory-with-coordinates']);

/**
 * A structure's frames as its row and the inspector use them (model/scene.ts, `SceneTrajectory`): read from the shared
 * trajectory source (its own by ref), acting through the API. A seek in flight takes only the latest index after it (a
 * scrubber's drag); a step counts from the frame last asked for (pending or in flight), so quick presses add up.
 */
function sceneTrajectory(
  source: TrajectorySource,
  api: ViewerApi,
  structure: string,
  unit: SceneTrajectory['unit'],
): SceneTrajectory {
  const info = derive(source.byRef, (m) => m.get(structure));
  const count = derive(info, (i) => i?.count ?? 0);
  const index = derive(source.indexByRef, (m) => m.get(structure) ?? 0);
  // Every structure's (a play of all) or this one's.
  const playing = derive(source.playingRefs, (r) => r.has(structure));
  const mode = observableState<'loop' | 'palindrome' | 'once'>('loop');
  const seconds = observableState(5);
  /** The frame to show next (not sent yet), and the one being shown (sent, not done). */
  let wanted: number | undefined;
  let sending: number | undefined;
  const flush = async () => {
    if (sending !== undefined) return;
    try {
      while (wanted !== undefined) {
        sending = wanted;
        wanted = undefined;
        await api.run('trajectory.frame', { structure, index: sending });
      }
    } finally {
      sending = undefined;
    }
  };
  return {
    unit,
    count,
    index,
    playing,
    mode,
    seconds,
    seek(i) {
      wanted = i;
      void flush();
    },
    step(by) {
      const c = count.value;
      if (c < 2) return;
      // From the frame last asked for (pending, else in flight), not the frame shown: quick presses add up.
      wanted = ((((wanted ?? sending ?? index.value) + by) % c) + c) % c;
      void flush();
    },
    togglePlay() {
      void (playing.value
        ? api.run('trajectory.stop')
        : api.run('trajectory.play', { structure, mode: mode.value, seconds: seconds.value }));
    },
    markShown(on) {
      source.markRowShown(structure, on);
    },
  };
}

function describeStructure(s: Structure): string {
  return `${s.elementCount.toLocaleString('en')} atoms`;
}

function reprLabel(plugin: PluginContext, r: StructureRepresentationRef): string {
  const name = (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name;
  if (!name) return r.cell.obj?.label ?? 'Representation';
  return plugin.representation.structure.registry.get(name)?.label ?? name;
}

function colorLabel(plugin: PluginContext, r: StructureRepresentationRef): string | undefined {
  const name = (r.cell.transform.params as { colorTheme?: { name: string } } | undefined)?.colorTheme?.name;
  if (!name) return undefined;
  return plugin.representation.structure.themes.colorThemeRegistry.get(name)?.label ?? name;
}

/** Mol*'s component menu: add representation, color, modify by the selection, select this. */
/** A component changed by the current selection (v5's Modify by Selection): include, subtract or intersect it. */
function modifyBySelection(
  api: ViewerApi,
  ref: string,
  structureRef: string,
  noSelection: Observable<boolean>,
): SceneMenuItem {
  return {
    id: 'modify',
    label: 'Modify by selection',
    icon: 'squares-unite',
    // Three icons on its row, there at all times, unavailable with nothing selected (David, 2026-10-15).
    inline: true,
    items: (['add', 'remove', 'intersect'] as const).map((mode) => ({
      disabled: noSelection,
      id: mode,
      label:
        mode === 'add'
          ? 'Include the selection'
          : mode === 'remove'
            ? 'Subtract the selection'
            : 'Intersect with the selection',
      icon: mode === 'add' ? 'squares-unite' : mode === 'remove' ? 'squares-subtract' : 'squares-intersect',
      run: () =>
        void api.run('component.modify', { component: ref, structure: structureRef, target: 'selection', mode }),
    })),
  };
}

/** The preset components' labels as target keywords (targets.md §1). */
const ComponentKeywords: Readonly<Record<string, string>> = {
  polymer: 'polymer',
  protein: 'protein',
  'nucleic acid': 'nucleic',
  ligand: 'ligand',
  ion: 'ion',
  water: 'water',
  carbohydrate: 'branched',
  branched: 'branched',
};

function componentMenu(
  plugin: PluginContext,
  api: ViewerApi,
  ref: string,
  structureRef: string,
  loci: () => StructureElement.Loci | undefined,
  openCommand: OpenCommand | undefined,
  label: string,
  noSelection: Observable<boolean>,
): SceneMenuItem[] {
  // With the command menu: each opens its command for this component (its fields filled in, the rest to choose there),
  // so the long lists (types, themes) are the command menu's to search, not a fold in a panel.
  if (openCommand) {
    // A preset component (Polymer, Ligand, Water, …) is also a target by its keyword; another is typed there.
    const keyword = ComponentKeywords[label.toLowerCase()];
    // Named as the person reads them (the command's fields show their values), when clicked: the structure's name and
    // the component's label, each its state ref only where the name is not unique.
    const named = () => {
      const structures = plugin.managers.structure.hierarchy.current.structures;
      // The tree's names (the entry, "1CBS"); one two structures share reads "1CBS · Assembly 1": then the ref.
      const name = structureNames(structures).get(structureRef);
      const structure = name && !name.includes(' · ') ? name : structureRef;
      const siblings = structures.find((s) => s.cell.transform.ref === structureRef)?.components ?? [];
      const own = siblings.find((c) => c.cell.transform.ref === ref)?.cell.obj?.label ?? label;
      const sameLabel = siblings.filter((c) => (c.cell.obj?.label ?? '') === own).length;
      return { structure, component: own && sameLabel === 1 ? own : ref };
    };
    return [
      {
        id: 'add-repr',
        label: 'Add representation…',
        icon: 'plus',
        run: () => openCommand('representation.set', named()),
      },
      {
        id: 'coloring',
        label: 'Set coloring…',
        icon: 'palette',
        run: () =>
          openCommand('representation.color-theme', {
            structure: named().structure,
            ...(keyword ? { target: keyword } : {}),
          }),
      },
      // This component changed by the current selection: three choices, at once (no form to fill).
      modifyBySelection(api, ref, structureRef, noSelection),
      {
        id: 'select',
        label: 'Select this',
        icon: 'mouse-pointer-2',
        run: () => {
          const l = loci();
          if (l) void api.run('select.set', { target: { loci: l }, mode: 'set' });
        },
      },
    ];
  }
  const types = plugin.representation.structure.registry.types.map(([id, label]) => ({ id, label }));
  const themes = plugin.representation.structure.themes.colorThemeRegistry.types
    .slice(0, 24)
    .map(([id, label]) => ({ id, label }));
  return [
    {
      id: 'add-repr',
      label: 'Add representation',
      icon: 'plus',
      items: types.map((t) => ({
        id: t.id,
        label: t.label,
        run: () => void api.run('representation.set', { component: ref, structure: structureRef, type: t.id }),
      })),
    },
    {
      id: 'coloring',
      label: 'Set coloring',
      icon: 'palette',
      items: themes.map((t) => ({
        id: t.id,
        label: t.label,
        run: () => {
          const l = loci();
          if (l) void api.run('representation.color', { target: { loci: l }, theme: t.id });
        },
      })),
    },
    modifyBySelection(api, ref, structureRef, noSelection),
    {
      id: 'select',
      label: 'Select this',
      icon: 'mouse-pointer-2',
      run: () => {
        const l = loci();
        if (l) void api.run('select.set', { target: { loci: l }, mode: 'set' });
      },
    },
    // Remove comes from the tree (the row is removable), set apart in the danger color.
  ];
}

/** What a structure is made of (the Contents inspector): its counts and its entities with their chains and atoms. */
export function contentsSummary(s: Structure, reveal: (entity: string) => void): ContentsSummary {
  const { entities } = s.model;
  const byEntity = new Map<string, { chains: Set<string>; atoms: number }>();
  const chainKeys = new Set<string>();
  const loc = StructureElement.Location.create(s);
  for (const unit of s.units) {
    loc.unit = unit as Unit;
    loc.element = unit.elements[0];
    const entity = StructureProperties.chain.label_entity_id(loc);
    let e = byEntity.get(entity);
    if (!e) byEntity.set(entity, (e = { chains: new Set(), atoms: 0 }));
    e.chains.add(StructureProperties.chain.auth_asym_id(loc));
    e.atoms += unit.elements.length;
    chainKeys.add(`${unit.conformation.operator.name}:${StructureProperties.chain.label_asym_id(loc)}`);
  }
  const list: ContentsEntity[] = [];
  for (let i = 0; i < entities.data._rowCount; i++) {
    const id = entities.data.id.value(i);
    const e = byEntity.get(id);
    if (!e) continue;
    list.push({
      id,
      name: entities.data.pdbx_description.value(i).join(', ') || `Entity ${id}`,
      type: entities.data.type.value(i) || 'other',
      chains: [...e.chains],
      atoms: e.atoms,
    });
  }
  const n = (x: number) => x.toLocaleString('en');
  const counts = [
    { label: list.length === 1 ? 'Entity' : 'Entities', value: n(list.length) },
    { label: chainKeys.size === 1 ? 'Chain' : 'Chains', value: n(chainKeys.size) },
    { label: 'Residues', value: n(s.polymerResidueCount) },
    { label: 'Atoms', value: n(s.elementCount) },
  ];
  if (s.models.length > 1) counts.push({ label: 'Models', value: n(s.models.length) });
  return { counts, entities: list, reveal };
}

/** A structure's entities, each with its chains (lazy contents, decision 18: a row's eye hides it). */
function entities(
  plugin: PluginContext,
  s: Structure,
  structureRef: string,
  nodes: Map<string, Node>,
  api: ViewerApi,
): SceneNode[] {
  const shown = hiddenContents(plugin, structureRef);
  // As sets (an `includes` a chain was chains × hidden, at every sync, in a structure of thousands of chains).
  const hidden = { entities: new Set(shown.entities), chains: new Set(shown.chains) };
  const { entities } = s.model;
  const byEntity = new Map<string, Map<string, number>>();
  const loc = StructureElement.Location.create(s);
  for (const unit of s.units) {
    loc.unit = unit as Unit;
    loc.element = unit.elements[0];
    const entity = StructureProperties.chain.label_entity_id(loc);
    const chain = StructureProperties.chain.auth_asym_id(loc);
    let chains = byEntity.get(entity);
    if (!chains) byEntity.set(entity, (chains = new Map()));
    chains.set(chain, (chains.get(chain) ?? 0) + unit.elements.length);
  }
  const out: SceneNode[] = [];
  for (let i = 0; i < entities.data._rowCount; i++) {
    const id = entities.data.id.value(i);
    const chains = byEntity.get(id);
    if (!chains) continue;
    const key = `${structureRef}#entity:${id}`;
    let n = nodes.get(key);
    const label = entities.data.pdbx_description.value(i).join(', ') || `Entity ${id}`;
    if (!n) {
      n = new Node(key, 'entity', 'shapes', label, false, {
        visibility: (hide) =>
          void api.run('contents.visibility', { structure: structureRef, entity: id, visible: !hide }),
        highlight: (on) =>
          on
            ? void api.run('highlight.set', { target: { loci: contentsLoci(dataOf(plugin, structureRef) ?? s, id) } })
            : void api.run('highlight.clear'),
      });
      n.removable = false;
      nodes.set(key, n);
    }
    const entityHidden = hidden.entities.has(id);
    n.set(label, `${entities.data.type.value(i)} · ${chains.size} chain${chains.size === 1 ? '' : 's'}`, entityHidden);
    n.menu = contentsMenu(api, () => contentsLoci(s, id), label);
    const chainNodes: SceneNode[] = [];
    for (const [chain, atoms] of chains) {
      const ck = `${structureRef}#chain:${id}:${chain}`;
      let c = nodes.get(ck);
      if (!c) {
        c = new Node(ck, 'chain', 'sequence', `Chain ${chain}`, false, {
          // That entity's chain only (a ligand or waters with the same author chain stay).
          visibility: (hide) =>
            void api.run('contents.visibility', { structure: structureRef, entity: id, chain, visible: !hide }),
          highlight: (on) =>
            on
              ? void api.run('highlight.set', {
                  target: { loci: contentsLoci(dataOf(plugin, structureRef) ?? s, id, chain) },
                })
              : void api.run('highlight.clear'),
        });
        c.removable = false;
        nodes.set(ck, c);
      }
      // Hidden itself, or with its entity (its own eye then shows the entity's state).
      c.set(
        `Chain ${chain}`,
        `${atoms.toLocaleString('en')} atoms`,
        entityHidden || hidden.chains.has(`${id}:${chain}`) || hidden.chains.has(chain),
      );
      c.menu = contentsMenu(api, () => contentsLoci(s, id, chain), `${label} · chain ${chain}`);
      chainNodes.push(c);
    }
    n.setChildren(chainNodes);
    out.push(n);
  }
  return out;
}

/** "1 distance, 2 angles", or "none" (the group's description). */
const FocusPrefix = /^\[Focus\] /;

/** A component the focus behavior made (the target, the surroundings). */
function isFocusPart(c: StructureComponentRef): boolean {
  const tags = c.cell.transform.tags;
  return (
    !!tags &&
    (tags.includes(StructureFocusRepresentationTags.TargetSel) ||
      tags.includes(StructureFocusRepresentationTags.SurrSel))
  );
}

function focusParts(plugin: PluginContext, structureRef: string): StructureComponentRef[] {
  const s = plugin.managers.structure.hierarchy.current.refs.get(structureRef) as StructureRef | undefined;
  return s?.components?.filter(isFocusPart) ?? [];
}

function measurementCounts(list: readonly MeasuredCell[]): string {
  if (!list.length) return 'none';
  const counts = new Map<string, number>();
  for (const m of list) counts.set(m.kind, (counts.get(m.kind) ?? 0) + 1);
  return [...counts]
    .map(
      ([kind, n]) => `${n} ${MeasureKinds[kind as keyof typeof MeasureKinds].title.toLowerCase()}${n === 1 ? '' : 's'}`,
    )
    .join(', ');
}

/** Whether a structure has density to stream (v5's Volume Streaming action applies to it). */
function canStream(plugin: PluginContext, s: StructureRef): boolean {
  const cell = s.cell;
  return !!cell.obj && !!InitVolumeStreaming.definition.isApplicable?.(cell.obj, cell.transform, plugin);
}

/** A structure's data as it is now, by ref (a row's actions outlive the structure they were made with). */
function dataOf(plugin: PluginContext, ref: string): Structure | undefined {
  return plugin.state.data.cells.get(ref)?.obj?.data as Structure | undefined;
}

/** The atoms of an entity, or of one of its chains, in a structure. */
function contentsLoci(s: Structure, entity: string, chain?: string): StructureElement.Loci {
  const isEntity = MS.core.rel.eq([MS.ammp('label_entity_id'), entity]);
  const test = chain ? MS.core.logic.and([isEntity, MS.core.rel.eq([MS.ammp('auth_asym_id'), chain])]) : isEntity;
  const query = MS.struct.generator.atomGroups({ 'chain-test': test });
  return StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(query, s));
}

/** An entity's or a chain's row menu (parity.md §2: select, focus, a component of it). */
function contentsMenu(api: ViewerApi, loci: () => StructureElement.Loci, label: string): SceneMenuItem[] {
  return [
    {
      id: 'select',
      label: 'Select',
      icon: 'mouse-pointer-2',
      run: () => void api.run('select.set', { target: { loci: loci() }, mode: 'set' }),
    },
    { id: 'focus', label: 'Focus', icon: 'focus', run: () => void api.run('focus.set', { target: { loci: loci() } }) },
    {
      id: 'component',
      label: 'Make a component',
      icon: 'plus',
      run: () => void api.run('component.add', { target: { loci: loci() }, label }),
    },
  ];
}

/** A form's edits, applied a moment after the last (sliders change values at pointer rate). */
function debounced<T>(apply: (values: T) => void, ms = 150): (values: T) => void {
  let timer: ReturnType<typeof setTimeout> | undefined;
  return (values) => {
    clearTimeout(timer);
    timer = setTimeout(() => apply(values), ms);
  };
}

// ---- Volume streaming's inspector (parity.md §2; v5's Volume Streaming controls) ----

const AroundFocus = (v: Record<string, unknown>) => v.view !== 'selection-box' && v.view !== 'auto';

function streamingForm(on: StreamingState): PD.Params {
  const levels = on.entry.header.availablePrecisions.map((_, i) => [i + 1, `${i + 1}`] as [number, string]);
  const params: PD.Params = {
    view: PD.Select<string>('selection-box', VolumeStreaming.ViewTypeOptions, {
      description: 'What of the map is streamed: around the focus, around the camera, the whole structure, …',
    }),
    radius: PD.Numeric(5, { min: 0.5, max: 50, step: 0.5 }, { label: 'Radius (Å)', hideIf: AroundFocus }),
    cameraRadius: PD.Numeric(
      0.5,
      { min: 0.05, max: 1, step: 0.05 },
      { label: 'Radius (of the view)', hideIf: (v) => v.view !== 'camera-target' },
    ),
    detail: PD.Select<number>(1, levels, { label: 'Detail', hideIf: (v) => v.view === 'off' }),
  };
  return params;
}

function streamingValues(on: StreamingState): PD.Values<PD.Params> {
  const p = on.params.entry.params;
  const view = p.view.params as { radius?: number; dynamicDetailLevel?: number };
  const camera = p.view.name === 'camera-target';
  const values: Record<string, unknown> = {
    view: p.view.name,
    radius: !camera && view.radius !== undefined ? view.radius : 5,
    cameraRadius: camera && view.radius !== undefined ? view.radius : 0.5,
    detail: (camera ? (view.dynamicDetailLevel ?? 0) : p.detailLevel) + 1,
  };
  return values;
}

/**
 * The form's edits as `volume.streaming-view` calls: the fields changed in a moment, as one call (sliders). The
 * channels' levels are the inspector's level controls (`streamingLevels`).
 */
function streamingChanges(api: ViewerApi, structure: string): (values: PD.Values<PD.Params>, name: string) => void {
  let timer: ReturnType<typeof setTimeout> | undefined;
  let pending: { view?: string; radius?: number; detail?: number } = {};
  return (values, name) => {
    const v = values as Record<string, any>;
    if (name === 'view') pending.view = v.view;
    else if (name === 'radius') pending.radius = v.radius;
    else if (name === 'cameraRadius') pending.radius = v.cameraRadius;
    else if (name === 'detail') pending.detail = v.detail;
    clearTimeout(timer);
    timer = setTimeout(() => {
      const change = pending;
      pending = {};
      void api.run('volume.streaming-view', { structure, ...change });
    }, 150);
  };
}

// parity.md §2, the Conformation Transform: v5's form (translation, then a rotation about the structure's center) under
// the structure's facts, folded; a superposition's matrix shows as the same parts.
const TransformForm = {
  transform: PD.Group(
    {
      translation: PD.Vec3(Vec3.create(0, 0, 0), {}, { label: 'Translation', description: 'In Å' }),
      axis: PD.Vec3(Vec3.create(0, 0, 1), {}, { label: 'Rotation axis' }),
      angle: PD.Numeric(0, { min: -360, max: 360, step: 1 }, { label: 'Angle', description: 'In degrees' }),
    },
    { label: 'Transform', isExpanded: false, description: 'Move and rotate the structure (about its center)' },
  ),
};

/** A rigid transform as the form's parts: a rotation about `data`'s center, then a translation. */
function transformValues(m: Mat4, data: Structure): PD.Values<typeof TransformForm>['transform'] {
  const q = Mat4.getRotation(Quat(), m);
  const axis = Vec3();
  let angle = (Quat.getAxisAngle(axis, q) * 180) / Math.PI;
  if (angle > 180) angle -= 360;
  if (Math.abs(angle) < 1e-4) {
    angle = 0;
    Vec3.set(axis, 0, 0, 1);
  }
  // M = T(t)·T(c)·R·T(−c), so t = M's translation + R·c − c.
  const c = data.boundary.sphere.center;
  const t = Vec3.add(Vec3(), Mat4.getTranslation(Vec3(), m), Vec3.sub(Vec3(), Vec3.transformQuat(Vec3(), c, q), c));
  const round = (x: number) => Math.round(x * 1000) / 1000 || 0;
  return {
    translation: Vec3.create(round(t[0]), round(t[1]), round(t[2])),
    axis: Vec3.create(round(axis[0]), round(axis[1]), round(axis[2])),
    angle: round(angle),
  };
}

// parity.md §2, the structure's source and crystal symmetry: v5's Source panel (the model of a trajectory, the structure
// type: the model, an assembly, symmetry mates or unit cells) as fields over `structure.model`, `.assembly` and
// `.symmetry`, above the Transform. A field with one choice is left out.

interface SourceValues {
  type: string;
  radius: number;
  cells: number;
  bonds: BondProvider.Choice;
}

type StructureType = {
  name: string;
  params?: { id?: string; radius?: number; ijkMax?: Vec3; bonds?: BondProvider.Choice };
};

const DefaultBondChoice: BondProvider.Choice = { name: BondProvider.DefaultName, params: {} };

const typeOf = (s: StructureRef) => (s.cell.transform.params as { type: StructureType } | undefined)?.type;

/** What the structure can be built as: the model, its assemblies, and for a crystal its mates or cells. */
function structureTypes(s: StructureRef, data: Structure): [string, string][] {
  const symmetry = ModelSymmetry.Provider.get(data.model);
  const assemblies = symmetry?.assemblies ?? [];
  const crystal = !!symmetry && !SpacegroupCell.isZero(symmetry.spacegroup.cell);
  const options: [string, string][] = [['model', 'Model (as deposited)']];
  for (const a of assemblies) {
    options.push([`assembly:${a.id}`, `Assembly ${a.id}${a.details ? ` · ${a.details.replace(/_/g, ' ')}` : ''}`]);
  }
  if (crystal) options.push(['mates', 'Symmetry mates'], ['cells', 'Unit cells']);
  // Types no field here sets (the loader's auto, a symmetry assembly from `state.update`): shown while they hold.
  const now = typeOf(s)?.name;
  if (now === 'auto') options.unshift(['auto', assemblies.length ? 'Auto (the first assembly)' : 'Auto (the model)']);
  if (now === 'symmetry-assembly') options.push(['symmetry-assembly', 'Symmetry assembly']);
  return options;
}

function typeValue(s: StructureRef, data: Structure): string {
  const t = typeOf(s);
  switch (t?.name) {
    case 'assembly': {
      const id = t.params?.id || ModelSymmetry.Provider.get(data.model)?.assemblies[0]?.id;
      return id ? `assembly:${id}` : 'model';
    }
    case 'symmetry-mates':
      return 'mates';
    case 'symmetry':
      return 'cells';
    case 'auto':
    case 'symmetry-assembly':
      return t.name;
    default:
      return 'model';
  }
}

function sourceForm(s: StructureRef, data: Structure, bondProviders: BondProviderRegistry): PD.Params {
  const types = structureTypes(s, data);
  // bond-providers plan §4.1: hidden with the built-in bonds as the only choice, so such a viewer looks as before
  const model = s.model?.cell.obj?.data;
  const bondOptions = bondProviders.options(model).length;
  // The model shown is not here: a structure with several has a row of its own under it (its frames, a scrubber and
  // play; views/structure/trajectory-row.tsx), its controls in the inspector when picked.
  return {
    type: PD.Select<string>('model', types, {
      label: 'Type',
      description: 'What is built from the model: itself, an assembly, or its crystal symmetry',
      isHidden: types.length < 2,
    }),
    radius: PD.Numeric(
      5,
      { min: 0, max: 50, step: 1 },
      { label: 'Radius (Å)', description: 'The mates within this distance', hideIf: (v) => v.type !== 'mates' },
    ),
    cells: PD.Numeric(
      1,
      { min: 0, max: 3, step: 1 },
      { label: 'Cells around', description: '1 is the cells −1 to 1 on each axis', hideIf: (v) => v.type !== 'cells' },
    ),
    bonds: bondProviders.getParam(model, { isHidden: bondOptions < 2 || typeOf(s) === undefined }),
  };
}

/** The structure's form values: its source (model, type, radius, cells) and its transform. */
function structureValues(s: StructureRef, data: Structure): SourceValues & PD.Values<typeof TransformForm> {
  const t = typeOf(s);
  return {
    type: typeValue(s, data),
    radius: t?.name === 'symmetry-mates' ? (t.params?.radius ?? 5) : 5,
    cells: t?.name === 'symmetry' && t.params?.ijkMax ? t.params.ijkMax[0] : 1,
    bonds: t?.params?.bonds ?? DefaultBondChoice,
    transform: transformValues(currentTransform(s, data), data),
  };
}

/**
 * The source fields' edits as operations: a model or a type at once, the radius and cells a moment after the last
 * (sliders). `reload` puts the structure's values back after each (a refused change leaves the form as it was).
 */
function sourceChanges(
  api: ViewerApi,
  s: StructureRef,
  reload: () => void,
): (values: SourceValues, name: string) => void {
  const structure = s.cell.transform.ref;
  let timer: ReturnType<typeof setTimeout> | undefined;
  const apply = async (v: SourceValues, name: string) => {
    if (name === 'bonds') {
      return api.run('structure.bonds', { structure, provider: v.bonds.name, params: v.bonds.params });
    }
    if (v.type.startsWith('assembly:')) {
      return api.run('structure.assembly', { structure, id: v.type.slice('assembly:'.length) });
    }
    if (v.type === 'mates') return api.run('structure.symmetry', { structure, kind: 'mates', radius: v.radius });
    if (v.type === 'cells') return api.run('structure.symmetry', { structure, kind: 'cells', cells: v.cells });
    // The model: from an assembly through `assembly` (its own words), from symmetry through `symmetry` off.
    if (v.type === 'model') {
      const now = typeOf(s)?.name;
      return now === 'assembly' || now === 'auto'
        ? api.run('structure.assembly', { structure, id: '' })
        : api.run('structure.symmetry', { structure, kind: 'off' });
    }
    return undefined;
  };
  return (values, name) => {
    clearTimeout(timer);
    const run = () => void apply(values, name).finally(reload);
    // a provider's own params (sliders) a moment after the last; a change of provider at once
    const bondParams = name === 'bonds' && values.bonds.name === (typeOf(s)?.params?.bonds ?? DefaultBondChoice).name;
    if (name === 'radius' || name === 'cells' || bondParams) timer = setTimeout(run, 150);
    else run();
  };
}

/**
 * A structure representation's row menu: its type and color (what people change most; the current ones checked), then
 * the whole form in the inspector. Remove follows from the tree.
 */
function representationMenu(plugin: PluginContext, api: ViewerApi, ref: string, edit: () => void): SceneMenuItem[] {
  const choices = representationChoices(plugin, ref);
  // One change, one undo step: through the quick settings' own path (the params of another type keep the opacity).
  const change = (c: RepresentationQuickChange) => {
    const quick = createRepresentationQuick(plugin, api, ref);
    quick?.change(c);
    quick?.settle();
  };
  return [
    ...(choices
      ? [
          {
            id: 'type',
            label: 'Type',
            icon: 'shapes' as const,
            items: choices.types.map((t) => ({
              id: t.id,
              label: t.label,
              checked: t.id === choices.type,
              run: () => change({ type: t.id }),
            })),
          },
          {
            id: 'color',
            label: 'Color',
            icon: 'palette' as const,
            items: choices.themes.map((t) => ({
              id: t.id,
              label: t.label,
              checked: t.id === choices.theme,
              run: () => change(t.id === 'uniform' ? { theme: 'uniform' } : { theme: t.id, color: undefined }),
            })),
          },
        ]
      : []),
    { id: 'edit', label: 'Edit…', icon: 'sliders-horizontal' as const, run: edit },
  ];
}

/** The structure row's menu: presets, assemblies, crystal symmetry, density streaming, and the transform's reset. */
function structureMenu(plugin: PluginContext, api: ViewerApi, s: StructureRef): SceneMenuItem[] {
  const structure = s.cell.transform.ref;
  const data = s.cell.obj?.data;
  const items: SceneMenuItem[] = [];
  const presets = plugin.builders.structure.representation.providers.filter(
    (p) => !p.isApplicable || !s.cell.obj || p.isApplicable(s.cell.obj, plugin),
  );
  if (presets.length) {
    items.push({
      id: 'preset',
      label: 'Preset',
      icon: 'shapes',
      items: presets.map((p) => ({
        id: p.id,
        label: p.display.name,
        run: () => void api.run('structure.preset', { structure, preset: p.id }),
      })),
    });
  }
  if (data) {
    const types = structureTypes(s, data).filter(([id]) => id === 'model' || id.startsWith('assembly:'));
    if (types.length > 1) {
      items.push({
        id: 'assembly',
        label: 'Assembly',
        icon: 'copy',
        items: types.map(([id, label]) => ({
          id,
          label,
          run: () =>
            void api.run('structure.assembly', { structure, id: id === 'model' ? '' : id.slice('assembly:'.length) }),
        })),
      });
    }
    const symmetry = ModelSymmetry.Provider.get(data.model);
    if (symmetry && !SpacegroupCell.isZero(symmetry.spacegroup.cell)) {
      items.push({
        id: 'symmetry',
        label: 'Crystal symmetry',
        icon: 'box',
        items: (
          [
            ['mates', 'Symmetry mates'],
            ['cells', 'Unit cells'],
            ['off', 'Off'],
          ] as const
        ).map(([kind, label]) => ({
          id: kind,
          label,
          run: () => void api.run('structure.symmetry', { structure, kind }),
        })),
      });
    }
  }
  // parity.md §2, volume streaming: turned on from here (its row, once on, turns it off).
  if (api.describe('volume.streaming') && !streamingOf(plugin, s)) {
    items.push({
      id: 'stream-density',
      label: 'Stream density',
      icon: 'volume',
      run: () => void api.run('volume.streaming', { structure, on: true }),
    });
  }
  // What an operation offers for a structure (its `structureMenu`: model-archive's PAE fetch for an AlphaFold DB entry).
  for (const def of api.definitions.value.values()) {
    const offer = def.structureMenu;
    if (!offer || api.describe(def.id)?.supported !== true || !offer.applies(plugin, structure)) continue;
    items.push({
      id: def.id,
      label: offer.label,
      ...(offer.icon ? { icon: offer.icon } : {}),
      run: () => void api.run(def.id, { structure }),
    });
  }
  // parity.md §2: the Conformation Transform's reset (moving is the inspector's Transform).
  items.push({
    id: 'reset-transform',
    label: 'Reset transform',
    icon: 'rotate-ccw',
    run: () => void api.run('structure.transform-reset', { structure }),
  });
  return items;
}

// parity.md §2, the unit cell: its color and opacity (the operation's; the line size stays the representation's).
const UnitCellForm = {
  color: PD.Color(Color(0xffa500), { label: 'Color', isEssential: true }),
  opacity: PD.Numeric(1, { min: 0.05, max: 1, step: 0.05 }, { label: 'Opacity', isEssential: true }),
};

function unitCellValues(cell: StateObjectCell): PD.Values<typeof UnitCellForm> {
  const p = cell.transform.params as { cellColor?: Color; alpha?: number } | undefined;
  return { color: p?.cellColor ?? Color(0xffa500), opacity: p?.alpha ?? 1 };
}

/** What the inspector shows for a node: facts, and the parameters of a representation as a form. */
function inspect(
  plugin: PluginContext,
  api: ViewerApi,
  node: SceneNode,
  owner: string | undefined,
  unit: ObservableState<VolumeLevelUnit>,
  pick?: (node: SceneNode) => void,
): SceneInspection {
  const cell = plugin.state.data.cells.get(node.id);
  const slot = `inspector:${node.kind}`;
  if (node.kind === 'unit-cell' && cell && isParticleList(plugin.state.data.cells.get(cell.transform.parent)?.obj)) {
    const list = plugin.state.data.cells.get(cell.transform.parent)!.obj!.data as ParticleList;
    const form = new ParamsFormModel(UnitCellForm, unitCellValues(cell), {
      onChange: debounced((values) => {
        const now = plugin.state.data.cells.get(node.id);
        void api.run('particles.unit-cell', {
          particles: cell.transform.parent,
          visible: !now?.state.isHidden,
          color: Color.toHexStyle(values.color),
          opacity: values.opacity,
        });
      }),
    });
    return { node, form, slot, facts: list.cell ? cellFacts(list.cell) : [] };
  }
  if (node.kind === 'unit-cell' && cell) {
    const parent = plugin.state.data.cells.get(cell.transform.parent)?.obj as PluginStateObject.Molecule.Model;
    const symmetry = parent?.data && ModelSymmetry.Provider.get(parent.data as Model);
    const facts: [string, string][] = [];
    if (symmetry) {
      const { size, anglesInRadians } = symmetry.spacegroup.cell;
      const deg = (r: number) => `${Math.round(((r * 180) / Math.PI) * 100) / 100}°`;
      facts.push(['Space group', symmetry.spacegroup.name]);
      facts.push(['Cell', `${[...size].map((x) => Math.round(x * 100) / 100).join(' × ')} Å`]);
      facts.push(['Angles', [...anglesInRadians].map(deg).join(', ')]);
    }
    const form = new ParamsFormModel(UnitCellForm, unitCellValues(cell), {
      onChange: debounced((values) => {
        const now = plugin.state.data.cells.get(node.id);
        void api.run('structure.unit-cell', {
          structure: owner,
          visible: !now?.state.isHidden,
          color: Color.toHexStyle(values.color),
          opacity: values.opacity,
        });
      }),
    });
    return { node, form, slot, facts };
  }
  if (node.kind === 'extension' && owner) {
    const s = plugin.managers.structure.hierarchy.current.refs.get(owner);
    const on = s?.kind === 'structure' ? streamingOf(plugin, s) : undefined;
    if (on) {
      const state = () => {
        const now = plugin.managers.structure.hierarchy.current.refs.get(owner);
        return now?.kind === 'structure' ? streamingOf(plugin, now) : undefined;
      };
      return {
        node,
        form: new ParamsFormModel(streamingForm(on), streamingValues(on), { onChange: streamingChanges(api, owner) }),
        levels: streamingLevels(plugin, owner, state).map((port) => createVolumeLevels(api, port, unit)),
        slot,
        facts: [
          ['Entry', `${on.entry.dataId} (${on.entry.kind === 'em' ? 'EM' : 'X-ray'})`],
          ['Detail levels', String(on.entry.header.availablePrecisions.length)],
        ],
      };
    }
  }
  if (node.kind === 'structure' && cell?.obj?.data) {
    const s = cell.obj.data as Structure;
    const model = s.model;
    const chains = new Set(s.units.map((u) => u.chainGroupId)).size;
    const ref = plugin.managers.structure.hierarchy.current.refs.get(node.id);
    const structureRef = ref?.kind === 'structure' ? ref : undefined;
    const now = structureRef ? currentTransform(structureRef, s) : Mat4.identity();
    const transform = debounced(
      ({ transform }: PD.Values<typeof TransformForm>) =>
        void api.run('structure.transform', {
          structure: node.id,
          translation: transform.translation,
          axis: transform.axis,
          angle: transform.angle,
          center: 'centroid',
        }),
    );
    const params = { ...(structureRef ? sourceForm(structureRef, s, plugin.bondProviders) : {}), ...TransformForm };
    const values = structureRef ? structureValues(structureRef, s) : { transform: transformValues(now, s) };
    // After a source change that did nothing (a failure), the form shows what the structure is again.
    const reload = () => {
      const fresh = plugin.managers.structure.hierarchy.current.refs.get(node.id);
      const data = fresh?.kind === 'structure' ? fresh.cell.obj?.data : undefined;
      if (fresh?.kind === 'structure' && data) form.load(structureValues(fresh, data) as never);
    };
    const source = structureRef ? sourceChanges(api, structureRef, reload) : undefined;
    const form: ParamsFormModel<PD.Params> = new ParamsFormModel(params as PD.Params, values as never, {
      onChange: (v, name) => {
        if (name === 'transform') transform(v as PD.Values<typeof TransformForm>);
        else source?.(v as SourceValues, name);
      },
    });
    return {
      node,
      form,
      slot,
      facts: [
        ['Entry', model.entryId],
        ['Atoms', s.elementCount.toLocaleString('en')],
        ['Residues', s.polymerResidueCount.toLocaleString('en')],
        ['Chains', String(chains)],
        ['Model', model.label],
      ],
    };
  }
  if (node.kind === 'volume' && cell?.obj?.data) {
    const data = cell.obj.data as Volume;
    const stats = volumeStats(data);
    const signed = (x: number) => shortNumber(x).replace('-', '−');
    const levels = createVolumeLevels(
      api,
      mapLevels(plugin, () => volumeRefOf(plugin, node.id)),
      unit,
    );
    // Its direct volumes' transfer functions (named when there are several).
    const direct = (volumeRefOf(plugin, node.id)?.representations ?? []).filter(
      (r) => (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name === 'direct-volume',
    );
    const transfer = direct.flatMap((r, i) => {
      const t = createVolumeTransfer(plugin, api, () => volumeRefOf(plugin, node.id), r.cell.transform.ref, unit, {
        ...(direct.length > 1 ? { title: `Direct volume ${i + 1}` } : {}),
      });
      return t ? [t] : [];
    });
    return {
      node,
      form: undefined,
      slot,
      ...(levels.levels.value.length ? { levels: [levels] } : {}),
      ...(transfer.length ? { transfer } : {}),
      facts: [
        ...(data.sourceData.kind ? [['Source', data.sourceData.kind.toUpperCase()] as const] : []),
        ['Grid', gridText(stats)],
        ['Voxel', voxelText(stats)],
        ['Mean · σ', `${signed(stats.mean)} · ${signed(stats.sigma)}`],
        ['Range', `${signed(stats.min)} … ${signed(stats.max)}`],
      ],
    };
  }
  if (node.kind === 'contents') {
    // A structure's contents (its id is the structure's ref + '#contents'): what it is made of, at a glance; an entity
    // opens in the tree (its rows computed on first look).
    const structureRef = node.id.slice(0, node.id.indexOf('#'));
    const s = dataOf(plugin, structureRef);
    const reveal = (entity: string) => {
      node.load?.();
      if (!node.isExpanded.value) node.toggleExpanded();
      const row = node.children.value.find((c) => c.id === `${structureRef}#entity:${entity}`);
      if (row) pick?.(row);
    };
    if (s) return { node, form: undefined, slot, contents: contentsSummary(s, reveal) };
  }
  if (node.kind === 'particles' && cell?.obj?.data) {
    const data = particleListData(plugin.state.data, node.id) ?? (cell.obj.data as ParticleList);
    // How it is read from its file (the format's transform: tomograms, pixel size, radius, a frame): `particles.source`
    const parent = plugin.state.data.cells.get(cell.transform.parent);
    const def = cell.transform.transformer.definition.params?.(parent?.obj as never, plugin) as PD.Params | undefined;
    const form =
      def && Object.keys(def).length
        ? new ParamsFormModel(def, cell.transform.params as PD.Values<PD.Params>, {
            onChange: debounced((values) => void api.run('particles.source', { particles: node.id, params: values })),
          })
        : undefined;
    return { node, form, slot, facts: particlesFacts(data) };
  }
  const updateOperation = UpdateOperations[node.kind];
  if (updateOperation && cell) {
    const parent = plugin.state.data.cells.get(cell.transform.parent);
    const def = cell.transform.transformer.definition.params?.(parent?.obj as never, plugin) as PD.Params | undefined;
    if (def) {
      let timer: ReturnType<typeof setTimeout> | undefined;
      const form = new ParamsFormModel(def, cell.transform.params as PD.Values<PD.Params>, {
        // Sliders change values at pointer rate: one update a moment after the last.
        onChange: (values) => {
          clearTimeout(timer);
          timer = setTimeout(() => void api.run(updateOperation, { node: node.id, params: values }), 150);
        },
      });
      // An isosurface: its level above the whole form (the iso value, the color, the opacity at hand).
      const volume =
        node.kind === 'volume-representation'
          ? plugin.managers.volume.hierarchy.current.volumes.find((v) =>
              v.representations.some((r) => r.cell.transform.ref === node.id),
            )?.cell.transform.ref
          : undefined;
      const levels = volume
        ? createVolumeLevels(
            api,
            mapLevels(plugin, () => volumeRefOf(plugin, volume), node.id),
            unit,
          )
        : undefined;
      // A direct volume: its transfer function above the whole form.
      const direct = (cell.transform.params as { type?: { name: string } } | undefined)?.type?.name === 'direct-volume';
      const transfer =
        volume && direct
          ? createVolumeTransfer(plugin, api, () => volumeRefOf(plugin, volume), node.id, unit)
          : undefined;
      // A structure representation: its quick settings above the whole form.
      const quick = node.kind === 'representation' ? createRepresentationQuick(plugin, api, node.id) : undefined;
      return {
        node,
        form,
        slot,
        ...(levels?.levels.value.length ? { levels: [levels] } : {}),
        ...(transfer ? { transfer: [transfer] } : {}),
        ...(quick ? { quick } : {}),
      };
    }
  }
  // A component: its name editable in the inspector's title (the mocks', v5's label).
  if (node.kind === 'component') {
    return {
      node,
      form: undefined,
      slot,
      rename: (label) => void api.run('component.rename', { component: node.id, label }),
    };
  }
  return { node, form: undefined, slot };
}

// Particles (the lazy particles feature's operations, by id: this module imports none of its code).

const ParticleListUnitcell = 'ms-plugin.particle-list-unitcell-3d';

/** A particle list's unit cell cell (`ParticleListUnitcell3D` below it), if any */
function particleUnitCellOf(state: State, listRef: string): StateObjectCell | undefined {
  return state.select(
    StateSelection.Generators.byRef(listRef)
      .subtree()
      .filter((c) => c.transform.transformer.id === ParticleListUnitcell),
  )[0];
}

function isParticleList(obj: unknown): obj is PluginStateObject.Particle.List {
  return PluginStateObject.Particle.List.is(obj as never);
}

/** A list's row description: its count, and its targets when it has some ("1,200 particles · 3 targets"). */
function particlesText(data: ParticleList): string {
  const targets = data.targetInfo.size;
  return `${data.count.toLocaleString('en')} particles${targets > 1 ? ` · ${targets} targets` : ''}`;
}

/** What stands for a target: a structure's entry, a shape's or a volume's name */
function targetLabel(
  target: NonNullable<ParticleList['targetMapping']> extends ReadonlyMap<number, infer T> ? T : never,
) {
  switch (target.kind) {
    case 'structure':
      return target.structure.model.entryId || target.structure.label;
    case 'shape':
      return target.shape.name;
    case 'volume':
      return target.volume.sourceData.name || 'a volume';
  }
}

/** A list as its representations see it: the data at the end of its decorators (targets attached) */
function particleListData(state: State, ref: string): ParticleList | undefined {
  const root = StateTree.getDecoratorRoot(state.tree, ref);
  return (state.cells.get(root)?.obj?.data ?? state.cells.get(ref)?.obj?.data) as ParticleList | undefined;
}

const listOf = (names: readonly string[], max = 3) =>
  names.length <= max ? names.join(', ') : `${names.slice(0, max).join(', ')} and ${names.length - max} more`;

/** What a particle list is, at a glance: count, source, poses, targets, entities, compartments, attributes, cell. */
function particlesFacts(data: ParticleList): [string, string][] {
  const facts: [string, string][] = [['Particles', data.count.toLocaleString('en')]];
  if (data.sourceData?.kind) facts.push(['Source', data.sourceData.kind]);
  facts.push(['Orientations', data.rotations ? 'yes' : 'none']);
  // Each target by name (its entity's), and what stands for it once one is attached ("ribosome: 1crn")
  const targets = [...data.targetInfo].map(([id, t]) => {
    const name =
      t.name ?? (t.entity !== undefined ? data.entityInfo?.get(t.entity)?.name : undefined) ?? `target ${id}`;
    const object = data.targetMapping?.get(id);
    return object ? `${name}: ${targetLabel(object)}` : name;
  });
  if (targets.length) facts.push(['Targets', `${targets.length}: ${listOf(targets)}`]);
  if (data.entityInfo?.size) {
    facts.push(['Entities', `${data.entityInfo.size}: ${listOf([...data.entityInfo.values()].map((e) => e.name))}`]);
  }
  if (data.compartmentInfo?.size) facts.push(['Compartments', String(data.compartmentInfo.size)]);
  if (data.fibers?.count) facts.push(['Fibers', data.fibers.count.toLocaleString('en')]);
  for (const a of data.attributes?.values() ?? [])
    facts.push([a.label, `${shortNumber(a.min)} … ${shortNumber(a.max)}`]);
  if (data.cell) facts.push(...cellFacts(data.cell));
  return facts;
}

function cellFacts(cell: { size: ArrayLike<number>; anglesInRadians: ArrayLike<number> }): [string, string][] {
  const deg = (r: number) => `${Math.round(((r * 180) / Math.PI) * 100) / 100}°`;
  return [
    ['Cell', `${[...Array.from(cell.size)].map((x) => Math.round(x * 100) / 100).join(' × ')} Å`],
    ['Angles', Array.from(cell.anglesInRadians).map(deg).join(', ')],
  ];
}

/** A particle list's row menu: a representation added, all colored by a theme, the camera to it, its unit cell. */
function particlesMenu(
  plugin: PluginContext,
  api: ViewerApi,
  ref: string,
  openCommand: OpenCommand | undefined,
): SceneMenuItem[] | undefined {
  const l = plugin.managers.particles.hierarchy.current.lists.find((x) => x.cell.transform.ref === ref);
  const data = l && particleListData(plugin.state.data, ref);
  if (!l || !data) return undefined;
  const { registry, themes } = plugin.representation.particles;
  const present = new Set(
    l.representations.map((r) => (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name),
  );
  const types = registry.getApplicableList(data).filter((t) => !present.has(t.name));
  const colors = themes.colorThemeRegistry.getApplicableList({ particles: data });
  const items: SceneMenuItem[] = [];
  if (types.length) {
    items.push({
      id: 'add',
      label: 'Add representation',
      icon: 'plus',
      items: types.map((t) => ({
        id: t.name,
        label: t.provider.label,
        run: () => void api.run('particles.representation', { particles: ref, type: t.name }),
      })),
    });
  }
  if (l.representations.length && colors.length) {
    items.push({
      id: 'color',
      label: 'Color by',
      icon: 'palette',
      items: colors.map((t) => ({
        id: t.name,
        label: t.provider.label,
        run: () => void api.run('particles.color', { particles: ref, theme: t.name }),
      })),
    });
  }
  items.push({
    id: 'focus',
    label: 'Focus',
    icon: 'focus',
    run: () => void api.run('particles.focus', { particles: ref }),
  });
  if (data.cell) {
    const cell = particleUnitCellOf(plugin.state.data, ref);
    const shown = !!cell && !cell.state.isHidden;
    items.push({
      id: 'unit-cell',
      label: shown ? 'Hide unit cell' : 'Show unit cell',
      icon: 'box',
      run: () => void api.run('particles.unit-cell', { particles: ref, visible: !shown }),
    });
  }
  return items;
}

/** A particle representation's row menu: its color (the current one checked), then the whole form in the inspector. */
function particlesRepresentationMenu(
  plugin: PluginContext,
  api: ViewerApi,
  ref: string,
  edit: () => void,
): SceneMenuItem[] {
  const cell = plugin.state.data.cells.get(ref);
  const list = plugin.managers.particles.hierarchy.current.lists.find((l) =>
    l.representations.some((r) => r.cell.transform.ref === ref),
  );
  const data = list && particleListData(plugin.state.data, list.cell.transform.ref);
  const current = (cell?.transform.params as { colorTheme?: { name: string } } | undefined)?.colorTheme?.name;
  const colors = data
    ? plugin.representation.particles.themes.colorThemeRegistry.getApplicableList({ particles: data })
    : [];
  return [
    ...(colors.length
      ? [
          {
            id: 'color',
            label: 'Color',
            icon: 'palette' as const,
            items: colors.map((t) => ({
              id: t.name,
              label: t.provider.label,
              checked: t.name === current,
              run: () => void api.run('particles.color', { node: ref, theme: t.name }),
            })),
          },
        ]
      : []),
    { id: 'edit', label: 'Edit…', icon: 'sliders-horizontal' as const, run: edit },
  ];
}

const ParticleListFromTrajectory = 'ms-plugin.particle-list-from-trajectory';
const ParticleTrajectoryAnimation = 'built-in.animate-particle-trajectory';

/** A list read from a trajectory's frame: its frame count and the frame shown; none with one frame */
function particleFrames(plugin: PluginContext, cell: StateObjectCell): { count: number; index: number } | undefined {
  if (cell.transform.transformer.id !== ParticleListFromTrajectory) return undefined;
  const parent = plugin.state.data.cells.get(cell.transform.parent)?.obj;
  const count = (parent?.data as { frameCount?: number } | undefined)?.frameCount ?? 0;
  if (count < 2) return undefined;
  return { count, index: (cell.transform.params as { frameIndex?: number }).frameIndex ?? 0 };
}

function isPlayingParticleFrames(plugin: PluginContext, ref: string): boolean {
  const animation = plugin.managers.animation;
  if (!animation.isAnimating || animation.current?.anim.name !== ParticleTrajectoryAnimation) return false;
  const refs = (animation.current.paramValues as { refs?: readonly string[] }).refs;
  return !refs?.length || refs.includes(ref);
}

interface ParticleSceneTrajectory extends SceneTrajectory {
  refresh(frames: { count: number; index: number }): void;
  refreshPlaying(): void;
}

/** A particle trajectory's frames as the tree's frame row shows them, acting through `particles.frame`, `.play`, `.stop` */
function particleTrajectory(plugin: PluginContext, api: ViewerApi, ref: string): ParticleSceneTrajectory {
  const count = observableState(0);
  const index = observableState(0);
  const playing = observableState(false);
  const mode = observableState<'loop' | 'palindrome' | 'once'>('loop');
  const seconds = observableState(5);
  /** The frame to show next (not sent yet), and the one being shown (sent, not done): the latest seek wins */
  let wanted: number | undefined;
  let sending: number | undefined;
  const flush = async () => {
    if (sending !== undefined) return;
    try {
      while (wanted !== undefined) {
        sending = wanted;
        wanted = undefined;
        await api.run('particles.frame', { particles: ref, index: sending });
      }
    } finally {
      sending = undefined;
    }
  };
  return {
    unit: 'frames',
    count: count.asReadonly,
    index: index.asReadonly,
    playing: playing.asReadonly,
    mode,
    seconds,
    seek(i) {
      wanted = i;
      void flush();
    },
    step(by) {
      const c = count.value;
      if (c < 2) return;
      wanted = ((((wanted ?? sending ?? index.value) + by) % c) + c) % c;
      void flush();
    },
    togglePlay() {
      void (playing.value
        ? api.run('particles.stop')
        : api.run('particles.play', { particles: ref, mode: mode.value, seconds: seconds.value }));
    },
    // The canvas' frame bar steps back while the row is on screen, as for a structure's (the trajectory source)
    markShown(on) {
      trajectorySourceOf(plugin).markRowShown(ref, on);
    },
    refresh(frames) {
      if (count.value !== frames.count) count.set(frames.count);
      if (index.value !== frames.index) index.set(frames.index);
      this.refreshPlaying();
    },
    refreshPlaying() {
      const now = isPlayingParticleFrames(plugin, ref);
      if (playing.value !== now) playing.set(now);
    },
  };
}
