/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { StateObjectCell } from '@molstar/core/state';
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
import type { PluginStateObject } from '@molstar/plugin/state/objects';
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
import type { SceneInspection, SceneMenuItem, SceneModel, SceneNode, SceneNodeKind } from '../model/scene.js';
import { namesOf } from '../api/scope.js';
import { currentTransform } from '../api/structure.js';
import { VolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/behavior';
import {
  CreateVolumeStreamingBehavior,
  CreateVolumeStreamingInfo,
} from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import {
  gridText,
  type StreamingChannelsChange,
  type StreamingState,
  shortNumber,
  streamingChannels,
  streamingOf,
  streamingSummary,
  voxelText,
  volumeStats,
} from '../api/volume.js';
import { type MeasuredCell, MeasureKinds, readMeasurements } from './measurements.js';
import { structureNames } from './naming.js';
import type { SourceScope } from './state.js';

// The scene as the Structure tab shows it (state.md §1, row 2; view-models.md §3): one record per node, keyed by state
// ref, so a node keeps its identity across hierarchy rebuilds, and each field its own state, so a row updates only where
// it changed. A structure's contents (entities, then chains) are computed when first expanded. Every action goes through
// the API, bound to the person (README §3.3): the tree changes state only by operations.

const ReprIcons: Readonly<Record<string, IconName>> = {
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

class Node implements SceneNode {
  readonly #label: ObservableState<string>;
  readonly #description = observableState<string | undefined>(undefined);
  readonly #hidden = observableState(false);
  readonly #expanded: ObservableState<boolean>;
  readonly #children = observableState<readonly SceneNode[]>([]);
  icon: IconName;
  menu?: readonly SceneMenuItem[];
  removable?: boolean;
  lazy?: boolean;
  load?: () => void;

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

  set(label: string, description: string | undefined, hidden: boolean): void {
    this.#label.set(label);
    this.#description.set(description);
    this.#hidden.set(hidden);
  }

  setChildren(children: readonly SceneNode[]): void {
    this.#children.set(children);
  }

  toggleVisibility(): void {
    this.actions.visibility?.(!this.#hidden.value);
  }
  toggleExpanded(): void {
    this.#expanded.set(!this.#expanded.value);
    if (this.#expanded.value && this.lazy) this.load?.();
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
export function createSceneSource(scope: SourceScope, plugin: PluginContext, api: ViewerApi): SceneSource {
  const nodes = new Map<string, Node>();
  const roots = observableState<readonly SceneNode[]>([]);
  const picked = observableState<SceneNode | undefined>(undefined);
  const inspections = new Map<string, SceneInspection>();
  /** The structure a unit cell row is under (its operations name the structure). */
  const unitCellOwners = new Map<string, string>();
  /** The structure a volume streaming row is under, and the streaming's behavior cell (its params change). */
  const streamingOwners = new Map<string, string>();
  const streamingBehaviors = new Map<string, string>();

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
      // parity.md §2: the Conformation Transform's reset (moving is the inspector's Transform).
      n.menu = [
        {
          id: 'reset-transform',
          label: 'Reset transform',
          icon: 'rotate-ccw',
          run: () => void api.run('structure.transform-reset', { structure: ref }),
        },
      ];
      nodes.set(ref, n);
    }
    const data = s.cell.obj?.data;
    n.set(name, data ? describeStructure(data) : undefined, !!s.cell.state.isHidden);
    // An open inspector follows the transform (Undo, superposition), without reporting it as an edit.
    if (data) inspections.get(ref)?.form?.load({ transform: transformValues(currentTransform(s, data), data) });
    return n;
  };

  /** The unit cell of a structure's model (v5's model-level Unit Cell row; parity.md §2). */
  const unitCellNode = (u: ModelUnitcellRef, structureRef: string): Node => {
    const ref = u.cell.transform.ref;
    unitCellOwners.set(ref, structureRef);
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'unit-cell', 'box', 'Unit cell', false, {
        visibility: (hide) =>
          void api.run('structure.unit-cell', { structure: unitCellOwners.get(ref), visible: !hide }),
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
          const owner = streamingOwners.get(ref);
          const view = hide ? 'off' : on.entry.kind === 'em' ? 'auto' : 'selection-box';
          void api.run('volume.streaming-view', { structure: owner, view });
        },
        remove: () => void api.run('volume.streaming', { structure: streamingOwners.get(ref), enabled: false }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    n.set('Volume streaming', streamingSummary(on), on.params.entry.params.view.name === 'off');
    inspections.get(ref)?.form?.load(streamingValues(on));
    return n;
  };

  const componentNode = (c: StructureComponentRef, structureRef: string): Node => {
    const ref = c.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
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
      n.menu = componentMenu(plugin, api, ref, structureRef, loci);
      nodes.set(ref, n);
    }
    const reprs = c.representations.map((r) => reprLabel(plugin, r)).join(', ');
    const atoms = c.cell.obj?.data.elementCount;
    n.set(
      c.cell.obj?.label ?? 'Component',
      [reprs, atoms !== undefined ? `${atoms.toLocaleString('en')} atoms` : undefined].filter(Boolean).join(' · ') ||
        undefined,
      !!c.cell.state.isHidden,
    );
    return n;
  };

  const reprNode = (r: StructureRepresentationRef, componentRef: string, structureRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    const type = (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name ?? '';
    if (!n) {
      n = new Node(ref, 'representation', ReprIcons[type] ?? 'cartoon', reprLabel(plugin, r), false, {
        visibility: (hide) => void api.run('representation.visibility', { representation: ref, visible: !hide }),
        remove: () => void api.run('representation.remove', { component: componentRef, type, structure: structureRef }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    n.icon = ReprIcons[type] ?? 'cartoon';
    n.set(reprLabel(plugin, r), colorLabel(plugin, r), !!r.cell.state.isHidden);
    // An open inspector follows the cell (Undo, another surface), without reporting it as an edit.
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
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
      node.load = () => {
        node.lazy = false;
        const data = s.cell.obj?.data;
        if (!data) return;
        node.setChildren(entities(data, structureRef, nodes, api));
      };
      n = node;
      nodes.set(id, n);
    }
    const data = s.cell.obj?.data;
    n.set('Contents', data ? `${data.model.entities.data._rowCount} entities` : undefined, false);
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
      const locis = () => readMeasurements(plugin).find((x) => x.ref === ref)?.locis ?? [];
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

  // Volumes and particles (plan 2.5): top-level nodes like structures, their representations as children. Eyes and
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
    return n;
  };

  const volumeReprNode = (r: VolumeRepresentationRef, volumeRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'volume-representation', 'volume', '', false, {
        visibility: (hide) =>
          void api.run('volume.visibility', { volume: volumeRef, representation: ref, visible: !hide }),
        remove: () => void api.run('volume.remove', { volume: volumeRef, representation: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const iso = (r.cell.transform.params as { type?: { params?: { isoValue?: Volume.IsoValue } } } | undefined)?.type
      ?.params?.isoValue;
    n.set(
      r.cell.obj?.label ?? 'Representation',
      iso ? Volume.IsoValue.toString(iso) : undefined,
      !!r.cell.state.isHidden,
    );
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
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
    const count = l.cell.obj?.data.count;
    n.set(name, count !== undefined ? `${count.toLocaleString('en')} particles` : undefined, !!l.cell.state.isHidden);
    return n;
  };

  const particlesReprNode = (r: ParticleRepresentationRef, listRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'particles-representation', 'shapes', '', false, {
        visibility: (hide) =>
          void api.run('particles.visibility', { particles: listRef, representation: ref, visible: !hide }),
        remove: () => void api.run('particles.remove', { particles: listRef, representation: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    n.set(r.cell.obj?.label ?? 'Representation', undefined, !!r.cell.state.isHidden);
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    return n;
  };

  const sync = () => {
    const list = plugin.managers.structure.hierarchy.current.structures;
    const names = structureNames(list);
    const seen = new Set<string>();
    const top: SceneNode[] = [];
    for (const s of list) {
      const ref = s.cell.transform.ref;
      const sn = structureNode(s, names.get(ref) ?? ref);
      seen.add(ref);
      const contents = contentsNode(s, ref);
      seen.add(contents.id);
      const kids: SceneNode[] = [contents];
      for (const c of s.components) {
        const cn = componentNode(c, ref);
        seen.add(cn.id);
        const reprs = c.representations.map((r) => {
          const rn = reprNode(r, cn.id, ref);
          seen.add(rn.id);
          return rn;
        });
        cn.setChildren(reprs);
        kids.push(cn);
      }
      // The model's unit cell, under the first of its structures (another structure of the same model shares it).
      const cell = s.model?.unitcell;
      if (cell?.cell.obj?.data?.repr && !seen.has(cell.cell.transform.ref)) {
        const un = unitCellNode(cell, ref);
        seen.add(un.id);
        kids.push(un);
      }
      const streaming = streamingOf(plugin, s);
      if (streaming) {
        const vn = streamingNode(streaming, ref);
        seen.add(vn.id);
        kids.push(vn);
      }
      sn.setChildren(kids);
      top.push(sn);
    }
    const volumes = plugin.managers.volume.hierarchy.current.volumes;
    const volumeLabels = namesOf(volumes);
    for (const v of volumes) {
      const ref = v.cell.transform.ref;
      const vn = volumeNode(v, volumeLabels.get(ref) ?? ref);
      seen.add(ref);
      vn.setChildren(
        v.representations.map((r) => {
          const rn = volumeReprNode(r, ref);
          seen.add(rn.id);
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
      seen.add(ref);
      ln.setChildren(
        l.representations.map((r) => {
          const rn = particlesReprNode(r, ref);
          seen.add(rn.id);
          return rn;
        }),
      );
      top.push(ln);
    }
    // Measurements (they may span structures): one group after the data, a row each; their eyes and removes go
    // through `measure.*` (the measurements feature), their hover and focus through the measured atoms.
    const measured = readMeasurements(plugin);
    if (measured.length) {
      const gn = measurementsNode();
      seen.add(gn.id);
      gn.setChildren(
        measured.map((m) => {
          const mn = measurementNode(m);
          seen.add(mn.id);
          return mn;
        }),
      );
      top.push(gn);
    }
    // Lazy contents' entity and chain nodes stay with their structure.
    for (const [id, n] of nodes) {
      const owner = id.includes('#') ? id.slice(0, id.indexOf('#')) : id;
      if (!seen.has(id) && !seen.has(owner)) {
        nodes.delete(id);
        inspections.delete(id);
        unitCellOwners.delete(id);
        streamingOwners.delete(id);
        if (picked.value === n) picked.set(undefined);
      }
    }
    roots.set(top);
  };

  scope.add(plugin.managers.structure.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.volume.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.particles.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.structure.measurement.behaviors.state.subscribe(() => sync()));
  // Streaming's cells are not in the hierarchies: their coming and going rebuilds, a params change refreshes the row.
  const isStreaming = (t: { transformer: unknown } | undefined) =>
    t?.transformer === CreateVolumeStreamingInfo || t?.transformer === CreateVolumeStreamingBehavior;
  scope.add(plugin.state.data.events.cell.created.subscribe((e) => isStreaming(e.cell.transform) && sync()));
  scope.add(
    plugin.state.data.events.cell.removed.subscribe((e) => {
      if (!streamingOwners.has(e.ref) && !streamingBehaviors.has(e.ref)) return;
      streamingBehaviors.delete(e.ref);
      sync();
    }),
  );
  scope.add(
    plugin.state.data.events.object.updated.subscribe((e) => {
      const root = streamingBehaviors.get(e.ref);
      const owner = root ? streamingOwners.get(root) : undefined;
      const s = owner ? plugin.managers.structure.hierarchy.current.refs.get(owner) : undefined;
      const on = s?.kind === 'structure' ? streamingOf(plugin, s) : undefined;
      if (on) streamingNode(on, owner!);
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
    inspect: (n) => {
      let inspection = inspections.get(n.id);
      if (!inspection) {
        inspection = inspect(plugin, api, n, unitCellOwners.get(n.id) ?? streamingOwners.get(n.id));
        inspections.set(n.id, inspection);
      }
      return inspection;
    },
    sync,
  };
}

function Structure_toLoci(s: Structure): StructureElement.Loci {
  return StructureElement.Loci.all(s);
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
function componentMenu(
  plugin: PluginContext,
  api: ViewerApi,
  ref: string,
  structureRef: string,
  loci: () => StructureElement.Loci | undefined,
): SceneMenuItem[] {
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
    {
      id: 'modify',
      label: 'Modify by selection',
      icon: 'squares-unite',
      items: (['union', 'subtract', 'intersect'] as const).map((mode) => ({
        id: mode,
        label: mode === 'union' ? 'Include' : mode === 'subtract' ? 'Subtract' : 'Intersect',
        icon: mode === 'union' ? 'squares-unite' : mode === 'subtract' ? 'squares-subtract' : 'squares-intersect',
        run: () =>
          void api.run('component.modify', { component: ref, structure: structureRef, target: 'selection', mode }),
      })),
    },
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

/** A structure's entities, each with its chains (lazy contents, decision 18: a row's eye hides it). */
function entities(s: Structure, structureRef: string, nodes: Map<string, Node>, api: ViewerApi): SceneNode[] {
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
      });
      n.removable = false;
      nodes.set(key, n);
    }
    n.set(label, `${entities.data.type.value(i)} · ${chains.size} chain${chains.size === 1 ? '' : 's'}`, false);
    n.menu = contentsMenu(api, () => contentsLoci(s, id), label);
    const chainNodes: SceneNode[] = [];
    for (const [chain, atoms] of chains) {
      const ck = `${structureRef}#chain:${id}:${chain}`;
      let c = nodes.get(ck);
      if (!c) {
        c = new Node(ck, 'chain', 'sequence', `Chain ${chain}`, false, {
          visibility: (hide) => void api.run('contents.visibility', { structure: structureRef, chain, visible: !hide }),
        });
        c.removable = false;
        nodes.set(ck, c);
      }
      c.set(`Chain ${chain}`, `${atoms.toLocaleString('en')} atoms`, false);
      c.menu = contentsMenu(api, () => contentsLoci(s, id, chain), `${label} · chain ${chain}`);
      chainNodes.push(c);
    }
    n.setChildren(chainNodes);
    out.push(n);
  }
  return out;
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

/** v5's iso ranges per channel (σ), as its panel's sliders have them. */
const ChannelBounds: Readonly<Record<VolumeStreaming.ChannelType, readonly [number, number]>> = {
  em: [-5, 5],
  '2fo-fc': [0, 3],
  'fo-fc(+ve)': [1, 5],
  'fo-fc(-ve)': [-5, -1],
};

const ChannelLabels: Readonly<Record<VolumeStreaming.ChannelType, string>> = {
  em: 'EM',
  '2fo-fc': '2Fo-Fc',
  'fo-fc(+ve)': 'Fo-Fc(+ve)',
  'fo-fc(-ve)': 'Fo-Fc(-ve)',
};

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
  for (const c of streamingChannels(on.entry)) {
    const [min, max] = ChannelBounds[c];
    params[`iso:${c}`] = PD.Numeric(min, { min, max, step: 0.05 }, { label: `${ChannelLabels[c]} (σ)` });
    params[`color:${c}`] = PD.Color(Color(0), { label: `${ChannelLabels[c]} color` });
  }
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
  const sampling = on.entry.header.sampling[0];
  for (const c of streamingChannels(on.entry)) {
    const ch = (p.channels as Partial<Record<string, VolumeStreaming.ChannelParams>>)[c];
    if (!ch) continue;
    const stats = sampling?.valuesInfo[c.startsWith('fo-fc') ? 1 : 0];
    const iso =
      ch.isoValue.kind === 'relative' || !stats
        ? (ch.isoValue as Volume.IsoValue.Relative).relativeValue
        : Volume.IsoValue.toRelative(ch.isoValue, stats).relativeValue;
    values[`iso:${c}`] = Math.round(iso * 100) / 100;
    values[`color:${c}`] = ch.color;
  }
  return values;
}

/** The form's edits as `volume.streaming-view` calls: the fields changed in a moment, as one call (sliders). */
function streamingChanges(api: ViewerApi, structure: string): (values: PD.Values<PD.Params>, name: string) => void {
  let timer: ReturnType<typeof setTimeout> | undefined;
  let pending: { view?: string; radius?: number; detail?: number; channels: StreamingChannelsChange } = {
    channels: {},
  };
  return (values, name) => {
    const v = values as Record<string, any>;
    if (name === 'view') pending.view = v.view;
    else if (name === 'radius') pending.radius = v.radius;
    else if (name === 'cameraRadius') pending.radius = v.cameraRadius;
    else if (name === 'detail') pending.detail = v.detail;
    else {
      const [what, channel] = name.split(':');
      const key = what === 'iso' ? 'isoValue' : 'color';
      pending.channels = { ...pending.channels, [channel]: { ...pending.channels[channel], [key]: v[name] } };
    }
    clearTimeout(timer);
    timer = setTimeout(() => {
      const change = pending;
      pending = { channels: {} };
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
function inspect(plugin: PluginContext, api: ViewerApi, node: SceneNode, owner?: string): SceneInspection {
  const cell = plugin.state.data.cells.get(node.id);
  const slot = `inspector:${node.kind}`;
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
      return {
        node,
        form: new ParamsFormModel(streamingForm(on), streamingValues(on), { onChange: streamingChanges(api, owner) }),
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
    const now = ref?.kind === 'structure' ? currentTransform(ref, s) : Mat4.identity();
    const form = new ParamsFormModel(
      TransformForm,
      { transform: transformValues(now, s) },
      {
        onChange: debounced(
          ({ transform }) =>
            void api.run('structure.transform', {
              structure: node.id,
              translation: transform.translation,
              axis: transform.axis,
              angle: transform.angle,
              center: 'centroid',
            }),
        ),
      },
    );
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
    return {
      node,
      form: undefined,
      slot,
      facts: [
        ...(data.sourceData.kind ? [['Source', data.sourceData.kind.toUpperCase()] as const] : []),
        ['Grid', gridText(stats)],
        ['Voxel', voxelText(stats)],
        ['Mean · σ', `${signed(stats.mean)} · ${signed(stats.sigma)}`],
        ['Range', `${signed(stats.min)} … ${signed(stats.max)}`],
      ],
    };
  }
  if (node.kind === 'particles' && cell?.obj?.data) {
    const data = cell.obj.data as ParticleList;
    return {
      node,
      form: undefined,
      slot,
      facts: [
        ['Particles', data.count.toLocaleString('en')],
        ['Targets', String(data.targetInfo.size)],
      ],
    };
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
          timer = setTimeout(() => void api.run(updateOperation, { representation: node.id, params: values }), 150);
        },
      });
      return { node, form, slot };
    }
  }
  return { node, form: undefined, slot };
}
