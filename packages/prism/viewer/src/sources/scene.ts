/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { emptiedComponents, hiddenContents } from '../api/contents.js';
import { StructureFocusRepresentationTags } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { InitVolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { SpacegroupCell } from '@molstar/core/math/geometry';
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
import { type MeasuredCell, MeasureKinds, readMeasurements } from './measurements.js';
import { structureNames } from './naming.js';
import type { SourceScope } from './state.js';
import { createVolumeLevels, mapLevels, streamingLevels, type VolumeLevelsSource } from './volume-levels.js';
import { createVolumeTransfer, type VolumeTransferSource } from './volume-transfer.js';
import { createRepresentationQuick, type RepresentationQuickSource } from './representation-quick.js';

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
  #swatch?: ObservableState<number | undefined>;
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
/** Opens a command in the command menu with values (a component's own row: the component filled in). */
export type OpenCommand = (id: string, values: Readonly<Record<string, unknown>>) => void;

export function createSceneSource(
  scope: SourceScope,
  plugin: PluginContext,
  api: ViewerApi,
  openCommand?: OpenCommand,
): SceneSource {
  const nodes = new Map<string, Node>();
  const roots = observableState<readonly SceneNode[]>([]);
  const picked = observableState<SceneNode | undefined>(undefined);
  const inspections = new Map<string, SceneInspection>();
  /** The structure a unit cell row is under (its operations name the structure). */
  const unitCellOwners = new Map<string, string>();
  /** The structure a volume streaming row is under, and the streaming's behavior cell (its params change). */
  const streamingOwners = new Map<string, string>();
  const streamingBehaviors = new Map<string, string>();
  /** How map levels are shown and typed (σ or map units), one choice for every map. */
  const levelUnit = observableState<VolumeLevelUnit>('sigma');
  const refreshLevels = (id: string) => {
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
    n.set(name, data ? describeStructure(data) : undefined, !!s.cell.state.isHidden);
    // Its menu follows what the structure offers (its assemblies, a crystal's symmetry, streaming on or off).
    n.menu = structureMenu(plugin, api, s);
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
    refreshLevels(ref);
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
          if (!hide) void api.run('volume.streaming', { structure: structureRef, enabled: true });
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
      n.menu = componentMenu(plugin, api, ref, structureRef, loci, openCommand, c.cell.obj?.label ?? '');
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
    // hidden, a changed structure).
    const node = n;
    node.load = () => {
      node.lazy = false;
      const data = s.cell.obj?.data;
      if (data) node.setChildren(entities(plugin, data, structureRef, nodes, api));
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
        visibility: (hide) => {
          for (const { component, structure } of focusPartsNow) {
            void api.run('component.visibility', {
              component: component.cell.transform.ref,
              structure,
              visible: !hide,
            });
          }
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
    refreshLevels(ref);
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
    const allFocused: SceneNode[] = [];
    const allFocusParts: { component: StructureComponentRef; structure: string }[] = [];
    let focusedOn: string | undefined;
    for (const s of list) {
      const ref = s.cell.transform.ref;
      const sn = structureNode(s, names.get(ref) ?? ref);
      seen.add(ref);
      const contents = contentsNode(s, ref);
      seen.add(contents.id);
      const kids: SceneNode[] = [contents];
      const focused: SceneNode[] = [];
      for (const c of s.components) {
        const isFocus = isFocusPart(c);
        // An emptied focus component (the focus cleared) has no row.
        if (isFocus && !c.cell.obj?.data.elementCount) continue;
        const cn = componentNode(c, ref);
        if (isFocus) cn.removable = false;
        seen.add(cn.id);
        const reprs = c.representations.map((r) => {
          const rn = reprNode(r, cn.id, ref);
          seen.add(rn.id);
          return rn;
        });
        cn.setChildren(reprs);
        (isFocus ? focused : kids).push(cn);
      }
      // Components emptied by hidden contents, where they stand among the others (the state tree's order).
      const emptied = emptiedComponents(plugin, ref);
      if (emptied.length) {
        const order = [...plugin.state.data.tree.children.get(emptied[0].parent).values()];
        for (const e of emptied) {
          const en = emptiedNode(e.ref, e.label);
          seen.add(en.id);
          kids.push(en);
        }
        // The contents row first; the components by their place in the tree.
        const [first, ...rest] = kids;
        rest.sort((a, b) => order.indexOf(a.id) - order.indexOf(b.id));
        kids.splice(0, kids.length, first, ...rest);
      }
      const parts = focusParts(plugin, ref).filter((c) => !!c.cell.obj?.data.elementCount);
      if (parts.length) focusedOn = names.get(ref) ?? ref;
      for (const component of parts) allFocusParts.push({ component, structure: ref });
      allFocused.push(...focused);
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
      } else if (canStream(plugin, s)) {
        const vn = streamingOffNode(ref);
        seen.add(vn.id);
        kids.push(vn);
      }
      sn.setChildren(kids);
      top.push(sn);
    }
    if (list.length) {
      // Named by its structure when there are several.
      const fn = focusNode(allFocusParts, list.length > 1 ? focusedOn : undefined);
      seen.add(fn.id);
      fn.setChildren(allFocused);
      top.push(fn);
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
    // The group is there once something is loaded, empty or not: its inspector says how to measure.
    const measured = readMeasurements(plugin);
    if (measured.length || top.length) {
      const gn = measurementsNode();
      seen.add(gn.id);
      gn.set('Measurements', measurementCounts(measured), false);
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
        inspection = inspect(plugin, api, n, unitCellOwners.get(n.id) ?? streamingOwners.get(n.id), levelUnit);
        inspections.set(n.id, inspection);
      }
      return inspection;
    },
    sync,
  };
}

/** A volume of the hierarchy by its ref. */
function volumeRefOf(plugin: PluginContext, ref: string): VolumeRef | undefined {
  return plugin.managers.volume.hierarchy.current.volumes.find((v) => v.cell.transform.ref === ref);
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
): SceneMenuItem[] {
  // With the command menu: each opens its command for this component (its fields filled in, the rest to choose there),
  // so the long lists (types, themes) are the command menu's to search, not a fold in a panel.
  if (openCommand) {
    // A preset component (Polymer, Ligand, Water, …) is also a target by its keyword; another is typed there.
    const keyword = ComponentKeywords[label.toLowerCase()];
    return [
      {
        id: 'add-repr',
        label: 'Add representation…',
        icon: 'plus',
        run: () => openCommand('representation.set', { component: ref, structure: structureRef }),
      },
      {
        id: 'coloring',
        label: 'Set coloring…',
        icon: 'palette',
        run: () =>
          openCommand('representation.color-theme', {
            structure: structureRef,
            ...(keyword ? { target: keyword } : {}),
          }),
      },
      {
        id: 'modify',
        label: 'Modify by selection…',
        icon: 'squares-unite',
        run: () => openCommand('component.modify', { component: ref, structure: structureRef, target: 'selection' }),
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

/**
 * A structure's contents in a few facts: its entities by type (polymer, branched, non-polymer, water), its chains, its
 * polymer residues, atoms, and models when there are several.
 */
export function contentsFacts(s: Structure): [string, string][] {
  const { entities } = s.model;
  const present = new Set<string>();
  const chains = new Set<string>();
  const loc = StructureElement.Location.create(s);
  for (const unit of s.units) {
    loc.unit = unit as Unit;
    loc.element = unit.elements[0];
    present.add(StructureProperties.chain.label_entity_id(loc));
    chains.add(`${unit.conformation.operator.name}:${StructureProperties.chain.label_asym_id(loc)}`);
  }
  const byType = new Map<string, number>();
  for (let i = 0; i < entities.data._rowCount; i++) {
    if (!present.has(entities.data.id.value(i))) continue;
    const type = entities.data.type.value(i) || 'other';
    byType.set(type, (byType.get(type) ?? 0) + 1);
  }
  const order = ['polymer', 'branched', 'non-polymer', 'water', 'macrolide', 'other'];
  const types = [...byType].sort((a, b) => order.indexOf(a[0]) - order.indexOf(b[0]));
  const n = (x: number) => x.toLocaleString('en');
  const facts: [string, string][] = [
    ['Entities', `${n(present.size)}${types.length ? ` (${types.map(([t, c]) => `${c} ${t}`).join(', ')})` : ''}`],
    ['Chains', n(chains.size)],
    ['Residues', n(s.polymerResidueCount)],
    ['Atoms', n(s.elementCount)],
  ];
  if (s.models.length > 1) facts.push(['Models', n(s.models.length)]);
  return facts;
}

/** A structure's entities, each with its chains (lazy contents, decision 18: a row's eye hides it). */
function entities(
  plugin: PluginContext,
  s: Structure,
  structureRef: string,
  nodes: Map<string, Node>,
  api: ViewerApi,
): SceneNode[] {
  const hidden = hiddenContents(plugin, structureRef);
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
    const entityHidden = hidden.entities.includes(id);
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
        entityHidden || hidden.chains.includes(`${id}:${chain}`) || hidden.chains.includes(chain),
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
  model: number;
  type: string;
  radius: number;
  cells: number;
}

type StructureType = { name: string; params?: { id?: string; radius?: number; ijkMax?: Vec3 } };

const typeOf = (s: StructureRef) => (s.cell.transform.params as { type: StructureType } | undefined)?.type;
const modelCount = (s: StructureRef) => s.model?.trajectory?.cell.obj?.data.frameCount ?? 1;
const modelIndex = (s: StructureRef) =>
  (s.model?.cell.transform.params as { modelIndex?: number } | undefined)?.modelIndex ?? 0;

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

function sourceForm(s: StructureRef, data: Structure): PD.Params {
  const count = modelCount(s);
  const types = structureTypes(s, data);
  return {
    // A short trajectory as a list, a long one as a number.
    model:
      count <= 50
        ? PD.Select<number>(
            1,
            Array.from({ length: count }, (_, i) => [i + 1, `${i + 1}`] as [number, string]),
            { label: 'Model', description: `The model shown, of ${count}`, isHidden: count < 2 },
          )
        : PD.Numeric(
            1,
            { min: 1, max: count, step: 1 },
            { label: 'Model', description: `The model shown, of ${count}` },
          ),
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
  };
}

/** The structure's form values: its source (model, type, radius, cells) and its transform. */
function structureValues(s: StructureRef, data: Structure): SourceValues & PD.Values<typeof TransformForm> {
  const t = typeOf(s);
  return {
    model: modelIndex(s) + 1,
    type: typeValue(s, data),
    radius: t?.name === 'symmetry-mates' ? (t.params?.radius ?? 5) : 5,
    cells: t?.name === 'symmetry' && t.params?.ijkMax ? t.params.ijkMax[0] : 1,
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
    if (name === 'model') return api.run('structure.model', { structure, index: v.model - 1 });
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
    if (name === 'radius' || name === 'cells') timer = setTimeout(run, 150);
    else run();
  };
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
      run: () => void api.run('volume.streaming', { structure, enabled: true }),
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
): SceneInspection {
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
    const params = { ...(structureRef ? sourceForm(structureRef, s) : {}), ...TransformForm };
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
    // A structure's contents (its id is the structure's ref + '#contents'): what it is made of, at a glance.
    const s = dataOf(plugin, node.id.slice(0, node.id.indexOf('#')));
    if (s) return { node, form: undefined, slot, facts: contentsFacts(s) };
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
  return { node, form: undefined, slot };
}
