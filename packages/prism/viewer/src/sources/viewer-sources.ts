/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { StateTransform } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DParams, type Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { merge } from 'rxjs';
import { createInteractionSources, type InteractionSources } from './interaction.js';
import { createMarksSource, type MarksSource, type MarksView, marksView } from './marks.js';
import { structureNames } from './naming.js';
import { createPointerSources, type PointerSources } from './pointer.js';
import { createSettingsSources, type SettingsSources } from './settings.js';
import { fromBehavior, fromChanged, SourceScope, sameItems, scanList } from './state.js';
import { readStyleLayers, sameStyleLayers, type StyleLayerInfo } from './style-layers.js';
import { watchedState } from './watched.js';

// The viewer's sources (state.md §4): read-only states over the plugin, the only code besides `api/` that reads it.
// What views, view models, the API's results and the agent's account read. Built per area as each is built (plan).

/** The canvas as the viewer reads it; undefined headless (no canvas). */
export interface CanvasState {
  /** The background, `#rrggbb`. */
  readonly background: string;
  readonly transparentBackground: boolean;
  readonly illumination: boolean;
  /** The props themselves, for readers of other parts. */
  readonly props: Canvas3DProps;
}

/** A loaded structure, for the scene outline (`structure.get`) and the agent's account. */
export interface StructureInfo {
  /** Its state ref. */
  readonly ref: string;
  readonly label: string;
  /** The model's label (the entry id). */
  readonly model: string;
  readonly atoms: number;
  readonly residues: number;
  readonly chains: number;
  /** Its components, by label. */
  readonly components: readonly string[];
}

/** A running task (the status column). */
export interface TaskInfo {
  readonly id: number;
  readonly label: string;
  /** 0..1, or undefined when indeterminate. */
  readonly progress?: number;
  /** It holds the plugin (rendering an animation): an overlay with Cancel. */
  readonly overlay: boolean;
  readonly abort?: () => void;
}

export interface ViewerSources {
  readonly canvas: ReadonlyState<CanvasState | undefined>;
  readonly structures: ReadonlyState<readonly StructureInfo[]>;
  /** Nothing in the data tree. */
  readonly empty: ReadonlyState<boolean>;
  /** The plugin is updating or animating for a while (its busy flag). */
  readonly busy: ReadonlyState<boolean>;
  /** The plugin updates its state now (not debounced, as `busy` is): what `ViewerCore.state.working` follows. */
  readonly updating: ReadonlyState<boolean>;
  readonly tasks: ReadonlyState<readonly TaskInfo[]>;
  /** What is under the pointer: the hover label's lines, hover and click events. */
  readonly pointer: PointerSources;
  /** The selection mode, the picking level, the selection and the focus. */
  readonly interaction: InteractionSources;
  /** The file formats the plugin reads, by label (the drop area's hint): read when asked, so a lazy extension's are in. */
  readonly formats: readonly string[];
  /** The Settings tab's groups with their current values. */
  readonly settings: SettingsSources;
  /**
   * The style layers applied (color, transparency, glow, material, clipping; `style-layers.ts`): followed only while
   * subscribed, re-read once per animation frame, and only walked again when a layer cell changed.
   */
  readonly layers: ReadonlyState<readonly StyleLayerInfo[]>;
  /** A view's marks and hover on the scene (marks.md): in the 3D view and the sequence, never a history step. */
  readonly marks: MarksView;
}

/** The sources as the API's areas have them: the marks with their writers (`highlight`). */
export interface ApiSources extends ViewerSources {
  readonly marks: MarksSource;
}

/**
 * The sources as `viewer.sources` gives them: every one a state to read (api.md §6), the marks without their writers.
 * Getters stay getters (`formats` is read when asked).
 */
export function publicSources(sources: ApiSources): ViewerSources {
  return Object.create(sources, { marks: { value: marksView(sources.marks), enumerable: true } }) as ViewerSources;
}

/** The sources and their disposer, kept by whoever made them. */
export interface ViewerSourcesHost {
  readonly sources: ApiSources;
  dispose(): void;
}

export function createViewerSources(plugin: PluginContext): ViewerSourcesHost {
  const scope = new SourceScope();
  const canvas = fromChanged(
    scope,
    () => canvasState(plugin),
    merge(plugin.events.canvas3d.settingsUpdated, plugin.behaviors.canvas3d.initialized),
    (a, b) =>
      a === b ||
      (!!a &&
        !!b &&
        a.background === b.background &&
        a.transparentBackground === b.transparentBackground &&
        a.illumination === b.illumination &&
        // By value: `canvas3d.props` is a new object at each read, and `settingsUpdated` fires for a camera focus too.
        PD.areEqual(Canvas3DParams, a.props, b.props)),
  );
  const structures = fromChanged(
    scope,
    () => structureInfos(plugin),
    plugin.managers.structure.hierarchy.behaviors.selection,
    (a, b) => sameItems(a, b, sameStructure),
  );
  const empty = fromChanged(scope, () => isEmpty(plugin), plugin.state.data.events.changed);
  const busy = fromBehavior(scope, plugin.behaviors.state.isBusy);
  const updating = fromBehavior(scope, plugin.behaviors.state.isUpdating);
  const tasks = scanList<TaskEvent, TaskInfo>(
    scope,
    merge<TaskEvent[]>(plugin.managers.task.events.progress, plugin.managers.task.events.finished),
    (list, e: TaskEvent) => {
      if (!('progress' in e)) return list.filter((t) => t.id !== e.id);
      if (e.level === 'none') return list;
      const root = e.progress.root.progress;
      const info: TaskInfo = {
        id: e.id,
        label: root.message,
        ...(root.isIndeterminate ? {} : { progress: root.max > 0 ? root.current / root.max : undefined }),
        overlay: !!e.useOverlay,
        ...(root.canAbort ? { abort: () => plugin.managers.task.requestAbort(e.progress, 'Cancelled') } : {}),
      };
      const i = list.findIndex((t) => t.id === e.id);
      return i < 0 ? [...list, info] : list.map((t, j) => (j === i ? info : t));
    },
  );
  const pointer = createPointerSources(scope, plugin);
  const interaction = createInteractionSources(scope, plugin);
  const settings = createSettingsSources(plugin);
  const memo = {};
  const layers = watchedState(
    () => readStyleLayers(plugin, memo),
    () => [plugin.state.data.events.changed, plugin.managers.structure.hierarchy.behaviors.selection],
    { compare: sameStyleLayers },
  );
  return {
    sources: {
      canvas,
      structures,
      empty,
      busy,
      updating,
      tasks,
      pointer,
      interaction,
      get formats() {
        return plugin.dataFormats.options.map(([, label]) => label);
      },
      settings,
      layers,
      marks: createMarksSource(scope, plugin),
    },
    dispose: () => scope.dispose(),
  };
}

type TaskEvent =
  | { id: number; level: 'none' | 'background'; useOverlay?: boolean; progress: import('@molstar/core/task').Progress }
  | { id: number };

function canvasState(plugin: PluginContext): CanvasState | undefined {
  const props = plugin.canvas3d?.props;
  if (!props) return undefined;
  return {
    background: Color.toHexStyle(props.renderer.backgroundColor),
    transparentBackground: !!props.transparentBackground,
    illumination: !!props.illumination?.enabled,
    props,
  };
}

function isEmpty(plugin: PluginContext): boolean {
  const tree = plugin.state.data.tree;
  return tree.children.get(StateTransform.RootRef).size === 0;
}

function structureInfos(plugin: PluginContext): StructureInfo[] {
  const list = plugin.managers.structure.hierarchy.current.structures;
  const names = structureNames(list);
  return list.map((s) => {
    const data = s.cell.obj?.data;
    const model = s.model?.cell.obj?.data;
    return {
      ref: s.cell.transform.ref,
      label: names.get(s.cell.transform.ref) ?? s.cell.transform.ref,
      model: model?.entryId ?? model?.label ?? '',
      atoms: data?.elementCount ?? 0,
      residues: data?.polymerResidueCount ?? 0,
      chains: data ? new Set(data.units.map((u) => u.chainGroupId)).size : 0,
      components: s.components.map((c) => c.cell.obj?.label ?? c.key ?? ''),
    };
  });
}

function sameStructure(a: StructureInfo, b: StructureInfo): boolean {
  return (
    a.ref === b.ref &&
    a.label === b.label &&
    a.atoms === b.atoms &&
    a.residues === b.residues &&
    a.chains === b.chains &&
    sameItems(a.components, b.components)
  );
}
