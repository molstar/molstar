/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { StateTransform } from '@molstar/core/state';
import type { Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui';
import { merge } from 'rxjs';
import { fromBehavior, fromChanged, SourceScope, sameItems, scanList } from './state.js';

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
  readonly tasks: ReadonlyState<readonly TaskInfo[]>;
  dispose(): void;
}

export function createViewerSources(plugin: PluginContext): ViewerSources {
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
        a.props === b.props),
  );
  const structures = fromChanged(
    scope,
    () => structureInfos(plugin),
    plugin.managers.structure.hierarchy.behaviors.selection,
    (a, b) => sameItems(a, b, sameStructure),
  );
  const empty = fromChanged(scope, () => isEmpty(plugin), plugin.state.data.events.changed);
  const busy = fromBehavior(scope, plugin.behaviors.state.isBusy);
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
  return { canvas, structures, empty, busy, tasks, dispose: () => scope.dispose() };
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
  return plugin.managers.structure.hierarchy.current.structures.map((s) => {
    const data = s.cell.obj?.data;
    const model = s.model?.cell.obj?.data;
    return {
      ref: s.cell.transform.ref,
      label: s.cell.obj?.label ?? s.cell.transform.ref,
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
