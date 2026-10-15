/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell, StateTransform } from '@molstar/core/state';
import { StructureFocusRepresentationTags } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginState } from '@molstar/plugin/state';
import type { HistoryManager, HistoryPort } from '@molstar/prism-ui';
import { merge } from 'rxjs';

// The history's checkpoints (sequence.md §7.2): the plugin's state directly, never its snapshot manager (whose
// `getStateSnapshot` adds entries and whose restore clears; review F1). What undo restores: the data tree, behavior
// settings, canvas props, interactivity, the selection mode, and every owner's custom state (the sequence's, through
// `registerCustomState`). What it never touches: the camera, the selection and the focus (view state, not steps).

export interface ViewerCheckpoint {
  readonly plugin: PluginState.Snapshot;
  readonly selectionMode: boolean;
}

/** The snapshot's parts a checkpoint holds; `camera: false` keeps it from reading the canvas (headless). */
const CheckpointParams: PluginState.SnapshotParams = {
  data: true,
  behavior: true,
  componentManager: true,
  canvas3d: true,
  interactivity: true,
  structureSelection: false,
  camera: false,
  canvas3dContext: false,
  animation: false,
  startAnimation: false,
};

export function createCheckpointPort(plugin: PluginContext): HistoryPort<ViewerCheckpoint> {
  return {
    checkpoint: () => ({ plugin: plugin.state.getSnapshot(CheckpointParams), selectionMode: plugin.selectionMode }),
    restore: async (c) => {
      // The focus is view state: a restore leaves it (and the camera, which the checkpoint does not hold).
      await plugin.state.setSnapshot({ ...c.plugin, structureFocus: undefined });
      if (plugin.selectionMode !== c.selectionMode) plugin.selectionMode = c.selectionMode;
    },
    size: (c) => checkpointSize(c),
  };
}

/** The checkpoint's JSON size without data payloads (a file's text in a raw-data cell): what the budget measures. */
export function checkpointSize(c: ViewerCheckpoint): number {
  let size = 0;
  try {
    size = JSON.stringify(c.plugin, (_key, value) => {
      if (typeof value === 'string' && value.length > 1024) return '';
      if (ArrayBuffer.isView(value)) return '';
      return value;
    }).length;
  } catch {
    size = 0;
  }
  return size;
}

/** Cells owned by behaviors: they change outside the API by design and do not make checkpoints unsafe (§7.3). */
const BehaviorTags = new Set<string>(Object.values(StructureFocusRepresentationTags));

export function isBehaviorOwned(transform: StateTransform): boolean {
  if (transform.tags?.some((t) => BehaviorTags.has(t))) return true;
  return transform.transformer.id.includes('volume-streaming');
}

/**
 * Watches the data tree for changes outside the API (review F5): a cell created, removed or updated while the history
 * is not applying a run or a restore, and not owned by a behavior, makes every earlier checkpoint unsafe.
 */
export function watchOutsideChanges(plugin: PluginContext, history: HistoryManager<ViewerCheckpoint>): () => void {
  const state = plugin.state.data;
  const check = (transform: StateTransform | undefined) => {
    if (history.applying || !transform || isBehaviorOwned(transform)) return;
    history.invalidate('Changed outside the viewer: this step can no longer be undone.');
  };
  const transformOf = (ref: string): StateTransform | undefined =>
    (state.cells.get(ref) as StateObjectCell | undefined)?.transform ?? state.tree.transforms.get(ref);
  const removedTransforms = new Map<string, StateTransform>();
  const sub = merge(state.events.cell.created.pipe(), state.events.object.updated.pipe()).subscribe((e) =>
    check(transformOf(e.ref)),
  );
  // A removed cell is gone from the tree by the time the event comes: remember transforms as they are created.
  const created = state.events.cell.created.subscribe((e) => removedTransforms.set(e.ref, e.cell.transform));
  const removed = state.events.cell.removed.subscribe((e) => {
    check(removedTransforms.get(e.ref));
    removedTransforms.delete(e.ref);
  });
  for (const [ref, cell] of state.cells) removedTransforms.set(ref, cell.transform);
  return () => {
    sub.unsubscribe();
    created.unsubscribe();
    removed.unsubscribe();
  };
}
