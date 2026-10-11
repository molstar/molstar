/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

/** A saved snapshot (Mol*'s `PluginStateSnapshotManager.Entry`): what the list shows and the editor edits. */
export interface SnapshotEntry {
  id: string;
  name: string;
  /** Activates it from descriptions and labels (`[title](#key)`). */
  key?: string;
  /** Markdown, shown with the snapshot (the playback card). */
  description?: string;
  /** When it was taken ("2 min ago"). */
  timestamp?: string;
  /** A thumbnail when the snapshot was saved with an image (a CSS background for now). */
  image?: string;
}

/**
 * The plugin's snapshots and how they are saved (implemented over the snapshot manager and `PluginCommands.State`):
 * add the current state, apply one, edit its name, key and description, reorder, replace with the current state,
 * remove; download or open sessions and states.
 */
export interface SnapshotsModel {
  entries: ReadonlyState<readonly SnapshotEntry[]>;
  /** The one applied last (marked in the list, edited below it). */
  current: ReadonlyState<string | undefined>;
  /** What a new snapshot includes (Mol*'s `PluginState.SnapshotParams`: data, camera, transition, image, …). */
  saveOptions: ParamsFormModel<any>;
  /** A snapshot's own playback settings (how long it shows, how the camera gets there). */
  entryOptions(id: string): ParamsFormModel<any> | undefined;
  add(name?: string): void;
  apply(id: string): void;
  update(id: string, patch: Partial<Pick<SnapshotEntry, 'name' | 'key' | 'description'>>): void;
  /** Takes the current state into this snapshot (keeps its name, key and description). */
  replace(id: string): void;
  move(id: string, by: -1 | 1): void;
  /** Moves it to a position in the list (dragged by its handle). */
  moveTo(id: string, index: number): void;
  remove(id: string): void;
  clear(): void;
  /** The state with its data (.molx), or the state description alone (.molj; data loads from its sources). */
  download(kind: 'session' | 'state'): void;
  open(files: readonly File[]): void;
}
