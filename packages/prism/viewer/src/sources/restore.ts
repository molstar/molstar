/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { PluginState } from '@molstar/plugin/state';
import type { HistoryManager } from '@molstar/prism-ui/core';
import type { PrismStore } from '../spec.js';
import type { SourceScope } from './state.js';

// The light restore (journal.md §6): after each step, debounced, this tab's scene without its data (the state tree's
// transforms and params, the camera, the canvas, the sequence's state) when every data cell comes from an address, so a
// reload can fetch it again; a scene with a local file or data given by value is not kept, and the viewer says once
// what will not come back. It runs only when the host asks (a reload of the page, the URL asking for nothing to load:
// the app decides), after the viewer is ready; the history starts empty.

/** The key in the `restore` store. */
const RestoreKey = 'restore';
/** Data cells that hold their bytes (a file, data given by value): such a scene does not come back by itself. */
const HeldData = /(?:^|\.)(?:read-file|raw-data|import-string|import-json)$/;

const RestoreParams: PluginState.SnapshotParams = {
  data: true,
  behavior: true,
  componentManager: true,
  canvas3d: true,
  camera: true,
  interactivity: false,
  structureSelection: false,
  canvas3dContext: false,
  animation: false,
  startAnimation: false,
};

export interface LightRestore {
  /** A restore is kept for this tab (by this version of the app). */
  readonly available: boolean;
  /** Brings the kept scene back; false when there was none or it failed (it is dropped then). */
  run(): Promise<boolean>;
  /** Writes now (a test, or before the page goes). */
  flush(): void;
}

export interface LightRestoreOptions {
  /** After a step, the write waits this long for the next (ms; 1000). */
  readonly debounceMs?: number;
  /** Said once when the scene is not kept (a local file). */
  readonly notify?: (message: string) => void;
}

export function createLightRestore(
  scope: SourceScope,
  plugin: PluginContext,
  history: HistoryManager<unknown>,
  store: PrismStore,
  o: LightRestoreOptions = {},
): LightRestore {
  const debounceMs = o.debounceMs ?? 1000;
  let timer: ReturnType<typeof setTimeout> | undefined;
  let told = false;
  // A write that took long makes the next wait longer (journal.md §6: skipped while its last write took long).
  let wait = debounceMs;
  let restoring = false;

  const write = () => {
    timer = undefined;
    if (restoring) return;
    const t0 = performance.now();
    const cells = [...plugin.state.data.cells.values()];
    if (!cells.some((c) => c.transform.parent !== c.transform.ref && c.obj)) {
      store.remove(RestoreKey);
      return;
    }
    if (cells.some((c) => HeldData.test(c.transform.transformer.id))) {
      store.remove(RestoreKey);
      if (!told) o.notify?.('Local files are not kept across a reload; save a session to keep them.');
      told = true;
      return;
    }
    // The camera where there is one (a headless viewer has none).
    store.set(RestoreKey, plugin.state.getSnapshot({ ...RestoreParams, camera: !!plugin.canvas3d }));
    const ms = performance.now() - t0;
    wait = ms > 50 ? Math.min(wait * 2, 30_000) : debounceMs;
  };
  const schedule = () => {
    if (timer !== undefined) clearTimeout(timer);
    timer = setTimeout(write, wait);
  };
  scope.add(history.steps.subscribe(() => schedule()));
  scope.add(() => {
    if (timer !== undefined) clearTimeout(timer);
  });

  return {
    get available() {
      return store.get(RestoreKey) !== undefined;
    },
    async run() {
      const snapshot = store.get(RestoreKey) as PluginState.Snapshot | undefined;
      if (!snapshot) return false;
      restoring = true;
      try {
        await plugin.state.setSnapshot(snapshot);
        // Not a step: a fresh start that looks like the last one (journal.md §6, the history starts empty).
        history.clear();
        return true;
      } catch (e) {
        console.warn('[prism-viewer] the light restore failed:', e);
        store.remove(RestoreKey);
        return false;
      } finally {
        restoring = false;
      }
    },
    flush() {
      if (timer !== undefined) clearTimeout(timer);
      write();
    },
  };
}
