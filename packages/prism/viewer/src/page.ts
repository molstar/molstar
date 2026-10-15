/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import type { PrismViewerSpec } from './spec.js';
import { mountShell } from './ui/shell.js';
import { createViewerCore, type PrismViewer } from './viewer.js';

// `createPrismViewer` (viewer-spec.md §6): the viewer in a page. The target gets a shadow root with the shell and the
// canvas; `ready` is set once the canvas is mounted, and only then do URL calls and the light restore run.

export interface CreatePrismViewerOptions {
  /** The element the viewer fills (it gets a shadow root); it should have a size. */
  readonly target: HTMLElement;
  readonly spec: PrismViewerSpec;
}

/** Live viewers by their persist key: two on a page with the same key are refused (README decision 29). */
const liveKeys = new Set<string>();

export async function createPrismViewer(o: CreatePrismViewerOptions): Promise<PrismViewer> {
  const key = o.spec.persist?.key;
  if (key && liveKeys.has(key)) {
    throw new Error(
      `Another viewer on this page keeps its state under '${key}': give each viewer its own persist key.`,
    );
  }
  const core = await createViewerCore({
    spec: o.spec,
    createPlugin: (spec) => new PluginContext(spec),
    modules: ['canvas', 'files', 'network', 'storage', 'page'],
  });
  if (key) {
    liveKeys.add(key);
    core.onDispose(() => {
      liveKeys.delete(key);
    });
  }
  const mounted = await mountShell(core, o.target);
  core.onDispose(() => mounted.dispose());
  (core.views as { shell?: unknown }).shell = mounted.shell;
  core.ready$.set(true);
  return core;
}
