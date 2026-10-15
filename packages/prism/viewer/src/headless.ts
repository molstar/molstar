/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { PluginSpec } from '@molstar/plugin/spec';
import type { ApiModuleKey } from './api/area.js';
import { createDefaultPlugin } from './sources/create-plugin.js';
import type { PrismViewerSpec } from './spec.js';
import { createViewerCore, type PrismViewer } from './viewer.js';

// `@molstar/prism-viewer/headless` (viewer-spec.md §6, README §3.3): the same viewer without a DOM: the same API, tools
// and state account. With `plugin`, a `HeadlessPluginContext` (images) is made by the caller; without, a plain plugin
// with no canvas, and canvas operations say why they cannot run.

export interface HeadlessViewerOptions {
  readonly spec: PrismViewerSpec;
  /** Makes the plugin (a `HeadlessPluginContext` from `@molstar/plugin-headless` for images); plain by default. */
  readonly plugin?: (spec: PluginSpec) => PluginContext;
  /** What the host offers; by default files and the network, with a canvas when `plugin` is given. */
  readonly modules?: readonly ApiModuleKey[];
}

export async function createPrismViewerHeadless(o: HeadlessViewerOptions): Promise<PrismViewer> {
  const core = await createViewerCore({
    spec: o.spec,
    createPlugin: o.plugin ?? createDefaultPlugin,
    modules: o.modules ?? (o.plugin ? ['canvas', 'files', 'network'] : ['files', 'network']),
  });
  core.ready$.set(true);
  return core;
}

export type { PrismViewer } from './viewer.js';
