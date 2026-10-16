/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import type { PluginSpec } from '@molstar/plugin/spec';

/** The viewer's plugin when the host gives none: a plain `PluginContext` (no canvas until the page mounts one). */
export function createDefaultPlugin(spec: PluginSpec): PluginContext {
  return new PluginContext(spec);
}
