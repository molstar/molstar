/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell } from '@molstar/core/state';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject } from '../objects.js';

export function parseArray(input?: string): string[] {
  return (
    input
      ?.split(',')
      .map((s) => s.trim())
      .filter((s) => s.length > 0) ?? []
  );
}

export function findRepresentations(plugin: PluginContext, cells: StateObjectCell[]): StateObjectCell[] {
  if (!cells.length) return [];
  return plugin.state.data.selectQ((q) =>
    q
      .byValue(...cells)
      .subtree()
      .filter((c) => PluginStateObject.isRepresentation3D(c.obj)),
  );
}
