/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from '@molstar/core/state';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import type { PluginContext } from '@molstar/plugin/context';
import type { PrismViewerPluginParts } from '../spec.js';

// A lazy extension's plugin part on the running plugin (README §3.4, "Lazy extensions"): its registry entries
// registered (their undo kept), its config set (later wins, as at creation), its behaviors applied to the behavior tree
// as `initBehaviors()` does, custom properties first. A conflict (a registration, a behavior that is there) throws with
// nothing left changed: the behaviors are checked first, and what was done is undone when the apply fails. The undo
// puts back the config it set (unless something set it after).

export async function applyExtensionPlugin(
  plugin: PluginContext,
  parts: Partial<PrismViewerPluginParts>,
): Promise<() => Promise<void>> {
  const behaviors = parts.behaviors ?? [];
  for (const b of behaviors) {
    if (plugin.state.behaviors.transforms.has(b.transformer.id)) {
      throw new Error(`The behavior '${b.transformer.id}' is there already`);
    }
  }
  const undoRegistry = parts.registry?.length ? plugin.register(parts.registry) : () => {};
  const config = (parts.config ?? []).map(([item, value]) => ({ item, value, before: plugin.config.get(item) }));
  for (const c of config) plugin.config.set(c.item, c.value);
  const undoConfig = () => {
    for (let i = config.length - 1; i >= 0; i--) {
      const c = config[i];
      if (plugin.config.get(c.item) === c.value) plugin.config.set(c.item, c.before);
    }
  };
  const added: string[] = [];
  const apply = async (customProps: boolean) => {
    const tree = plugin.state.behaviors.build();
    let any = false;
    for (const b of behaviors) {
      const category = PluginBehavior.getCategoryId(b.transformer);
      if ((category === 'custom-props') !== customProps) continue;
      tree.to(category).apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
      added.push(b.transformer.id);
      any = true;
    }
    if (any) await plugin.runTask(plugin.state.behaviors.updateTree(tree, { doNotUpdateCurrent: true }));
  };
  const removeAdded = async () => {
    const present = added.filter((ref) => plugin.state.behaviors.transforms.has(ref));
    added.length = 0;
    if (!present.length) return;
    const tree = plugin.state.behaviors.build();
    for (const ref of present) tree.delete(ref);
    await plugin.runTask(plugin.state.behaviors.updateTree(tree));
  };
  try {
    await apply(true);
    await apply(false);
  } catch (e) {
    try {
      await removeAdded();
    } catch {
      // The failure to report is the first.
    }
    undoConfig();
    undoRegistry();
    throw e;
  }
  if (parts.canvas3d && plugin.canvas3d) plugin.setCanvas3DProps(parts.canvas3d);
  return async () => {
    await removeAdded();
    undoConfig();
    undoRegistry();
  };
}

/** Which of these transformer ids, format names or file extensions the plugin does not know (a lazy extension's). */
export function unknownToPlugin(
  plugin: PluginContext,
  need: { readonly transformers?: readonly string[]; readonly formats?: readonly string[] },
): { transformers: string[]; formats: string[] } {
  const transformers = (need.transformers ?? []).filter((id) => !hasTransformer(id));
  const formats = (need.formats ?? []).filter((name) => !plugin.dataFormats.get(name));
  return { transformers, formats };
}

function hasTransformer(id: string): boolean {
  return StateTransformer.has(id);
}
