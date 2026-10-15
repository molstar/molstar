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
// as `initBehaviors()` does, custom properties first. A registration conflict throws before anything changes.

export async function applyExtensionPlugin(
  plugin: PluginContext,
  parts: Partial<PrismViewerPluginParts>,
): Promise<() => Promise<void>> {
  const undoRegistry = parts.registry?.length ? plugin.register(parts.registry) : () => {};
  for (const [item, value] of parts.config ?? []) plugin.config.set(item, value);
  const behaviors = parts.behaviors ?? [];
  const added: string[] = [];
  const apply = async (customProps: boolean) => {
    const tree = plugin.state.behaviors.build();
    let any = false;
    for (const b of behaviors) {
      const category = PluginBehavior.getCategoryId(b.transformer);
      if ((category === 'custom-props') !== customProps) continue;
      if (plugin.state.behaviors.transforms.has(b.transformer.id)) {
        throw new Error(`The behavior '${b.transformer.id}' is there already`);
      }
      tree.to(category).apply(b.transformer, b.defaultParams, { ref: b.transformer.id });
      added.push(b.transformer.id);
      any = true;
    }
    if (any) await plugin.runTask(plugin.state.behaviors.updateTree(tree, { doNotUpdateCurrent: true }));
  };
  try {
    await apply(true);
    await apply(false);
  } catch (e) {
    undoRegistry();
    throw e;
  }
  if (parts.canvas3d && plugin.canvas3d) plugin.setCanvas3DProps(parts.canvas3d);
  return async () => {
    if (added.length) {
      const tree = plugin.state.behaviors.build();
      for (const ref of added) tree.delete(ref);
      await plugin.runTask(plugin.state.behaviors.updateTree(tree));
    }
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
