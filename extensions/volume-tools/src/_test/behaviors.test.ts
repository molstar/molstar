/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { StateTransformer } from '@molstar/core/state';
import { VolumeMaskBehavior } from '@molstar/volume-tools-extension/mask/behavior';
import { VolumeSegmentorBehavior } from '@molstar/volume-tools-extension/segmentor/behavior';
import { MaskSelectionColorThemeProvider } from '@molstar/volume-tools-extension/mask/theme';
import { BodyLabelColorThemeProvider } from '@molstar/volume-tools-extension/segmentor/theme';
import { VolumeSegmentorManager } from '@molstar/volume-tools-extension/segmentor/manager';

async function createPlugin() {
  const defaults = DefaultPluginSpec();
  const behaviors = [VolumeMaskBehavior, VolumeSegmentorBehavior].map((transformer) => ({
    transformer,
    defaultParams: {},
  }));
  const plugin = new PluginContext({ ...defaults, behaviors: [...defaults.behaviors, ...behaviors] } as any);
  await plugin.init();
  return plugin;
}

const colorThemes = (plugin: PluginContext) =>
  plugin.representation.volume.themes.colorThemeRegistry.list.map((t) => t.name);

async function toggle(plugin: PluginContext, behavior: StateTransformer, enabled: boolean) {
  if (enabled) await plugin.state.updateBehavior(behavior, () => {});
  else await plugin.runTask(plugin.state.behaviors.updateTree(plugin.state.behaviors.build().delete(behavior.id)));
}

describe('volume tools behaviors', () => {
  it('register their volume color themes through one entry and remove them when toggled off', async () => {
    const plugin = await createPlugin();
    const initial = colorThemes(plugin);
    expect(initial).toContain(MaskSelectionColorThemeProvider.name);
    expect(initial).toContain(BodyLabelColorThemeProvider.name);
    expect(VolumeSegmentorManager.get(plugin)).toBeDefined();

    await toggle(plugin, VolumeMaskBehavior, false);
    expect(colorThemes(plugin)).not.toContain(MaskSelectionColorThemeProvider.name);
    expect(colorThemes(plugin)).toContain(BodyLabelColorThemeProvider.name);

    await toggle(plugin, VolumeSegmentorBehavior, false);
    expect(colorThemes(plugin)).toEqual(
      initial.filter((n) => n !== MaskSelectionColorThemeProvider.name && n !== BodyLabelColorThemeProvider.name),
    );
    expect(VolumeSegmentorManager.get(plugin)).toBeUndefined();

    await toggle(plugin, VolumeMaskBehavior, true);
    await toggle(plugin, VolumeSegmentorBehavior, true);
    expect([...colorThemes(plugin)].sort()).toEqual([...initial].sort());
    expect(VolumeSegmentorManager.get(plugin)).toBeDefined();
    plugin.dispose();
  });
});
