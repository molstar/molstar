/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PresetTrajectoryHierarchy } from '../hierarchy-presets/catalog.js';
import { PresetStructureRepresentations } from '../representation-presets/catalog.js';

describe('built-in preset short keys in a plugin with the default spec', () => {
  it('resolve through their aliases to the preset with the matching id', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    const { hierarchy, representation } = plugin.builders.structure;

    for (const [key, preset] of Object.entries(PresetTrajectoryHierarchy)) {
      expect(hierarchy.has(key)).toBe(true);
      expect(hierarchy.has(preset.id)).toBe(true);
      expect(hierarchy.resolveProvider(key)).toBe(preset);
      expect(hierarchy.resolveProvider(preset.id)).toBe(preset);
    }
    for (const [key, preset] of Object.entries(PresetStructureRepresentations)) {
      expect(representation.has(key)).toBe(true);
      expect(representation.has(preset.id)).toBe(true);
      expect(representation.resolveProvider(key)).toBe(preset);
      expect(representation.resolveProvider(preset.id)).toBe(preset);
    }
    plugin.dispose();
  });

  it('are the only names the default spec registers besides the ids', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    const { hierarchy, representation } = plugin.builders.structure;
    expect(hierarchy.providers.map((p) => p.alias).sort()).toEqual(Object.keys(PresetTrajectoryHierarchy).sort());
    expect(representation.providers.map((p) => p.alias).sort()).toEqual(
      Object.keys(PresetStructureRepresentations).sort(),
    );
    plugin.dispose();
  });

  it('give a clear error for an unknown id or alias', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    const { hierarchy, representation } = plugin.builders.structure;
    expect(hierarchy.has('no-such-preset')).toBe(false);
    expect(hierarchy.resolveProvider('no-such-preset')).toBeUndefined();
    expect(representation.resolveProvider('no-such-preset')).toBeUndefined();
    const message = "Preset 'no-such-preset' is not registered in this plugin";
    expect(() => hierarchy.applyPreset({} as any, 'no-such-preset')).toThrow(message);
    expect(() => representation.applyPreset({} as any, 'no-such-preset')).toThrow(message);
    plugin.dispose();
  });
});
