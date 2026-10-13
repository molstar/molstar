/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '@molstar/model/model/structure';
import { StructureFocusRepresentation } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { StructureFocusRepresentationId } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation/id';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import { StructureRepresentationPresetProvider } from '../representation-presets/types.js';

const { updateFocusRepr } = StructureRepresentationPresetProvider;

async function createPlugin(withFocus = true) {
  const spec = DefaultPluginSpec();
  const plugin = new PluginContext({
    registry: spec.registry,
    behaviors: withFocus ? [PluginSpec.Behavior(StructureFocusRepresentation)] : [],
  });
  await plugin.init();
  return plugin;
}

/** Providers are registered by identity counts; remove until the provider is gone. */
function purge(registry: { has(p: any): boolean; remove(p: any): void }, provider: any) {
  while (registry.has(provider)) registry.remove(provider);
}

function focusParams(plugin: PluginContext) {
  return plugin.state.behaviors.cells.get(StructureFocusRepresentationId)!.params!.values;
}

describe('updateFocusRepr', () => {
  it('sets the requested theme on the target and surroundings', async () => {
    const plugin = await createPlugin();
    await updateFocusRepr(plugin, Structure.Empty, 'chain-id', undefined);
    expect(focusParams(plugin).targetParams.colorTheme.name).toBe('chain-id');
    expect(focusParams(plugin).surroundingsParams.colorTheme.name).toBe('chain-id');
    plugin.dispose();
  });

  it('uses the default color theme of the representation type without a requested theme', async () => {
    const plugin = await createPlugin();
    await updateFocusRepr(plugin, Structure.Empty, 'chain-id', undefined);
    await updateFocusRepr(plugin, Structure.Empty, undefined, undefined);
    const defaultTheme = plugin.representation.structure.registry.get('ball-and-stick').defaultColorTheme.name;
    expect(defaultTheme).toBe('element-symbol');
    expect(focusParams(plugin).targetParams.colorTheme.name).toBe(defaultTheme);
    expect(focusParams(plugin).surroundingsParams.colorTheme.name).toBe(defaultTheme);
    plugin.dispose();
  });

  it('does nothing and does not insert the behavior when it is absent', async () => {
    const plugin = await createPlugin(false);
    expect(updateFocusRepr(plugin, Structure.Empty, 'chain-id', undefined)).toBeUndefined();
    expect(plugin.state.hasBehavior(StructureFocusRepresentationId)).toBe(false);
    plugin.dispose();
  });

  it('skips an unregistered theme without a warning', async () => {
    const plugin = await createPlugin();
    const warn = jest.spyOn(plugin.log, 'warn');
    const before = focusParams(plugin);
    expect(updateFocusRepr(plugin, Structure.Empty, 'no-such-theme' as any, undefined)).toBeUndefined();
    expect(focusParams(plugin)).toBe(before);
    expect(warn).not.toHaveBeenCalled();
    plugin.dispose();
  });

  it('skips representation types that are not registered', async () => {
    const plugin = await createPlugin();
    const { registry } = plugin.representation.structure;
    purge(registry, registry.get('ball-and-stick'));
    expect(registry.has('ball-and-stick')).toBe(false);

    const warn = jest.spyOn(plugin.log, 'warn');
    const before = focusParams(plugin);
    expect(updateFocusRepr(plugin, Structure.Empty, 'chain-id', undefined)).toBeUndefined();
    expect(focusParams(plugin)).toBe(before);
    expect(warn).not.toHaveBeenCalled();
    plugin.dispose();
  });

  it('skips only the representation whose type is not registered', async () => {
    const plugin = await createPlugin();
    const { registry } = plugin.representation.structure;
    const spacefill = registry.get('spacefill');
    await plugin.state.updateBehavior(StructureFocusRepresentationId, (p: any) => {
      p.targetParams.type.name = 'spacefill';
    });
    purge(registry, spacefill);
    await updateFocusRepr(plugin, Structure.Empty, 'chain-id', undefined);
    expect(focusParams(plugin).surroundingsParams.colorTheme.name).toBe('chain-id');
    expect(focusParams(plugin).targetParams.colorTheme.name).not.toBe('chain-id');
    plugin.dispose();
  });

  it('does not throw when nothing is registered', async () => {
    const plugin = await createPlugin();
    const { registry, themes } = plugin.representation.structure;
    registry.clear();
    themes.colorThemeRegistry.clear();
    expect(() => updateFocusRepr(plugin, Structure.Empty, undefined, undefined)).not.toThrow();
    expect(updateFocusRepr(plugin, Structure.Empty, undefined, undefined)).toBeUndefined();
    plugin.dispose();
  });
});

describe('structure focus representation behavior without its default providers', () => {
  it('builds default params that name unregistered representations and themes', async () => {
    const plugin = await createPlugin(false);
    const { registry, themes } = plugin.representation.structure;
    purge(registry, registry.get('ball-and-stick'));
    purge(themes.sizeThemeRegistry, themes.sizeThemeRegistry.get('physical'));

    const warn = jest.spyOn(plugin.log, 'warn');
    let params: any;
    expect(() => {
      params = StructureFocusRepresentation.createDefaultParams(void 0 as any, plugin);
    }).not.toThrow();
    expect(params.targetParams).toBeDefined();
    expect(params.surroundingsParams).toBeDefined();
    // lenient path: unregistered names are reported, then the registry defaults are used
    expect(warn).toHaveBeenCalled();
    plugin.dispose();
  });

  it('builds default params for a plugin without any structure representation', async () => {
    const plugin = await createPlugin(false);
    const { registry, themes } = plugin.representation.structure;
    registry.clear();
    themes.colorThemeRegistry.clear();
    themes.sizeThemeRegistry.clear();
    expect(() => StructureFocusRepresentation.createDefaultParams(void 0 as any, plugin)).not.toThrow();
    plugin.dispose();
  });

  it('builds the interactions defaults without the Interactions providers registered', async () => {
    const plugin = await createPlugin(false);
    // the default spec registry does not register the Interactions representation, only the behavior does
    const params = StructureFocusRepresentation.createDefaultParams(void 0 as any, plugin);
    expect(params.nciParams.type.name).toBe('interactions');
    expect(plugin.representation.structure.registry.has('interactions')).toBe(false);
    plugin.dispose();
  });
});
