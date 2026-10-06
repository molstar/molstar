/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginBehaviors } from '@molstar/plugin/behavior';
import { StructureFocusRepresentation } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { StructureFocusRepresentationId } from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation/id';
import { PluginContext } from '@molstar/plugin/context';
import { PluginSpec } from '@molstar/plugin/spec';
import type { PluginState } from '@molstar/plugin/state';
import { current, all } from '@molstar/plugin/state/queries/structure/basic';
import { StructureComponentManager } from '../structure/component.js';

async function createPlugin(behaviors: PluginSpec['behaviors'] = []) {
  const plugin = new PluginContext({ behaviors });
  await plugin.init();
  return plugin;
}

function changedInteractions(plugin: PluginContext) {
  const options = plugin.managers.structure.component.state.options;
  const interactions = PD.getDefaultValues(StructureComponentManager.OptionsParams.interactions.params);
  return {
    ...options,
    interactions: {
      ...interactions,
      contacts: { ...interactions.contacts, lineOfSightDistFactor: 2.5 },
    },
  };
}

describe('StructureComponentManager option handlers', () => {
  it('calls a registered handler when its option changes and not when it is unchanged', async () => {
    const plugin = await createPlugin();
    const manager = plugin.managers.structure.component;
    const handler = jest.fn();
    manager.registerOptionHandler('interactions', handler);

    const next = changedInteractions(plugin);
    await manager.setOptions(next);
    expect(handler).toHaveBeenCalledTimes(1);
    expect(handler).toHaveBeenCalledWith(next.interactions, manager);

    await manager.setOptions({ ...next, hydrogens: 'hide-all' });
    expect(handler).toHaveBeenCalledTimes(1);
    plugin.dispose();
  });

  it('removes the handler with the returned undo', async () => {
    const plugin = await createPlugin();
    const manager = plugin.managers.structure.component;
    const handler = jest.fn();
    const undo = manager.registerOptionHandler('interactions', handler);
    undo();
    undo();

    await manager.setOptions(changedInteractions(plugin));
    expect(handler).not.toHaveBeenCalled();
    plugin.dispose();
  });

  it('keeps other handlers when one is removed', async () => {
    const plugin = await createPlugin();
    const manager = plugin.managers.structure.component;
    const a = jest.fn();
    const b = jest.fn();
    const undoA = manager.registerOptionHandler('interactions', a);
    manager.registerOptionHandler('interactions', b);
    undoA();

    await manager.setOptions(changedInteractions(plugin));
    expect(a).not.toHaveBeenCalled();
    expect(b).toHaveBeenCalledTimes(1);
    plugin.dispose();
  });

  it('preserves the option and ignores it when no handler is registered', async () => {
    const plugin = await createPlugin();
    const manager = plugin.managers.structure.component;
    const next = changedInteractions(plugin);
    await expect(manager.setOptions(next)).resolves.toBeUndefined();
    expect(manager.state.options.interactions).toBe(next.interactions);
    expect(manager.state.options.interactions.contacts.lineOfSightDistFactor).toBe(2.5);
    plugin.dispose();
  });

  it('is registered by the Interactions behavior and removed when the behavior is unregistered', async () => {
    const plugin = await createPlugin([PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions)]);
    const handlers = () => (plugin.managers.structure.component as any).optionHandlers.get('interactions');
    expect(handlers()).toHaveLength(1);

    await plugin.state.behaviors.build().delete(PluginBehaviors.CustomProps.Interactions.id).commit();
    expect(handlers()).toBeUndefined();
    plugin.dispose();
  });

  it('does not register a handler without the Interactions behavior', async () => {
    const plugin = await createPlugin();
    expect((plugin.managers.structure.component as any).optionHandlers.get('interactions')).toBeUndefined();
    plugin.dispose();
  });
});

describe('StructureComponentManager options snapshot', () => {
  it('round-trips options.interactions through snapshot JSON without the Interactions behavior', async () => {
    const source = await createPlugin();
    const next = changedInteractions(source);
    await source.managers.structure.component.setOptions(next);

    const snapshot = source.state.getSnapshot({
      data: false,
      behavior: false,
      animation: false,
      startAnimation: false,
      camera: false,
      canvas3d: false,
      interactivity: false,
      structureSelection: false,
      canvas3dContext: false,
      componentManager: true,
    } as unknown as PluginState.SnapshotParams);
    const json = JSON.parse(JSON.stringify(snapshot));
    expect(json.structureComponentManager.options.interactions).toEqual(next.interactions);

    const target = await createPlugin();
    await target.state.setSnapshot(json);
    expect(target.managers.structure.component.state.options.interactions).toEqual(next.interactions);
    source.dispose();
    target.dispose();
  });
});

describe('StructureComponentManager.setOptions and the focus representation behavior', () => {
  it('does not insert the focus representation behavior when it is absent', async () => {
    const plugin = await createPlugin();
    expect(plugin.state.hasBehavior(StructureFocusRepresentationId)).toBe(false);

    await plugin.managers.structure.component.setOptions({
      ...plugin.managers.structure.component.state.options,
      ignoreLight: true,
    });
    expect(plugin.state.hasBehavior(StructureFocusRepresentationId)).toBe(false);
    plugin.dispose();
  });

  it('the id module holds the transformer id of the behavior', () => {
    expect(StructureFocusRepresentation.id).toBe(StructureFocusRepresentationId);
  });

  it('updates the focus representation behavior when it is present', async () => {
    const plugin = await createPlugin([PluginSpec.Behavior(StructureFocusRepresentation)]);
    expect(plugin.state.hasBehavior(StructureFocusRepresentationId)).toBe(true);

    await plugin.managers.structure.component.setOptions({
      ...plugin.managers.structure.component.state.options,
      ignoreLight: true,
      hydrogens: 'hide-all',
    });
    const params = plugin.state.behaviors.cells.get(StructureFocusRepresentationId)!.params!.values;
    expect(params.ignoreLight).toBe(true);
    expect(params.ignoreHydrogens).toBe(true);
    plugin.dispose();
  });
});

describe('StructureComponentManager default selection query', () => {
  it('picks current by identity, not by position', async () => {
    const plugin = await createPlugin();
    const { options } = plugin.query.structure.registry;
    const index = options.findIndex((o) => o[0] === current);
    expect(index).toBeGreaterThanOrEqual(0);
    expect(StructureComponentManager.getDefaultSelectionQuery(options)).toBe(current);

    // current is no longer the second option
    const reordered = [options[index], ...options.filter((_, i) => i !== index)];
    expect(reordered[1][0]).not.toBe(current);
    expect(StructureComponentManager.getDefaultSelectionQuery(reordered)).toBe(current);

    const params = StructureComponentManager.getAddParams(plugin);
    expect(PD.getDefaultValues(params).selection).toBe(current);
    plugin.dispose();
  });

  it('picks the first option when current is not registered', async () => {
    const plugin = await createPlugin();
    plugin.query.structure.registry.remove(current);
    const { options } = plugin.query.structure.registry;
    expect(options[0][0]).toBe(all);

    expect(StructureComponentManager.getDefaultSelectionQuery(options)).toBe(all);
    expect(PD.getDefaultValues(StructureComponentManager.getAddParams(plugin)).selection).toBe(all);
    expect(PD.getDefaultValues(StructureComponentManager.getThemeParams(plugin, undefined)).selection).toBe(all);
    plugin.dispose();
  });

  it('returns an empty select for an empty registry', async () => {
    const plugin = await createPlugin();
    const registry = plugin.query.structure.registry;
    for (const q of [...registry.list]) registry.remove(q);
    expect(registry.options).toHaveLength(0);

    expect(StructureComponentManager.getDefaultSelectionQuery(registry.options)).toBeUndefined();
    const add = StructureComponentManager.getAddParams(plugin);
    expect(add.selection.options).toEqual([]);
    expect(() => PD.getDefaultValues(add)).not.toThrow();
    const theme = StructureComponentManager.getThemeParams(plugin, undefined);
    expect(theme.selection.options).toEqual([]);
    plugin.dispose();
  });
});
