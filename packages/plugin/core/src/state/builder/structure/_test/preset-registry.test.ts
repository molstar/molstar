/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { PluginConfig, PluginConfigManager } from '../../../../config.js';
import { TrajectoryHierarchyBuilder } from '../hierarchy.js';
import { StructureRepresentationBuilder } from '../representation.js';
import { PresetTrajectoryHierarchy } from '../hierarchy-presets/catalog.js';
import { PresetStructureRepresentations } from '../representation-presets/catalog.js';
import { TrajectoryHierarchyPresetProvider } from '../hierarchy-presets/types.js';
import { StructureRepresentationPresetProvider } from '../representation-presets/types.js';

function createPlugin() {
  return { config: new PluginConfigManager() } as unknown as PluginContext;
}

function hierarchyPreset(id: string, alias?: string) {
  return TrajectoryHierarchyPresetProvider({
    id,
    alias,
    display: { name: id },
    apply: async () => ({}),
  });
}

function representationPreset(id: string, alias?: string) {
  return StructureRepresentationPresetProvider({
    id,
    alias,
    display: { name: id },
    apply: async () => ({}),
  });
}

describe('TrajectoryHierarchyBuilder presets', () => {
  it('preloads the built-in catalog in order', () => {
    const builder = new TrajectoryHierarchyBuilder(createPlugin());
    expect(builder.providers).toEqual(Object.values(PresetTrajectoryHierarchy));
  });

  it('resolves by id, then alias', () => {
    const builder = new TrajectoryHierarchyBuilder(createPlugin());
    for (const p of Object.values(PresetTrajectoryHierarchy)) {
      expect(builder.resolveProvider(p.id)).toBe(p);
      expect(builder.resolveProvider(p.alias)).toBe(p);
      expect(builder.has(p.id)).toBe(true);
      expect(builder.has(p.alias)).toBe(true);
    }
    expect(builder.resolveProvider('default')).toBe(PresetTrajectoryHierarchy.default);
    expect(builder.resolveProvider('preset-trajectory-default')).toBe(PresetTrajectoryHierarchy.default);
    expect(builder.resolveProvider('nope')).toBeUndefined();
    expect(builder.has('nope')).toBe(false);

    const custom = hierarchyPreset('custom-id', 'custom-alias');
    expect(builder.resolveProvider(custom)).toBe(custom);
    builder.registerPreset(custom);
    expect(builder.resolveProvider('custom-id')).toBe(custom);
    expect(builder.resolveProvider('custom-alias')).toBe(custom);
  });

  it('throws for an unresolved string', () => {
    const builder = new TrajectoryHierarchyBuilder(createPlugin());
    expect(() => builder.applyPreset({} as any, 'nope')).toThrow("Preset 'nope' is not registered in this plugin");
  });

  it('counts registrations of the same object and removes at zero', () => {
    const builder = new TrajectoryHierarchyBuilder(createPlugin());
    const p = hierarchyPreset('a', 'a-alias');
    builder.registerPreset(p);
    builder.registerPreset(p);
    expect(builder.providers.filter((x) => x === p).length).toBe(1);

    builder.unregisterPreset(p);
    expect(builder.has('a')).toBe(true);
    expect(builder.has('a-alias')).toBe(true);

    builder.unregisterPreset('a');
    expect(builder.has('a')).toBe(false);
    expect(builder.has('a-alias')).toBe(false);
    expect(builder.providers.includes(p)).toBe(false);

    // unknown is a no-op
    builder.unregisterPreset(p);
    builder.unregisterPreset('a');
    builder.unregisterPreset('never-registered');
  });

  it('throws on a different object under an existing id or alias', () => {
    const builder = new TrajectoryHierarchyBuilder(createPlugin());
    const a = hierarchyPreset('a', 'a-alias');
    builder.registerPreset(a);

    expect(() => builder.registerPreset(hierarchyPreset('a'))).toThrow(/id 'a'/);
    expect(() => builder.registerPreset(hierarchyPreset('b', 'a-alias'))).toThrow(/alias 'a-alias'/);
    // alias equal to another preset's id
    expect(() => builder.registerPreset(hierarchyPreset('c', 'a'))).toThrow(/alias 'a'/);
    // id equal to another preset's alias
    expect(() => builder.registerPreset(hierarchyPreset('a-alias'))).toThrow(/id 'a-alias'/);
    // built-in alias
    expect(() => builder.registerPreset(hierarchyPreset('x', 'default'))).toThrow(/alias 'default'/);
    // failed registrations changed nothing
    expect(builder.has('b')).toBe(false);
    expect(builder.has('c')).toBe(false);
    expect(builder.has('x')).toBe(false);
  });

  it('findConflict reports without changing anything', () => {
    const builder = new TrajectoryHierarchyBuilder(createPlugin());
    const a = hierarchyPreset('a', 'a-alias');
    const count = builder.providers.length;
    expect(builder.findConflict(a)).toBeUndefined();
    expect(builder.providers.length).toBe(count);
    builder.registerPreset(a);
    expect(builder.findConflict(a)).toBeUndefined();

    const other = hierarchyPreset('b', 'a-alias');
    const message = builder.findConflict(other);
    expect(message).toBeDefined();
    expect(() => builder.registerPreset(other)).toThrow(message!);
    expect(builder.providers.length).toBe(count + 1);
  });

  it('getPresetSelect defaults to the configured preset when registered, otherwise the first option', () => {
    const plugin = createPlugin();
    const builder = new TrajectoryHierarchyBuilder(plugin);
    expect(builder.getPresetSelect().defaultValue).toBe('preset-trajectory-default');

    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, 'preset-trajectory-unitcell');
    expect(builder.getPresetSelect().defaultValue).toBe('preset-trajectory-unitcell');

    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, 'not-registered');
    expect(builder.getPresetSelect().defaultValue).toBe('preset-trajectory-default');
    expect(builder.getPresetSelect().options.map((o) => o[0])).toContain(builder.getPresetSelect().defaultValue);
  });
});

describe('StructureRepresentationBuilder presets', () => {
  it('preloads the built-in catalog in order', () => {
    const builder = new StructureRepresentationBuilder(createPlugin());
    expect(builder.providers).toEqual(Object.values(PresetStructureRepresentations));
  });

  it('resolves by id, then alias', () => {
    const builder = new StructureRepresentationBuilder(createPlugin());
    for (const p of Object.values(PresetStructureRepresentations)) {
      expect(builder.resolveProvider(p.id)).toBe(p);
      expect(builder.resolveProvider(p.alias)).toBe(p);
      expect(builder.has(p.id)).toBe(true);
      expect(builder.has(p.alias)).toBe(true);
    }
    expect(builder.resolveProvider('auto')).toBe(PresetStructureRepresentations.auto);
    expect(builder.resolveProvider('preset-structure-representation-auto')).toBe(PresetStructureRepresentations.auto);
    expect(builder.resolveProvider('nope')).toBeUndefined();
    const custom = representationPreset('custom');
    expect(builder.resolveProvider(custom)).toBe(custom);
  });

  it('throws for an unresolved string', () => {
    const builder = new StructureRepresentationBuilder(createPlugin());
    expect(() => builder.applyPreset({} as any, 'nope')).toThrow("Preset 'nope' is not registered in this plugin");
  });

  it('no longer resolves unregistered built-ins', () => {
    const builder = new StructureRepresentationBuilder(createPlugin());
    builder.unregisterPreset(PresetStructureRepresentations.auto);
    expect(builder.resolveProvider('auto')).toBeUndefined();
    expect(() => builder.applyPreset({} as any, 'auto')).toThrow("Preset 'auto' is not registered in this plugin");
    expect(builder.providers).not.toContain(PresetStructureRepresentations.auto);
  });

  it('counts registrations of the same object and removes at zero', () => {
    const builder = new StructureRepresentationBuilder(createPlugin());
    const p = representationPreset('a', 'a-alias');
    builder.registerPreset(p);
    builder.registerPreset(p);
    builder.unregisterPreset(p);
    expect(builder.resolveProvider('a-alias')).toBe(p);
    builder.unregisterPreset('a');
    expect(builder.resolveProvider('a')).toBeUndefined();
    expect(builder.resolveProvider('a-alias')).toBeUndefined();
    builder.unregisterPreset(p);
  });

  it('throws on id and alias conflicts and findConflict agrees', () => {
    const builder = new StructureRepresentationBuilder(createPlugin());
    const a = representationPreset('a', 'a-alias');
    builder.registerPreset(a);

    for (const other of [
      representationPreset('a'),
      representationPreset('b', 'a-alias'),
      representationPreset('c', 'a'),
      representationPreset('a-alias'),
      representationPreset('x', 'auto'),
    ]) {
      const message = builder.findConflict(other);
      expect(message).toBeDefined();
      expect(() => builder.registerPreset(other)).toThrow(message!);
    }
    expect(builder.findConflict(a)).toBeUndefined();
    expect(builder.has('b')).toBe(false);
    expect(builder.has('x')).toBe(false);
  });

  it('getPresetSelect defaults to the configured preset when registered, otherwise the first option', () => {
    const plugin = createPlugin();
    const builder = new StructureRepresentationBuilder(plugin);
    expect(builder.getPresetSelect().defaultValue).toBe('preset-structure-representation-auto');

    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, 'preset-structure-representation-mesoscale');
    expect(builder.getPresetSelect().defaultValue).toBe('preset-structure-representation-mesoscale');

    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, 'not-registered');
    expect(builder.getPresetSelect().defaultValue).toBe('preset-structure-representation-empty');
  });
});
