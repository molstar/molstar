/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import {
  DefaultActions,
  DefaultAnimations,
  DefaultDragAndDrop,
  DefaultFormats,
  DefaultMarkdownExtensions,
  DefaultPresets,
  DefaultRegistry,
  DefaultSelectionQueries,
  DefaultStructureRepresentations,
  DefaultThemes,
  DefaultVolumeRepresentations,
} from '@molstar/plugin/default-registry';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ExternalColorThemes } from '@molstar/plugin/themes/external';

describe('DefaultRegistry', () => {
  it('lists the ten named entries in order', () => {
    expect(DefaultRegistry).toEqual([
      DefaultThemes,
      DefaultStructureRepresentations,
      DefaultVolumeRepresentations,
      DefaultActions,
      DefaultFormats,
      DefaultPresets,
      DefaultSelectionQueries,
      DefaultMarkdownExtensions,
      DefaultDragAndDrop,
      DefaultAnimations,
    ]);
    expect(DefaultRegistry.length).toBe(new Set(DefaultRegistry).size);
  });

  it('keeps the format order and registers the open-files drag and drop fallback', () => {
    const names = DefaultFormats.formats!.map((f) => f.name);
    // volume, topology, coordinates, shape, trajectory (particles: `ParticlesRegistry`)
    expect(names.slice(0, 2)).toEqual(['ccp4', 'dsn6']);
    expect(names.indexOf('psf')).toBeGreaterThan(names.indexOf('mtz'));
    expect(names.indexOf('ply')).toBeGreaterThan(names.indexOf('psf'));
    expect(names.indexOf('mmcif')).toBeGreaterThan(names.indexOf('ply'));
    expect(names.filter((n) => n.includes('particles'))).toEqual([]);
    expect(new Set(names).size).toBe(names.length);
    expect(DefaultDragAndDrop.dragAndDrop!.map((e) => [e.name, e.fallback])).toEqual([['open-files', true]]);
  });

  it('lists the default preset counts', () => {
    expect(DefaultPresets.structure!.presets!.hierarchy!.length).toBe(5);
    expect(DefaultPresets.structure!.presets!.representation!.length).toBe(11);
  });

  it('DefaultThemes includes the external color themes in the structure and volume scopes, no particle themes', () => {
    expect(DefaultThemes.particles).toBeUndefined();
    for (const scope of ['structure', 'volume'] as const) {
      const names = DefaultThemes[scope]!.themes!.color!.map((t) => t.name);
      expect(names).toEqual(expect.arrayContaining(['external-structure', 'external-volume']));
      expect(names.filter((n) => n.startsWith('particle-'))).toEqual([]);
      const sizes = DefaultThemes[scope]!.themes!.size!.map((t) => t.name);
      expect(sizes.length).toBeGreaterThan(0);
      expect(sizes).not.toContain('particle-size');
    }
  });

  it('DefaultPluginSpec lists the registry and no removed spec fields', () => {
    const spec = DefaultPluginSpec();
    expect(spec.registry).toBe(DefaultRegistry);
    expect(Object.keys(spec).sort()).toEqual(['behaviors', 'registry']);
  });

  it('registers into a plugin with the default spec', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();

    expect(plugin.representation.structure.registry.default?.name).toBe('cartoon');
    expect(plugin.representation.volume.registry.default?.name).toBe('direct-volume');
    // particles are opt-in (`ParticlesRegistry`)
    expect(plugin.representation.particles.registry.list).toEqual([]);
    expect(plugin.dataFormats.has('relion_star_particles')).toBe(false);
    expect(plugin.managers.animation.animations[0].name).toBe('built-in.animate-model-index');
    for (const scope of ['structure', 'volume'] as const) {
      const themes = plugin.representation[scope].themes;
      expect(themes.colorThemeRegistry.has('external-structure')).toBe(true);
      expect(themes.colorThemeRegistry.has('external-volume')).toBe(true);
    }

    // Registering the default registry again only counts: it neither throws nor duplicates anything.
    const formats = plugin.dataFormats.list.length;
    const queries = plugin.query.structure.registry.list.length;
    const undo = plugin.register(DefaultRegistry);
    expect(plugin.dataFormats.list.length).toBe(formats);
    expect(plugin.query.structure.registry.list.length).toBe(queries);
    undo();
    expect(plugin.dataFormats.list.length).toBe(formats);
    expect(plugin.query.structure.registry.list.length).toBe(queries);
    plugin.dispose();
  });
});

describe('ExternalColorThemes', () => {
  it('lists both providers under the color themes of all three scopes', () => {
    for (const scope of ['structure', 'volume', 'particles'] as const) {
      expect(ExternalColorThemes[scope]!.themes!.color!.map((t) => t.name)).toEqual([
        'external-structure',
        'external-volume',
      ]);
    }
  });

  it('registers the themes in a plugin that lists the entry', async () => {
    const plugin = new PluginContext({ behaviors: [], registry: [ExternalColorThemes] });
    await plugin.init();
    for (const scope of ['structure', 'volume', 'particles'] as const) {
      expect(plugin.representation[scope].themes.colorThemeRegistry.has('external-volume')).toBe(true);
    }
    plugin.dispose();
  });
});
