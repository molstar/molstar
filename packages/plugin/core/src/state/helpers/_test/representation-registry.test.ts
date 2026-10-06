/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '@molstar/model/model/structure';
import { Volume } from '@molstar/model/model/volume';
import { createFakePlugin } from '../../_test/fake-plugin.js';
import {
  createStructureColorThemeParams,
  createStructureRepresentationParams,
  createStructureSizeThemeParams,
} from '../structure-representation-params.js';
import {
  createVolumeColorThemeParams,
  createVolumeRepresentationParams,
  createVolumeSizeThemeParams,
} from '../volume-representation-params.js';

const empty = { name: '', params: {} };

describe('representation param helpers with unregistered names', () => {
  it('warn for an unregistered structure representation and still return params', () => {
    const { plugin, warnings } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, { type: 'no-such-repr' as any });
    expect(params.type.name).toBeDefined();
    expect(warnings).toEqual([
      "Structure representation 'no-such-repr' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('warn for unregistered structure color and size theme names', () => {
    const { plugin, warnings } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, {
      type: 'cartoon',
      color: 'no-such-color' as any,
      size: 'no-such-size' as any,
    });
    expect(params.type.name).toBe('cartoon');
    expect(warnings).toEqual([
      "Structure color theme 'no-such-color' is not registered in this plugin; the registry default is used",
      "Structure size theme 'no-such-size' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('warn for an unregistered default theme of the representation', () => {
    const { plugin, warnings } = createFakePlugin();
    const { colorThemeRegistry } = plugin.representation.structure.themes;
    const defaultColor = plugin.representation.structure.registry.get('cartoon').defaultColorTheme.name;
    const provider = colorThemeRegistry.get(defaultColor);
    colorThemeRegistry.remove(provider);
    createStructureRepresentationParams(plugin, Structure.Empty, { type: 'cartoon' });
    expect(warnings).toEqual([
      `Structure color theme '${defaultColor}' is not registered in this plugin; the registry default is used`,
    ]);
  });

  it('do not warn for registered names', () => {
    const { plugin, warnings } = createFakePlugin();
    createStructureRepresentationParams(plugin, Structure.Empty, { type: 'cartoon', color: 'chain-id' });
    createStructureRepresentationParams(plugin, Structure.Empty);
    createStructureColorThemeParams(plugin, Structure.Empty, 'cartoon', 'chain-id');
    createStructureSizeThemeParams(plugin, Structure.Empty, 'cartoon', 'uniform');
    createVolumeRepresentationParams(plugin, Volume.One, { type: 'isosurface' });
    expect(warnings).toEqual([]);
  });

  it('use provider objects directly without warnings', () => {
    const { plugin, warnings } = createFakePlugin();
    const provider = plugin.representation.structure.registry.get('cartoon');
    const params = createStructureRepresentationParams(plugin, Structure.Empty, { type: provider });
    expect(params.type.name).toBe('cartoon');
    expect(warnings).toEqual([]);
  });

  it('warn in the structure theme helpers', () => {
    const { plugin, warnings } = createFakePlugin();
    const color = createStructureColorThemeParams(plugin, Structure.Empty, 'no-such-repr', 'no-such-color');
    const size = createStructureSizeThemeParams(plugin, Structure.Empty, 'cartoon', 'no-such-size');
    expect(color).toBeDefined();
    expect(size).toBeDefined();
    expect(warnings).toEqual([
      "Structure representation 'no-such-repr' is not registered in this plugin; the registry default is used",
      "Structure color theme 'no-such-color' is not registered in this plugin; the registry default is used",
      "Structure size theme 'no-such-size' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('warn in the volume helpers', () => {
    const { plugin, warnings } = createFakePlugin();
    createVolumeRepresentationParams(plugin, Volume.One, { type: 'no-such-repr' as any });
    createVolumeColorThemeParams(plugin, Volume.One, 'isosurface', 'no-such-color');
    createVolumeSizeThemeParams(plugin, Volume.One, 'isosurface', 'no-such-size');
    expect(warnings).toEqual([
      "Volume representation 'no-such-repr' is not registered in this plugin; the registry default is used",
      "Volume color theme 'no-such-color' is not registered in this plugin; the registry default is used",
      "Volume size theme 'no-such-size' is not registered in this plugin; the registry default is used",
    ]);
  });
});

describe('representation param helpers with empty registries', () => {
  it('throw with data when representations are missing', () => {
    const { plugin } = createFakePlugin(['representations']);
    const message = 'No structure representations are registered in this plugin';
    expect(() => createStructureRepresentationParams(plugin, Structure.Empty)).toThrow(message);
    expect(() => createStructureColorThemeParams(plugin, Structure.Empty, undefined)).toThrow(message);
    expect(() => createStructureSizeThemeParams(plugin, Structure.Empty, undefined)).toThrow(message);
    const volumeMessage = 'No volume representations are registered in this plugin';
    expect(() => createVolumeRepresentationParams(plugin, Volume.One)).toThrow(volumeMessage);
    expect(() => createVolumeColorThemeParams(plugin, Volume.One, undefined)).toThrow(volumeMessage);
    expect(() => createVolumeSizeThemeParams(plugin, Volume.One, undefined)).toThrow(volumeMessage);
  });

  it('throw with data when themes are missing', () => {
    const colors = createFakePlugin(['color-themes']).plugin;
    expect(() => createStructureRepresentationParams(colors, Structure.Empty)).toThrow(
      'No structure color themes are registered in this plugin',
    );
    expect(() => createVolumeRepresentationParams(colors, Volume.One)).toThrow(
      'No volume color themes are registered in this plugin',
    );
    expect(() => createStructureColorThemeParams(colors, Structure.Empty, undefined)).toThrow(
      'No structure color themes are registered in this plugin',
    );
    // the size theme helper does not need color themes
    expect(() => createStructureSizeThemeParams(colors, Structure.Empty, undefined)).not.toThrow();

    const sizes = createFakePlugin(['size-themes']).plugin;
    expect(() => createStructureRepresentationParams(sizes, Structure.Empty)).toThrow(
      'No structure size themes are registered in this plugin',
    );
    expect(() => createVolumeSizeThemeParams(sizes, Volume.One, undefined)).toThrow(
      'No volume size themes are registered in this plugin',
    );
  });

  it('return empty mapped params without data', () => {
    for (const kind of ['representations', 'color-themes', 'size-themes'] as const) {
      const { plugin, warnings } = createFakePlugin([kind]);
      expect(createStructureRepresentationParams(plugin)).toEqual({
        type: empty,
        colorTheme: empty,
        sizeTheme: empty,
      });
      expect(createVolumeRepresentationParams(plugin)).toEqual({ type: empty, colorTheme: empty, sizeTheme: empty });
      expect(warnings).toEqual([]);
    }

    const { plugin } = createFakePlugin(['representations']);
    expect(createStructureColorThemeParams(plugin, undefined, undefined)).toEqual(empty);
    expect(createStructureSizeThemeParams(plugin, undefined, undefined)).toEqual(empty);
    expect(createVolumeColorThemeParams(plugin, undefined, undefined)).toEqual(empty);
    expect(createVolumeSizeThemeParams(plugin, undefined, undefined)).toEqual(empty);
  });

  it('return regular params without data when the registries are populated', () => {
    const { plugin } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin);
    expect(params.type.name).toBe('cartoon');
    expect(params.colorTheme.name).not.toBe('');
  });
});
