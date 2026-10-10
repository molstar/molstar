/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '@molstar/model/model/structure';
import { Volume } from '@molstar/model/model/volume';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { ElementSymbolColorThemeProvider } from '@molstar/graphics/theme/color/element-symbol';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';
import { createFakePlugin } from '../../_test/fake-plugin.js';
import { createStructureRepresentationParams } from '../structure-representation-params.js';
import { createVolumeRepresentationParams } from '../volume-representation-params.js';

const NotRegistered = /not registered/;

describe('mixed name and provider props', () => {
  it('resolve a provider type with a theme name', () => {
    const { plugin, warnings } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, {
      type: CartoonRepresentationProvider,
      color: 'element-symbol',
    });
    expect(params.type.name).toBe('cartoon');
    expect(params.colorTheme.name).toBe('element-symbol');
    // the size theme is the representation's default one
    expect(params.sizeTheme.name).toBe(CartoonRepresentationProvider.defaultSizeTheme.name);
    expect(warnings).toEqual([]);
  });

  it('resolve a type name with provider themes', () => {
    const { plugin, warnings } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, {
      type: 'ball-and-stick',
      color: ChainIdColorThemeProvider,
      size: UniformSizeThemeProvider,
    });
    expect(params.type.name).toBe('ball-and-stick');
    expect(params.colorTheme.name).toBe('chain-id');
    expect(params.sizeTheme.name).toBe('uniform');
    expect(warnings).toEqual([]);
  });

  it('resolve each field independently', () => {
    const { plugin, warnings } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, {
      type: BallAndStickRepresentationProvider,
      color: ChainIdColorThemeProvider,
      size: 'uniform',
    });
    expect(params.type.name).toBe('ball-and-stick');
    expect(params.colorTheme.name).toBe('chain-id');
    expect(params.sizeTheme.name).toBe('uniform');
    expect(warnings).toEqual([]);
  });

  it('use the registry default for fields that are not given', () => {
    const { plugin } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, {
      color: ElementSymbolColorThemeProvider,
    });
    expect(params.type.name).toBe('cartoon');
    expect(params.colorTheme.name).toBe('element-symbol');
    expect(params.sizeTheme.name).toBe(CartoonRepresentationProvider.defaultSizeTheme.name);
  });

  it('give a theme name and the theme provider the same params', () => {
    const { plugin } = createFakePlugin();
    const type = BallAndStickRepresentationProvider;
    const byName = createStructureRepresentationParams(plugin, Structure.Empty, {
      type,
      color: type.defaultColorTheme.name,
    });
    const byProvider = createStructureRepresentationParams(plugin, Structure.Empty, {
      type,
      color: ElementSymbolColorThemeProvider,
    });
    expect(byName).toEqual(byProvider);
    expect(byName).toEqual(createStructureRepresentationParams(plugin, Structure.Empty, { type }));
  });

  it('merge the params of each field over the defaults of that field', () => {
    const { plugin } = createFakePlugin();
    const params = createStructureRepresentationParams(plugin, Structure.Empty, {
      type: CartoonRepresentationProvider,
      typeParams: { alpha: 0.5 },
      color: 'uniform',
      colorParams: { value: UniformColorThemeProvider.defaultValues.value },
      size: PhysicalSizeThemeProvider,
      sizeParams: { scale: 2 },
    });
    expect(params.type.params.alpha).toBe(0.5);
    expect(params.colorTheme.name).toBe('uniform');
    expect(params.sizeTheme).toEqual({ name: 'physical', params: expect.objectContaining({ scale: 2 }) });
  });

  it('warn only for an unregistered name, not for the provider fields', () => {
    const { plugin, warnings } = createFakePlugin();
    createStructureRepresentationParams(plugin, Structure.Empty, {
      type: CartoonRepresentationProvider,
      color: 'no-such-color',
      size: UniformSizeThemeProvider,
    });
    expect(warnings).toEqual([
      "Structure color theme 'no-such-color' is not registered in this plugin; the registry default is used",
    ]);

    warnings.length = 0;
    createStructureRepresentationParams(plugin, Structure.Empty, {
      type: 'no-such-repr',
      color: ChainIdColorThemeProvider,
    });
    expect(warnings).toEqual([
      "Structure representation 'no-such-repr' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('warn for an unregistered default theme of a provider type', () => {
    const { plugin, warnings } = createFakePlugin();
    const { colorThemeRegistry } = plugin.representation.structure.themes;
    colorThemeRegistry.remove(colorThemeRegistry.get('chain-id'));
    createStructureRepresentationParams(plugin, Structure.Empty, {
      type: CartoonRepresentationProvider,
      size: 'uniform',
    });
    expect(warnings.filter((w) => NotRegistered.test(w))).toEqual([
      "Structure color theme 'chain-id' is not registered in this plugin; the registry default is used",
    ]);
  });

  it('resolve volume fields independently', () => {
    const { plugin, warnings } = createFakePlugin();
    const params = createVolumeRepresentationParams(plugin, Volume.One, {
      type: IsosurfaceRepresentationProvider,
      color: 'uniform',
      size: UniformSizeThemeProvider,
    });
    expect(params.type.name).toBe('isosurface');
    expect(params.colorTheme.name).toBe('uniform');
    expect(params.sizeTheme.name).toBe('uniform');
    expect(warnings).toEqual([]);

    const byName = createVolumeRepresentationParams(plugin, Volume.One, {
      type: 'isosurface',
      color: UniformColorThemeProvider,
    });
    expect(byName.type.name).toBe('isosurface');
    expect(byName.colorTheme.name).toBe('uniform');
    expect(warnings).toEqual([]);
  });

  it('keep the empty registry errors for mixed props', () => {
    const { plugin } = createFakePlugin(['color-themes']);
    expect(() =>
      createStructureRepresentationParams(plugin, Structure.Empty, {
        type: CartoonRepresentationProvider,
        color: 'chain-id',
      }),
    ).toThrow('No structure color themes are registered in this plugin');
  });
});

/** Compile-time checks of the params types derived per field. */
function paramTypes() {
  const plugin = createFakePlugin().plugin;
  // a provider field takes the params of that provider
  createStructureRepresentationParams(plugin, undefined, {
    type: CartoonRepresentationProvider,
    typeParams: { alpha: 1 },
  });
  createStructureRepresentationParams(plugin, undefined, {
    type: CartoonRepresentationProvider,
    // @ts-expect-error not a param of the cartoon representation
    typeParams: { nope: 1 },
  });
  // a built-in name takes the params of the built-in provider, whatever the other fields are
  createStructureRepresentationParams(plugin, undefined, {
    type: CartoonRepresentationProvider,
    color: 'uniform',
    colorParams: { saturation: 1 },
  });
  createStructureRepresentationParams(plugin, undefined, {
    color: 'uniform',
    // @ts-expect-error not a param of the uniform color theme
    colorParams: { nope: 1 },
  });
  // any other name takes anything
  const name: string = 'custom';
  createStructureRepresentationParams(plugin, undefined, { type: name, typeParams: { anything: 1 } });
}
void paramTypes;
