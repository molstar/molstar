/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { VolumeRepresentationRegistry } from '@molstar/graphics/repr/volume/registry';
import { Volume } from '@molstar/model/model/volume';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RepresentationProvider } from '@molstar/graphics/repr/representation';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import type { StateTransformer } from '@molstar/core/state';
import type { VolumeRepresentation3D } from './representation.js';

/** A name is looked up in the registry, a provider object is used as is. */
function resolveType(ctx: PluginContext, type: string | RepresentationProvider<Volume>) {
  return typeof type === 'string' ? ctx.representation.volume.registry.get(type) : type;
}

function resolveColor(ctx: PluginContext, color: string | ColorTheme.Provider) {
  return typeof color === 'string' ? ctx.representation.volume.themes.colorThemeRegistry.get(color) : color;
}

function resolveSize(ctx: PluginContext, size: string | SizeTheme.Provider) {
  return typeof size === 'string' ? ctx.representation.volume.themes.sizeThemeRegistry.get(size) : size;
}

export namespace VolumeRepresentation3DHelpers {
  export function getDefaultParams(
    ctx: PluginContext,
    name: VolumeRepresentationRegistry.BuiltIn | RepresentationProvider<Volume>,
    volume: Volume,
    volumeParams?: Partial<PD.Values<PD.Params>>,
    colorName?: ColorTheme.BuiltIn | ColorTheme.Provider,
    colorParams?: Partial<ColorTheme.Props>,
    sizeName?: SizeTheme.BuiltIn | SizeTheme.Provider,
    sizeParams?: Partial<SizeTheme.Props>,
  ): StateTransformer.Params<VolumeRepresentation3D> {
    const type = resolveType(ctx, name);
    const colorType = resolveColor(ctx, colorName || type.defaultColorTheme.name);
    const sizeType = resolveSize(ctx, sizeName || type.defaultSizeTheme.name);
    const volumeDefaultParams = PD.getDefaultValues(type.getParams(ctx.representation.volume.themes, volume));
    return {
      type: {
        name: type.name,
        params: volumeParams ? { ...volumeDefaultParams, ...volumeParams } : volumeDefaultParams,
      },
      colorTheme: {
        name: colorType.name,
        params: colorParams ? { ...colorType.defaultValues, ...colorParams } : colorType.defaultValues,
      },
      sizeTheme: {
        name: sizeType.name,
        params: sizeParams ? { ...sizeType.defaultValues, ...sizeParams } : sizeType.defaultValues,
      },
    };
  }

  export function getDefaultParamsStatic(
    ctx: PluginContext,
    name: VolumeRepresentationRegistry.BuiltIn | RepresentationProvider<Volume>,
    volumeParams?: Partial<PD.Values<PD.Params>>,
    colorName?: ColorTheme.BuiltIn | ColorTheme.Provider,
    colorParams?: Partial<ColorTheme.Props>,
    sizeName?: SizeTheme.BuiltIn | SizeTheme.Provider,
    sizeParams?: Partial<SizeTheme.Props>,
  ): StateTransformer.Params<VolumeRepresentation3D> {
    const type = resolveType(ctx, name);
    const colorType = resolveColor(ctx, colorName || type.defaultColorTheme.name);
    const sizeType = resolveSize(ctx, sizeName || type.defaultSizeTheme.name);
    return {
      type: { name: type.name, params: volumeParams ? { ...type.defaultValues, ...volumeParams } : type.defaultValues },
      colorTheme: {
        name: type.defaultColorTheme.name,
        params: colorParams ? { ...colorType.defaultValues, ...colorParams } : colorType.defaultValues,
      },
      sizeTheme: {
        name: type.defaultSizeTheme.name,
        params: sizeParams ? { ...sizeType.defaultValues, ...sizeParams } : sizeType.defaultValues,
      },
    };
  }

  export function getDescription(props: any) {
    if (props.isoValue) {
      return Volume.IsoValue.toString(props.isoValue);
    } else if (props.renderMode?.params?.isoValue) {
      return Volume.IsoValue.toString(props.renderMode?.params?.isoValue);
    }
  }
}
