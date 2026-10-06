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
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import type { StateTransformer } from '@molstar/core/state';
import type { VolumeRepresentation3D } from './representation.js';

export namespace VolumeRepresentation3DHelpers {
  export function getDefaultParams(
    ctx: PluginContext,
    name: VolumeRepresentationRegistry.BuiltIn,
    volume: Volume,
    volumeParams?: Partial<PD.Values<PD.Params>>,
    colorName?: ColorTheme.BuiltIn,
    colorParams?: Partial<ColorTheme.Props>,
    sizeName?: SizeTheme.BuiltIn,
    sizeParams?: Partial<SizeTheme.Props>,
  ): StateTransformer.Params<VolumeRepresentation3D> {
    const type = ctx.representation.volume.registry.get(name);

    const colorType = ctx.representation.volume.themes.colorThemeRegistry.get(colorName || type.defaultColorTheme.name);
    const sizeType = ctx.representation.volume.themes.sizeThemeRegistry.get(sizeName || type.defaultSizeTheme.name);
    const volumeDefaultParams = PD.getDefaultValues(type.getParams(ctx.representation.volume.themes, volume));
    return {
      type: { name, params: volumeParams ? { ...volumeDefaultParams, ...volumeParams } : volumeDefaultParams },
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
    name: VolumeRepresentationRegistry.BuiltIn,
    volumeParams?: Partial<PD.Values<PD.Params>>,
    colorName?: ColorTheme.BuiltIn,
    colorParams?: Partial<ColorTheme.Props>,
    sizeName?: SizeTheme.BuiltIn,
    sizeParams?: Partial<SizeTheme.Props>,
  ): StateTransformer.Params<VolumeRepresentation3D> {
    const type = ctx.representation.volume.registry.get(name);
    const colorType = ctx.representation.volume.themes.colorThemeRegistry.get(colorName || type.defaultColorTheme.name);
    const sizeType = ctx.representation.volume.themes.sizeThemeRegistry.get(sizeName || type.defaultSizeTheme.name);
    return {
      type: { name, params: volumeParams ? { ...type.defaultValues, ...volumeParams } : type.defaultValues },
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
