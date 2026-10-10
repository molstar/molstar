/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Volume } from '@molstar/model/model/volume';
import type { PluginContext } from '@molstar/plugin/context';
import type { RepresentationProvider } from '@molstar/graphics/repr/representation';
import type { VolumeRepresentationRegistry } from '@molstar/graphics/repr/volume/registry';
import type { StateTransformer } from '@molstar/core/state';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  assertColorThemes,
  assertRepresentations,
  assertRepresentationScope,
  assertSizeThemes,
  emptyMappedValue,
  hasColorThemes,
  hasRepresentations,
  hasRepresentationScope,
  hasSizeThemes,
  warnUnregisteredNames,
} from './representation-registry.js';
import type { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';

/** A representation of the scope, given as a provider or as a registered name. */
export type VolumeRepresentationRef = RepresentationProvider<Volume> | string;
/** A color theme, given as a provider or as a registered name. */
export type VolumeColorThemeRef = ColorTheme.Provider | string;
/** A size theme, given as a provider or as a registered name. */
export type VolumeSizeThemeRef = SizeTheme.Provider | string;

/** Params of a representation reference: a provider's params, the params of a built-in name, or `{}` for other names. */
export type VolumeRepresentationParamsOf<R extends VolumeRepresentationRef> = [R] extends [
  RepresentationProvider<Volume>,
]
  ? Partial<RepresentationProvider.ParamValues<R>>
  : [R] extends [VolumeRepresentationRegistry.BuiltIn]
    ? VolumeRepresentationRegistry.BuiltInParams<R>
    : {};
/** Params of a color theme reference: a provider's params, the params of a built-in name, or `{}` for other names. */
export type VolumeColorThemeParamsOf<C extends VolumeColorThemeRef> = [C] extends [ColorTheme.Provider]
  ? Partial<ColorTheme.ParamValues<C>>
  : [C] extends [ColorTheme.BuiltIn]
    ? ColorTheme.BuiltInParams<C>
    : {};
/** Params of a size theme reference: a provider's params, the params of a built-in name, or `{}` for other names. */
export type VolumeSizeThemeParamsOf<S extends VolumeSizeThemeRef> = [S] extends [SizeTheme.Provider]
  ? Partial<SizeTheme.ParamValues<S>>
  : [S] extends [SizeTheme.BuiltIn]
    ? SizeTheme.BuiltInParams<S>
    : {};

/**
 * Each of `type`, `color`, and `size` is a provider or a name and is resolved independently, so a provider can be
 * combined with a theme name. The params types follow the field: a provider's params, the params of a built-in name, or
 * `{}` for any other name.
 */
export interface VolumeRepresentationProps<
  R extends VolumeRepresentationRef = VolumeRepresentationRef,
  C extends VolumeColorThemeRef = VolumeColorThemeRef,
  S extends VolumeSizeThemeRef = VolumeSizeThemeRef,
> {
  type?: R;
  typeParams?: VolumeRepresentationParamsOf<R>;
  color?: C;
  colorParams?: VolumeColorThemeParamsOf<C>;
  size?: S;
  sizeParams?: VolumeSizeThemeParamsOf<S>;
}

/** Props with representation, color theme, and size theme names. */
export type VolumeRepresentationBuiltInProps<
  R extends VolumeRepresentationRegistry.BuiltIn = VolumeRepresentationRegistry.BuiltIn,
  C extends ColorTheme.BuiltIn = ColorTheme.BuiltIn,
  S extends SizeTheme.BuiltIn = SizeTheme.BuiltIn,
> = VolumeRepresentationProps<R, C, S>;

export function createVolumeRepresentationParams<
  R extends VolumeRepresentationRef = VolumeRepresentationRef,
  C extends VolumeColorThemeRef = VolumeColorThemeRef,
  S extends VolumeSizeThemeRef = VolumeSizeThemeRef,
>(
  ctx: PluginContext,
  volume?: Volume,
  props: VolumeRepresentationProps<R, C, S> = {},
): StateTransformer.Params<VolumeRepresentation3D> {
  const scope = ctx.representation.volume;
  if (volume) assertRepresentationScope('volume', scope);
  else if (!hasRepresentationScope(scope)) {
    // nothing to build params from; applying them fails in the transformer
    return { type: emptyMappedValue(), colorTheme: emptyMappedValue(), sizeTheme: emptyMappedValue() };
  }

  warnUnregisteredNames(ctx, 'volume', scope, {
    type: props.type || undefined,
    color: props.color || undefined,
    size: props.size || undefined,
    checkColor: true,
    checkSize: true,
  });
  return createParams(ctx, volume || Volume.One, props);
}

export function getVolumeThemeTypes(ctx: PluginContext, volume?: Volume) {
  const { themes: themeCtx } = ctx.representation.volume;
  if (!volume) return themeCtx.colorThemeRegistry.types;
  return themeCtx.colorThemeRegistry.getApplicableTypes({ volume });
}

export function createVolumeColorThemeParams<T extends ColorTheme.BuiltIn>(
  ctx: PluginContext,
  volume: Volume | undefined,
  typeName: string | undefined,
  themeName: T,
  params?: ColorTheme.BuiltInParams<T>,
): StateTransformer.Params<VolumeRepresentation3D>['colorTheme'];
export function createVolumeColorThemeParams(
  ctx: PluginContext,
  volume: Volume | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<VolumeRepresentation3D>['colorTheme'];
export function createVolumeColorThemeParams(
  ctx: PluginContext,
  volume: Volume | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<VolumeRepresentation3D>['colorTheme'] {
  const scope = ctx.representation.volume;
  const { registry, themes } = scope;
  if (volume) {
    assertRepresentations('volume', scope);
    assertColorThemes('volume', scope);
  } else if (!(hasRepresentations(scope) && hasColorThemes(scope))) {
    // nothing to build params from; applying them fails in the transformer
    return emptyMappedValue();
  }
  warnUnregisteredNames(ctx, 'volume', scope, {
    type: typeName || undefined,
    color: themeName || undefined,
    checkColor: true,
    checkSize: false,
  });
  const repr = registry.get(typeName || registry.default!.name);
  const color = themes.colorThemeRegistry.get(themeName || repr.defaultColorTheme.name);
  const colorDefaultParams = PD.getDefaultValues(color.getParams({ volume: volume || Volume.One }));
  if (color.name === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
  return { name: color.name, params: Object.assign(colorDefaultParams, params) };
}

export function createVolumeSizeThemeParams<T extends SizeTheme.BuiltIn>(
  ctx: PluginContext,
  volume: Volume | undefined,
  typeName: string | undefined,
  themeName: T,
  params?: SizeTheme.BuiltInParams<T>,
): StateTransformer.Params<VolumeRepresentation3D>['sizeTheme'];
export function createVolumeSizeThemeParams(
  ctx: PluginContext,
  volume: Volume | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<VolumeRepresentation3D>['sizeTheme'];
export function createVolumeSizeThemeParams(
  ctx: PluginContext,
  volume: Volume | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<VolumeRepresentation3D>['sizeTheme'] {
  const scope = ctx.representation.volume;
  const { registry, themes } = scope;
  if (volume) {
    assertRepresentations('volume', scope);
    assertSizeThemes('volume', scope);
  } else if (!(hasRepresentations(scope) && hasSizeThemes(scope))) {
    // nothing to build params from; applying them fails in the transformer
    return emptyMappedValue();
  }
  warnUnregisteredNames(ctx, 'volume', scope, {
    type: typeName || undefined,
    size: themeName || undefined,
    checkSize: true,
    checkColor: false,
  });
  const repr = registry.get(typeName || registry.default!.name);
  const size = themes.sizeThemeRegistry.get(themeName || repr.defaultSizeTheme.name);
  const sizeDefaultParams = PD.getDefaultValues(size.getParams({ volume: volume || Volume.One }));
  if (size.name === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
  return { name: size.name, params: Object.assign(sizeDefaultParams, params) };
}

function createParams(
  ctx: PluginContext,
  volume: Volume,
  props: VolumeRepresentationProps = {},
): StateTransformer.Params<VolumeRepresentation3D> {
  const { registry, themes: themeCtx } = ctx.representation.volume;
  const themeDataCtx = { volume };

  // Each field is a provider or a name and is resolved on its own.
  const repr =
    (typeof props.type === 'string' ? props.type && registry.get(props.type) : props.type) ||
    registry.get(registry.default!.name);
  const reprDefaultParams = PD.getDefaultValues(repr.getParams(themeCtx, volume));
  const reprParams = Object.assign(reprDefaultParams, props.typeParams);

  const color =
    (typeof props.color === 'string' ? props.color && themeCtx.colorThemeRegistry.get(props.color) : props.color) ||
    themeCtx.colorThemeRegistry.get(repr.defaultColorTheme.name);
  const colorDefaultParams = PD.getDefaultValues(color.getParams(themeDataCtx));
  if (color.name === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
  const colorParams = Object.assign(colorDefaultParams, props.colorParams);

  const size =
    (typeof props.size === 'string' ? props.size && themeCtx.sizeThemeRegistry.get(props.size) : props.size) ||
    themeCtx.sizeThemeRegistry.get(repr.defaultSizeTheme.name);
  const sizeDefaultParams = PD.getDefaultValues(size.getParams(themeDataCtx));
  if (size.name === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
  const sizeParams = Object.assign(sizeDefaultParams, props.sizeParams);

  return {
    type: { name: repr.name, params: reprParams },
    colorTheme: { name: color.name, params: colorParams },
    sizeTheme: { name: size.name, params: sizeParams },
  };
}
