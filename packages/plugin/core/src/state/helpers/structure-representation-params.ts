/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { RepresentationProvider } from '@molstar/graphics/repr/representation';
import type { StructureRepresentationRegistry } from '@molstar/graphics/repr/structure/registry';
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
import type { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';

export function isSurfaceRepresentationType(name: string) {
  return name.endsWith('-surface');
}

/** A representation of the scope, given as a provider or as a registered name. */
export type StructureRepresentationRef = RepresentationProvider<Structure> | string;
/** A color theme, given as a provider or as a registered name. */
export type StructureColorThemeRef = ColorTheme.Provider | string;
/** A size theme, given as a provider or as a registered name. */
export type StructureSizeThemeRef = SizeTheme.Provider | string;

/** Params of a representation reference: a provider's params, the params of a built-in name, or `{}` for other names. */
export type StructureRepresentationParamsOf<R extends StructureRepresentationRef> = [R] extends [
  RepresentationProvider<Structure>,
]
  ? Partial<RepresentationProvider.ParamValues<R>>
  : [R] extends [StructureRepresentationRegistry.BuiltIn]
    ? StructureRepresentationRegistry.BuiltInParams<R>
    : {};
/** Params of a color theme reference: a provider's params, the params of a built-in name, or `{}` for other names. */
export type StructureColorThemeParamsOf<C extends StructureColorThemeRef> = [C] extends [ColorTheme.Provider]
  ? Partial<ColorTheme.ParamValues<C>>
  : [C] extends [ColorTheme.BuiltIn]
    ? ColorTheme.BuiltInParams<C>
    : {};
/** Params of a size theme reference: a provider's params, the params of a built-in name, or `{}` for other names. */
export type StructureSizeThemeParamsOf<S extends StructureSizeThemeRef> = [S] extends [SizeTheme.Provider]
  ? Partial<SizeTheme.ParamValues<S>>
  : [S] extends [SizeTheme.BuiltIn]
    ? SizeTheme.BuiltInParams<S>
    : {};

/**
 * Each of `type`, `color`, and `size` is a provider or a name and is resolved independently, so a provider can be
 * combined with a theme name. The params types follow the field: a provider's params, the params of a built-in name, or
 * `{}` for any other name.
 */
export interface StructureRepresentationProps<
  R extends StructureRepresentationRef = StructureRepresentationRef,
  C extends StructureColorThemeRef = StructureColorThemeRef,
  S extends StructureSizeThemeRef = StructureSizeThemeRef,
> {
  type?: R;
  typeParams?: StructureRepresentationParamsOf<R>;
  color?: C;
  colorParams?: StructureColorThemeParamsOf<C>;
  size?: S;
  sizeParams?: StructureSizeThemeParamsOf<S>;
}

/** Props with representation, color theme, and size theme names. */
export type StructureRepresentationBuiltInProps<
  R extends StructureRepresentationRegistry.BuiltIn = StructureRepresentationRegistry.BuiltIn,
  C extends ColorTheme.BuiltIn = ColorTheme.BuiltIn,
  S extends SizeTheme.BuiltIn = SizeTheme.BuiltIn,
> = StructureRepresentationProps<R, C, S>;

export function createStructureRepresentationParams<
  R extends StructureRepresentationRef = StructureRepresentationRef,
  C extends StructureColorThemeRef = StructureColorThemeRef,
  S extends StructureSizeThemeRef = StructureSizeThemeRef,
>(
  ctx: PluginContext,
  structure?: Structure,
  props: StructureRepresentationProps<R, C, S> = {},
): StateTransformer.Params<StructureRepresentation3D> {
  const scope = ctx.representation.structure;
  if (structure) assertRepresentationScope('structure', scope);
  else if (!hasRepresentationScope(scope)) {
    // nothing to build params from; applying them fails in the transformer
    return { type: emptyMappedValue(), colorTheme: emptyMappedValue(), sizeTheme: emptyMappedValue() };
  }

  warnUnregisteredNames(ctx, 'structure', scope, {
    type: props.type || undefined,
    color: props.color || undefined,
    size: props.size || undefined,
    checkColor: true,
    checkSize: true,
  });
  return createParams(ctx, structure || Structure.Empty, props);
}

export function getStructureThemeTypes(ctx: PluginContext, structure?: Structure) {
  const { themes: themeCtx } = ctx.representation.structure;
  if (!structure) return themeCtx.colorThemeRegistry.types;
  return themeCtx.colorThemeRegistry.getApplicableTypes({ structure });
}

export function createStructureColorThemeParams<T extends ColorTheme.BuiltIn>(
  ctx: PluginContext,
  structure: Structure | undefined,
  typeName: string | undefined,
  themeName: T,
  params?: ColorTheme.BuiltInParams<T>,
): StateTransformer.Params<StructureRepresentation3D>['colorTheme'];
export function createStructureColorThemeParams(
  ctx: PluginContext,
  structure: Structure | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<StructureRepresentation3D>['colorTheme'];
export function createStructureColorThemeParams(
  ctx: PluginContext,
  structure: Structure | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<StructureRepresentation3D>['colorTheme'] {
  const scope = ctx.representation.structure;
  const { registry, themes } = scope;
  if (structure) {
    assertRepresentations('structure', scope);
    assertColorThemes('structure', scope);
  } else if (!(hasRepresentations(scope) && hasColorThemes(scope))) {
    // nothing to build params from; applying them fails in the transformer
    return emptyMappedValue();
  }
  warnUnregisteredNames(ctx, 'structure', scope, {
    type: typeName || undefined,
    color: themeName || undefined,
    checkColor: true,
    checkSize: false,
  });
  const repr = registry.get(typeName || registry.default!.name);
  const color = themes.colorThemeRegistry.get(themeName || repr.defaultColorTheme.name);
  const colorDefaultParams = PD.getDefaultValues(color.getParams({ structure: structure || Structure.Empty }));
  if (color.name === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
  return { name: color.name, params: Object.assign(colorDefaultParams, params) };
}

export function createStructureSizeThemeParams<T extends SizeTheme.BuiltIn>(
  ctx: PluginContext,
  structure: Structure | undefined,
  typeName: string | undefined,
  themeName: T,
  params?: SizeTheme.BuiltInParams<T>,
): StateTransformer.Params<StructureRepresentation3D>['sizeTheme'];
export function createStructureSizeThemeParams(
  ctx: PluginContext,
  structure: Structure | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<StructureRepresentation3D>['sizeTheme'];
export function createStructureSizeThemeParams(
  ctx: PluginContext,
  structure: Structure | undefined,
  typeName: string | undefined,
  themeName?: string,
  params?: any,
): StateTransformer.Params<StructureRepresentation3D>['sizeTheme'] {
  const scope = ctx.representation.structure;
  const { registry, themes } = scope;
  if (structure) {
    assertRepresentations('structure', scope);
    assertSizeThemes('structure', scope);
  } else if (!(hasRepresentations(scope) && hasSizeThemes(scope))) {
    // nothing to build params from; applying them fails in the transformer
    return emptyMappedValue();
  }
  warnUnregisteredNames(ctx, 'structure', scope, {
    type: typeName || undefined,
    size: themeName || undefined,
    checkSize: true,
    checkColor: false,
  });
  const repr = registry.get(typeName || registry.default!.name);
  const size = themes.sizeThemeRegistry.get(themeName || repr.defaultSizeTheme.name);
  const sizeDefaultParams = PD.getDefaultValues(size.getParams({ structure: structure || Structure.Empty }));
  if (size.name === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
  return { name: size.name, params: Object.assign(sizeDefaultParams, params) };
}

function createParams(
  ctx: PluginContext,
  structure: Structure,
  props: StructureRepresentationProps = {},
): StateTransformer.Params<StructureRepresentation3D> {
  const { registry, themes: themeCtx } = ctx.representation.structure;
  const themeDataCtx = { structure };

  // Each field is a provider or a name and is resolved on its own.
  const repr =
    (typeof props.type === 'string' ? props.type && registry.get(props.type) : props.type) ||
    registry.get(registry.default!.name);
  const reprDefaultParams = PD.getDefaultValues(repr.getParams(themeCtx, structure));
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
