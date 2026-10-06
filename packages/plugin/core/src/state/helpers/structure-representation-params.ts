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

export interface StructureRepresentationBuiltInProps<
  R extends StructureRepresentationRegistry.BuiltIn = StructureRepresentationRegistry.BuiltIn,
  C extends ColorTheme.BuiltIn = ColorTheme.BuiltIn,
  S extends SizeTheme.BuiltIn = SizeTheme.BuiltIn,
> {
  /** Using any registered name will work, but code completion will break */
  type?: R;
  typeParams?: StructureRepresentationRegistry.BuiltInParams<R>;
  /** Using any registered name will work, but code completion will break */
  color?: C;
  colorParams?: ColorTheme.BuiltInParams<C>;
  /** Using any registered name will work, but code completion will break */
  size?: S;
  sizeParams?: SizeTheme.BuiltInParams<S>;
}

export interface StructureRepresentationProps<
  R extends RepresentationProvider<Structure> = RepresentationProvider<Structure>,
  C extends ColorTheme.Provider = ColorTheme.Provider,
  S extends SizeTheme.Provider = SizeTheme.Provider,
> {
  type?: R;
  typeParams?: Partial<RepresentationProvider.ParamValues<R>>;
  color?: C;
  colorParams?: Partial<ColorTheme.ParamValues<C>>;
  size?: S;
  sizeParams?: Partial<SizeTheme.ParamValues<S>>;
}

export function createStructureRepresentationParams<
  R extends StructureRepresentationRegistry.BuiltIn,
  C extends ColorTheme.BuiltIn,
  S extends SizeTheme.BuiltIn,
>(
  ctx: PluginContext,
  structure?: Structure,
  props?: StructureRepresentationBuiltInProps<R, C, S>,
): StateTransformer.Params<StructureRepresentation3D>;
export function createStructureRepresentationParams<
  R extends RepresentationProvider<Structure>,
  C extends ColorTheme.Provider,
  S extends SizeTheme.Provider,
>(
  ctx: PluginContext,
  structure?: Structure,
  props?: StructureRepresentationProps<R, C, S>,
): StateTransformer.Params<StructureRepresentation3D>;
export function createStructureRepresentationParams(
  ctx: PluginContext,
  structure?: Structure,
  props: any = {},
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
  const p = props as StructureRepresentationBuiltInProps;
  if (typeof p.type === 'string' || typeof p.color === 'string' || typeof p.size === 'string')
    return createParamsByName(ctx, structure || Structure.Empty, props);
  return createParamsProvider(ctx, structure || Structure.Empty, props);
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

function createParamsByName(
  ctx: PluginContext,
  structure: Structure,
  props: StructureRepresentationBuiltInProps,
): StateTransformer.Params<StructureRepresentation3D> {
  const typeProvider =
    (props.type && ctx.representation.structure.registry.get(props.type)) ||
    ctx.representation.structure.registry.get(ctx.representation.structure.registry.default!.name);
  const colorProvider =
    (props.color && ctx.representation.structure.themes.colorThemeRegistry.get(props.color)) ||
    ctx.representation.structure.themes.colorThemeRegistry.get(typeProvider.defaultColorTheme.name);
  const sizeProvider =
    (props.size && ctx.representation.structure.themes.sizeThemeRegistry.get(props.size)) ||
    ctx.representation.structure.themes.sizeThemeRegistry.get(typeProvider.defaultSizeTheme.name);

  return createParamsProvider(ctx, structure, {
    type: typeProvider,
    typeParams: props.typeParams,
    color: colorProvider,
    colorParams: props.colorParams,
    size: sizeProvider,
    sizeParams: props.sizeParams,
  });
}

function createParamsProvider(
  ctx: PluginContext,
  structure: Structure,
  props: StructureRepresentationProps = {},
): StateTransformer.Params<StructureRepresentation3D> {
  const { themes: themeCtx } = ctx.representation.structure;
  const themeDataCtx = { structure };

  const repr =
    props.type || ctx.representation.structure.registry.get(ctx.representation.structure.registry.default!.name);
  const reprDefaultParams = PD.getDefaultValues(repr.getParams(themeCtx, structure));
  const reprParams = Object.assign(reprDefaultParams, props.typeParams);

  const color = props.color || themeCtx.colorThemeRegistry.get(repr.defaultColorTheme.name);
  const colorDefaultParams = PD.getDefaultValues(color.getParams(themeDataCtx));
  if (color.name === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
  const colorParams = Object.assign(colorDefaultParams, props.colorParams);

  const size = props.size || themeCtx.sizeThemeRegistry.get(repr.defaultSizeTheme.name);
  const sizeDefaultParams = PD.getDefaultValues(size.getParams(themeDataCtx));
  if (size.name === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
  const sizeParams = Object.assign(sizeDefaultParams, props.sizeParams);

  return {
    type: { name: repr.name, params: reprParams },
    colorTheme: { name: color.name, params: colorParams },
    sizeTheme: { name: size.name, params: sizeParams },
  };
}
