/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { RepresentationProvider } from '../../mol-repr/representation';
import { StructureRepresentationRegistry } from '../../mol-repr/structure/registry';
import { StateTransformer } from '../../mol-state';
import { ColorTheme } from '../../mol-theme/color';
import { SizeTheme } from '../../mol-theme/size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentation3D } from '../transforms/representation';

export interface StructureRepresentationBuiltInProps<
    R extends StructureRepresentationRegistry.BuiltIn = StructureRepresentationRegistry.BuiltIn,
    C extends ColorTheme.BuiltIn = ColorTheme.BuiltIn,
    S extends SizeTheme.BuiltIn = SizeTheme.BuiltIn> {
    /** Using any registered name will work, but code completion will break */
    type?: R,
    typeParams?: StructureRepresentationRegistry.BuiltInParams<R>,
    /** Using any registered name will work, but code completion will break */
    color?: C,
    colorParams?: ColorTheme.BuiltInParams<C>,
    /** Using any registered name will work, but code completion will break */
    size?: S,
    sizeParams?: SizeTheme.BuiltInParams<S>
}

export interface StructureRepresentationProps<
    R extends RepresentationProvider<Structure> = RepresentationProvider<Structure>,
    C extends ColorTheme.Provider = ColorTheme.Provider,
    S extends SizeTheme.Provider = SizeTheme.Provider> {
    type?: R,
    typeParams?: Partial<RepresentationProvider.ParamValues<R>>,
    color?: C,
    colorParams?: Partial<ColorTheme.ParamValues<C>>,
    size?: S,
    sizeParams?: Partial<SizeTheme.ParamValues<S>>
}

export function createStructureRepresentationParams<R extends StructureRepresentationRegistry.BuiltIn, C extends ColorTheme.BuiltIn, S extends SizeTheme.BuiltIn>(ctx: PluginContext, structure?: Structure, props?: StructureRepresentationBuiltInProps<R, C, S>): StateTransformer.Params<StructureRepresentation3D>
export function createStructureRepresentationParams<R extends RepresentationProvider<Structure>, C extends ColorTheme.Provider, S extends SizeTheme.Provider>(ctx: PluginContext, structure?: Structure, props?: StructureRepresentationProps<R, C, S>): StateTransformer.Params<StructureRepresentation3D>
export function createStructureRepresentationParams(ctx: PluginContext, structure?: Structure, props: any = {}): StateTransformer.Params<StructureRepresentation3D>  {
    const p = props as StructureRepresentationBuiltInProps;
    if (typeof p.type === 'string' || typeof p.color === 'string' || typeof p.size === 'string') return createParamsByName(ctx, structure || Structure.Empty, props);
    return createParamsProvider(ctx, structure || Structure.Empty, props);
}

export function getStructureThemeTypes(ctx: PluginContext, structure?: Structure) {
    const { themes: themeCtx } = ctx.representation.structure;
    if (!structure) return themeCtx.colorThemeRegistry.types;
    return themeCtx.colorThemeRegistry.getApplicableTypes({ structure });
}

export function createStructureColorThemeParams<T extends ColorTheme.BuiltIn>(ctx: PluginContext, structure: Structure | undefined, typeName: string | undefined, themeName: T, params?: ColorTheme.BuiltInParams<T>): StateTransformer.Params<StructureRepresentation3D>['colorTheme']
export function createStructureColorThemeParams(ctx: PluginContext, structure: Structure | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<StructureRepresentation3D>['colorTheme']
export function createStructureColorThemeParams(ctx: PluginContext, structure: Structure | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<StructureRepresentation3D>['colorTheme'] {
    const { registry, themes } = ctx.representation.structure;
    const repr = registry.get(typeName || registry.default.name);
    const color = themes.colorThemeRegistry.get(themeName || repr.defaultColorTheme.name);
    const colorDefaultParams = PD.getDefaultValues(color.getParams({ structure: structure || Structure.Empty }));
    if (color.name === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
    return { name: color.name, params: Object.assign(colorDefaultParams, params) };
}

export function createStructureSizeThemeParams<T extends SizeTheme.BuiltIn>(ctx: PluginContext, structure: Structure | undefined, typeName: string | undefined, themeName: T, params?: SizeTheme.BuiltInParams<T>): StateTransformer.Params<StructureRepresentation3D>['sizeTheme']
export function createStructureSizeThemeParams(ctx: PluginContext, structure: Structure | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<StructureRepresentation3D>['sizeTheme']
export function createStructureSizeThemeParams(ctx: PluginContext, structure: Structure | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<StructureRepresentation3D>['sizeTheme'] {
    const { registry, themes } = ctx.representation.structure;
    const repr = registry.get(typeName || registry.default.name);
    const size = themes.sizeThemeRegistry.get(themeName || repr.defaultSizeTheme.name);
    const sizeDefaultParams = PD.getDefaultValues(size.getParams({ structure: structure || Structure.Empty }));
    if (size.name === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
    return { name: size.name, params: Object.assign(sizeDefaultParams, params) };
}

function createParamsByName(ctx: PluginContext, structure: Structure, props: StructureRepresentationBuiltInProps): StateTransformer.Params<StructureRepresentation3D> {
    const typeProvider = (props.type && ctx.representation.structure.registry.get(props.type))
        || ctx.representation.structure.registry.default.provider;
    const colorProvider = (props.color && ctx.representation.structure.themes.colorThemeRegistry.get(props.color))
        || ctx.representation.structure.themes.colorThemeRegistry.get(typeProvider.defaultColorTheme.name);
    const sizeProvider = (props.size && ctx.representation.structure.themes.sizeThemeRegistry.get(props.size))
        || ctx.representation.structure.themes.sizeThemeRegistry.get(typeProvider.defaultSizeTheme.name);

    return createParamsProvider(ctx, structure, {
        type: typeProvider,
        typeParams: props.typeParams,
        color: colorProvider,
        colorParams: props.colorParams,
        size: sizeProvider,
        sizeParams: props.sizeParams
    });
}

function createParamsProvider(ctx: PluginContext, structure: Structure, props: StructureRepresentationProps = {}): StateTransformer.Params<StructureRepresentation3D> {
    const { themes: themeCtx } = ctx.representation.structure;
    const themeDataCtx = { structure };

    const repr = props.type || ctx.representation.structure.registry.default.provider;
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

    return ({
        type: { name: repr.name, params: reprParams },
        colorTheme: { name: color.name, params: colorParams },
        sizeTheme: { name: size.name, params: sizeParams }
    });
}