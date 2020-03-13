/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../../mol-model/structure';
import { PluginContext } from '../../mol-plugin/context';
import { RepresentationProvider } from '../../mol-repr/representation';
import { BuiltInStructureRepresentations, BuiltInStructureRepresentationsName } from '../../mol-repr/structure/registry'
import { StateTransformer } from '../../mol-state';
import { BuiltInColorThemeName, BuiltInColorThemes, ColorTheme } from '../../mol-theme/color';
import { BuiltInSizeThemeName, BuiltInSizeThemes, SizeTheme } from '../../mol-theme/size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { StructureRepresentation3D } from '../transforms/representation';

export interface StructureRepresentationBuiltInProps<
    R extends BuiltInStructureRepresentationsName = BuiltInStructureRepresentationsName,
    C extends BuiltInColorThemeName = BuiltInColorThemeName,
    S extends BuiltInSizeThemeName = BuiltInSizeThemeName> {
    /** Using any registered name will work, but code completion will break */
    type?: R,
    typeParams?: Partial<RepresentationProvider.ParamValues<BuiltInStructureRepresentations[R]>>,
    /** Using any registered name will work, but code completion will break */
    color?: C,
    colorParams?: Partial<ColorTheme.ParamValues<BuiltInColorThemes[C]>>,
    /** Using any registered name will work, but code completion will break */
    size?: S,
    sizeParams?: Partial<SizeTheme.ParamValues<BuiltInSizeThemes[S]>>
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

export function createStructureRepresentationParams
    <R extends BuiltInStructureRepresentationsName, C extends BuiltInColorThemeName, S extends BuiltInSizeThemeName>
    (ctx: PluginContext, structure?: Structure, props?: StructureRepresentationBuiltInProps<R, C, S>): StateTransformer.Params<StructureRepresentation3D>
export function createStructureRepresentationParams
    <R extends RepresentationProvider<Structure>, C extends ColorTheme.Provider, S extends SizeTheme.Provider>
    (ctx: PluginContext, structure?: Structure, props?: StructureRepresentationProps<R, C, S>): StateTransformer.Params<StructureRepresentation3D> 
export function createStructureRepresentationParams(ctx: PluginContext, structure?: Structure, props: any = {}): StateTransformer.Params<StructureRepresentation3D>  {
    const p = props as StructureRepresentationBuiltInProps;
    if (typeof p.type === 'string' || typeof p.color === 'string' || typeof p.size === 'string') return createParamsByName(ctx, structure || Structure.Empty, props);
    return createParamsProvider(ctx, structure || Structure.Empty, props);
}

function createParamsByName(ctx: PluginContext, structure: Structure, props: StructureRepresentationBuiltInProps): StateTransformer.Params<StructureRepresentation3D> {
    const typeProvider = (props.type && ctx.structureRepresentation.registry.get(props.type))
        || ctx.structureRepresentation.registry.default.provider;
    const colorProvider = (props.color && ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(props.color)) 
        || ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(typeProvider.defaultColorTheme.name);
    const sizeProvider = (props.size && ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(props.size))
        || ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(typeProvider.defaultSizeTheme.name);

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
    const { themeCtx } = ctx.structureRepresentation
    const themeDataCtx = { structure };
    
    const repr = props.type || ctx.structureRepresentation.registry.default.provider;
    const reprDefaultParams = PD.getDefaultValues(repr.getParams(themeCtx, structure));
    const reprParams = Object.assign(reprDefaultParams, props.typeParams);
    
    const color = props.color || themeCtx.colorThemeRegistry.get(repr.defaultColorTheme.name);
    const colorName = themeCtx.colorThemeRegistry.getName(color);
    const colorDefaultParams = PD.getDefaultValues(color.getParams(themeDataCtx));
    if (colorName === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
    const colorParams = Object.assign(colorDefaultParams, props.colorParams);

    const size = props.size || themeCtx.sizeThemeRegistry.get(repr.defaultSizeTheme.name);
    const sizeName = themeCtx.sizeThemeRegistry.getName(size);
    const sizeDefaultParams = PD.getDefaultValues(size.getParams(themeDataCtx));
    if (sizeName === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
    const sizeParams = Object.assign(sizeDefaultParams, props.sizeParams);

    return ({
        type: { name: ctx.structureRepresentation.registry.getName(repr), params: reprParams },
        colorTheme: { name: colorName, params: colorParams },
        sizeTheme: { name: sizeName, params: sizeParams }
    });
}