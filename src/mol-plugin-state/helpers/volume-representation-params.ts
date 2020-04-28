/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Volume } from '../../mol-model/volume';
import { PluginContext } from '../../mol-plugin/context';
import { RepresentationProvider } from '../../mol-repr/representation';
import { VolumeRepresentationRegistry } from '../../mol-repr/volume/registry';
import { StateTransformer } from '../../mol-state';
import { ColorTheme } from '../../mol-theme/color';
import { SizeTheme } from '../../mol-theme/size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { VolumeRepresentation3D } from '../transforms/representation';

export interface VolumeRepresentationBuiltInProps<
    R extends VolumeRepresentationRegistry.BuiltIn = VolumeRepresentationRegistry.BuiltIn,
    C extends ColorTheme.BuiltIn = ColorTheme.BuiltIn,
    S extends SizeTheme.BuiltIn = SizeTheme.BuiltIn> {
    /** Using any registered name will work, but code completion will break */
    type?: R,
    typeParams?: VolumeRepresentationRegistry.BuiltInParams<R>,
    /** Using any registered name will work, but code completion will break */
    color?: C,
    colorParams?: ColorTheme.BuiltInParams<C>,
    /** Using any registered name will work, but code completion will break */
    size?: S,
    sizeParams?: SizeTheme.BuiltInParams<S>
}

export interface VolumeRepresentationProps<
    R extends RepresentationProvider<Volume> = RepresentationProvider<Volume>,
    C extends ColorTheme.Provider = ColorTheme.Provider,
    S extends SizeTheme.Provider = SizeTheme.Provider> {
    type?: R,
    typeParams?: Partial<RepresentationProvider.ParamValues<R>>,
    color?: C,
    colorParams?: Partial<ColorTheme.ParamValues<C>>,
    size?: S,
    sizeParams?: Partial<SizeTheme.ParamValues<S>>
}

export function createVolumeRepresentationParams<R extends VolumeRepresentationRegistry.BuiltIn, C extends ColorTheme.BuiltIn, S extends SizeTheme.BuiltIn>(ctx: PluginContext, volume?: Volume, props?: VolumeRepresentationBuiltInProps<R, C, S>): StateTransformer.Params<VolumeRepresentation3D>
export function createVolumeRepresentationParams<R extends RepresentationProvider<Volume>, C extends ColorTheme.Provider, S extends SizeTheme.Provider>(ctx: PluginContext, volume?: Volume, props?: VolumeRepresentationProps<R, C, S>): StateTransformer.Params<VolumeRepresentation3D>
export function createVolumeRepresentationParams(ctx: PluginContext, volume?: Volume, props: any = {}): StateTransformer.Params<VolumeRepresentation3D>  {
    const p = props as VolumeRepresentationBuiltInProps;
    if (typeof p.type === 'string' || typeof p.color === 'string' || typeof p.size === 'string') return createParamsByName(ctx, volume || Volume.One, props);
    return createParamsProvider(ctx, volume || Volume.One, props);
}

export function getVolumeThemeTypes(ctx: PluginContext, volume?: Volume) {
    const { themes: themeCtx } = ctx.representation.volume;
    if (!volume) return themeCtx.colorThemeRegistry.types;
    return themeCtx.colorThemeRegistry.getApplicableTypes({ volume });
}

export function createVolumeColorThemeParams<T extends ColorTheme.BuiltIn>(ctx: PluginContext, volume: Volume | undefined, typeName: string | undefined, themeName: T, params?: ColorTheme.BuiltInParams<T>): StateTransformer.Params<VolumeRepresentation3D>['colorTheme']
export function createVolumeColorThemeParams(ctx: PluginContext, volume: Volume | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<VolumeRepresentation3D>['colorTheme']
export function createVolumeColorThemeParams(ctx: PluginContext, volume: Volume | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<VolumeRepresentation3D>['colorTheme'] {
    const { registry, themes } = ctx.representation.volume;
    const repr = registry.get(typeName || registry.default.name);
    const color = themes.colorThemeRegistry.get(themeName || repr.defaultColorTheme.name);
    const colorDefaultParams = PD.getDefaultValues(color.getParams({ volume: volume || Volume.One }));
    if (color.name === repr.defaultColorTheme.name) Object.assign(colorDefaultParams, repr.defaultColorTheme.props);
    return { name: color.name, params: Object.assign(colorDefaultParams, params) };
}

export function createVolumeSizeThemeParams<T extends SizeTheme.BuiltIn>(ctx: PluginContext, volume: Volume | undefined, typeName: string | undefined, themeName: T, params?: SizeTheme.BuiltInParams<T>): StateTransformer.Params<VolumeRepresentation3D>['sizeTheme']
export function createVolumeSizeThemeParams(ctx: PluginContext, volume: Volume | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<VolumeRepresentation3D>['sizeTheme']
export function createVolumeSizeThemeParams(ctx: PluginContext, volume: Volume | undefined, typeName: string | undefined, themeName?: string, params?: any): StateTransformer.Params<VolumeRepresentation3D>['sizeTheme'] {
    const { registry, themes } = ctx.representation.volume;
    const repr = registry.get(typeName || registry.default.name);
    const size = themes.sizeThemeRegistry.get(themeName || repr.defaultSizeTheme.name);
    const sizeDefaultParams = PD.getDefaultValues(size.getParams({ volume: volume || Volume.One }));
    if (size.name === repr.defaultSizeTheme.name) Object.assign(sizeDefaultParams, repr.defaultSizeTheme.props);
    return { name: size.name, params: Object.assign(sizeDefaultParams, params) };
}

function createParamsByName(ctx: PluginContext, volume: Volume, props: VolumeRepresentationBuiltInProps): StateTransformer.Params<VolumeRepresentation3D> {
    const typeProvider = (props.type && ctx.representation.volume.registry.get(props.type))
        || ctx.representation.volume.registry.default.provider;
    const colorProvider = (props.color && ctx.representation.volume.themes.colorThemeRegistry.get(props.color))
        || ctx.representation.volume.themes.colorThemeRegistry.get(typeProvider.defaultColorTheme.name);
    const sizeProvider = (props.size && ctx.representation.volume.themes.sizeThemeRegistry.get(props.size))
        || ctx.representation.volume.themes.sizeThemeRegistry.get(typeProvider.defaultSizeTheme.name);

    return createParamsProvider(ctx, volume, {
        type: typeProvider,
        typeParams: props.typeParams,
        color: colorProvider,
        colorParams: props.colorParams,
        size: sizeProvider,
        sizeParams: props.sizeParams
    });
}

function createParamsProvider(ctx: PluginContext, volume: Volume, props: VolumeRepresentationProps = {}): StateTransformer.Params<VolumeRepresentation3D> {
    const { themes: themeCtx } = ctx.representation.volume;
    const themeDataCtx = { volume };

    const repr = props.type || ctx.representation.volume.registry.default.provider;
    const reprDefaultParams = PD.getDefaultValues(repr.getParams(themeCtx, volume));
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