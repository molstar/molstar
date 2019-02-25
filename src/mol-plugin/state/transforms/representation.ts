/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateTransformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { createTheme } from 'mol-theme/theme';
import { BuiltInStructureRepresentationsName } from 'mol-repr/structure/registry';
import { Structure } from 'mol-model/structure';
import { StructureParams } from 'mol-repr/structure/representation';
import { ExplodeRepresentation3D } from 'mol-plugin/behavior/dynamic/representation';
import { VolumeData, VolumeIsoValue } from 'mol-model/volume';
import { BuiltInVolumeRepresentationsName } from 'mol-repr/volume/registry';
import { VolumeParams } from 'mol-repr/volume/representation';
import { BuiltInColorThemeName, ColorTheme } from 'mol-theme/color';
import { BuiltInSizeThemeName, SizeTheme } from 'mol-theme/size';

export namespace StructureRepresentation3DHelpers {
    export function getDefaultParams(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structure: Structure, structureParams?: Partial<PD.Values<StructureParams>>): StateTransformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(name);

        const themeDataCtx = { structure };
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).getParams(themeDataCtx);
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        const structureDefaultParams = PD.getDefaultValues(type.getParams(ctx.structureRepresentation.themeCtx, structure))
        return ({
            type: { name, params: structureParams ? { ...structureDefaultParams, ...structureParams } : structureDefaultParams },
            colorTheme: { name: type.defaultColorTheme, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function getDefaultParamsWithTheme(ctx: PluginContext, reprName: BuiltInStructureRepresentationsName, colorName: BuiltInColorThemeName | undefined, structure: Structure, structureParams?: Partial<PD.Values<StructureParams>>): StateTransformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(reprName);

        const themeDataCtx = { structure };
        const color = colorName || type.defaultColorTheme;
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(color).getParams(themeDataCtx);
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        const structureDefaultParams = PD.getDefaultValues(type.getParams(ctx.structureRepresentation.themeCtx, structure))
        return ({
            type: { name: reprName, params: structureParams ? { ...structureDefaultParams, ...structureParams } : structureDefaultParams },
            colorTheme: { name: color, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function getDefaultParamsStatic(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structureParams?: Partial<PD.Values<StructureParams>>): StateTransformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(name);
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).defaultValues;
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).defaultValues
        return ({
            type: { name, params: structureParams ? { ...type.defaultValues, ...structureParams } : type.defaultValues },
            colorTheme: { name: type.defaultColorTheme, params: colorParams },
            sizeTheme: { name: type.defaultSizeTheme, params: sizeParams }
        })
    }
}
export { StructureRepresentation3D }
type StructureRepresentation3D = typeof StructureRepresentation3D
const StructureRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'structure-representation-3d',
    display: '3D Representation',
    from: SO.Molecule.Structure,
    to: SO.Molecule.Representation3D,
    params: (a, ctx: PluginContext) => {
        const { registry, themeCtx } = ctx.structureRepresentation
        const type = registry.get(registry.default.name);

        if (!a) {
            return {
                type: PD.Mapped<any>(
                    registry.default.name,
                    registry.types,
                    name => PD.Group<any>(registry.get(name).getParams(themeCtx, Structure.Empty))),
                colorTheme: PD.Mapped<any>(
                    type.defaultColorTheme,
                    themeCtx.colorThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({ structure: Structure.Empty }))
                ),
                sizeTheme: PD.Mapped<any>(
                    type.defaultSizeTheme,
                    themeCtx.sizeThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({ structure: Structure.Empty }))
                )
            }
        }

        const dataCtx = { structure: a.data }
        return ({
            type: PD.Mapped<any>(
                registry.default.name,
                registry.types,
                name => PD.Group<any>(registry.get(name).getParams(themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme,
                themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx))
            ),
            sizeTheme: PD.Mapped<any>(
                type.defaultSizeTheme,
                themeCtx.sizeThemeRegistry.types,
                name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx))
            )
        })
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        // TODO: allow for small molecules
        return oldParams.type.name === newParams.type.name;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const provider = plugin.structureRepresentation.registry.get(params.type.name)
            const props = params.type.params || {}
            const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.structureRepresentation.themeCtx }, provider.getParams)
            repr.setTheme(createTheme(plugin.structureRepresentation.themeCtx, { structure: a.data }, params))
            // TODO set initial state, repr.setState({})
            await repr.createOrUpdate(props, a.data).runInContext(ctx);
            return new SO.Molecule.Representation3D(repr, { label: provider.label });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;
            const props = { ...b.data.props, ...newParams.type.params }
            b.data.setTheme(createTheme(plugin.structureRepresentation.themeCtx, { structure: a.data }, newParams))
            await b.data.createOrUpdate(props, a.data).runInContext(ctx);
            return StateTransformer.UpdateResult.Updated;
        });
    }
});

export { ExplodeStructureRepresentation3D }
type ExplodeStructureRepresentation3D = typeof ExplodeStructureRepresentation3D
const ExplodeStructureRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'explode-structure-representation-3d',
    display: 'Explode 3D Representation',
    from: SO.Molecule.Representation3D,
    to: ExplodeRepresentation3D.Obj,
    params: ExplodeRepresentation3D.Params
})({
    canAutoUpdate() {
        return true;
    },
    apply({ params }, plugin: PluginContext) {
        return new ExplodeRepresentation3D.Obj(new ExplodeRepresentation3D.Behavior(plugin, params), { label: `Explosion T = ${params.t.toFixed(2)}` });
    },
    update({ b, newParams }) {
        return Task.create('Update Explosion', async () => {
            const updated = await b.data.update(newParams);
            b.label = `Explosion T = ${newParams.t.toFixed(2)}`;
            return updated ? StateTransformer.UpdateResult.Updated : StateTransformer.UpdateResult.Unchanged;
        });
    }
});

//

export namespace VolumeRepresentation3DHelpers {
    export function getDefaultParams(ctx: PluginContext, name: BuiltInVolumeRepresentationsName, volume: VolumeData, volumeParams?: Partial<PD.Values<VolumeParams>>): StateTransformer.Params<VolumeRepresentation3D> {
        const type = ctx.volumeRepresentation.registry.get(name);

        const themeDataCtx = { volume };
        const colorParams = ctx.volumeRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).getParams(themeDataCtx);
        const sizeParams = ctx.volumeRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        const volumeDefaultParams = PD.getDefaultValues(type.getParams(ctx.volumeRepresentation.themeCtx, volume))
        return ({
            type: { name, params: volumeParams ? { ...volumeDefaultParams, ...volumeParams } : volumeDefaultParams },
            colorTheme: { name: type.defaultColorTheme, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function getDefaultParamsStatic(ctx: PluginContext, name: BuiltInVolumeRepresentationsName, volumeParams?: Partial<PD.Values<PD.Params>>, colorName?: BuiltInColorThemeName, colorParams?: Partial<ColorTheme.Props>, sizeName?: BuiltInSizeThemeName, sizeParams?: Partial<SizeTheme.Props>): StateTransformer.Params<VolumeRepresentation3D> {
        const type = ctx.volumeRepresentation.registry.get(name);
        const colorType = ctx.volumeRepresentation.themeCtx.colorThemeRegistry.get(colorName || type.defaultColorTheme);
        const sizeType = ctx.volumeRepresentation.themeCtx.sizeThemeRegistry.get(sizeName || type.defaultSizeTheme);
        return ({
            type: { name, params: volumeParams ? { ...type.defaultValues, ...volumeParams } : type.defaultValues },
            colorTheme: { name: type.defaultColorTheme, params: colorParams ? { ...colorType.defaultValues, ...colorParams } : colorType.defaultValues },
            sizeTheme: { name: type.defaultSizeTheme, params: sizeParams ? { ...sizeType.defaultValues, ...sizeParams } : sizeType.defaultValues }
        })
    }

    export function getDescription(props: any) {
        return props.isoValue && VolumeIsoValue.toString(props.isoValue)
    }
}
export { VolumeRepresentation3D }
type VolumeRepresentation3D = typeof VolumeRepresentation3D
const VolumeRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'volume-representation-3d',
    display: '3D Representation',
    from: SO.Volume.Data,
    to: SO.Volume.Representation3D,
    params: (a, ctx: PluginContext) => {
        const { registry, themeCtx } = ctx.volumeRepresentation
        const type = registry.get(registry.default.name);

        if (!a) {
            return {
                type: PD.Mapped<any>(
                    registry.default.name,
                    registry.types,
                    name => PD.Group<any>(registry.get(name).getParams(themeCtx, VolumeData.One ))),
                colorTheme: PD.Mapped<any>(
                    type.defaultColorTheme,
                    themeCtx.colorThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams({ volume: VolumeData.One }))
                ),
                sizeTheme: PD.Mapped<any>(
                    type.defaultSizeTheme,
                    themeCtx.sizeThemeRegistry.types,
                    name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams({ volume: VolumeData.One }))
                )
            }
        }

        const dataCtx = { volume: a.data }
        return ({
            type: PD.Mapped<any>(
                registry.default.name,
                registry.types,
                name => PD.Group<any>(registry.get(name).getParams(themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme,
                themeCtx.colorThemeRegistry.getApplicableTypes(dataCtx),
                name => PD.Group<any>(themeCtx.colorThemeRegistry.get(name).getParams(dataCtx))
            ),
            sizeTheme: PD.Mapped<any>(
                type.defaultSizeTheme,
                themeCtx.sizeThemeRegistry.types,
                name => PD.Group<any>(themeCtx.sizeThemeRegistry.get(name).getParams(dataCtx))
            )
        })
    }
})({
    canAutoUpdate({ oldParams, newParams }) {
        // TODO: allow for small molecules
        return oldParams.type.name === newParams.type.name;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Volume Representation', async ctx => {
            const provider = plugin.volumeRepresentation.registry.get(params.type.name)
            const props = params.type.params || {}
            const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.volumeRepresentation.themeCtx }, provider.getParams)
            repr.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: a.data }, params))
            // TODO set initial state, repr.setState({})
            await repr.createOrUpdate(props, a.data).runInContext(ctx);
            return new SO.Volume.Representation3D(repr, { label: provider.label, description: VolumeRepresentation3DHelpers.getDescription(props) });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Volume Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return StateTransformer.UpdateResult.Recreate;
            const props = { ...b.data.props, ...newParams.type.params }
            b.data.setTheme(createTheme(plugin.volumeRepresentation.themeCtx, { volume: a.data }, newParams))
            await b.data.createOrUpdate(props, a.data).runInContext(ctx);
            b.description = VolumeRepresentation3DHelpers.getDescription(props)
            return StateTransformer.UpdateResult.Updated;
        });
    }
});