/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { createTheme } from 'mol-theme/theme';
import { BuiltInStructureRepresentationsName } from 'mol-repr/structure/registry';
import { Structure } from 'mol-model/structure';
import { StructureParams } from 'mol-repr/structure/representation';
import { ShapeRepresentation } from 'mol-repr/shape/representation';

export namespace StructureRepresentation3DHelpers {
    export function getDefaultParams(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structure: Structure, structureParams?: Partial<PD.Values<StructureParams>>): Transformer.Params<StructureRepresentation3D> {
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

    export function getDefaultParamsStatic(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structureParams?: Partial<PD.Values<StructureParams>>): Transformer.Params<StructureRepresentation3D> {
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
            if (newParams.type.name !== oldParams.type.name) return Transformer.UpdateResult.Recreate;
            const props = { ...b.data.props, ...newParams.type.params }
            b.data.setTheme(createTheme(plugin.structureRepresentation.themeCtx, { structure: a.data }, newParams))
            await b.data.createOrUpdate(props, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});

export { ShapeRepresentation3D }
type ShapeRepresentation3D = typeof ShapeRepresentation3D
const ShapeRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'shape-representation-3d',
    display: '3D Representation',
    from: SO.Shape.Provider,
    to: SO.Shape.Representation3D,
    params: (a, ctx: PluginContext) => {
        return { }
    }
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Shape Representation', async ctx => {
            const props = { ...PD.getDefaultValues(a.data.geometryUtils.Params), params }
            const repr = ShapeRepresentation(a.data.getShape, a.data.geometryUtils)
            // TODO set initial state, repr.setState({})
            await repr.createOrUpdate(props, a.data.data).runInContext(ctx);
            return new SO.Shape.Representation3D(repr, { label: a.data.label });
        });
    }
});