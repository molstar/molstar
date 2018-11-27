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
import { UnitsMeshParams } from 'mol-repr/structure/units-visual';

export namespace StructureRepresentation3DHelpers {
    export function getDefaultParams(ctx: PluginContext, name: BuiltInStructureRepresentationsName, structure: Structure, meshParams?: Partial<PD.Values<UnitsMeshParams>>): Transformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(name);

        const themeDataCtx = { structure };
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).getParams(themeDataCtx);
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        return ({
            type: { name, params: meshParams ? { ...type.defaultValues, ...meshParams } : type.defaultValues },
            colorTheme: { name: type.defaultColorTheme, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
        })
    }

    export function getDefaultParamsStatic(ctx: PluginContext, name: BuiltInStructureRepresentationsName, meshParams?: Partial<PD.Values<UnitsMeshParams>>): Transformer.Params<StructureRepresentation3D> {
        const type = ctx.structureRepresentation.registry.get(name);

        // TODO: there should be "static default properties" for the themes.
        const themeDataCtx = { };
        const colorParams = ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(type.defaultColorTheme).getParams(themeDataCtx);
        const sizeParams = ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(type.defaultSizeTheme).getParams(themeDataCtx)
        return ({
            type: { name, params: meshParams ? { ...type.defaultValues, ...meshParams } : type.defaultValues },
            colorTheme: { name: type.defaultColorTheme, params: PD.getDefaultValues(colorParams) },
            sizeTheme: { name: type.defaultSizeTheme, params: PD.getDefaultValues(sizeParams) }
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
        const type = ctx.structureRepresentation.registry.get(ctx.structureRepresentation.registry.default.name);
        return ({
            type: PD.Mapped<any>(
                ctx.structureRepresentation.registry.default.name,
                ctx.structureRepresentation.registry.types,
                name => PD.Group<any>(ctx.structureRepresentation.registry.get(name).getParams(ctx.structureRepresentation.themeCtx, a.data))),
            colorTheme: PD.Mapped<any>(
                type.defaultColorTheme,
                ctx.structureRepresentation.themeCtx.colorThemeRegistry.types,
                name => PD.Group<any>(ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(name).getParams({ structure: a.data }))
            ),
            sizeTheme: PD.Mapped<any>(
                type.defaultSizeTheme,
                ctx.structureRepresentation.themeCtx.sizeThemeRegistry.types,
                name => PD.Group<any>(ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(name).getParams({ structure: a.data }))
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