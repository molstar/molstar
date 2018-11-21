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

export { StructureRepresentation3D }
namespace StructureRepresentation3D {
    export interface Params {
        type: { name: string, params: any /** TODO is there "common type" */ },
        colorTheme: { name: string, params: any /** TODO is there "common type" */ },
        sizeTheme: { name: string, params: any /** TODO is there "common type" */ },
    }
}
const StructureRepresentation3D = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Representation3D, StructureRepresentation3D.Params>({
    name: 'structure-representation-3d',
    display: { name: '3D Representation' },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Representation3D],
    params: (a, ctx: PluginContext) => ({
        type: PD.Mapped(
            ctx.structureRepresentation.registry.default.name,
            ctx.structureRepresentation.registry.types,
            name => PD.Group<any>(ctx.structureRepresentation.registry.get(name).getParams(ctx.structureRepresentation.themeCtx, a.data))),
        colorTheme: PD.Mapped(
            // TODO how to get a default color theme dependent on the repr type?
            ctx.structureRepresentation.themeCtx.colorThemeRegistry.default.name,
            ctx.structureRepresentation.themeCtx.colorThemeRegistry.types,
            name => PD.Group<any>(ctx.structureRepresentation.themeCtx.colorThemeRegistry.get(name).getParams({ structure: a.data }))
        ),
        sizeTheme: PD.Mapped(
            // TODO how to get a default size theme dependent on the repr type?
            ctx.structureRepresentation.themeCtx.sizeThemeRegistry.default.name,
            ctx.structureRepresentation.themeCtx.sizeThemeRegistry.types,
            name => PD.Group<any>(ctx.structureRepresentation.themeCtx.sizeThemeRegistry.get(name).getParams({ structure: a.data }))
        ),
    }),
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