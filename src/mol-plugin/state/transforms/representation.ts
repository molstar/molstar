/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export { StructureRepresentation3D }
namespace StructureRepresentation3D {
    export interface Params {
        type: { name: string, params: any /** todo is there "common type" */ },
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
            name => PD.Group<any>(ctx.structureRepresentation.registry.get(name).getParams(ctx.structureRepresentation.themeCtx, a.data)))
    }),
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const provider = plugin.structureRepresentation.registry.get(params.type.name)
            const repr = provider.factory({ webgl: plugin.canvas3d.webgl, ...plugin.structureRepresentation.themeCtx }, provider.getParams)
            await repr.createOrUpdate(params.type.params || {}, a.data).runInContext(ctx);
            return new SO.Molecule.Representation3D(repr, { label: provider.label });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return Transformer.UpdateResult.Recreate;

            await b.data.createOrUpdate({ ...b.data.props, ...newParams.type.params }, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});