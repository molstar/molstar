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

export { CreateStructureRepresentation }
namespace CreateStructureRepresentation {
    export interface Params {
        type: { name: string, params: any /** todo is there "common type" */ }
    }
}
const CreateStructureRepresentation = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Representation3D, CreateStructureRepresentation.Params>({
    name: 'create-structure-representation',
    display: { name: 'Create 3D Representation' },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Representation3D],
    params: {
        default: (a, ctx: PluginContext) => ({
            type: {
                name: ctx.structureReprensentation.registry.default.name,
                params: PD.getDefaultValues(ctx.structureReprensentation.registry.default.provider.getParams(ctx.structureReprensentation.themeCtx, a.data))
            }
        }),
        definition: (a, ctx: PluginContext) => ({
            type: PD.Mapped('Type', '',
                ctx.structureReprensentation.registry.default.name,
                ctx.structureReprensentation.registry.types,
                name => ctx.structureReprensentation.registry.get(name)!.getParams(ctx.structureReprensentation.themeCtx, a.data))
        })
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const repr = plugin.structureReprensentation.registry.create(params.type.name, plugin.structureReprensentation.themeCtx, a.data)
            await repr.createOrUpdate({ webgl: plugin.canvas3d.webgl, ...plugin.structureReprensentation.themeCtx }, params.type.params || {}, {}, a.data).runInContext(ctx);
            return new SO.Molecule.Representation3D(repr, { label: params.type.name });
        });
    },
    update({ a, b, oldParams, newParams }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            if (newParams.type.name !== oldParams.type.name) return Transformer.UpdateResult.Recreate;

            await b.data.createOrUpdate({ webgl: plugin.canvas3d.webgl, ...plugin.structureReprensentation.themeCtx }, { ...b.data.props, ...newParams.type.params }, {}, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});