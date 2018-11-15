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
        type: { name: string, params: any /** todo is there "common type" */ },
        // TODO
        // colorTheme: { name: string, params: any /** todo is there "common type" */ }
        // sizeTheme: { name: string, params: any /** todo is there "common type" */ }
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
                params: ctx.structureReprensentation.registry.default.provider.defaultValues
            }
        }),
        definition: (a, ctx: PluginContext) => ({
            type: PD.Mapped(
                ctx.structureReprensentation.registry.default.name,
                ctx.structureReprensentation.registry.types,
                name => PD.Group(
                    ctx.structureReprensentation.registry.get(name).getParams(ctx.structureReprensentation.themeCtx, a.data),
                    { label: 'Params' }
                )
            )
        })
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const provider = plugin.structureReprensentation.registry.get(params.type.name)
            const repr = provider.factory(provider.getParams)
            await repr.createOrUpdate({ webgl: plugin.canvas3d.webgl, ...plugin.structureReprensentation.themeCtx }, params.type.params || {}, {}, a.data).runInContext(ctx);
            return new SO.Molecule.Representation3D(repr, { label: provider.label });
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