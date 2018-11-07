/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';
//import { CartoonRepresentation, DefaultCartoonProps } from 'mol-repr/structure/representation/cartoon';
import { BallAndStickRepresentation, DefaultBallAndStickProps } from 'mol-repr/structure/representation/ball-and-stick';

export { CreateStructureRepresentation }
namespace CreateStructureRepresentation { export interface Params { } }
const CreateStructureRepresentation = PluginStateTransform.Create<SO.Structure, SO.StructureRepresentation3D, CreateStructureRepresentation.Params>({
    name: 'create-structure-representation',
    display: { name: 'Create 3D Representation' },
    from: [SO.Structure],
    to: [SO.StructureRepresentation3D],
    apply({ a, params }) {
        return Task.create('Structure Representation', async ctx => {
            const repr = BallAndStickRepresentation(); // CartoonRepresentation();
            await repr.createOrUpdate({ /* TODO add `webgl: WebGLContext` */ }, DefaultBallAndStickProps, a.data).runInContext(ctx);
            return new SO.StructureRepresentation3D({ label: 'Cartoon' }, { repr });
        });
    },
    update({ a, b }) {
        return Task.create('Structure Representation', async ctx => {
            await b.data.repr.createOrUpdate({ /* TODO add `webgl: WebGLContext` */ }, b.data.repr.props, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});