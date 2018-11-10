/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../objects';
import { PluginStateObject as SO } from '../objects';
import { BallAndStickRepresentation, DefaultBallAndStickProps } from 'mol-repr/structure/representation/ball-and-stick';
import { PluginContext } from 'mol-plugin/context';

export { CreateStructureRepresentation }
namespace CreateStructureRepresentation { export interface Params { } }
const CreateStructureRepresentation = PluginStateTransform.Create<SO.Molecule.Structure, SO.Molecule.Representation3D, CreateStructureRepresentation.Params>({
    name: 'create-structure-representation',
    display: { name: 'Create 3D Representation' },
    from: [SO.Molecule.Structure],
    to: [SO.Molecule.Representation3D],
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            const repr = BallAndStickRepresentation(); // CartoonRepresentation();
            await repr.createOrUpdate({ webgl: plugin.canvas3d.webgl }, DefaultBallAndStickProps, a.data).runInContext(ctx);
            return new SO.Molecule.Representation3D(repr);
        });
    },
    update({ a, b }, plugin: PluginContext) {
        return Task.create('Structure Representation', async ctx => {
            await b.data.createOrUpdate({ webgl: plugin.canvas3d.webgl }, b.data.props, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});