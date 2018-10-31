/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CartoonRepresentation, DefaultCartoonProps } from 'mol-geo/representation/structure/representation/cartoon';
import { Transformer } from 'mol-state';
import { Task } from 'mol-task';
import { PluginStateTransform } from '../base';
import { PluginStateObjects as SO } from '../objects';

export const CreateStructureRepresentation = PluginStateTransform.Create<SO.Structure, SO.StructureRepresentation3D, { }>({
    name: 'create-structure-representation',
    from: [SO.Structure],
    to: [SO.StructureRepresentation3D],
    defaultParams: () => ({ modelIndex: 0 }),
    apply({ a, params }) {
        return Task.create('Structure Representation', async ctx => {
            const repr = CartoonRepresentation();
            await repr.createOrUpdate({ ...DefaultCartoonProps }, a.data).runInContext(ctx);
            return new SO.StructureRepresentation3D({ label: 'Cartoon' }, { repr });
        });
    },
    update({ a, b }) {
        return Task.create('Structure Representation', async ctx => {
            await b.data.repr.createOrUpdate(b.data.repr.props, a.data).runInContext(ctx);
            return Transformer.UpdateResult.Updated;
        });
    }
});