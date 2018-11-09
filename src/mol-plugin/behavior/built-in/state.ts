/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '../behavior';
import { PluginCommands } from '../../command';

export const SetCurrentObject = PluginBehavior.create({
    name: 'set-current-data-object-behavior',
    ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.SetCurrentObject, ({ state, ref }, ctx) => state.setCurrent(ref)),
    display: { name: 'Set Current Handler', group: 'Data' }
});

export const Update = PluginBehavior.create({
    name: 'update-data-behavior',
    ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.Update, ({ state, tree }, ctx) => ctx.runTask(state.update(tree))),
    display: { name: 'Update Data Handler', group: 'Data' }
});

export const RemoveObject = PluginBehavior.create({
    name: 'remove-object-data-behavior',
    ctor: PluginBehavior.simpleCommandHandler(PluginCommands.State.RemoveObject, ({ state, ref }, ctx) => {
        const tree = state.tree.build().delete(ref).getTree();
        return ctx.runTask(state.update(tree));
    }),
    display: { name: 'Remove Object Handler', group: 'Data' }
});