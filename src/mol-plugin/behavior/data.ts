
/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from './behavior';
import { PluginCommands } from 'mol-plugin/command';

// export class SetCurrentObject implements PluginBehavior<undefined> {
//     private sub: PluginCommand.Subscription | undefined = void 0;

//     register(): void {
//         this.sub = PluginCommands.Data.SetCurrentObject.subscribe(this.ctx, ({ ref }) => this.ctx.state.data.setCurrent(ref));
//     }
//     unregister(): void {
//         if (this.sub) this.sub.unsubscribe();
//         this.sub = void 0;
//     }

//     constructor(private ctx: PluginContext) { }
// }

export const SetCurrentObject = PluginBehavior.create({
    name: 'set-current-data-object-behavior',
    ctor: PluginBehavior.commandHandler(PluginCommands.Data.SetCurrentObject, ({ ref }, ctx) => ctx.state.data.setCurrent(ref)),
    display: { name: 'Set Current Handler' }
});

export const Update = PluginBehavior.create({
    name: 'update-data-behavior',
    ctor: PluginBehavior.commandHandler(PluginCommands.Data.Update, ({ tree }, ctx) => ctx.runTask(ctx.state.data.update(tree))),
    display: { name: 'Update Data Handler' }
});