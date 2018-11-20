/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from 'mol-plugin/context';
import { PluginCommands } from 'mol-plugin/command';

export function registerDefault(ctx: PluginContext) {
    Canvas3DSetSettings(ctx);
}

export function Canvas3DSetSettings(ctx: PluginContext) {
    PluginCommands.Canvas3D.SetSettings.subscribe(ctx, e => {
        ctx.canvas3d.setProps(e.settings);
        ctx.events.canvad3d.settingsUpdated.next();
    })
}
