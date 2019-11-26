/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginCommands } from '../../../mol-plugin/command';

export function registerDefault(ctx: PluginContext) {
    Canvas3DSetSettings(ctx);
    InteractivitySetProps(ctx);
}

export function Canvas3DSetSettings(ctx: PluginContext) {
    PluginCommands.Canvas3D.SetSettings.subscribe(ctx, e => {
        ctx.canvas3d?.setProps(e.settings);
        ctx.events.canvas3d.settingsUpdated.next();
    })
}

export function InteractivitySetProps(ctx: PluginContext) {
    PluginCommands.Interactivity.SetProps.subscribe(ctx, e => {
        ctx.interactivity.setProps(e.props);
        ctx.events.interactivity.propsUpdated.next();
    })
}
