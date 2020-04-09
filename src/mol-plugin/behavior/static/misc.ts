/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginCommands } from '../../commands';
import { DefaultCanvas3DParams } from '../../../mol-canvas3d/canvas3d';

export function registerDefault(ctx: PluginContext) {
    Canvas3DSetSettings(ctx);
}

export function Canvas3DSetSettings(ctx: PluginContext) {
    PluginCommands.Canvas3D.ResetSettings.subscribe(ctx, () => {
        ctx.canvas3d?.setProps(DefaultCanvas3DParams);
    });

    PluginCommands.Canvas3D.SetSettings.subscribe(ctx, e => {
        if (!ctx.canvas3d) return;

        ctx.canvas3d?.setProps(e.settings);
        ctx.events.canvas3d.settingsUpdated.next();
    });
}
