/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { PluginCommands } from '@molstar/plugin/commands';
import { DefaultCanvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';

export function registerDefault(ctx: PluginContext) {
    Canvas3DSetSettings(ctx);
}

export function Canvas3DSetSettings(ctx: PluginContext) {
    PluginCommands.Canvas3D.ResetSettings.subscribe(ctx, () => {
        ctx.canvas3d?.setProps(DefaultCanvas3DParams);
        ctx.events.canvas3d.settingsUpdated.next(void 0);
    });

    PluginCommands.Canvas3D.SetSettings.subscribe(ctx, e => {
        if (!ctx.canvas3d) return;

        ctx.canvas3d?.setProps(e.settings);
        ctx.events.canvas3d.settingsUpdated.next(void 0);
    });
}
