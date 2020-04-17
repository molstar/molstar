/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginCommands } from '../../commands';

export function registerDefault(ctx: PluginContext) {
    Reset(ctx);
    Focus(ctx);
    SetSnapshot(ctx);
}

export function Reset(ctx: PluginContext) {
    PluginCommands.Camera.Reset.subscribe(ctx, options => {
        ctx.managers.camera.reset(options?.snapshot, options?.durationMs);
    });
}

export function SetSnapshot(ctx: PluginContext) {
    PluginCommands.Camera.SetSnapshot.subscribe(ctx, ({ snapshot, durationMs }) => {
        ctx.managers.camera.setSnapshot(snapshot, durationMs);
    });
}

export function Focus(ctx: PluginContext) {
    PluginCommands.Camera.Focus.subscribe(ctx, ({ center, radius, durationMs }) => {
        ctx.managers.camera.focusSphere({ center, radius }, { durationMs });
        ctx.events.canvas3d.settingsUpdated.next();
    });
}