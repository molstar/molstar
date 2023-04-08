/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginCommands } from '../../commands';

export function registerDefault(ctx: PluginContext) {
    Reset(ctx);
    Focus(ctx);
    SetSnapshot(ctx);
    OrientAxes(ctx);
    ResetAxes(ctx);
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
        ctx.events.canvas3d.settingsUpdated.next(void 0);
    });
}

export function OrientAxes(ctx: PluginContext) {
    PluginCommands.Camera.OrientAxes.subscribe(ctx, ({ structures, durationMs }) => {
        ctx.managers.camera.orientAxes(structures, durationMs);
    });
}

export function ResetAxes(ctx: PluginContext) {
    PluginCommands.Camera.ResetAxes.subscribe(ctx, ({ durationMs }) => {
        ctx.managers.camera.resetAxes(durationMs);
    });
}
