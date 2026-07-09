/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginCommands } from '../../commands';

export function registerDefault(ctx: PluginContext) {
    Reset(ctx);
    Focus(ctx);
    FocusObject(ctx);
    SetSnapshot(ctx);
    OrientAxes(ctx);
    ResetAxes(ctx);
}

export function Reset(ctx: PluginContext) {
    PluginCommands.Camera.Reset.subscribe(ctx, ({ snapshot, durationMs, easing, shape }) => {
        ctx.managers.camera.reset(snapshot, durationMs, { easing, shape });
    });
}

export function SetSnapshot(ctx: PluginContext) {
    PluginCommands.Camera.SetSnapshot.subscribe(ctx, ({ snapshot, durationMs, easing, shape }) => {
        ctx.managers.camera.setSnapshot(snapshot, durationMs, { easing, shape });
    });
}

export function Focus(ctx: PluginContext) {
    PluginCommands.Camera.Focus.subscribe(ctx, ({ center, radius, durationMs, easing, shape }) => {
        ctx.managers.camera.focusSphere({ center, radius }, { durationMs, easing, shape });
        ctx.events.canvas3d.settingsUpdated.next(undefined);
    });
}

export function FocusObject(ctx: PluginContext) {
    PluginCommands.Camera.FocusObject.subscribe(ctx, options => {
        ctx.managers.camera.focusObject(options);
    });
}

export function OrientAxes(ctx: PluginContext) {
    PluginCommands.Camera.OrientAxes.subscribe(ctx, ({ structures, durationMs, easing, shape }) => {
        ctx.managers.camera.orientAxes(structures, durationMs, { easing, shape });
    });
}

export function ResetAxes(ctx: PluginContext) {
    PluginCommands.Camera.ResetAxes.subscribe(ctx, ({ durationMs, easing, shape }) => {
        ctx.managers.camera.resetAxes(durationMs, { easing, shape });
    });
}
