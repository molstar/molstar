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
    PluginCommands.Camera.Reset.subscribe(ctx, ({ snapshot, durationMs, easing, trajectory }) => {
        ctx.managers.camera.reset(snapshot, durationMs, { easing, trajectory });
    });
}

export function SetSnapshot(ctx: PluginContext) {
    PluginCommands.Camera.SetSnapshot.subscribe(ctx, ({ snapshot, durationMs, easing, trajectory }) => {
        ctx.managers.camera.setSnapshot(snapshot, durationMs, { easing, trajectory });
    });
}

export function Focus(ctx: PluginContext) {
    PluginCommands.Camera.Focus.subscribe(ctx, ({ center, radius, durationMs, easing, trajectory }) => {
        ctx.managers.camera.focusSphere({ center, radius }, { durationMs, easing, trajectory });
        ctx.events.canvas3d.settingsUpdated.next(undefined);
    });
}

export function FocusObject(ctx: PluginContext) {
    PluginCommands.Camera.FocusObject.subscribe(ctx, options => {
        ctx.managers.camera.focusObject(options);
    });
}

export function OrientAxes(ctx: PluginContext) {
    PluginCommands.Camera.OrientAxes.subscribe(ctx, ({ structures, durationMs, easing, trajectory }) => {
        ctx.managers.camera.orientAxes(structures, durationMs, { easing, trajectory });
    });
}

export function ResetAxes(ctx: PluginContext) {
    PluginCommands.Camera.ResetAxes.subscribe(ctx, ({ durationMs, easing, trajectory }) => {
        ctx.managers.camera.resetAxes(durationMs, { easing, trajectory });
    });
}
