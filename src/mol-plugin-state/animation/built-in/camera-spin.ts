/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { clamp } from '../../../mol-math/interpolate';
import { Quat, Vec3 } from '../../../mol-math/linear-algebra/3d';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateAnimation } from '../model';

const _dir = Vec3(), _axis = Vec3(), _rot = Quat();

export const AnimateCameraSpin = PluginStateAnimation.create({
    name: 'built-in.animate-camera-spin',
    display: { name: 'Camera Spin' },
    params: () => ({
        durationInMs: PD.Numeric(4000, { min: 100, max: 20000, step: 100 }),
        direction: PD.Select<'cw' | 'ccw'>('cw', [['cw', 'Clockwise'], ['ccw', 'Counter Clockwise']], { cycle: true }),
        skipLastFrame: PD.Boolean(true)
    }),
    initialState: () => ({ }),
    async apply(animState: { snapshot: Camera.Snapshot }, t, ctx) {
        if (t.current === 0) {
            return { kind: 'next', state: animState };
        } else if (ctx.params.skipLastFrame && t.current >= ctx.params.durationInMs) {
            return { kind: 'finished' };
        }

        const camera = ctx.plugin.canvas3d?.camera!;
        if (camera.state.radiusMax < 0.0001) {
            return { kind: 'finished' };
        }

        const delta = clamp((t.current - t.lastApplied) / ctx.params.durationInMs, 0, 1);
        const angle = 2 * Math.PI * delta * (ctx.params.direction === 'ccw' ? -1 : 1);

        Vec3.sub(_dir, camera.position, camera.target);
        Vec3.normalize(_axis, camera.up);
        Quat.setAxisAngle(_rot, _axis, angle);
        Vec3.transformQuat(_dir, _dir, _rot);
        const position = Vec3.add(Vec3(), camera.target, _dir);
        ctx.plugin.canvas3d?.requestCameraReset({ snapshot: { position }, durationMs: 0 });

        if (t.current >= ctx.params.durationInMs) {
            return { kind: 'finished' };
        }

        return { kind: 'next', state: animState };
    }
});