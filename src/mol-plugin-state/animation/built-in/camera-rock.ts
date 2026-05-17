/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { clamp } from '../../../mol-math/interpolate';
import { Quat } from '../../../mol-math/linear-algebra/3d/quat';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { degToRad } from '../../../mol-math/misc';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateAnimation } from '../model';

const _dir = Vec3(), _axis = Vec3(), _rot = Quat(), _up = Vec3(), _side = Vec3();

type State = { snapshot: Camera.Snapshot };

export const AnimateCameraRock = PluginStateAnimation.create({
    name: 'built-in.animate-camera-rock',
    display: { name: 'Camera Rock', description: 'Rock the 3D scene around the x-axis in view space' },
    isExportable: true,
    params: () => ({
        durationInMs: PD.Numeric(4000, { min: 100, max: 20000, step: 100 }),
        speed: PD.Numeric(1, { min: 1, max: 10, step: 1 }, { description: 'How many times to rock from side to side.' }),
        angle: PD.Numeric(10, { min: 0, max: 180, step: 1 }, { description: 'How many degrees to rotate in each direction.' }),
        axis: PD.Vec3(Vec3.create(0, -1, 0), {}, { description: 'Axis of rotation in camera space' }),
    }),
    initialState: (p, ctx) => ({ snapshot: ctx.canvas3d!.camera.getSnapshot() }) as State,
    getDuration: p => ({ kind: 'fixed', durationMs: p.durationInMs }),
    teardown: (_, state: State, ctx) => {
        ctx.canvas3d?.requestCameraReset({ snapshot: state.snapshot, durationMs: 0 });
    },

    async apply(animState: State, t, ctx) {
        if (t.current === 0) {
            return { kind: 'next', state: animState };
        }

        const snapshot = animState.snapshot;
        if (snapshot.radiusMax < 0.0001) {
            return { kind: 'finished' };
        }

        const phase = t.animation
            ? t.animation?.currentFrame / (t.animation.frameCount + 1)
            : clamp(t.current / ctx.params.durationInMs, 0, 1);
        const angle = Math.sin(phase * ctx.params.speed * Math.PI * 2) * degToRad(ctx.params.angle);

        Vec3.sub(_dir, snapshot.position, snapshot.target);

        // Transform axis from camera space to world space
        Vec3.normalize(_axis, _dir); // Z = view direction
        Vec3.normalize(_up, snapshot.up); // Y = up
        Vec3.cross(_side, _up, _axis); // X = right
        Vec3.normalize(_side, _side);
        const a = ctx.params.axis;
        Vec3.set(_axis,
            a[0] * _side[0] + a[1] * _up[0] + a[2] * _axis[0],
            a[0] * _side[1] + a[1] * _up[1] + a[2] * _axis[1],
            a[0] * _side[2] + a[1] * _up[2] + a[2] * _axis[2]
        );
        Vec3.normalize(_axis, _axis);

        Quat.setAxisAngle(_rot, _axis, angle);
        Vec3.transformQuat(_dir, _dir, _rot);
        Vec3.transformQuat(_up, snapshot.up, _rot);
        const position = Vec3.add(Vec3(), snapshot.target, _dir);
        ctx.plugin.canvas3d?.requestCameraReset({ snapshot: { ...snapshot, position, up: _up }, durationMs: 0 });

        if (phase >= 0.99999) {
            return { kind: 'finished' };
        }

        return { kind: 'next', state: animState };
    }
});