/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

const _dir = Vec3(), _axis = Vec3(), _rot = Quat();

type State = { snapshot: Camera.Snapshot };

export const AnimateCameraRock = PluginStateAnimation.create({
    name: 'built-in.animate-camera-rock',
    display: { name: 'Camera Rock', description: 'Rock the 3D scene around the x-axis in view space' },
    isExportable: true,
    params: () => ({
        durationInMs: PD.Numeric(4000, { min: 100, max: 20000, step: 100 }),
        speed: PD.Numeric(1, { min: 1, max: 10, step: 1 }, { description: 'How many times to rock from side to side.' }),
        angle: PD.Numeric(10, { min: 0, max: 180, step: 1 }, { description: 'How many degrees to rotate in each direction.' }),
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
        Vec3.normalize(_axis, snapshot.up);
        Quat.setAxisAngle(_rot, _axis, angle);
        Vec3.transformQuat(_dir, _dir, _rot);
        const position = Vec3.add(Vec3(), snapshot.target, _dir);
        ctx.plugin.canvas3d?.requestCameraReset({ snapshot: { ...snapshot, position }, durationMs: 0 });

        if (phase >= 0.99999) {
            return { kind: 'finished' };
        }

        return { kind: 'next', state: animState };
    }
});