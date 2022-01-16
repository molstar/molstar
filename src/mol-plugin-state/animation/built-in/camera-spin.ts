/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { clamp } from '../../../mol-math/interpolate';
import { Quat } from '../../../mol-math/linear-algebra/3d/quat';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateAnimation } from '../model';

const _dir = Vec3(), _axis = Vec3(), _rot = Quat();

type State = { snapshot: Camera.Snapshot };

export const AnimateCameraSpin = PluginStateAnimation.create({
    name: 'built-in.animate-camera-spin',
    display: { name: 'Camera Spin', description: 'Spin the 3D scene around the x-axis in view space' },
    isExportable: true,
    params: () => ({
        durationInMs: PD.Numeric(4000, { min: 100, max: 20000, step: 100 }),
        speed: PD.Numeric(1, { min: 1, max: 10, step: 1 }, { description: 'How many times to spin in the specified duration.' }),
        direction: PD.Select<'cw' | 'ccw'>('cw', [['cw', 'Clockwise'], ['ccw', 'Counter Clockwise']], { cycle: true })
    }),
    initialState: (_, ctx) => ({ snapshot: ctx.canvas3d?.camera.getSnapshot()! }) as State,
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
        const angle = 2 * Math.PI * phase * ctx.params.speed * (ctx.params.direction === 'ccw' ? -1 : 1);

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