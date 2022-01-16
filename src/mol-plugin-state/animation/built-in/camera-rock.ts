/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { clamp, smoothstep } from '../../../mol-math/interpolate';
import { Quat } from '../../../mol-math/linear-algebra/3d/quat';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { degToRad } from '../../../mol-math/misc';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateAnimation } from '../model';

const _dir = Vec3(), _axis = Vec3(), _rot = Quat();

type State = {
    snapshot: Camera.Snapshot,
    angles: number[]
};

export const AnimateCameraRock = PluginStateAnimation.create({
    name: 'built-in.animate-camera-rock',
    display: { name: 'Camera Rock', description: 'Rock the 3D scene around the x-axis in view space' },
    isExportable: true,
    params: () => ({
        durationInMs: PD.Numeric(4000, { min: 100, max: 20000, step: 100 }),
        angle: PD.Numeric(10, { min: 0, max: 90, step: 1 }, { description: 'How many degrees to rotate in each direction.' }),
    }),
    initialState: (p, ctx) => {
        const angles: number[] = [];
        const frameSpeed = 1 / 1000;
        const deltaT = 1000 / 30;

        // TODO get rid of the 3.3 factor (compensates for using `smoothstep`)
        const maxAngle = degToRad(p.angle * 3.3);

        let angleSum = 0;
        let direction = 1;
        let zeroAngleCount = 0;
        let prevAngle = 0;
        while (true) {
            const alpha = smoothstep(0, 1, Math.abs(angleSum) / maxAngle);
            const rockSpeed = 60 * Math.min(Math.abs(deltaT), 1000 / 8) / 1000 * frameSpeed;
            angleSum += Math.abs(rockSpeed);
            const angle = prevAngle + rockSpeed * direction * (1.1 - alpha);
            angles.push(angle);
            if (Math.sign(prevAngle) !== Math.sign(angle)) {
                zeroAngleCount += 1;
                if (zeroAngleCount === 3) break;
            }
            prevAngle = angle;
            if (angleSum >= maxAngle) {
                direction *= -1;
                angleSum = -maxAngle;
            }
        }

        return {
            snapshot: ctx.canvas3d!.camera.getSnapshot(),
            angles
        } as State;
    },
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
        const angle = animState.angles[Math.round(phase * (animState.angles.length - 1))];

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