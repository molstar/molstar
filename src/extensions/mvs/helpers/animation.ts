/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { produce } from 'immer';
import { Snapshot } from '../mvs-data';
import { Tree } from '../tree/generic/tree-schema';
import { clamp, lerp } from '../../../mol-math/interpolate';
import { MVSAnimationEasing, MVSAnimationNode, MVSAnimationSchema } from '../tree/mvs/mvs-animation';
import { MVSTree } from '../tree/mvs/mvs-tree';
import * as EasingFns from '../../../mol-math/easing';
import { addDefaults } from '../tree/generic/tree-utils';
import { RuntimeContext } from '../../../mol-task';
import { EPSILON, Mat3, Mat4, Quat, Vec3 } from '../../../mol-math/linear-algebra';


export async function generateStateTransition(ctx: RuntimeContext, snapshot: Snapshot) {
    if (!snapshot.animation) return undefined;

    const tree = addDefaults(snapshot.animation, MVSAnimationSchema);
    const transitions = tree.children?.filter(child => child.kind === 'interpolate');
    if (!transitions?.length) return undefined;

    const duration = Math.max(...transitions.map(t => (t.params.start_ms ?? 0) + t.params.duration_ms));

    const frames: MVSTree[] = [];
    const dt = tree.params?.frame_time_ms ?? (1000 / 60);
    const N = Math.ceil(duration / dt);

    for (let i = 0; i <= N; i++) {
        const t = i * dt;
        const root = createSnapshot(snapshot.root, transitions, t);
        frames.push(root);

        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'Generating transition...' });
        }
    }

    return { tree, frametimeMs: dt, frames };
}

const EasingFnMap: Record<MVSAnimationEasing, (t: number) => number> = {
    'linear': t => t,
    'bounce-in': EasingFns.bounceIn,
    'bounce-out': EasingFns.bounceOut,
    'bounce-in-out': EasingFns.bounceInOut,
    'circle-in': EasingFns.circleIn,
    'circle-out': EasingFns.circleOut,
    'circle-in-out': EasingFns.circleInOut,
    'cubic-in': EasingFns.cubicIn,
    'cubic-out': EasingFns.cubicOut,
    'cubic-in-out': EasingFns.cubicInOut,
    'exp-in': EasingFns.expIn,
    'exp-out': EasingFns.expOut,
    'exp-in-out': EasingFns.expInOut,
    'quad-in': EasingFns.quadIn,
    'quad-out': EasingFns.quadOut,
    'quad-in-out': EasingFns.quadInOut,
    'sin-in': EasingFns.sinIn,
    'sin-out': EasingFns.sinOut,
    'sin-in-out': EasingFns.sinInOut,
};

function createSnapshot(tree: MVSTree, transitions: MVSAnimationNode<'interpolate'>[], time: number) {
    return produce(tree, (draft) => {
        for (const transition of transitions) {
            const node = findNode(draft, transition.params.target_ref);
            if (!node) continue;

            const target = transition.params.property[0] === 'custom' ? node?.custom : node?.params;
            if (!target) continue;

            const offset = transition.params.property[0] === 'custom' ? 1 : 0;
            const startTime = transition.params.start_ms ?? 0;
            const startValue: any = transition.params.from ?? select(target, transition.params.property, offset);
            const endTime = startTime + transition.params.duration_ms;
            const endValue: any = transition.params.to;

            if (time <= startTime) {
                assign(target, transition.params.property, startValue, offset);
                continue;
            }

            if (time >= endTime - EPSILON) {
                assign(target, transition.params.property, endValue, offset);
                continue;
            }

            const deltaT = endTime - startTime;
            const easing = EasingFnMap[transition.params.easing ?? 'linear'] ?? EasingFnMap['linear'];
            const t = easing(clamp((time - startTime) / deltaT, 0, 1));

            let next: any;
            if (transition.params.type === 'scalar') {
                next = interpolateScalar(startValue, endValue, t, transition.params.noise_magnitude ?? 0);
            } else if (transition.params.type === 'vec3') {
                next = interpolateVec3(startValue, endValue, t, transition.params.noise_magnitude ?? 0, !!transition.params.spherical);
            } else if (transition.params.type === 'rotation_matrix') {
                next = interpolateRotation(startValue, endValue, t, transition.params.noise_magnitude ?? 0);
            }

            assign(target, transition.params.property, next, offset);
        }
    });
}

function interpolateScalar(start: number, end: number, t: number, noise: number) {
    let v = lerp(start, end, t);
    if (noise) {
        v += (Math.random() - 0.5) * noise;
    }
    return v;
}

const Vec3Noise = Vec3();
function interpolateVec3(start: Vec3, end: Vec3, t: number, noise: number, isSpherical: boolean) {
    const v = isSpherical
        ? Vec3.slerp(Vec3(), start, end, t)
        : Vec3.lerp(Vec3(), start, end, t);

    if (noise) {
        Vec3.random(Vec3Noise, noise);
        Vec3.add(v, v, Vec3Noise);
    }
    return v;
}

const RotationState = {
    start: Quat(),
    end: Quat(),
    v: Quat(),
    noise: Quat(),
    axis: Vec3(),
    temp: Mat4(),
};
function interpolateRotation(start: Mat3, end: Mat3, t: number, noise: number) {
    Quat.fromMat3(RotationState.start, start);
    Quat.fromMat3(RotationState.end, end);
    Quat.slerp(RotationState.v, RotationState.start, RotationState.end, t);
    if (noise) {
        Vec3.random(RotationState.axis, 1);
        Quat.setAxisAngle(RotationState.noise, RotationState.axis, 2 * Math.PI * noise * (Math.random() - 0.5));
        Quat.multiply(RotationState.v, RotationState.noise, RotationState.v);
    }
    Mat4.fromQuat(RotationState.temp, RotationState.v);
    return Mat3.fromMat4(Mat3(), RotationState.temp);
}

function select(params: any, path: string | (string | number)[], offset: number) {
    if (typeof path === 'string') {
        return params?.[path];
    }

    let f = params;
    for (let i = offset; i < path.length; i++) {
        if (!f) break;
        f = f[path[i]];
    }
    return f;
}

function assign(params: any, path: string | (string | number)[], value: any, offset: number) {
    if (!params) return;

    if (typeof path === 'string') {
        params[path] = value;
        return;
    }

    let f = params;
    for (let i = offset; i < path.length; i++) {
        if (!f) break;
        if (i === path.length - 1) {
            f[path[i]] = value;
        } else {
            f = f[path[i]];
        }
    }
}

function findNode(tree: Tree, ref: string): Tree | undefined {
    if (tree.ref === ref) return tree;
    if (!tree.children) return undefined;
    for (const child of tree.children) {
        const result = findNode(child, ref);
        if (result) return result;
    }
    return undefined;
}