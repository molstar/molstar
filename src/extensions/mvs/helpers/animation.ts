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
import { Color } from '../../../mol-util/color';
import { makeContinuousPaletteCheckpoints, MVSContinuousPaletteProps, MVSDiscretePaletteProps } from '../components/annotation-color-theme';
import { palettePropsFromMVSPalette } from '../load-helpers';
import { SortedArray } from '../../../mol-data/int';


export async function generateStateTransition(ctx: RuntimeContext, snapshot: Snapshot) {
    if (!snapshot.animation) return undefined;

    const tree = addDefaults(snapshot.animation, MVSAnimationSchema);
    const transitions = tree.children?.filter(child => child.kind === 'interpolate');
    if (!transitions?.length) return undefined;

    const duration = Math.max(...transitions.map(t => (t.params.start_ms ?? 0) + t.params.duration_ms));

    const frames: MVSTree[] = [];
    const dt = tree.params?.frame_time_ms ?? (1000 / 60);
    const N = Math.ceil(duration / dt);

    const cache = new Map<any, TransitionCacheEntry>();

    for (let i = 0; i <= N; i++) {
        const t = i * dt;
        const root = createSnapshot(snapshot.root, transitions, t, cache);
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

interface TransitionCacheEntry {
    paletteFn?: (value: number) => Color
}

function createSnapshot(tree: MVSTree, transitions: MVSAnimationNode<'interpolate'>[], time: number, cache: Map<any, TransitionCacheEntry>) {
    return produce(tree, (draft) => {
        for (const transition of transitions) {
            const node = findNode(draft, transition.params.target_ref);
            if (!node) continue;

            const target = transition.params.property[0] === 'custom' ? node?.custom : node?.params;
            if (!target) continue;

            if (transition.params.kind === 'color') {
                if (!cache.has(transition)) {
                    cache.set(transition, {
                        paletteFn: makePaletteFunction(transition)
                    });
                }
            }

            const paletteFn = cache.get(transition)?.paletteFn;

            const transformTarget = transition.params.kind === 'rotation_matrix' || transition.params.kind === 'vec3'
                ? transition.params.transform_matrix_target
                : undefined;

            const offset = transition.params.property[0] === 'custom' ? 1 : 0;
            const startTime = transition.params.start_ms ?? 0;
            const startValue: any = transition.params.kind === 'color'
                ? Color.toHexStyle(paletteFn!(0))
                : transition.params.from ?? select(target, transition.params.property, offset, transformTarget);
            const endTime = startTime + transition.params.duration_ms;
            const endValue: any = transition.params.kind === 'color'
                ? Color.toHexStyle(paletteFn!(1))
                : transition.params.to;

            if (time <= startTime) {
                assign(target, transition.params.property, startValue, offset, transformTarget);
                continue;
            }

            if (time >= endTime - EPSILON) {
                assign(target, transition.params.property, endValue, offset, transformTarget);
                continue;
            }

            const deltaT = endTime - startTime;
            const easing = EasingFnMap[transition.params.easing ?? 'linear'] ?? EasingFnMap['linear'];
            const t = easing(clamp((time - startTime) / deltaT, 0, 1));

            let next: any;
            if (transition.params.kind === 'scalar') {
                next = interpolateScalar(startValue, endValue, t, transition.params.noise_magnitude ?? 0);
            } else if (transition.params.kind === 'vec3') {
                next = interpolateVec3(startValue, endValue, t, transition.params.noise_magnitude ?? 0, !!transition.params.spherical);
            } else if (transition.params.kind === 'rotation_matrix') {
                next = interpolateRotation(startValue, endValue, t, transition.params.noise_magnitude ?? 0);
            } else if (transition.params.kind === 'color') {
                const color = paletteFn!(t);
                next = Color.toHexStyle(color);
            }

            assign(target, transition.params.property, next, offset, transformTarget);
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

type TransformTarget = 'rotation' | 'translation' | null | undefined

function getTransformTarget(value: any, target: TransformTarget) {
    if (!value || !target) return value;

    switch (target) {
        case 'rotation':
            return Mat3.fromMat4(Mat3(), value);
        case 'translation':
            return Mat4.getTranslation(Vec3(), value);
        default:
            return value;
    }
}

function select(params: any, path: string | (string | number)[], offset: number, transformTarget: TransformTarget) {
    if (typeof path === 'string') {
        return getTransformTarget(params?.[path], transformTarget);
    }

    let f = params;
    for (let i = offset; i < path.length; i++) {
        if (!f) break;
        f = f[path[i]];
    }

    return getTransformTarget(f, transformTarget);
}

function assignTransformTarget(target: any, xform: TransformTarget, value: any) {
    if (!xform) return;

    switch (xform) {
        case 'rotation':
            Mat4.fromMat3(target, value);
        case 'translation':
            Mat4.setTranslation(target, value);
    }
}

function assign(params: any, path: string | (string | number)[], value: any, offset: number, transformTarget: TransformTarget) {
    if (!params) return;

    if (typeof path === 'string') {
        if (transformTarget) assignTransformTarget(params[path], transformTarget, value);
        else params[path] = value;
        return;
    }

    let f = params;
    for (let i = offset; i < path.length; i++) {
        if (!f) break;
        if (i === path.length - 1) {
            if (transformTarget) assignTransformTarget(f[path[i]], transformTarget, value);
            else f[path[i]] = value;
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

function makePaletteFunction(props: MVSAnimationNode<'interpolate'>): ((value: number) => Color) | undefined {
    if (props.params.kind !== 'color') return undefined;

    const params = palettePropsFromMVSPalette(props.params.palette);
    if (params.name === 'discrete') return makePaletteFunctionDiscrete(params.params);
    if (params.name === 'continuous') return makePaletteFunctionContinuous(params.params);
    throw new Error(`NotImplementedError: makePaletteFunction for ${(props as any).name}`);
}


function makePaletteFunctionDiscrete(props: MVSDiscretePaletteProps): (value: number) => Color {
    const defaultColor = Color(0x0);
    if (props.colors.length === 0) return () => defaultColor;

    return (value: number) => {
        const x = clamp(value, 0, 1);
        for (let i = props.colors.length - 1; i >= 0; i--) {
            const { color, fromValue, toValue } = props.colors[i];
            if (fromValue <= x && x <= toValue) return color;
        }
        return defaultColor;
    };
}

function makePaletteFunctionContinuous(props: MVSContinuousPaletteProps): (value: number) => Color {
    const defaultColor = Color(0x0);
    const { colors, checkpoints } = makeContinuousPaletteCheckpoints(props);
    if (colors.length === 0) return () => defaultColor;

    const underflowColor = props.setUnderflowColor ? props.underflowColor : defaultColor;
    const overflowColor = props.setOverflowColor ? props.overflowColor : defaultColor;

    return (value: number) => {
        const x = clamp(value, 0, 1);
        const gteIdx = SortedArray.findPredecessorIndex(checkpoints, x); // Index of the first greater or equal checkpoint
        if (gteIdx === 0) {
            if (x === checkpoints[0]) return colors[0];
            else return underflowColor;
        }
        if (gteIdx === checkpoints.length) {
            return overflowColor;
        }
        const q = (x - checkpoints[gteIdx - 1]) / (checkpoints[gteIdx] - checkpoints[gteIdx - 1]);
        return Color.interpolate(colors[gteIdx - 1], colors[gteIdx], q);
    };
}