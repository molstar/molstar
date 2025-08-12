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

            const startTime = transition.params.start_ms ?? 0;
            const endTime = startTime + transition.params.duration_ms;

            const deltaT = endTime - startTime;
            let t = clamp((time - startTime) / deltaT, 0, 1);

            if (transition.params.kind === 'transform_matrix') {
                processTransformMatrix(transition, target, t);
                continue;
            }

            if (transition.params.kind === 'color') {
                if (!cache.has(transition)) {
                    cache.set(transition, {
                        paletteFn: makePaletteFunction(transition)
                    });
                }
            }

            const paletteFn = cache.get(transition)?.paletteFn;

            const offset = transition.params.property[0] === 'custom' ? 1 : 0;
            const startValue: any = transition.params.kind === 'color'
                ? Color.toHexStyle(paletteFn!(0))
                : transition.params.start ?? select(target, transition.params.property, offset);
            const endValue: any = transition.params.kind === 'color'
                ? Color.toHexStyle(paletteFn!(1))
                : transition.params.end;

            if (time <= startTime) {
                assign(target, transition.params.property, startValue, offset);
                continue;
            }

            if (time >= endTime - EPSILON) {
                assign(target, transition.params.property, endValue, offset);
                continue;
            }

            const easing = EasingFnMap[transition.params.easing ?? 'linear'] ?? EasingFnMap['linear'];
            t = easing(t);

            let next: any;
            if (transition.params.kind === 'scalar') {
                next = interpolateScalar(startValue, endValue, t, transition.params.noise_magnitude ?? 0);
            } else if (transition.params.kind === 'vec3') {
                next = interpolateVectors(startValue, endValue, t, transition.params.noise_magnitude ?? 0, !!transition.params.spherical);
            } else if (transition.params.kind === 'rotation_matrix') {
                next = interpolateRotation(startValue, endValue, t, transition.params.noise_magnitude ?? 0);
            } else if (transition.params.kind === 'color') {
                const color = paletteFn!(t);
                next = Color.toHexStyle(color);
            }

            assign(target, transition.params.property, next, offset);
        }
    });
}

const TransformState = {
    pivotTranslation: Mat4(),
    pivotTranslationInv: Mat4(),
    rotation: Mat4(),
    scale: Mat4(),
    translation: Mat4(),
    pivotNeg: Vec3(),
    temp: Mat4(),
};
function processTransformMatrix(transition: MVSAnimationNode<'interpolate'>, target: any, t: number) {
    if (transition.params.kind !== 'transform_matrix') return;

    const offset = transition.params.property[0] === 'custom' ? 1 : 0;
    const transform = select(target, transition.params.property, offset) ?? Mat4.identity();

    const startRotation = transition.params.rotation_start ?? Mat3.fromMat4(Mat3(), transform);
    const startTranslation = transition.params.translation_start ?? Mat4.getTranslation(Vec3(), transform);
    const startScale = transition.params.scale_start ?? Mat4.getScaling(Vec3(), transform);

    const endRotation = transition.params.rotation_end ?? startRotation;
    const endTranslation = transition.params.translation_end ?? startTranslation;
    const endScale = transition.params.scale_end ?? startScale;

    let easing = EasingFnMap[transition.params.rotation_easing ?? 'linear'] ?? EasingFnMap['linear'];
    const rotation = interpolateRotation(startRotation as Mat3, endRotation as Mat3, easing(t), transition.params.rotation_noise_magnitude ?? 0);
    easing = EasingFnMap[transition.params.translation_easing ?? 'linear'] ?? EasingFnMap['linear'];
    const translation = interpolateVec3(startTranslation as Vec3, endTranslation as Vec3, easing(t), transition.params.translation_noise_magnitude ?? 0, false);
    easing = EasingFnMap[transition.params.scale_easing ?? 'linear'] ?? EasingFnMap['linear'];
    const scale = interpolateVec3(startScale as Vec3, endScale as Vec3, easing(t), transition.params.scale_noise_magnitude ?? 0, false);

    const pivot = transition.params.pivot ?? Vec3.zero();

    Mat4.fromTranslation(TransformState.translation, translation);
    Mat4.fromScaling(TransformState.scale, scale);
    Mat4.setIdentity(TransformState.rotation);
    Mat4.fromMat3(TransformState.rotation, rotation);
    Mat4.fromTranslation(TransformState.pivotTranslation, pivot as Vec3);
    Mat4.fromTranslation(TransformState.pivotTranslationInv, Vec3.negate(TransformState.pivotNeg, pivot as Vec3));

    // translation . pivot . scale . rotation . pivotInv
    const result = Mat4();
    Mat4.mul(result, TransformState.rotation, TransformState.pivotTranslationInv);
    Mat4.mul(result, TransformState.scale, result);
    Mat4.mul(result, TransformState.translation, result);

    assign(target, transition.params.property, result, offset);
}

function interpolateScalar(start: number, end: number, t: number, noise: number) {
    let v = lerp(start, end, t);
    if (noise) {
        v += (Math.random() - 0.5) * noise;
    }
    return v;
}

const InterpolateVectorsState = {
    start: Vec3(),
    end: Vec3(),
    v: Vec3(),
};
function interpolateVectors(start: number[], end: number[], t: number, noise: number, isSpherical: boolean) {
    if (start === end && !noise) return start;

    const ret: number[] = Array.from<number>({ length: start.length }).fill(0.1);

    for (let i = 0; i < start.length; i += 3) {
        const s = Vec3.fromArray(InterpolateVectorsState.start, start, i);
        const e = Vec3.fromArray(InterpolateVectorsState.end, end, i);

        const v = isSpherical
            ? Vec3.slerp(InterpolateVectorsState.v, s, e, t)
            : Vec3.lerp(InterpolateVectorsState.v, s, e, t);

        if (noise) {
            Vec3.random(Vec3Noise, noise);
            Vec3.add(v, v, Vec3Noise);
        }

        Vec3.toArray(v, ret, i);
    }

    return ret;
}

const Vec3Noise = Vec3();
function interpolateVec3(start: Vec3, end: Vec3, t: number, noise: number, isSpherical: boolean) {
    if (start === end && !noise) return start;

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
    if (start === end && !noise) return start;

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
        return Color.interpolateHcl(colors[gteIdx - 1], colors[gteIdx], q);
    };
}