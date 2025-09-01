/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { SortedArray } from '../../../mol-data/int';
import * as EasingFns from '../../../mol-math/easing';
import { clamp, lerp } from '../../../mol-math/interpolate';
import { EPSILON, Mat3, Mat4, Quat, Vec3 } from '../../../mol-math/linear-algebra';
import { RuntimeContext } from '../../../mol-task';
import { deepEqual } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { decodeColor } from '../../../mol-util/color/utils';
import { produce } from '../../../mol-util/produce';
import { makeContinuousPaletteCheckpoints, MVSContinuousPaletteProps, MVSDiscretePaletteProps } from '../components/annotation-color-theme';
import { palettePropsFromMVSPalette } from '../load-helpers';
import { Snapshot } from '../mvs-data';
import { MVSAnimationEasing, MVSAnimationNode, MVSAnimationSchema } from '../tree/animation/animation-tree';
import { Tree } from '../tree/generic/tree-schema';
import { addDefaults } from '../tree/generic/tree-utils';
import { MVSTree } from '../tree/mvs/mvs-tree';
import { ColorT } from '../tree/mvs/param-types';

export async function generateStateTransition(ctx: RuntimeContext, snapshot: Snapshot, snapshotIndex: number, snapshotCount: number) {
    if (!snapshot.animation) return undefined;

    const tree = addDefaults(snapshot.animation, MVSAnimationSchema);
    const transitions = tree.children?.filter(child => child.kind === 'interpolate');
    if (!transitions?.length) return undefined;

    const duration = Math.max(
        snapshot.animation.params?.duration_ms ?? 0,
        ...transitions.map(t => (t.params.start_ms ?? 0) + t.params.duration_ms)
    );

    const frames: [tree: MVSTree, time: number][] = [];
    const dt = tree.params?.frame_time_ms ?? (1000 / 60);
    const N = Math.ceil(duration / dt);

    const nodeMap = makeNodeMap(snapshot.root, new Map(), []);
    const cache = new Map<any, InterpolationCacheEntry>();

    const transitionGroups = groupTranstions(transitions);

    let prevRoot: MVSTree | undefined;
    for (let i = 0; i <= N; i++) {
        const t = i * dt;
        const root = createSnapshot(snapshot.root, transitionGroups, t, cache, nodeMap);

        if (root === prevRoot || (prevRoot && deepEqual(root, prevRoot))) {
            frames[frames.length - 1][1] += dt;
        } else {
            frames.push([root, dt]);
        }

        prevRoot = root;

        if (ctx.shouldUpdate) {
            await ctx.update({ message: `Generating transition for snapshot ${snapshotIndex + 1}/${snapshotCount}`, current: i + 1, max: N });
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

interface InterpolationCacheEntry {
    paletteFn?: (value: number) => Color,
    startColor?: Color | Record<number | string, Color>,
    endColor?: Color | Record<number | string, Color>,
    rotation?: { axis: Vec3, angle: number, start: Quat, end: Quat },
}

function getTransitionKey(transition: MVSAnimationNode<'interpolate'>) {
    const prop = transition.params.property;
    if (Array.isArray(prop)) {
        return `${transition.params.target_ref}:${prop.join('.')}`;
    }
    return `${transition.params.target_ref}:${prop}`;
}

function groupTranstions(transitions: MVSAnimationNode<'interpolate'>[]) {
    const map = new Map<string, MVSAnimationNode<'interpolate'>[]>();
    const groups: MVSAnimationNode<'interpolate'>[][] = [];
    for (const t of transitions) {
        const key = getTransitionKey(t);
        if (!map.has(key)) {
            const group: MVSAnimationNode<'interpolate'>[] = [];
            map.set(key, group);
            groups.push(group);
        }
        map.get(key)!.push(t);
    }
    for (const group of groups) {
        group.sort((a, b) => {
            const s = (a.params.start_ms ?? 0) - (b.params.start_ms ?? 0);
            if (s !== 0) return s;
            return a.params.duration_ms - b.params.duration_ms;
        });
    }
    return groups;
}

function createSnapshot(tree: MVSTree, transitionGroups: MVSAnimationNode<'interpolate'>[][], time: number, cache: Map<any, InterpolationCacheEntry>, nodeMap: Map<string, (string | number)[]>) {
    let modified = false;
    const ret = produce(tree, (draft) => {
        for (const transitionGroup of transitionGroups) {

            const pivot = transitionGroup[0];
            const nodePath = nodeMap.get(pivot.params.target_ref);
            if (!nodePath) continue;

            const node = select(draft, nodePath, 0);
            const target = pivot.params.property[0] === 'custom' ? node?.custom : node?.params;
            if (!target) continue;


            const offset = pivot.params.property[0] === 'custom' ? 1 : 0;

            let transition: MVSAnimationNode<'interpolate'> = pivot;
            let previous: MVSAnimationNode<'interpolate'> | undefined;

            for (let i = transitionGroup.length - 1; i > 0; i--) {
                const current = transitionGroup[i];
                const currentStart = current.params.start_ms ?? 0;
                if (time >= currentStart) {
                    transition = current;
                    previous = i > 0 ? transitionGroup[i - 1] : undefined;
                    break;
                }
            }

            if (!cache.has(transition)) {
                cache.set(transition, {});
            }

            const cacheEntry: InterpolationCacheEntry = cache.get(transition)!;

            const startTime: number = transition.params.start_ms ?? 0;
            const durationMs: number = transition.params.duration_ms ?? 0;
            const t = (time - startTime) / durationMs;

            let next: any;
            if (transition.params.kind === 'transform_matrix') {
                next = processTransformMatrix(transition, target, clamp(t, 0, 1), cacheEntry, offset, previous);
            } else {
                next = processScalarLike(transition, target, t, cacheEntry, offset, previous);
            }

            if (next === undefined) {
                continue;
            }

            modified = true;
            assign(target, transition.params.property, next, offset);
        }
    });
    return modified ? ret : tree;
}

function applyFrequency(t: number, frequency: number, alternate: boolean) {
    let v = (t * (frequency || 1));
    if (v < 1) return v;

    if (!alternate) {
        v = (v % 1);
        if (v === 0) return 1;
        return v;
    }

    if (Math.abs(v - 1) < EPSILON) return 1;
    v = v % 2;
    if (v > 1) return 2 - v;
    return v;
}

function getPreviousScalarEnd(previous: MVSAnimationNode<'interpolate'> | undefined) {
    if (!previous || previous.params.kind === 'transform_matrix') return undefined;
    return previous.params.end;
}

function processScalarLike(transition: MVSAnimationNode<'interpolate'>, target: any, time: number, cacheEntry: InterpolationCacheEntry, offset: number, previous: MVSAnimationNode<'interpolate'> | undefined) {
    if (transition.params.kind === 'transform_matrix') return;
    if (previous && previous.params.kind === 'transform_matrix') return;

    const startValue = transition.params.start ?? getPreviousScalarEnd(previous) ?? select(target, transition.params.property, offset);
    if (transition.params.kind === 'color' && !cacheEntry.paletteFn) {
        cacheEntry.paletteFn = makePaletteFunction(transition);
    }

    const endValue: any = transition.params.end;

    if (time <= 0) return startValue;
    else if (time >= 1 - EPSILON && !transition.params.alternate_direction && transition.params.kind !== 'color') return endValue;

    let t = clamp(time, 0, 1);
    t = applyFrequency(t, transition.params.frequency ?? 1, !!transition.params.alternate_direction);

    const easing = EasingFnMap[transition.params.easing ?? 'linear'] ?? EasingFnMap['linear'];
    t = easing(t);

    if (transition.params.kind === 'scalar') {
        return interpolateScalars(startValue, endValue, t, transition.params.noise_magnitude ?? 0, !!transition.params.discrete);
    } else if (transition.params.kind === 'vec3') {
        return interpolateVectors(startValue, endValue, t, transition.params.noise_magnitude ?? 0, !!transition.params.spherical);
    } else if (transition.params.kind === 'rotation_matrix') {
        return interpolateRotation(startValue, endValue, t, transition.params.noise_magnitude ?? 0, cacheEntry);
    } else if (transition.params.kind === 'color') {
        if (cacheEntry.paletteFn) {
            const color = cacheEntry.paletteFn(t);
            return Color.toHexStyle(color);
        }

        const baseColors = typeof startValue === 'object' ? select(target, transition.params.property, offset) : undefined;
        return interpolateColors(startValue, endValue, t, cacheEntry, baseColors);
    }
}

function getPreviousMatrixEnd(previous: MVSAnimationNode<'interpolate'> | undefined, prop: 'rotation_start' | 'translation_start' | 'scale_start') {
    if (!previous || previous.params.kind !== 'transform_matrix') return undefined;
    return previous.params[prop];
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
function processTransformMatrix(transition: MVSAnimationNode<'interpolate'>, target: any, time: number, cache: InterpolationCacheEntry, offset: number, previous: MVSAnimationNode<'interpolate'> | undefined) {
    if (transition.params.kind !== 'transform_matrix') return;
    if (previous && previous.params.kind !== 'transform_matrix') return;

    const transform = select(target, transition.params.property, offset) ?? Mat4.identity();

    const startRotation = transition.params.rotation_start ?? getPreviousMatrixEnd(previous, 'rotation_start') ?? Mat3.fromMat4(Mat3(), transform);
    const startTranslation = transition.params.translation_start ?? getPreviousMatrixEnd(previous, 'translation_start') ?? Mat4.getTranslation(Vec3(), transform);
    const startScale = transition.params.scale_start ?? getPreviousMatrixEnd(previous, 'scale_start') ?? Mat4.getScaling(Vec3(), transform);

    const endRotation = transition.params.rotation_end;
    const endTranslation = transition.params.translation_end;
    const endScale = transition.params.scale_end;

    let rotation, translation, scale;

    if (time <= 0) {
        rotation = startRotation as Mat3;
        translation = startTranslation as Vec3;
        scale = startScale as Vec3;
    } else {
        const clampedTime = clamp(time, 0, 1);

        let t = applyFrequency(clampedTime, transition.params.rotation_frequency ?? 1, !!transition.params.rotation_alternate_direction);
        let easing = EasingFnMap[transition.params.rotation_easing ?? 'linear'] ?? EasingFnMap['linear'];
        rotation = interpolateRotation(startRotation as Mat3, endRotation as Mat3, easing(t), transition.params.rotation_noise_magnitude ?? 0, cache);

        t = applyFrequency(clampedTime, transition.params.translation_frequency ?? 1, !!transition.params.translation_alternate_direction);
        easing = EasingFnMap[transition.params.translation_easing ?? 'linear'] ?? EasingFnMap['linear'];
        translation = interpolateVec3(startTranslation as Vec3, endTranslation as Vec3 | undefined, easing(t), transition.params.translation_noise_magnitude ?? 0, false);

        t = applyFrequency(clampedTime, transition.params.scale_frequency ?? 1, !!transition.params.scale_alternate_direction);
        easing = EasingFnMap[transition.params.scale_easing ?? 'linear'] ?? EasingFnMap['linear'];
        scale = interpolateVec3(startScale as Vec3, endScale as Vec3 | undefined, easing(t), transition.params.scale_noise_magnitude ?? 0, false);
    }

    const pivot = transition.params.pivot ?? Vec3.zero();

    Mat4.fromTranslation(TransformState.translation, translation);
    Mat4.fromScaling(TransformState.scale, scale);
    Mat4.setIdentity(TransformState.rotation);
    Mat4.fromMat3(TransformState.rotation, rotation);
    Mat4.fromTranslation(TransformState.pivotTranslation, pivot as Vec3);
    Mat4.fromTranslation(TransformState.pivotTranslationInv, Vec3.negate(TransformState.pivotNeg, pivot as Vec3));

    // translation . pivot . rotation . scale . pivotInv
    const result = Mat4();
    Mat4.mul(result, TransformState.scale, TransformState.pivotTranslationInv);
    Mat4.mul(result, TransformState.rotation, result);
    Mat4.mul(result, TransformState.translation, result);

    return result;
}

function interpolateScalars(start: number | number[], end: number | number[] | undefined, t: number, noise: number, discrete: boolean) {
    if (Array.isArray(start)) {
        const ret = Array.from<number>({ length: start.length }).fill(0.1);
        if (!end || !Array.isArray(end)) {
            for (let i = 0; i < start.length; i++) {
                ret[i] = interpolateScalar(start[i], end, t, noise, discrete);
            }
            return ret;
        }

        for (let i = 0; i < start.length; i++) {
            ret[i] = interpolateScalar(start[i], end[i], t, noise, discrete);
        }
        return ret;
    }

    if (Array.isArray(end)) {
        const ret = Array.from<number>({ length: end.length }).fill(0.1);
        for (let i = 0; i < end.length; i++) {
            ret[i] = interpolateScalar(start, end[i], t, noise, discrete);
        }
        return ret;
    }

    return interpolateScalar(start, end, t, noise, discrete);
}

function interpolateScalar(start: number, end: number | undefined, t: number, noise: number, discrete: boolean) {
    let v = typeof end === 'number' ? lerp(start, end, t) : start;
    if (noise) {
        v += (Math.random() - 0.5) * noise;
    }
    if (discrete) {
        v = Math.round(v);
    }
    return v;
}

const InterpolateVectorsState = {
    start: Vec3(),
    end: Vec3(),
    v: Vec3(),
};
function interpolateVectors(start: number[], end: number[] | undefined, t: number, noise: number, isSpherical: boolean) {
    if ((!end || start === end) && !noise) return start;

    const ret: number[] = Array.from<number>({ length: start.length }).fill(0.1);

    for (let i = 0; i < start.length; i += 3) {
        const s = Vec3.fromArray(InterpolateVectorsState.start, start, i);

        let v: Vec3;
        if (end) {
            const e = Vec3.fromArray(InterpolateVectorsState.end, end, i);
            v = isSpherical
                ? Vec3.slerp(InterpolateVectorsState.v, s, e, t)
                : Vec3.lerp(InterpolateVectorsState.v, s, e, t);
        } else {
            v = Vec3.clone(s);
        }

        if (noise && t <= 1 - EPSILON) {
            Vec3.random(Vec3Noise, noise);
            Vec3.add(v, v, Vec3Noise);
        }

        Vec3.toArray(v, ret, i);
    }

    return ret;
}

const Vec3Noise = Vec3();
function interpolateVec3(start: Vec3, end: Vec3 | undefined, t: number, noise: number, isSpherical: boolean) {
    if ((!end || start === end) && !noise) return start;

    let v: Vec3;

    if (end) {
        v = isSpherical
            ? Vec3.slerp(Vec3(), start, end, t)
            : Vec3.lerp(Vec3(), start, end, t);
    } else {
        v = Vec3.clone(start);
    }

    if (noise && t <= 1 - EPSILON) {
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
function interpolateRotation(start: Mat3, end: Mat3 | undefined, t: number, noise: number, cache: InterpolationCacheEntry) {
    if ((!end || start === end) && !noise) return start;

    if (end) {
        if (!cache.rotation) {
            cache.rotation = {
                ...relativeAxisAngle(start, end),
                start: Quat.fromMat3(Quat(), start),
                end: Quat.fromMat3(Quat(), end),
            };
        }

        const { axis, angle, start: startQ, end: endQ } = cache.rotation;

        if (angle < 1e-6) {
            // start ≈ end: make a clean spin about the detected (or default) axis
            Quat.setAxisAngle(RotationState.v, axis, t * 2 * Math.PI); // Make a full turn
        } else {
            // Normal case: stick with your existing slerp between start/end
            Quat.slerp(RotationState.v, startQ, endQ, t);
        }
    } else {
        Quat.fromMat3(RotationState.v, start);
    }

    if (noise && t <= 1 - EPSILON) {
        Vec3.random(RotationState.axis, 1);
        Quat.setAxisAngle(RotationState.noise, RotationState.axis, 2 * Math.PI * noise * (Math.random() - 0.5));
        Quat.multiply(RotationState.v, RotationState.noise, RotationState.v);
    }
    Mat4.fromQuat(RotationState.temp, RotationState.v);
    return Mat3.fromMat4(Mat3(), RotationState.temp);
}

function decodeColors(color: ColorT | Record<number | string, ColorT> | undefined, baseColors: Record<number | string, ColorT> | undefined) {
    if (color === undefined || color === null) return undefined;

    if (typeof color === 'object') {
        const ret: Record<number | string, Color> = {};
        if (baseColors) {
            for (const key of Object.keys(baseColors)) {
                const decoded = decodeColor(baseColors[key]);
                if (decoded !== undefined) {
                    ret[key] = decoded;
                }
            }
        }
        for (const key of Object.keys(color)) {
            const decoded = decodeColor(color[key]);
            if (decoded !== undefined) {
                ret[key] = decoded;
            }
        }
        return ret;
    }

    return decodeColor(color);
}

function interpolateColors(start: ColorT | Record<number, ColorT>, end: ColorT | Record<number, ColorT> | undefined, time: number, cacheEntry: InterpolationCacheEntry, baseColors: Record<number, ColorT> | undefined) {
    const t = clamp(time, 0, 1);

    if (cacheEntry.paletteFn) {
        const c = cacheEntry.paletteFn(t);
        return Color.toHexStyle(c);
    }

    if (cacheEntry.startColor === undefined) {
        cacheEntry.startColor = decodeColors(start, baseColors);
    }
    if (cacheEntry.endColor === undefined) {
        cacheEntry.endColor = decodeColors(end, undefined);
    }

    const { startColor, endColor } = cacheEntry;

    if (typeof startColor === 'object') {
        if (typeof baseColors !== 'object') {
            throw new Error('Cannot interpolate from scalar color to color mapping');
        }

        const ret = { ...baseColors as any, ...startColor as any };
        if (typeof endColor === 'object') {
            for (const key of Object.keys(endColor)) {
                ret[key] = Color.toHexStyle(Color.interpolate(startColor[key], endColor[key], t));
            }
        } else if (typeof endColor === 'number') {
            for (const key of Object.keys(startColor)) {
                ret[key] = Color.toHexStyle(Color.interpolate(startColor[key], endColor, t));
            }
        }
        return ret;
    }
    if (typeof endColor === 'object') {
        throw new Error('Cannot interpolate from scalar color to color mapping');
    }

    if (typeof endColor === 'number' && typeof startColor === 'number') {
        return Color.toHexStyle(Color.interpolate(startColor, endColor, t));
    }

    return start;
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

function makeNodeMap(tree: Tree, map: Map<string, (string | number)[]>, currentPath: (string | number)[]) {
    if (tree.ref) {
        map.set(tree.ref, [...currentPath]);
    }

    if (!tree.children) return map;

    currentPath.push('children');
    for (let i = 0; i < tree.children.length; i++) {
        const child = tree.children[i];
        currentPath.push(i);
        makeNodeMap(child, map, currentPath);
        currentPath.pop();
    }
    currentPath.pop();

    return map;
}

function makePaletteFunction(props: MVSAnimationNode<'interpolate'>): ((value: number) => Color) | undefined {
    if (props.params.kind !== 'color' || !props.params.palette) return undefined;

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

const RelativeAxisAngleState = {
    Rt: Mat3(),
    R: Mat3(),
};
function relativeAxisAngle(start: Mat3, end: Mat3): { axis: Vec3, angle: number } {
    // R_rel = end * start^T
    const R0 = start, R1 = end;
    const Rt = Mat3.transpose(RelativeAxisAngleState.Rt, R0);
    const R = Mat3.mul(RelativeAxisAngleState.R, R1, Rt);

    const tr = R[0] + R[4] + R[8]; // trace
    let angle = Math.acos(clamp((tr - 1) * 0.5, -1, 1)); // in [0, π]
    const axis = Vec3();

    const eps = 1e-6;
    const sinA = Math.sin(angle);

    if (angle < eps) {
        // Near identity: axis undefined; return any unit axis (choose something stable)
        Vec3.set(axis, 0, 0, 1);
        angle = 0.0;
        return { axis, angle };
    }

    if (Math.PI - angle > 1e-4) {
        // General case
        axis[0] = (R[5] - R[7]) / (2 * sinA); // (r32 - r23)
        axis[1] = (R[6] - R[2]) / (2 * sinA); // (r13 - r31)
        axis[2] = (R[1] - R[3]) / (2 * sinA); // (r21 - r12)
        Vec3.normalize(axis, axis);
        return { axis, angle };
    }

    // angle ~ π: use diagonal-based extraction for stability
    // Compute squared components then pick the largest to avoid precision loss
    const xx = Math.max(0, (R[0] + 1) * 0.5);
    const yy = Math.max(0, (R[4] + 1) * 0.5);
    const zz = Math.max(0, (R[8] + 1) * 0.5);

    let x = Math.sqrt(xx), y = Math.sqrt(yy), z = Math.sqrt(zz);

    if (x >= y && x >= z) {
        x = Math.max(x, 1e-8);
        y = (R[1] + R[3]) / (4 * x);
        z = (R[2] + R[6]) / (4 * x);
        Vec3.set(axis, x, y, z);
    } else if (y >= x && y >= z) {
        y = Math.max(y, 1e-8);
        x = (R[1] + R[3]) / (4 * y);
        z = (R[5] + R[7]) / (4 * y);
        Vec3.set(axis, x, y, z);
    } else {
        z = Math.max(z, 1e-8);
        x = (R[2] + R[6]) / (4 * z);
        y = (R[5] + R[7]) / (4 * z);
        Vec3.set(axis, x, y, z);
    }

    Vec3.normalize(axis, axis);
    return { axis, angle: Math.PI };
}