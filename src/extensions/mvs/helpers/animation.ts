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


export async function generateStateAnimation(ctx: RuntimeContext, snapshot: Snapshot) {
    if (!snapshot.animation) return undefined;

    const tree = addDefaults(snapshot.animation, MVSAnimationSchema);
    const transitions = tree.children?.filter(child => child.kind === 'transition');
    if (!transitions?.length) return undefined;

    const duration = Math.min(snapshot.metadata.linger_duration_ms, Math.max(...transitions.map(t => t.params.end_ms)));

    const frames: MVSTree[] = [];
    const dt = tree.params?.frame_time_ms ?? (1000 / 60);
    const N = Math.ceil(duration / dt);

    for (let i = 0; i <= N; i++) {
        const t = i * dt;
        const root = createSnapshot(snapshot.root, transitions, t);
        if (!root) break;
        frames.push(root);

        if (ctx.shouldUpdate) {
            await ctx.update({ message: 'Generating animation...' });
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

function createSnapshot(tree: MVSTree, transitions: MVSAnimationNode<'transition'>[], time: number) {
    return produce(tree, (draft) => {
        for (const transition of transitions) {
            const node = findNode(draft, transition.params.target_ref);
            if (!node) continue;

            const target = transition.params.property[0] === 'custom' ? node?.custom : node?.params;
            if (!target) continue;

            const offset = transition.params.property[0] === 'custom' ? 1 : 0;
            const startTime = transition.params.start_ms ?? 0;
            const startValue = transition.params.start_value ?? select(target, transition.params.property, offset);
            const endTime = transition.params.end_ms;

            if (time <= startTime) {
                assign(target, transition.params.property, startValue, offset);
                continue;
            }

            const endValue = transition.params.end_value;
            if (time >= endTime) {
                assign(target, transition.params.property, endValue, offset);
                continue;
            }

            const deltaT = endTime - startTime;
            const t = clamp((time - startTime) / deltaT, 0, 1);
            const easing = EasingFnMap[transition.params.easing ?? 'linear'] ?? EasingFnMap['linear'];
            const next = lerp(startValue, endValue, easing(t));

            assign(target, transition.params.property, next, offset);
        }
    });
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