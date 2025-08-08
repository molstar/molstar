import { produce } from 'immer';
import { Snapshot } from '../mvs-data';
import { Tree } from '../tree/generic/tree-schema';
import { MVSNode } from '../tree/mvs/mvs-tree';
import { clamp, lerp } from '../../../mol-math/interpolate';

export function generateTransitions(snapshot: Snapshot) {
    const transitions = snapshot.root.children?.filter(child => child.kind === 'transition');
    if (!transitions?.length) return undefined;

    const duration = Math.max(...transitions.map(t => t.params.duration_ms));

    const snapshots: Snapshot[] = [];

    const dt = 1000 / 60;
    for (let t = 0; t < duration + dt; t += dt) {
        const root = createSnapshot(snapshot.root, transitions, t);
        if (!root) break;
        snapshots.push({ metadata: snapshot.metadata, root: root as any });
    }

    return { duration, snapshots };
}

function createSnapshot(tree: Tree, transitions: MVSNode<'transition'>[], time: number) {
    return produce(tree, (draft) => {
        for (const transition of transitions) {
            const node = findNode(draft, transition.params.target_ref);
            const value = select(node?.params, transition.params.property);
            const t = clamp(time / transition.params.duration_ms, 0, 1);
            const next = lerp(value, transition.params.target_value, t);
            assign(node?.params, transition.params.property, next);
        }
    });
}

function select(params: any, path: (string | number)[]) {
    if (!params) return;
    let f = params[path[0]];
    for (let i = 1; i < path.length; i++) {
        if (!f) break;
        f = f[path[i]];
    }
    return f;
}

function assign(params: any, path: (string | number)[], value: any) {
    if (!params) return;
    let f = params[path[0]];
    for (let i = 1; i < path.length; i++) {
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