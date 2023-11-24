/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { canonicalJsonString } from '../../../../mol-util/json';
import { DefaultsForTree, Kind, SubTree, SubTreeOfKind, Tree, TreeFor, TreeSchema, TreeSchemaWithAllRequired, getParams } from './tree-schema';


/** Run DFS (depth-first search) algorithm on a rooted tree.
 * Runs `visit` function when a node is discovered (before visiting any descendants).
 * Runs `postVisit` function when leaving a node (after all descendants have been visited). */
export function dfs<TTree extends Tree>(root: TTree, visit?: (node: SubTree<TTree>, parent?: SubTree<TTree>) => any, postVisit?: (node: SubTree<TTree>, parent?: SubTree<TTree>) => any) {
    return _dfs<SubTree<TTree>>(root, undefined, visit, postVisit);
}
function _dfs<TTree extends Tree>(root: TTree, parent: SubTree<TTree> | undefined, visit?: (node: SubTree<TTree>, parent?: SubTree<TTree>) => any, postVisit?: (node: SubTree<TTree>, parent?: SubTree<TTree>) => any) {
    if (visit) visit(root, parent);
    for (const child of root.children ?? []) {
        _dfs<SubTree<TTree>>(child, root, visit, postVisit);
    }
    if (postVisit) postVisit(root, parent);
}

/** Convert a tree into a pretty-printed string. */
export function treeToString(tree: Tree) {
    let level = 0;
    const lines: string[] = [];
    dfs(tree, node => lines.push('  '.repeat(level++) + `- ${node.kind} ${formatObject(node.params ?? {})}`), node => level--);
    return lines.join('\n');
}

/** Convert object to a human-friendly string (similar to JSON.stringify but without quoting keys) */
function formatObject(obj: {} | undefined): string {
    if (!obj) return 'undefined';
    return JSON.stringify(obj).replace(/,("\w+":)/g, ', $1').replace(/"(\w+)":/g, '$1: ');
}


/** Create a copy of a tree node, ignoring children. */
export function copyNodeWithoutChildren<TTree extends Tree>(node: TTree): TTree {
    return {
        kind: node.kind,
        params: node.params ? { ...node.params } : undefined,
    } as TTree;
}
/** Create a copy of a tree node, including a shallow copy of children. */
export function copyNode<TTree extends Tree>(node: TTree): TTree {
    return {
        kind: node.kind,
        params: node.params ? { ...node.params } : undefined,
        children: node.children ? [...node.children] : undefined,
    } as TTree;
}

/** Create a deep copy of a tree. */
export function copyTree<T extends Tree>(root: T): T {
    return convertTree(root, {}) as T;
}

/** Set of rules for converting a tree of one schema into a different schema.
 * Each rule defines how to convert a node of a specific kind, e.g.
 * `{A: node => [], B: node => [{kind: 'X',...}], C: node => [{kind: 'Y',...}, {kind: 'Z',...}]}`:
 * nodes of kind `A` will be deleted (their children moved to parent),
 * nodes of kind `B` will be converted to kind `X`,
 * nodes of kind `C` will be converted to `Y` with a child `Z` (original children moved to `Z`),
 * nodes of other kinds will just be copied. */
export type ConversionRules<A extends Tree, B extends Tree> = {
    [kind in Kind<SubTree<A>>]?: (node: SubTreeOfKind<A, kind>, parent?: SubTree<A>) => SubTree<B>[]
};

/** Apply a set of conversion rules to a tree to change to a different schema. */
export function convertTree<A extends Tree, B extends Tree>(root: A, conversions: ConversionRules<A, B>): SubTree<B> {
    const mapping = new Map<SubTree<A>, SubTree<B>>();
    let convertedRoot: SubTree<B>;
    dfs<A>(root, (node, parent) => {
        const conversion = conversions[node.kind as (typeof node)['kind']] as ((n: typeof node, p?: SubTree<A>) => SubTree<B>[]) | undefined;
        if (conversion) {
            const convertidos = conversion(node, parent);
            if (!parent && convertidos.length === 0) throw new Error('Cannot convert root to empty path');
            let convParent = parent ? mapping.get(parent) : undefined;
            for (const conv of convertidos) {
                if (convParent) {
                    (convParent.children ??= []).push(conv);
                } else {
                    convertedRoot = conv;
                }
                convParent = conv;
            }
            mapping.set(node, convParent!);
        } else {
            const converted = copyNodeWithoutChildren(node);
            if (parent) {
                (mapping.get(parent)!.children ??= []).push(converted);
            } else {
                convertedRoot = converted;
            }
            mapping.set(node, converted);
        }
    });
    return convertedRoot!;
}

/** Create a copy of the tree where twins (siblings of the same kind with the same params) are merged into one node.
 * Applies only to the node kinds listed in `condenseNodes` (or all if undefined) except node kinds in `skipNodes`. */
export function condenseTree<T extends Tree>(root: T, condenseNodes?: Set<Kind<Tree>>, skipNodes?: Set<Kind<Tree>>): T {
    const map = new Map<string, SubTree<T>>();
    const result = copyTree(root);
    dfs<T>(result, node => {
        map.clear();
        const newChildren: SubTree<T>[] = [];
        for (const child of node.children ?? []) {
            let twin: SubTree<T> | undefined = undefined;
            const doApply = (!condenseNodes || condenseNodes.has(child.kind)) && !skipNodes?.has(child.kind);
            if (doApply) {
                const key = child.kind + canonicalJsonString(getParams(child));
                twin = map.get(key);
                if (!twin) map.set(key, child);
            }
            if (twin) {
                (twin.children ??= []).push(...child.children ?? []);
            } else {
                newChildren.push(child as SubTree<T>);
            }
        }
        node.children = newChildren;
    });
    return result;
}

/** Create a copy of the tree where missing optional params for each node are added based on `defaults`. */
export function addDefaults<S extends TreeSchema>(tree: TreeFor<S>, defaults: DefaultsForTree<S>): TreeFor<TreeSchemaWithAllRequired<S>> {
    const rules: ConversionRules<TreeFor<S>, TreeFor<S>> = {};
    for (const kind in defaults) {
        rules[kind] = node => [{ kind: node.kind, params: { ...defaults[kind], ...node.params } } as any];
    }
    return convertTree(tree, rules) as any;
}
