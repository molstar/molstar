/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { canonicalJsonString } from '../../../../mol-util/json';
import { addParamDefaults } from './params-schema';
import { CustomProps, Kind, Node, Subtree, SubtreeOfKind, Tree, TreeFor, TreeSchema, TreeSchemaWithAllRequired, getParams } from './tree-schema';


/** Run DFS (depth-first search) algorithm on a rooted tree.
 * Runs `visit` function when a node is discovered (before visiting any descendants).
 * Runs `postVisit` function when leaving a node (after all descendants have been visited). */
export function dfs<TTree extends Tree>(root: TTree, visit?: (node: Subtree<TTree>, parent?: Subtree<TTree>) => any, postVisit?: (node: Subtree<TTree>, parent?: Subtree<TTree>) => any) {
    return _dfs<Subtree<TTree>>(root, undefined, visit, postVisit);
}
function _dfs<TTree extends Tree>(root: TTree, parent: Subtree<TTree> | undefined, visit?: (node: Subtree<TTree>, parent?: Subtree<TTree>) => any, postVisit?: (node: Subtree<TTree>, parent?: Subtree<TTree>) => any) {
    if (visit) visit(root, parent);
    for (const child of root.children ?? []) {
        _dfs<Subtree<TTree>>(child, root, visit, postVisit);
    }
    if (postVisit) postVisit(root, parent);
}

/** Convert a tree into a pretty-printed string. */
export function treeToString(tree: Tree) {
    let level = 0;
    const lines: string[] = [];
    dfs(tree, node => lines.push('  '.repeat(level++) + nodeToString(node)), node => level--);
    return lines.join('\n');
}
function nodeToString(node: Node) {
    return `- ${node.kind} ${formatObject(node.params ?? {})}${formatCustomProps(node.custom)}${formatRef(node.ref)}`;
}

/** Convert object to a human-friendly string (similar to JSON.stringify but without quoting keys) */
export function formatObject(obj: {} | undefined): string {
    if (!obj) return 'undefined';
    return JSON.stringify(obj).replace(/,("\w+":)/g, ', $1').replace(/"(\w+)":/g, '$1: ');
}

/** Return human-friendly string with node custom properties, if any */
function formatCustomProps(customProps: CustomProps | undefined): string {
    if (!customProps || Object.keys(customProps).length === 0) return '';
    return `, custom: ${formatObject(customProps)}`;
}

/** Return human-friendly string with node ref, if any */
function formatRef(ref: string | undefined): string {
    if (ref === undefined) return '';
    return `, ref: "${ref}"`;
}


/** Create a copy of a tree node, ignoring children. */
export function copyNodeWithoutChildren<TTree extends Tree>(node: TTree): TTree {
    return {
        kind: node.kind,
        params: node.params ? { ...node.params } : undefined,
        custom: node.custom ? { ...node.custom } : undefined,
        ref: node.ref,
    } as TTree;
}
/** Create a copy of a tree node, including a shallow copy of children. */
export function copyNode<TTree extends Tree>(node: TTree): TTree {
    return {
        kind: node.kind,
        params: node.params ? { ...node.params } : undefined,
        custom: node.custom ? { ...node.custom } : undefined,
        ref: node.ref,
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
    [kind in Kind<Subtree<A>>]?: (node: SubtreeOfKind<A, kind>, parent?: Subtree<A>) => { subtree: Subtree<B>[] }
};

/** Apply a set of conversion rules to a tree to change to a different schema. */
export function convertTree<A extends Tree, B extends Tree>(root: A, conversions: ConversionRules<A, B>): Subtree<B> {
    const mapping = new Map<Subtree<A>, Subtree<B>>();
    let convertedRoot: Subtree<B>;
    dfs<A>(root, (node, parent) => {
        const conversion = conversions[node.kind as (typeof node)['kind']] as ((n: typeof node, p?: Subtree<A>) => { subtree: Subtree<B>[] }) | undefined;
        if (conversion) {
            const converted = conversion(node, parent);
            if (!parent && converted?.subtree.length === 0) throw new Error('Cannot convert root to empty path');
            let convParent = parent ? mapping.get(parent) : undefined;
            for (const conv of converted.subtree) {
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
    const map = new Map<string, Subtree<T>>();
    const result = copyTree(root);
    dfs<T>(result, node => {
        map.clear();
        const newChildren: Subtree<T>[] = [];
        for (const child of node.children ?? []) {
            let twin: Subtree<T> | undefined = undefined;
            const doApply = (!condenseNodes || condenseNodes.has(child.kind)) && !skipNodes?.has(child.kind);
            if (doApply) {
                const key = child.kind + canonicalJsonString(getParams(child));
                twin = map.get(key);
                if (!twin) map.set(key, child);
            }
            if (twin) {
                (twin.children ??= []).push(...child.children ?? []);
            } else {
                newChildren.push(child as Subtree<T>);
            }
        }
        node.children = newChildren;
    });
    return result;
}

/** Create a copy of the tree where missing optional params for each node are added based on `defaults`. */
export function addDefaults<S extends TreeSchema>(tree: TreeFor<S>, treeSchema: S): TreeFor<TreeSchemaWithAllRequired<S>> {
    type TTree = TreeFor<S>;
    const rules: ConversionRules<TTree, TTree> = {};
    for (const kind in treeSchema.nodes) {
        rules[kind as Kind<Subtree<TTree>>] = node => ({
            subtree: [{
                kind: node.kind,
                params: addParamDefaults(treeSchema.nodes[kind].params, node.params as any),
                custom: node.custom,
                ref: node.ref,
            } as Node as any]
        });
    }
    return convertTree(tree, rules) as any;
}

/** Resolve any URI params in a tree, in place. URI params are those listed in `uriParamNames`.
 * Relative URIs are treated as relative to `baseUri`, which can in turn be relative to the window URL (if available). */
export function resolveUris<T extends Tree>(tree: T, baseUri: string, uriParamNames: string[]): void {
    dfs(tree, node => {
        const params = node.params as Record<string, any> | undefined;
        if (!params) return;
        for (const name in params) {
            if (uriParamNames.includes(name)) {
                const uri = params[name];
                if (typeof uri === 'string') {
                    params[name] = resolveUri(uri, baseUri, windowUrl());
                }
            }
        }
    });
}

/** Resolve a sequence of URI references (relative URIs), where each reference is either absolute or relative to the next one
 * (i.e. the last one is the base URI). Skip any `undefined`.
 * E.g. `resolveUri('./unexpected.png', '/spanish/inquisition/expectations.html', 'https://example.org/spam/spam/spam')`
 * returns `'https://example.org/spanish/inquisition/unexpected.png'`. */
function resolveUri(...refs: (string | undefined)[]): string | undefined {
    let result: string | undefined = undefined;
    for (const ref of refs.reverse()) {
        if (ref !== undefined) {
            if (result === undefined) result = ref;
            else result = new URL(ref, result).href;
        }
    }
    return result;
}

/** Return URL of the current page when running in a browser; `undefined` when running in Node. */
function windowUrl(): string | undefined {
    return (typeof window !== 'undefined') ? window.location.href : undefined;
}
