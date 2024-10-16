/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObject, StateObjectSelector, StateTransform, StateTransformer } from '../../mol-state';
import { stringHash } from './helpers/utils';
import { Kind, Subtree, SubtreeOfKind, Tree } from './tree/generic/tree-schema';
import { dfs } from './tree/generic/tree-utils';


/** Function responsible for loading a tree node `node` into Mol*.
 * Should apply changes within `updateParent.update` but not commit them.
 * Should modify `context` accordingly, if it is needed for loading other nodes later.
 * `updateParent.selector` is the result of loading the node's parent into Mol* state hierarchy (or the hierarchy root in case of root node). */
export type LoadingAction<TNode extends Tree, TContext> = (updateParent: UpdateTarget, node: TNode, context: TContext) => UpdateTarget | undefined

/** Loading actions for loading a tree into Mol*, per node kind. */
export type LoadingActions<TTree extends Tree, TContext> = { [kind in Kind<Subtree<TTree>>]?: LoadingAction<SubtreeOfKind<TTree, kind>, TContext> }

/** Type for defining custom behavior when loading trees, usually based on node custom properties. */
export interface LoadingExtension<TTree extends Tree, TContext, TExtensionContext> {
    id: string,
    description: string,
    /** Runs before the tree is loaded */
    createExtensionContext: (tree: TTree, context: TContext) => TExtensionContext,
    /** Runs after the tree is loaded */
    disposeExtensionContext?: (extensionContext: TExtensionContext, tree: TTree, context: TContext) => void,
    /** Runs on every node of the tree */
    action: (updateTarget: UpdateTarget, node: Subtree<TTree>, context: TContext, extensionContext: TExtensionContext) => void,
}


/** Load a tree into Mol*, by applying loading actions in DFS order and then commiting at once.
 * If `options.replaceExisting`, remove all objects in the current Mol* state; otherwise add to the current state. */
export async function loadTree<TTree extends Tree, TContext>(
    plugin: PluginContext,
    tree: TTree,
    loadingActions: LoadingActions<TTree, TContext>,
    context: TContext,
    options?: { replaceExisting?: boolean, extensions?: LoadingExtension<TTree, TContext, any>[] }
) {
    const mapping = new Map<Subtree<TTree>, UpdateTarget | undefined>();
    const updateRoot: UpdateTarget = UpdateTarget.create(plugin, options?.replaceExisting ?? false);
    if (options?.replaceExisting) {
        UpdateTarget.deleteChildren(updateRoot);
    }
    const extensionContexts = (options?.extensions ?? []).map(ext => ({ ext, extCtx: ext.createExtensionContext(tree, context) }));
    const mvsRefMap = new Map<string, string>();
    dfs<TTree>(tree, (node, parent) => {
        const kind: Kind<typeof node> = node.kind;
        let msNode: UpdateTarget | undefined;
        const updateParent = parent ? mapping.get(parent) : updateRoot;
        const action = loadingActions[kind] as LoadingAction<typeof node, TContext> | undefined;
        if (action) {
            if (updateParent) {
                msNode = action(updateParent, node, context);
                if (msNode && node.ref) {
                    UpdateTarget.tag(msNode, mvsRefTags(node.ref));
                    mvsRefMap.set(node.ref, msNode.selector.ref);
                }
                mapping.set(node, msNode);
            } else {
                console.warn(`No target found for this "${node.kind}" node`);
                return;
            }
        }
        if (updateParent) {
            for (const { ext, extCtx } of extensionContexts) {
                ext.action(msNode ?? updateParent, node, context, extCtx);
            }
        }
    });

    for (const target of updateRoot.targetManager.allTargets) {
        UpdateTarget.dependsOn(target, mvsRefMap);
    }

    extensionContexts.forEach(e => e.ext.disposeExtensionContext?.(e.extCtx, tree, context));
    await UpdateTarget.commit(updateRoot);
}


/** A wrapper for updating Mol* state, while using deterministic transform refs.
 * ```
 * updateTarget = UpdateTarget.create(plugin); // like update = plugin.build();
 * UpdateTarget.apply(updateTarget, transformer, params); // like update.to(selector).apply(transformer, params);
 * await UpdateTarget.commit(updateTarget); // like await update.commit();
 * ```
 */
export interface UpdateTarget {
    readonly update: StateBuilder.Root,
    readonly selector: StateObjectSelector,
    readonly targetManager: TargetManager,
    readonly mvsDependencyRefs: Set<string>,
}
export const UpdateTarget = {
    /** Create a new update, with `selector` pointing to the root. */
    create(plugin: PluginContext, replaceExisting: boolean): UpdateTarget {
        const update = plugin.build();
        const msTarget = update.toRoot().selector;
        return { update, selector: msTarget, targetManager: new TargetManager(plugin, replaceExisting), mvsDependencyRefs: new Set() };
    },
    /** Add a child node to `target.selector`, return a new `UpdateTarget` pointing to the new child. */
    apply<A extends StateObject, B extends StateObject, P extends {}>(target: UpdateTarget, transformer: StateTransformer<A, B, P>, params?: Partial<P>, options?: Partial<StateTransform.Options>): UpdateTarget {
        let refSuffix: string = transformer.id;
        if (transformer.id === StructureRepresentation3D.id) {
            const reprType = (params as any)?.type?.name ?? '';
            refSuffix += `:${reprType}`;
        }
        const ref = target.targetManager.getChildRef(target.selector, refSuffix);
        const msResult = target.update.to(target.selector).apply(transformer, params, { ...options, ref }).selector;
        const result: UpdateTarget = { ...target, selector: msResult, mvsDependencyRefs: new Set() };
        target.targetManager.allTargets.push(result);
        return result;
    },
    setMvsDependencies(target: UpdateTarget, refs: string[] | Set<string>): UpdateTarget {
        refs.forEach(ref => target.mvsDependencyRefs.add(ref));
        return target;
    },
    dependsOn(target: UpdateTarget, mapping: Map<string, string>): UpdateTarget {
        if (!target.mvsDependencyRefs.size) return target;
        const dependsOn = Array.from(target.mvsDependencyRefs).map(d => mapping.get(d)!).filter(d => d);
        if (!dependsOn.length) return target;
        target.update.to(target.selector).dependsOn(dependsOn);
        return target;
    },
    /** Add tags to `target.selector` */
    tag(target: UpdateTarget, tags: string[]): UpdateTarget {
        if (tags.length > 0) {
            target.update.to(target.selector).tag(tags);
        }
        return target;
    },
    /** Delete all children of `target.selector`. */
    deleteChildren(target: UpdateTarget): UpdateTarget {
        const children = target.update.currentTree.children.get(target.selector.ref);
        children.forEach(child => target.update.delete(child));
        return target;
    },
    /** Commit all changes done in the current update. */
    commit(target: UpdateTarget): Promise<void> {
        return target.update.commit();
    },
};

/** Manages transform refs in a deterministic way. Uses refs like !mvs:3ce3664304d32c5d:0 */
class TargetManager {
    /** For each hash (e.g. 3ce3664304d32c5d), store the number of already used refs with that hash. */
    private _counter: Record<string, number> = {};
    constructor(plugin: PluginContext, replaceExisting: boolean) {
        if (!replaceExisting) {
            plugin.state.data.cells.forEach(cell => {
                const ref = cell.transform.ref;
                if (ref.startsWith('!mvs:')) {
                    const [_, hash, idNumber] = ref.split(':');
                    const nextIdNumber = parseInt(idNumber) + 1;
                    if (nextIdNumber > (this._counter[hash] ?? 0)) {
                        this._counter[hash] = nextIdNumber;
                    }
                }
            });
        }
    }
    /** Return ref for a new node with given `hash`; update the counter accordingly. */
    private nextRef(hash: string): string {
        this._counter[hash] ??= 0;
        const idNumber = this._counter[hash]++;
        return `!mvs:${hash}:${idNumber}`;
    }
    /** Return ref for a new node based on parent and desired suffix. */
    getChildRef(parent: StateObjectSelector, suffix: string): string {
        const hashBase = parent.ref.replace(/^!mvs:/, '') + ':' + suffix;
        const hash = stringHash(hashBase);
        const result = this.nextRef(hash);
        return result;
    }

    readonly allTargets: UpdateTarget[] = [];
}

/** Create node tags based of MVS node.ref */
export function mvsRefTags(mvsNodeRef: string | undefined): string[] {
    if (mvsNodeRef === undefined) return [];
    else return [`mvs-ref:${mvsNodeRef}`];
}
