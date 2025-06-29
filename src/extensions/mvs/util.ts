/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StateObjectSelector, StateTree } from '../../mol-state';


/**
 * Queries all MolViewSpec references in the current state of the plugin.
 */
export function queryMVSRef(plugin: PluginContext, ref: string) {
    return plugin.state.data.selectQ(q => q.root.subtree().withTag(`mvs-ref:${ref}`));
}

/**
 * Creates a mapping of all MolViewSpec references in the current state of the plugin.
 */
export function createMVSRefMap(plugin: PluginContext) {
    const tree = plugin.state.data.tree;
    const mapping = new Map<string, StateObjectSelector[]>();
    StateTree.doPreOrder(tree, tree.root, { mapping, plugin }, (n, _, s) => {
        if (!n.tags) return;
        for (const tag of n.tags) {
            if (!tag.startsWith('mvs-ref:')) continue;
            const mvsRef = tag.substring(8);
            const selector = new StateObjectSelector(n.ref, s.plugin.state.data);
            if (s.mapping.has(mvsRef)) s.mapping.get(mvsRef)!.push(selector);
            else s.mapping.set(mvsRef, [selector]);
        }
    });

    return mapping;
}