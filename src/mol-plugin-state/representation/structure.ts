/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayFind } from '../../mol-data/util';
import { Structure } from '../../mol-model/structure';
import { StateTransform, StateTree, StateSelection } from '../../mol-state';
import { Task } from '../../mol-task';
import { isProductionMode } from '../../mol-util/debug';
import { objectForEach } from '../../mol-util/object';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginContext } from '../../mol-plugin/context';
import { PresetStructureReprentations } from './structure/preset';
import { StructureRepresentationProvider, RepresentationProviderTags } from './structure/provider';
import { UniqueArray } from '../../mol-data/generic';

// TODO: support quality
// TODO: support ignore hydrogens

export class StructureRepresentationManager {
    private providers: StructureRepresentationProvider[] = [];
    private providerMap: Map<string, StructureRepresentationProvider> = new Map();

    readonly defaultProvider = PresetStructureReprentations.auto;

    hasProvider(s: Structure) {
        for (const p of this.providers) {
            if (!p.isApplicable || p.isApplicable(s, this.plugin)) return true;
        }
        return false;
    }

    getOptions(s: Structure) {
        const options: [string, string][] = [];
        const map: { [K in string]: PD.Any } = Object.create(null);
        for (const p of this.providers) {
            if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;

            options.push([p.id, p.display.name]);
            map[p.id] = p.params ? PD.Group(p.params(s, this.plugin)) : PD.EmptyGroup()
        }
        if (options.length === 0) return PD.MappedStatic('', { '': PD.EmptyGroup() });
        return PD.MappedStatic(options[0][0], map, { options });
    }

    hasManagedRepresentation(ref: StateTransform.Ref) {
        const tree = this.plugin.state.dataState.tree;
        // TODO: make this state selection function?
        return StateTree.doPreOrder(tree, tree.transforms.get(ref), { found: false, map: this.providerMap }, (n, _, s) => {
            if (!n.tags) return;
            for (const t of n.tags) {
                if (s.map.has(t)) {
                    s.found = true;
                    return false;
                }
            }
        }).found;
    }

    getManagedRepresentations(ref: StateTransform.Ref) {
        // TODO: check if Structure etc.
        const tree = this.plugin.state.dataState.tree;
        return StateTree.doPreOrder(tree, tree.transforms.get(ref), { found: UniqueArray.create<string, StructureRepresentationProvider>(), map: this.providerMap }, (n, _, s) => {
            if (!n.tags) return;
            for (const t of n.tags) {
                if (s.map.has(t)) UniqueArray.add(s.found, t, s.map.get(t)!);
            }
        }).found.array;
    }

    register(provider: StructureRepresentationProvider) {
        if (this.providerMap.has(provider.id)) {
            throw new Error(`Repr. provider with id '${provider.id}' already registered.`);
        }
        // TODO: sort by group
        this.providers.push(provider);
        this.providerMap.set(provider.id, provider);
    }

    remove(providerOrId: StructureRepresentationProvider | string, structureRoot?: StateTransform.Ref) {
        const root = structureRoot || StateTransform.RootRef;
        const id = typeof providerOrId === 'string' ? providerOrId : providerOrId.id;

        const state = this.plugin.state.dataState;
        const reprs = StateSelection.findWithAllTags(state.tree, root, new Set([id, RepresentationProviderTags.Representation]));

        const builder = state.build();
        for (const r of reprs) {
            builder.delete(r.ref);
        }

        const tree = builder.currentTree;
        const selections = StateSelection.findWithAllTags(tree, root, new Set([RepresentationProviderTags.Selection]));

        for (const s of selections) {
            if (!tree.children.has(s.ref) || tree.children.get(s.ref).size === 0) builder.delete(s.ref);
        }

        if (builder.editInfo.count === 0) return;
        return this.plugin.runTask(state.updateTree(builder));
    }

    apply<P = any, S = {}>(ref: StateTransform.Ref, providerOrId: StructureRepresentationProvider<P, S> | string, params?: P) {
        const provider = typeof providerOrId === 'string'
            ? arrayFind(this.providers, p => p.id === providerOrId)
            : providerOrId;
        if (!provider) return;

        const state = this.plugin.state.dataState;
        const cell = state.cells.get(ref);
        if (!cell || !cell.obj || cell.status !== 'ok') {
            if (!isProductionMode) console.warn(`Applying structure repr. provider to bad cell.`);
            return;
        }

        const prms = params || (provider.params
            ? PD.getDefaultValues(provider.params(cell.obj.data, this.plugin))
            : {})


        const task = Task.create<S>(`${provider.display.name}`, ctx => provider.apply(ctx, state, cell, prms, this.plugin) as Promise<S>);
        return this.plugin.runTask(task);
    }

    // init() {
    //     objectForEach(PresetStructureReprentations, r => this.register(r));
    // }

    constructor(public plugin: PluginContext) {
        objectForEach(PresetStructureReprentations, r => this.register(r));
    }
}