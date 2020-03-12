/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayFind } from '../../../mol-data/util';
import { Structure } from '../../../mol-model/structure';
import { StateTransform, StateTree, StateSelection, StateObjectRef, StateBuilder } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { isProductionMode } from '../../../mol-util/debug';
import { objectForEach } from '../../../mol-util/object';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginContext } from '../../../mol-plugin/context';
import { PresetStructureReprentations } from './preset';
import { StructureRepresentationProvider, RepresentationProviderTags } from './provider';
import { UniqueArray } from '../../../mol-data/generic';
import { PluginStateObject } from '../../objects';
import { StructureRepresentation3D, StructureRepresentation3DHelpers } from '../../transforms/representation';
import { RepresentationProvider } from '../../../mol-repr/representation';
import { SizeTheme, BuiltInSizeThemeName } from '../../../mol-theme/size';
import { ColorTheme, BuiltInColorThemeName } from '../../../mol-theme/color';
import { BuiltInStructureRepresentationsName } from '../../../mol-repr/structure/registry';

export type StructureRepresentationProviderRef = keyof PresetStructureReprentations | StructureRepresentationProvider | string

export class StructureRepresentationBuilder {
    private providers: StructureRepresentationProvider[] = [];
    private providerMap: Map<string, StructureRepresentationProvider> = new Map();
    private get dataState() { return this.plugin.state.dataState; }

    readonly defaultProvider = PresetStructureReprentations.auto;

    private resolveProvider(ref: StructureRepresentationProviderRef) {
        return typeof ref === 'string'
            ? PresetStructureReprentations[ref as keyof PresetStructureReprentations] ?? arrayFind(this.providers, p => p.id === ref)
            : ref;
    }

    hasPreset(s: Structure) {
        for (const p of this.providers) {
            if (!p.isApplicable || p.isApplicable(s, this.plugin)) return true;
        }
        return false;
    }

    get providerList(): ReadonlyArray<StructureRepresentationProvider> { return this.providers; }

    getPresets(s: Structure) {
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

    hasPresetRepresentation(ref: StateObjectRef) {
        // TODO: make this state selection function?
        const tree = this.plugin.state.dataState.tree;
        const root = StateObjectRef.resolve(this.plugin.state.dataState, ref);
        if (!root) return false;
        return StateTree.doPreOrder(tree, root.transform,  { found: false, map: this.providerMap }, (n, _, s) => {
            if (!n.tags) return;
            for (const t of n.tags) {
                if (s.map.has(t)) {
                    s.found = true;
                    return false;
                }
            }
        }).found;
    }

    getPresetRepresentations(ref: StateObjectRef) {
        const tree = this.plugin.state.dataState.tree;
        const root = StateObjectRef.resolve(this.plugin.state.dataState, ref);
        if (!root) return [];
        return StateTree.doPreOrder(tree, root.transform, { found: UniqueArray.create<string, StructureRepresentationProvider>(), map: this.providerMap }, (n, _, s) => {
            if (!n.tags) return;
            for (const t of n.tags) {
                if (s.map.has(t)) UniqueArray.add(s.found, t, s.map.get(t)!);
            }
        }).found.array;
    }

    registerPreset(provider: StructureRepresentationProvider) {
        if (this.providerMap.has(provider.id)) {
            throw new Error(`Repr. provider with id '${provider.id}' already registered.`);
        }
        // TODO: sort by group
        this.providers.push(provider);
        this.providerMap.set(provider.id, provider);
    }

    removePreset(providerRef: StructureRepresentationProviderRef, structureRoot?: StateObjectRef) {
        const id = this.resolveProvider(providerRef)?.id;
        if (!id) return;

        const state = this.plugin.state.dataState;
        const root = StateObjectRef.resolveRef(structureRoot) || StateTransform.RootRef;
        const reprs = StateSelection.findWithAllTags(state.tree, root, new Set([id, RepresentationProviderTags.Representation]));

        const builder = state.build();
        for (const r of reprs) {
            builder.delete(r.ref);
        }

        const tree = builder.currentTree;
        const selections = StateSelection.findWithAllTags(tree, root, new Set([RepresentationProviderTags.Component]));

        for (const s of selections) {
            if (!tree.children.has(s.ref) || tree.children.get(s.ref).size === 0) builder.delete(s.ref);
        }

        if (builder.editInfo.count === 0) return;
        return this.plugin.runTask(state.updateTree(builder));
    }

    structurePreset<K extends keyof PresetStructureReprentations>(parent: StateObjectRef, preset: K, params?: StructureRepresentationProvider.Params<PresetStructureReprentations[K]>): Promise<StructureRepresentationProvider.State<PresetStructureReprentations[K]>> | undefined
    structurePreset<P = any, S = {}>(parent: StateObjectRef, providers: StructureRepresentationProvider<P, S>, params?: P): Promise<S> | undefined
    structurePreset(parent: StateObjectRef, providerId: string, params?: any): Promise<any> | undefined
    structurePreset(parent: StateObjectRef, providerRef: string | StructureRepresentationProvider, params?: any): Promise<any> | undefined {
        const provider = this.resolveProvider(providerRef);
        if (!provider) return;

        const state = this.plugin.state.dataState;
        const cell = StateObjectRef.resolveAndCheck(state, parent);
        if (!cell) {
            if (!isProductionMode) console.warn(`Applying structure repr. provider to bad cell.`);
            return;
        }

        const prms = params || (provider.params
            ? PD.getDefaultValues(provider.params(cell.obj!.data, this.plugin) as PD.Params)
            : {})


        const task = Task.create(`${provider.display.name}`, ctx => provider.apply(ctx, state, cell, prms, this.plugin) as Promise<any>);
        return this.plugin.runTask(task);
    }

    async addRepresentation<R extends RepresentationProvider<Structure, any, any>, C extends ColorTheme.Provider<any>, S extends SizeTheme.Provider<any>>
        (structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props?: StructureRepresentation3DHelpers.Props<R, C, S>) {

        const data = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;
        if (!data) return;

        const params = StructureRepresentation3DHelpers.createParams(this.plugin, data, props);
        const repr = this.dataState.build()
            .to(structure)
            .apply(StructureRepresentation3D, params, { tags: RepresentationProviderTags.Representation });

        await this.plugin.runTask(this.dataState.updateTree(repr));
        return  repr.selector;
    }

    builtInRepresentation<R extends BuiltInStructureRepresentationsName, C extends BuiltInColorThemeName, S extends BuiltInSizeThemeName>
        (builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props?: StructureRepresentation3DHelpers.BuildInProps<R, C, S>) {

        if (!structure) return;
        const data = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;
        if (!data) return;

        const params = StructureRepresentation3DHelpers.createBuiltInParams(this.plugin, data, props);
        return builder
            .to(structure)
            .apply(StructureRepresentation3D, params, { tags: RepresentationProviderTags.Representation })
            .selector;
    }

    constructor(public plugin: PluginContext) {
        objectForEach(PresetStructureReprentations, r => this.registerPreset(r));
    }
}