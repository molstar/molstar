/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayFind } from '../../../mol-data/util';
import { Structure } from '../../../mol-model/structure';
import { PluginContext } from '../../../mol-plugin/context';
import { StateBuilder, StateObjectRef, StateObjectSelector } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { isProductionMode } from '../../../mol-util/debug';
import { objectForEach } from '../../../mol-util/object';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { createStructureRepresentationParams, StructureRepresentationBuiltInProps, StructureRepresentationProps } from '../../helpers/structure-representation-params';
import { PluginStateObject } from '../../objects';
import { StructureRepresentation3D } from '../../transforms/representation';
import { PresetStructureReprentations } from './preset';
import { RepresentationProviderTags, StructureRepresentationProvider } from './provider';

export type StructureRepresentationProviderRef = keyof PresetStructureReprentations | StructureRepresentationProvider | string

export class StructureRepresentationBuilder {
    private _providers: StructureRepresentationProvider[] = [];
    private providerMap: Map<string, StructureRepresentationProvider> = new Map();
    private get dataState() { return this.plugin.state.dataState; }

    readonly defaultProvider = PresetStructureReprentations.auto;

    private resolveProvider(ref: StructureRepresentationProviderRef) {
        return typeof ref === 'string'
            ? PresetStructureReprentations[ref as keyof PresetStructureReprentations] ?? arrayFind(this._providers, p => p.id === ref)
            : ref;
    }

    hasPreset(s: Structure) {
        for (const p of this._providers) {
            if (!p.isApplicable || p.isApplicable(s, this.plugin)) return true;
        }
        return false;
    }

    get providers(): ReadonlyArray<StructureRepresentationProvider> { return this._providers; }

    getPresets(s?: Structure) {
        if (!s) return this.providers;
        const ret = [];
        for (const p of this._providers) {
            if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;
            ret.push(p);
        }
        return ret;
    }

    getPresetsWithOptions(s: Structure) {
        const options: [string, string][] = [];
        const map: { [K in string]: PD.Any } = Object.create(null);
        for (const p of this._providers) {
            if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;

            options.push([p.id, p.display.name]);
            map[p.id] = p.params ? PD.Group(p.params(s, this.plugin)) : PD.EmptyGroup()
        }
        if (options.length === 0) return PD.MappedStatic('', { '': PD.EmptyGroup() });
        return PD.MappedStatic(options[0][0], map, { options });
    }

    registerPreset(provider: StructureRepresentationProvider) {
        if (this.providerMap.has(provider.id)) {
            throw new Error(`Repr. provider with id '${provider.id}' already registered.`);
        }
        this._providers.push(provider);
        this.providerMap.set(provider.id, provider);
    }

    applyPreset<K extends keyof PresetStructureReprentations>(parent: StateObjectRef, preset: K, params?: StructureRepresentationProvider.Params<PresetStructureReprentations[K]>): Promise<StructureRepresentationProvider.State<PresetStructureReprentations[K]>> | undefined
    applyPreset<P = any, S = {}>(parent: StateObjectRef, providers: StructureRepresentationProvider<P, S>, params?: P): Promise<S> | undefined
    applyPreset(parent: StateObjectRef, providerId: string, params?: any): Promise<any> | undefined
    applyPreset(parent: StateObjectRef, providerRef: string | StructureRepresentationProvider, params?: any): Promise<any> | undefined {
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

    async addRepresentation<P extends StructureRepresentationBuiltInProps>(structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props?: P): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>
    async addRepresentation<P extends StructureRepresentationProps>(structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props?: P): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>
    async addRepresentation(structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props?: any) {
        const data = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;
        if (!data) return;

        const params = createStructureRepresentationParams(this.plugin, data, props);
        const repr = this.dataState.build()
            .to(structure)
            .apply(StructureRepresentation3D, params, { tags: RepresentationProviderTags.Representation });

        await this.plugin.runTask(this.dataState.updateTree(repr));
        return  repr.selector;
    }

    async buildRepresentation<P extends StructureRepresentationBuiltInProps>(builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props?: P): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>
    async buildRepresentation<P extends StructureRepresentationProps>(builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props?: P): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>
    async buildRepresentation(builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props?: any) {
        if (!structure) return;
        const data = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;
        if (!data) return;

        const params = createStructureRepresentationParams(this.plugin, data, props);
        return builder
            .to(structure)
            .apply(StructureRepresentation3D, params, { tags: RepresentationProviderTags.Representation })
            .selector;
    }

    constructor(public plugin: PluginContext) {
        objectForEach(PresetStructureReprentations, r => this.registerPreset(r));
    }
}