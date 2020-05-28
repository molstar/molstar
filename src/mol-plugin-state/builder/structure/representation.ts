/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayFind } from '../../../mol-data/util';
import { PluginContext } from '../../../mol-plugin/context';
import { StateBuilder, StateObjectRef, StateObjectSelector, StateTransform } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { isProductionMode } from '../../../mol-util/debug';
import { objectForEach } from '../../../mol-util/object';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { createStructureRepresentationParams, StructureRepresentationBuiltInProps, StructureRepresentationProps } from '../../helpers/structure-representation-params';
import { PluginStateObject } from '../../objects';
import { StructureRepresentation3D } from '../../transforms/representation';
import { PresetStructureRepresentations, StructureRepresentationPresetProvider } from './representation-preset';
import { arrayRemoveInPlace } from '../../../mol-util/array';
import { PluginConfig } from '../../../mol-plugin/config';

// TODO factor out code shared with TrajectoryHierarchyBuilder?

export type StructureRepresentationPresetProviderRef = keyof PresetStructureRepresentations | StructureRepresentationPresetProvider | string

export class StructureRepresentationBuilder {
    private _providers: StructureRepresentationPresetProvider[] = [];
    private providerMap: Map<string, StructureRepresentationPresetProvider> = new Map();
    private get dataState() { return this.plugin.state.data; }

    readonly defaultProvider = PresetStructureRepresentations.auto;

    private resolveProvider(ref: StructureRepresentationPresetProviderRef) {
        return typeof ref === 'string'
            ? PresetStructureRepresentations[ref as keyof PresetStructureRepresentations] ?? arrayFind(this._providers, p => p.id === ref)
            : ref;
    }

    hasPreset(s: PluginStateObject.Molecule.Structure) {
        for (const p of this._providers) {
            if (!p.isApplicable || p.isApplicable(s, this.plugin)) return true;
        }
        return false;
    }

    get providers(): ReadonlyArray<StructureRepresentationPresetProvider> { return this._providers; }

    getPresets(s?: PluginStateObject.Molecule.Structure) {
        if (!s) return this.providers;
        const ret = [];
        for (const p of this._providers) {
            if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;
            ret.push(p);
        }
        return ret;
    }

    getPresetSelect(s?: PluginStateObject.Molecule.Structure): PD.Select<string> {
        const options: [string, string, string | undefined][] = [];
        for (const p of this._providers) {
            if (s && p.isApplicable && !p.isApplicable(s, this.plugin)) continue;
            options.push([p.id, p.display.name, p.display.group]);
        }
        return PD.Select('auto', options);
    }

    getPresetsWithOptions(s: PluginStateObject.Molecule.Structure) {
        const options: [string, string][] = [];
        const map: { [K in string]: PD.Any } = Object.create(null);
        for (const p of this._providers) {
            if (p.isApplicable && !p.isApplicable(s, this.plugin)) continue;

            options.push([p.id, p.display.name]);
            map[p.id] = p.params ? PD.Group(p.params(s, this.plugin)) : PD.EmptyGroup();
        }
        if (options.length === 0) return PD.MappedStatic('', { '': PD.EmptyGroup() });
        return PD.MappedStatic(options[0][0], map, { options });
    }

    registerPreset(provider: StructureRepresentationPresetProvider) {
        if (this.providerMap.has(provider.id)) {
            throw new Error(`Representation provider with id '${provider.id}' already registered.`);
        }
        this._providers.push(provider);
        this.providerMap.set(provider.id, provider);
    }

    unregisterPreset(provider: StructureRepresentationPresetProvider) {
        this.providerMap.delete(provider.id);
        arrayRemoveInPlace(this._providers, provider);
    }

    applyPreset<K extends keyof PresetStructureRepresentations>(parent: StateObjectRef<PluginStateObject.Molecule.Structure>, preset: K, params?: StructureRepresentationPresetProvider.Params<PresetStructureRepresentations[K]>): Promise<StructureRepresentationPresetProvider.State<PresetStructureRepresentations[K]>> | undefined
    applyPreset<P = any, S = {}>(parent: StateObjectRef<PluginStateObject.Molecule.Structure>, provider: StructureRepresentationPresetProvider<P, S>, params?: P): Promise<S> | undefined
    applyPreset(parent: StateObjectRef<PluginStateObject.Molecule.Structure>, providerId: string, params?: any): Promise<any> | undefined
    applyPreset(parent: StateObjectRef, providerRef: string | StructureRepresentationPresetProvider, params?: any): Promise<any> | undefined {
        const provider = this.resolveProvider(providerRef);
        if (!provider) return;

        const state = this.plugin.state.data;
        const cell = StateObjectRef.resolveAndCheck(state, parent);
        if (!cell) {
            if (!isProductionMode) console.warn(`Applying structure repr. provider to bad cell.`);
            return;
        }

        const pd = provider.params?.(cell.obj, this.plugin) as PD.Params || {};
        let prms = params || (provider.params
            ? PD.getDefaultValues(pd)
            : {});

        const defaults = this.plugin.config.get(PluginConfig.Structure.DefaultRepresentationPresetParams);
        prms = PD.merge(pd, defaults, prms);

        const task = Task.create(`${provider.display.name}`, () => provider.apply(cell, prms, this.plugin) as Promise<any>);
        return this.plugin.runTask(task);
    }

    async addRepresentation<P extends StructureRepresentationBuiltInProps>(structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props: P, options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>
    async addRepresentation<P extends StructureRepresentationProps>(structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props: P, options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>): Promise<StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>>
    async addRepresentation(structure: StateObjectRef<PluginStateObject.Molecule.Structure>, props: any, options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>) {
        const repr = this.dataState.build();
        const selector = this.buildRepresentation(repr, structure, props, options);
        if (!selector) return;

        await repr.commit();
        return selector;
    }

    buildRepresentation<P extends StructureRepresentationBuiltInProps>(builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props: P, options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>): StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>
    buildRepresentation<P extends StructureRepresentationProps>(builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props: P, options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>): StateObjectSelector<PluginStateObject.Molecule.Structure.Representation3D>
    buildRepresentation(builder: StateBuilder.Root, structure: StateObjectRef<PluginStateObject.Molecule.Structure> | undefined, props: any, options?: Partial<StructureRepresentationBuilder.AddRepresentationOptions>) {
        if (!structure) return;
        const data = StateObjectRef.resolveAndCheck(this.dataState, structure)?.obj?.data;
        if (!data) return;

        const params = createStructureRepresentationParams(this.plugin, data, props);
        return options?.tag
            ? builder.to(structure).applyOrUpdateTagged(options.tag, StructureRepresentation3D, params, { state: options?.initialState }).selector
            : builder.to(structure).apply(StructureRepresentation3D, params, { state: options?.initialState }).selector;
    }

    constructor(public plugin: PluginContext) {
        objectForEach(PresetStructureRepresentations, r => this.registerPreset(r));
    }
}

export namespace StructureRepresentationBuilder {
    export interface AddRepresentationOptions {
        initialState?: Partial<StateTransform.State>,
        tag?: string
    }
}