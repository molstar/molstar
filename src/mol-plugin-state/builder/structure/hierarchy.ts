/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { arrayFind } from '../../../mol-data/util';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectRef } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { isProductionMode } from '../../../mol-util/debug';
import { objectForEach } from '../../../mol-util/object';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';
import { PresetTrajectoryHierarchy, TrajectoryHierarchyPresetProvider } from './hierarchy-preset';
import { arrayRemoveInPlace } from '../../../mol-util/array';

// TODO factor out code shared with StructureRepresentationBuilder?

export type TrajectoryHierarchyPresetProviderRef = keyof PresetTrajectoryHierarchy | TrajectoryHierarchyPresetProvider | string

export class TrajectoryHierarchyBuilder {
    private _providers: TrajectoryHierarchyPresetProvider[] = [];
    private providerMap: Map<string, TrajectoryHierarchyPresetProvider> = new Map();

    readonly defaultProvider = PresetTrajectoryHierarchy.default;

    private resolveProvider(ref: TrajectoryHierarchyPresetProviderRef) {
        return typeof ref === 'string'
            ? PresetTrajectoryHierarchy[ref as keyof PresetTrajectoryHierarchy] ?? arrayFind(this._providers, p => p.id === ref)
            : ref;
    }

    hasPreset(t: PluginStateObject.Molecule.Trajectory) {
        for (const p of this._providers) {
            if (!p.isApplicable || p.isApplicable(t, this.plugin)) return true;
        }
        return false;
    }

    get providers(): ReadonlyArray<TrajectoryHierarchyPresetProvider> { return this._providers; }

    getPresets(t?: PluginStateObject.Molecule.Trajectory) {
        if (!t) return this.providers;
        const ret = [];
        for (const p of this._providers) {
            if (p.isApplicable && !p.isApplicable(t, this.plugin)) continue;
            ret.push(p);
        }
        return ret;
    }

    getPresetSelect(t?: PluginStateObject.Molecule.Trajectory): PD.Select<string> {
        const options: [string, string][] = [];
        for (const p of this._providers) {
            if (t && p.isApplicable && !p.isApplicable(t, this.plugin)) continue;
            options.push([p.id, p.display.name]);
        }
        return PD.Select('auto', options);
    }

    getPresetsWithOptions(t: PluginStateObject.Molecule.Trajectory) {
        const options: [string, string][] = [];
        const map: { [K in string]: PD.Any } = Object.create(null);
        for (const p of this._providers) {
            if (p.isApplicable && !p.isApplicable(t, this.plugin)) continue;

            options.push([p.id, p.display.name]);
            map[p.id] = p.params ? PD.Group(p.params(t, this.plugin)) : PD.EmptyGroup();
        }
        if (options.length === 0) return PD.MappedStatic('', { '': PD.EmptyGroup() });
        return PD.MappedStatic(options[0][0], map, { options });
    }

    registerPreset(provider: TrajectoryHierarchyPresetProvider) {
        if (this.providerMap.has(provider.id)) {
            throw new Error(`Hierarchy provider with id '${provider.id}' already registered.`);
        }
        this._providers.push(provider);
        this.providerMap.set(provider.id, provider);
    }

    unregisterPreset(provider: TrajectoryHierarchyPresetProvider) {
        this.providerMap.delete(provider.id);
        arrayRemoveInPlace(this._providers, provider);
    }

    applyPreset<K extends keyof PresetTrajectoryHierarchy>(parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>, preset: K, params?: Partial<TrajectoryHierarchyPresetProvider.Params<PresetTrajectoryHierarchy[K]>>): Promise<TrajectoryHierarchyPresetProvider.State<PresetTrajectoryHierarchy[K]>> | undefined
    applyPreset<P = any, S = {}>(parent: StateObjectRef<PluginStateObject.Molecule.Trajectory>, provider: TrajectoryHierarchyPresetProvider<P, S>, params?: P): Promise<S> | undefined
    applyPreset(parent: StateObjectRef, providerRef: string | TrajectoryHierarchyPresetProvider, params?: any): Promise<any> | undefined {
        const provider = this.resolveProvider(providerRef);
        if (!provider) return;

        const state = this.plugin.state.data;
        const cell = StateObjectRef.resolveAndCheck(state, parent);
        if (!cell) {
            if (!isProductionMode) console.warn(`Applying hierarchy preset provider to bad cell.`);
            return;
        }

        const prms = params || (provider.params
            ? PD.getDefaultValues(provider.params(cell.obj, this.plugin) as PD.Params)
            : {});

        const task = Task.create(`${provider.display.name}`, () => provider.apply(cell, prms, this.plugin) as Promise<any>);
        return this.plugin.runTask(task);
    }

    constructor(public plugin: PluginContext) {
        objectForEach(PresetTrajectoryHierarchy, r => this.registerPreset(r));
    }
}