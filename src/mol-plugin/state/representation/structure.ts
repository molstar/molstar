/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayFind } from '../../../mol-data/util';
import { Structure } from '../../../mol-model/structure';
import { StateTransform } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { isProductionMode } from '../../../mol-util/debug';
import { objectForEach } from '../../../mol-util/object';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginContext } from '../../context';
import { PresetStructureReprentations } from './structure/preset';
import { StructureRepresentationProvider } from './structure/providers';

export class StructureRepresentationManager {
    private providers: StructureRepresentationProvider[] = [];

    readonly defaultProvider = PresetStructureReprentations.default;

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

    register(provider: StructureRepresentationProvider) {
        // TODO: sort by group
        this.providers.push(provider);
    }

    apply<P>(ref: StateTransform.Ref, providerOrId: StructureRepresentationProvider<P> | string, params?: P) {
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

        const apply = provider.apply(state, cell, prms, this.plugin);

        if (Task.is(apply)) return this.plugin.runTask(apply);
        return apply;
    }

    // init() {
    //     objectForEach(PresetStructureReprentations, r => this.register(r));
    // }

    constructor(public plugin: PluginContext) {
        objectForEach(PresetStructureReprentations, r => this.register(r));
    }
}