/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { State, StateObject, StateSelection, StateObjectCell } from 'mol-state';
import { PluginContext } from '../context';
import { PluginStateObject } from '../state/objects';

export { SubstructureParentHelper };

class SubstructureParentHelper {
    private root = new Map<Structure, string>();
    private tracked = new Map<string, Structure>();

    get(s: Structure): StateObjectCell<PluginStateObject.Molecule.Structure> | undefined {
        const r = this.root.get(s);
        if (!r) return;
        return this.plugin.state.dataState.cells.get(r);
    }

    private addMapping(state: State, ref: string, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;
        const parent = state.select(StateSelection.Generators.byRef(ref).rootOfType([PluginStateObject.Molecule.Structure]))[0];
        this.tracked.set(ref, obj.data);
        if (!parent) {
            this.root.set(obj.data, ref);
        } else {
            this.root.set(obj.data, parent.transform.ref);
        }
    }

    private removeMapping(ref: string) {
        if (!this.tracked.has(ref)) return;
        const s = this.tracked.get(ref)!;
        this.tracked.delete(ref);
        this.root.get(s);
        this.root.delete(s);
    }

    private updateMapping(state: State, ref: string, oldObj: StateObject | undefined, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;

        this.removeMapping(ref);
        this.addMapping(state, ref, obj);
    }

    constructor(private plugin: PluginContext) {
        plugin.state.dataState.events.object.created.subscribe(e => {
            this.addMapping(e.state, e.ref, e.obj);
        });

        plugin.state.dataState.events.object.removed.subscribe(e => {
            this.removeMapping(e.ref);
        });

        plugin.state.dataState.events.object.updated.subscribe(e => {
            this.updateMapping(e.state, e.ref, e.oldObj, e.obj);
        });
    }
}