/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { State, StateObject, StateObjectCell, StateSelection } from '../../mol-state';
import { PluginContext } from '../context';

export { SubstructureParentHelper };

class SubstructureParentHelper {
    // private decorators = new Map<string, string[]>();
    private root = new Map<Structure, { ref: string, count: number }>();
    private tracked = new Map<string, Structure>();

    /** Returns the root node of given structure if existing, takes decorators into account */
    get(s: Structure, ignoreDecorators = false): StateObjectCell<PluginStateObject.Molecule.Structure> | undefined {
        const r = this.root.get(s);
        if (!r) return;
        return this.plugin.state.data.cells.get(r.ref);
    }

    private addMapping(state: State, ref: string, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;

        this.tracked.set(ref, obj.data);

        // if the structure is already present in the tree, do not rewrite the root.
        if (this.root.has(obj.data)) {
            const e = this.root.get(obj.data)!;
            e.count++;
        } else {
            const parent = state.select(StateSelection.Generators.byRef(ref).rootOfType([PluginStateObject.Molecule.Structure]))[0];

            if (!parent) {
                this.root.set(obj.data, { ref, count: 1 });
            } else {
                this.root.set(obj.data, { ref: parent.transform.ref, count: 1 });
            }
        }
    }

    private removeMapping(ref: string) {
        if (!this.tracked.has(ref)) return;

        const s = this.tracked.get(ref)!;
        this.tracked.delete(ref);

        const root = this.root.get(s)!;

        if (root.count > 1) {
            root.count--;
        } else {
            this.root.delete(s);
        }
    }

    private updateMapping(state: State, ref: string, oldObj: StateObject | undefined, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;

        this.removeMapping(ref);
        this.addMapping(state, ref, obj);
    }

    constructor(private plugin: PluginContext) {
        plugin.state.data.events.object.created.subscribe(e => {
            this.addMapping(e.state, e.ref, e.obj);
        });

        plugin.state.data.events.object.removed.subscribe(e => {
            this.removeMapping(e.ref);
        });

        plugin.state.data.events.object.updated.subscribe(e => {
            this.updateMapping(e.state, e.ref, e.oldObj, e.obj);
        });
    }
}