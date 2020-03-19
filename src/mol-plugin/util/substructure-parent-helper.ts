/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../../mol-model/structure';
import { State, StateObject, StateSelection, StateObjectCell, StateTransform } from '../../mol-state';
import { PluginContext } from '../context';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { arraySetAdd, arraySetRemove } from '../../mol-util/array';

export { SubstructureParentHelper };

class SubstructureParentHelper {
    private decorators = new Map<string, string[]>();
    private root = new Map<Structure, { ref: string, count: number }>();
    private tracked = new Map<string, Structure>();

    /** Returns the root node of given structure if existing, takes decorators into account */
    get(s: Structure, ignoreDecorators = false): StateObjectCell<PluginStateObject.Molecule.Structure> | undefined {
        const r = this.root.get(s);
        if (!r) return;
        const decorators = this.decorators.get(r.ref);
        if (ignoreDecorators || !decorators) return this.plugin.state.data.cells.get(r.ref);
        return this.plugin.state.data.cells.get(this.findDeepestDecorator(r.ref, decorators));
    }

    private findDeepestDecorator(ref: string, decorators: string[]) {
        if (decorators.length === 0) return ref;
        if (decorators.length === 1) return decorators[0];

        const cells = this.plugin.state.data.cells;
        let depth = 0, ret = ref;
        for (const dr of decorators) {
            let c = cells.get(dr);
            let d = 0;
            while (c && c.transform.ref !== StateTransform.RootRef) {
                d++;
                c = cells.get(c.transform.parent);
            }
            if (d > depth) {
                ret = dr;
                depth = d;
            }
        }
        return ret;
    }

    private addDecorator(root: string, ref: string) {
        if (this.decorators.has(root)) {
            arraySetAdd(this.decorators.get(root)!, ref);
        } else {
            this.decorators.set(root, [ref]);
        }
    }

    private tryRemoveDecorator(root: string, ref: string) {
        if (this.decorators.has(root)) {
            const xs = this.decorators.get(root)!;
            arraySetRemove(xs, ref);
            if (xs.length === 0) this.decorators.delete(root);
        }
    }

    private addMapping(state: State, ref: string, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;

        this.tracked.set(ref, obj.data);

        let parentRef;
        // if the structure is already present in the tree, do not rewrite the root.
        if (this.root.has(obj.data)) {
            const e = this.root.get(obj.data)!;
            parentRef = e.ref;
            e.count++;
        } else {
            const parent = state.select(StateSelection.Generators.byRef(ref).rootOfType([PluginStateObject.Molecule.Structure]))[0];

            if (!parent) {
                this.root.set(obj.data, { ref, count: 1 });
            } else {
                parentRef = parent.transform.ref;
                this.root.set(obj.data, { ref: parentRef, count: 1 });
            }
        }

        if (!parentRef) return;

        const cell = state.cells.get(ref);
        if (cell?.transform.isDecorator) {
            this.addDecorator(parentRef, ref);
        }
    }

    private removeMapping(ref: string) {
        if (!this.tracked.has(ref)) return;

        const s = this.tracked.get(ref)!;
        this.tracked.delete(ref);

        const root = this.root.get(s)!;

        this.tryRemoveDecorator(root.ref, ref);

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