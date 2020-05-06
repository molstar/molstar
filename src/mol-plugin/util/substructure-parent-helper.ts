/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { State, StateObject, StateObjectCell, StateSelection } from '../../mol-state';
import { PluginContext } from '../context';
import { RxEventHelper } from '../../mol-util/rx-event-helper';

export { SubstructureParentHelper };

class SubstructureParentHelper {
    private ev = RxEventHelper.create();

    readonly events = {
        updated: this.ev<{ ref: string, oldObj: PluginStateObject.Molecule.Structure | undefined, obj: PluginStateObject.Molecule.Structure }>(),
        removed: this.ev<{ ref: string, obj: PluginStateObject.Molecule.Structure | undefined }>(),
    }

    // private decorators = new Map<string, string[]>();
    private root = new Map<Structure, { ref: string, count: number }>();
    private tracked = new Map<string, Structure>();

    private getDecorator(root: string): string {
        const tree = this.plugin.state.data.tree;
        const children = tree.children.get(root);
        if (children.size !== 1) return root;
        const child = children.first();
        if (tree.transforms.get(child).transformer.definition.isDecorator) {
            return this.getDecorator(child);
        }
        return root;
    }

    /** Returns the root node of given structure if existing, takes decorators into account */
    get(s: Structure, ignoreDecorators = false): StateObjectCell<PluginStateObject.Molecule.Structure> | undefined {
        const r = this.root.get(s);
        if (!r) return;
        if (ignoreDecorators) return this.plugin.state.data.cells.get(r.ref);
        return this.plugin.state.data.cells.get(this.getDecorator(r.ref));
    }

    private addMapping(state: State, ref: string, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return false;

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
        return true;
    }

    private removeMapping(ref: string) {
        if (!this.tracked.has(ref)) return false;

        const s = this.tracked.get(ref)!;
        this.tracked.delete(ref);

        const root = this.root.get(s)!;

        if (root.count > 1) {
            root.count--;
        } else {
            this.root.delete(s);
        }
        return true;
    }

    private updateMapping(state: State, ref: string, oldObj: StateObject | undefined, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return false;

        this.removeMapping(ref);
        this.addMapping(state, ref, obj);
        return true;
    }

    dispose() {
        this.ev.dispose();
    }

    constructor(private plugin: PluginContext) {
        plugin.state.data.events.object.created.subscribe(e => {
            this.addMapping(e.state, e.ref, e.obj);
        });

        plugin.state.data.events.object.removed.subscribe(e => {
            if (this.removeMapping(e.ref)) {
                this.events.removed.next({ ref: e.ref, obj: e.obj });
            }
        });

        plugin.state.data.events.object.updated.subscribe(e => {
            if (this.updateMapping(e.state, e.ref, e.oldObj, e.obj)) {
                this.events.updated.next({ ref: e.ref, oldObj: e.oldObj, obj: e.obj });
            }
        });
    }
}