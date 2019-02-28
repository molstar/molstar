/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EmptyLoci, Loci } from 'mol-model/loci';
import { Structure, StructureElement } from 'mol-model/structure';
import { State, StateObject, StateSelection } from 'mol-state';
import { PluginContext } from '../context';
import { PluginStateObject } from '../state/objects';
import { OrderedSet } from 'mol-data/int';

export { StructureElementSelectionManager };

class StructureElementSelectionManager {
    private entries = new Map<string, SelectionEntry>();

    // maps structure to a parent StateObjectCell
    private mapping = {
        root: new Map<Structure, string>(),
        tracked: new Map<string, Structure>()
    };

    private getEntry(s: Structure) {
        if (!this.mapping.root.has(s)) return;
        const key = this.mapping.root.get(s)!;
        if (!this.entries.has(key)) {
            const entry = SelectionEntry(s);
            this.entries.set(key, entry);
            return entry;
        }

        return this.entries.get(key)!;
    }

    add(loci: StructureElement.Loci): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        entry.selection = StructureElement.Loci.union(entry.selection, loci);
        return entry.selection;
    }

    remove(loci: StructureElement.Loci): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        entry.selection = StructureElement.Loci.subtract(entry.selection, loci);
        return entry.selection.elements.length === 0 ? EmptyLoci : entry.selection;
    }

    set(loci: StructureElement.Loci): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        entry.selection = loci;
        return entry.selection.elements.length === 0 ? EmptyLoci : entry.selection;
    }

    get(structure: Structure) {
        const entry = this.getEntry(structure);
        if (!entry) return EmptyLoci;
        return entry.selection;
    }

    has(loci: StructureElement.Loci) {
        const entry = this.getEntry(loci.structure);
        if (!entry) return false;
        return StructureElement.Loci.areIntersecting(loci, entry.selection);
    }

    tryGetRange(loci: StructureElement.Loci): StructureElement.Loci | undefined {
        if (loci.elements.length !== 1) return;
        const entry = this.getEntry(loci.structure);
        if (!entry) return;

        let xs = loci.elements[0];
        let e: StructureElement.Loci['elements'][0] | undefined;
        for (const _e of entry.selection.elements) {
            if (xs.unit === _e.unit) {
                e = _e;
                break;
            }
        }
        if (!e) return;

        const predIdx = OrderedSet.findPredecessorIndex(e.indices, OrderedSet.min(xs.indices));
        if (predIdx === 0) return;

        const fst = predIdx < OrderedSet.size(e.indices)
            ? OrderedSet.getAt(e.indices, predIdx)
            : OrderedSet.getAt(e.indices, predIdx - 1) + 1 as StructureElement.UnitIndex;

        return StructureElement.Loci(entry.selection.structure, [{
            unit: e.unit,
            indices: OrderedSet.ofRange(fst, OrderedSet.max(xs.indices))
        }]);
    }

    private prevHighlight: StructureElement.Loci | undefined = void 0;

    accumulateInteractiveHighlight(loci: StructureElement.Loci) {
        if (this.prevHighlight) {
            this.prevHighlight = StructureElement.Loci.union(this.prevHighlight, loci);
        } else {
            this.prevHighlight = loci;
        }
        return this.prevHighlight;
    }

    clearInteractiveHighlight() {
        const ret = this.prevHighlight;
        this.prevHighlight = void 0;
        return ret || EmptyLoci;
    }

    private addMapping(state: State, ref: string, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;
        const parent = state.select(StateSelection.Generators.byRef(ref).rootOfType([PluginStateObject.Molecule.Structure]))[0];
        this.mapping.tracked.set(ref, obj.data);
        if (!parent) {
            this.mapping.root.set(obj.data, ref);
        } else {
            this.mapping.root.set(obj.data, parent.transform.ref);
        }
    }

    private removeMapping(ref: string) {
        if (!this.mapping.tracked.has(ref)) return;
        const s = this.mapping.tracked.get(ref)!;
        this.mapping.tracked.delete(ref);
        const root = this.mapping.root.get(s);
        this.mapping.root.delete(s);
        if (root === ref) this.entries.delete(ref);
    }

    private updateMapping(state: State, ref: string, oldObj: StateObject | undefined, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;

        if (this.entries.has(ref)) {
            if (!PluginStateObject.Molecule.Structure.is(oldObj) || oldObj === obj || oldObj.data === obj.data) return;

            // remap the old selection to be related to the new object if possible.
            if (Structure.areUnitAndIndicesEqual(oldObj.data, obj.data)) {
                this.entries.set(ref, remapSelectionEntry(this.entries.get(ref)!, obj.data));
                return;
            }

            // clear the selection
            this.entries.set(ref, SelectionEntry(obj.data));
        } else {
            this.removeMapping(ref);
            this.addMapping(state, ref, obj);
        }
    }

    constructor(plugin: PluginContext) {
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

interface SelectionEntry {
    selection: StructureElement.Loci
}

function SelectionEntry(s: Structure): SelectionEntry {
    return {
        selection: StructureElement.Loci(s, [])
    };
}

function remapSelectionEntry(e: SelectionEntry, s: Structure): SelectionEntry {
    return {
        selection: StructureElement.Loci.remap(e.selection, s)
    };
}