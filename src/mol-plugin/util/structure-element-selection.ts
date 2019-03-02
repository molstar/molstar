/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int';
import { EmptyLoci, Loci } from 'mol-model/loci';
import { Structure, StructureElement } from 'mol-model/structure';
import { StateObject } from 'mol-state';
import { PluginContext } from '../context';
import { PluginStateObject } from '../state/objects';

export { StructureElementSelectionManager };

class StructureElementSelectionManager {
    private entries = new Map<string, SelectionEntry>();

    private getEntry(s: Structure) {
        const cell = this.plugin.helpers.substructureParent.get(s);
        if (!cell) return;
        const ref = cell.transform.ref;
        if (!this.entries.has(ref)) {
            const entry = SelectionEntry(s);
            this.entries.set(ref, entry);
            return entry;
        }

        return this.entries.get(ref)!;
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

    clear() {
        const keys = this.entries.keys();
        const selections: StructureElement.Loci[] = [];
        while (true) {
            const k = keys.next();
            if (k.done) break;
            const s = this.entries.get(k.value)!;
            if (s.selection.elements.length > 0) selections.push(s.selection);
            s.selection = StructureElement.Loci(s.selection.structure, []);
        }
        return selections;
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

        let predIdx = OrderedSet.findPredecessorIndex(e.indices, OrderedSet.min(xs.indices));
        if (predIdx === 0) return;

        let fst;

        if (predIdx < OrderedSet.size(e.indices)) {
            fst = OrderedSet.getAt(e.indices, predIdx)
            if (fst > OrderedSet.min(xs.indices)) fst = OrderedSet.getAt(e.indices, predIdx - 1) + 1 as StructureElement.UnitIndex;
        } else {
            fst = OrderedSet.getAt(e.indices, predIdx - 1) + 1 as StructureElement.UnitIndex;
        }

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

    private onRemove(ref: string) {
        if (this.entries.has(ref)) this.entries.delete(ref);
    }

    private onUpdate(ref: string, oldObj: StateObject | undefined, obj: StateObject) {
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
        }
    }

    constructor(private plugin: PluginContext) {
        plugin.state.dataState.events.object.removed.subscribe(e => this.onRemove(e.ref));
        plugin.state.dataState.events.object.updated.subscribe(e => this.onUpdate(e.ref, e.oldObj, e.obj));
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