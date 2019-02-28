/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement, Structure } from 'mol-model/structure';
import { PluginContext } from '../../context';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { State, StateSelection, StateObject } from 'mol-state';
import { mapObjectMap } from 'mol-util/object';

export { StructureLociManager }

class StructureLociManager {
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

    add(loci: StructureElement.Loci | Structure.Loci, type: 'selection' | 'highlight'): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        const xs = entry.elements;
        xs[type] = Structure.isLoci(loci) ? StructureElement.Loci.all(loci.structure) : StructureElement.Loci.union(xs[type], loci);
        return xs[type];
    }

    remove(loci: StructureElement.Loci | Structure.Loci, type: 'selection' | 'highlight'): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        const xs = entry.elements;
        xs[type] = Structure.isLoci(loci) ? StructureElement.Loci(loci.structure, []) : StructureElement.Loci.subtract(xs[type], loci);
        return xs[type].elements.length === 0 ? EmptyLoci : xs[type];
    }

    set(loci: StructureElement.Loci | Structure.Loci, type: 'selection' | 'highlight'): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        const xs = entry.elements;
        xs[type] = Structure.isLoci(loci) ? StructureElement.Loci.all(loci.structure) : loci;
        return xs[type].elements.length === 0 ? EmptyLoci : xs[type];
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
    elements: { [category: string]: StructureElement.Loci }
}

function SelectionEntry(s: Structure): SelectionEntry {
    return {
        elements: { }
    };
}

function remapSelectionEntry(e: SelectionEntry, s: Structure): SelectionEntry {
    return {
        elements: mapObjectMap(e.elements, (l: StructureElement.Loci) => StructureElement.Loci.remap(l, s))
    };
}