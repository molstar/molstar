/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement, Structure } from 'mol-model/structure';
import { PluginContext } from '../../context';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { State } from 'mol-state';

export { StructureElementSelectionManager }

class StructureElementSelectionManager {
    private entries = new Map<string, Entry>();

    // maps structure to a parent StateObjectCell
    private structures = {
        root: new Map<Structure, string>(),
        tracked: new Map<string, Structure>()
    }

    private getEntry(s: Structure) {
        if (!this.structures.root.has(s)) return;
        const key = this.structures.root.get(s)!;
        if (!this.entries.has(key)) {
            const entry: Entry = {
                selection: StructureElement.Loci(s, []),
                highlight: StructureElement.Loci(s, []),
            };
            this.entries.set(key, entry);
            return entry;
        }

        return this.entries.get(key)!;
    }

    add(loci: StructureElement.Loci | Structure.Loci, type: 'selection' | 'highlight'): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        entry[type] = Structure.isLoci(loci) ? StructureElement.Loci.all(loci.structure) : StructureElement.Loci.union(entry[type], loci);
        return entry[type];
    }

    remove(loci: StructureElement.Loci | Structure.Loci, type: 'selection' | 'highlight'): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        entry[type] = Structure.isLoci(loci) ? StructureElement.Loci(loci.structure, []) : StructureElement.Loci.subtract(entry[type], loci);
        return entry[type].elements.length === 0 ? EmptyLoci : entry[type];
    }

    set(loci: StructureElement.Loci | Structure.Loci, type: 'selection' | 'highlight'): Loci {
        const entry = this.getEntry(loci.structure);
        if (!entry) return EmptyLoci;
        entry[type] = Structure.isLoci(loci) ? StructureElement.Loci.all(loci.structure) : loci;
        return entry[type].elements.length === 0 ? EmptyLoci : entry[type];
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

    private trackCell(state: State, ref: string) {
        const cell = state.cells.get(ref);
        if (!cell || !cell.obj || !PluginStateObject.Molecule.Structure.is(cell.obj)) return;
        const parent = state.selectQ(q => q.byRef(cell.transform.ref).rootOfType([PluginStateObject.Molecule.Structure]))[0];
        this.structures.tracked.set(cell.transform.ref, cell.obj.data);
        if (!parent) {
            this.structures.root.set(cell.obj.data, cell.transform.ref);
        } else {
            this.structures.root.set(cell.obj.data, parent.transform.ref);
        }
    }

    private untrackCell(state: State, ref: string) {
        if (!this.structures.tracked.has(ref)) return;
        const s = this.structures.tracked.get(ref)!;
        this.structures.tracked.delete(ref);
        this.structures.root.delete(s);
    }

    constructor(plugin: PluginContext) {
        plugin.state.dataState.events.object.created.subscribe(e => {
            this.trackCell(e.state, e.ref);
        });

        plugin.state.dataState.events.object.removed.subscribe(e => {
            this.untrackCell(e.state, e.ref);
        });

        plugin.state.dataState.events.object.updated.subscribe(e => {
            this.untrackCell(e.state, e.ref);
            this.trackCell(e.state, e.ref);
        });
    }
}

interface Entry {
    selection: StructureElement.Loci,
    highlight: StructureElement.Loci
}
