/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginComponent } from '../../component';
import { PluginContext } from '../../../mol-plugin/context';
import { StructureElement, Structure } from '../../../mol-model/structure';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Boundary } from '../../../mol-model/structure/structure/util/boundary';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { structureElementStatsLabel } from '../../../mol-theme/label';
import { OrderedSet } from '../../../mol-data/int';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';
import { arrayRemoveAtInPlace } from '../../../mol-util/array';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { StateObject, StateSelection } from '../../../mol-state';
import { PluginStateObject } from '../../objects';
import { StructureSelectionQuery } from '../../../mol-plugin/util/structure-selection-query';
import { Task } from '../../../mol-task';

interface StructureSelectionManagerState {
    entries: Map<string, SelectionEntry>,
    history: HistoryEntry[],
    stats?: SelectionStats
}

const boundaryHelper = new BoundaryHelper();
const HISTORY_CAPACITY = 8;

export type StructureSelectionModifier = 'add' | 'remove' | 'set'

export class StructureSelectionManager extends PluginComponent<StructureSelectionManagerState> {    
    readonly events = {
        changed: this.ev<undefined>()
    }

    private referenceLoci: Loci | undefined

    get entries() { return this.state.entries; }
    get history() { return this.state.history; }
    get stats() {
        if (this.state.stats) return this.state.stats;
        this.state.stats = this.calcStats();
        return this.state.stats;
    }

    private getEntry(s: Structure) {
        const cell = this.plugin.helpers.substructureParent.get(s);
        if (!cell) return;
        const ref = cell.transform.ref;
        if (!this.entries.has(ref)) {
            const entry = new SelectionEntry(StructureElement.Loci(s, []));
            this.entries.set(ref, entry);
            return entry;
        }

        return this.entries.get(ref)!;
    }

    private calcStats(): SelectionStats {
        let structureCount = 0
        let elementCount = 0
        const stats = StructureElement.Stats.create()

        this.entries.forEach(v => {
            const { elements } = v.selection
            if (elements.length) {
                structureCount += 1
                for (let i = 0, il = elements.length; i < il; ++i) {
                    elementCount += OrderedSet.size(elements[i].indices)
                }
                StructureElement.Stats.add(stats, stats, StructureElement.Stats.ofLoci(v.selection))
            }
        })

        const label = structureElementStatsLabel(stats, { countsOnly: true })

        return { structureCount, elementCount, label }
    }

    private add(loci: Loci): boolean {
        if (!StructureElement.Loci.is(loci)) return false;
        
        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = StructureElement.Loci.union(entry.selection, loci);
        this.addHistory(loci);
        this.referenceLoci = loci
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    private remove(loci: Loci) {
        if (!StructureElement.Loci.is(loci)) return false;

        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = StructureElement.Loci.subtract(entry.selection, loci);
        this.removeHistory(loci);
        this.referenceLoci = loci
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    private set(loci: Loci) {
        if (!StructureElement.Loci.is(loci)) return false;

        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = loci;
        this.referenceLoci = undefined;
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    private addHistory(loci: StructureElement.Loci) {
        if (Loci.isEmpty(loci)) return;

        let idx = 0, entry: HistoryEntry | undefined = void 0;
        for (const l of this.history) {
            if (Loci.areEqual(l.loci, loci)) {
                entry = l;
                break;
            }
            idx++;
        }

        if (entry) {
            arrayRemoveAtInPlace(this.history, idx);
            this.history.unshift(entry);
            return;
        }

        const stats = StructureElement.Stats.ofLoci(loci);
        const label = structureElementStatsLabel(stats)

        this.history.unshift({ loci, label });
        if (this.history.length > HISTORY_CAPACITY) this.history.pop();
    }

    private removeHistory(loci: Loci) {
        if (Loci.isEmpty(loci)) return;

        let idx = 0, found = false;
        for (const l of this.history) {
            if (Loci.areEqual(l.loci, loci)) {
                found = true;
                break;
            }
            idx++;
        }

        if (found) {
            arrayRemoveAtInPlace(this.history, idx);
        }
    }

    private onRemove(ref: string) {
        if (this.entries.has(ref)) {
            this.entries.delete(ref);
            // TODO: property update the latest loci
            this.state.history = [];
            this.referenceLoci = undefined
        }
    }

    private onUpdate(ref: string, oldObj: StateObject | undefined, obj: StateObject) {
        if (!PluginStateObject.Molecule.Structure.is(obj)) return;

        if (this.entries.has(ref)) {
            if (!PluginStateObject.Molecule.Structure.is(oldObj) || oldObj === obj || oldObj.data === obj.data) return;

            // TODO: property update the latest loci & reference loci
            this.state.history = [];
            this.referenceLoci = undefined

            // remap the old selection to be related to the new object if possible.
            if (Structure.areUnitAndIndicesEqual(oldObj.data, obj.data)) {
                this.entries.set(ref, remapSelectionEntry(this.entries.get(ref)!, obj.data));
                return;
            }

            // clear the selection
            this.entries.set(ref, new SelectionEntry(StructureElement.Loci(obj.data, [])));
        }
    }

    /** Removes all selections and returns them */
    clear() {
        const keys = this.entries.keys();
        const selections: StructureElement.Loci[] = [];
        while (true) {
            const k = keys.next();
            if (k.done) break;
            const s = this.entries.get(k.value)!;
            if (!StructureElement.Loci.isEmpty(s.selection)) selections.push(s.selection);
            s.selection = StructureElement.Loci(s.selection.structure, []);
        }
        this.referenceLoci = undefined
        this.state.stats = void 0;
        this.events.changed.next()
        return selections;
    }

    getLoci(structure: Structure) {
        const entry = this.getEntry(structure);
        if (!entry) return EmptyLoci;
        return entry.selection;
    }

    getStructure(structure: Structure) {
        const entry = this.getEntry(structure);
        if (!entry) return;
        return entry.structure;
    }

    has(loci: Loci) {
        if (StructureElement.Loci.is(loci)) {
            const entry = this.getEntry(loci.structure);
            if (entry) {
                return StructureElement.Loci.isSubset(entry.selection, loci);
            }
        }
        return false;
    }

    tryGetRange(loci: Loci): StructureElement.Loci | undefined {
        if (!StructureElement.Loci.is(loci)) return;
        if (loci.elements.length !== 1) return;
        const entry = this.getEntry(loci.structure);
        if (!entry) return;

        const xs = loci.elements[0];
        if (!xs) return;

        const ref = this.referenceLoci
        if (!ref || !StructureElement.Loci.is(ref) || ref.structure.root !== loci.structure.root) return;

        let e: StructureElement.Loci['elements'][0] | undefined;
        for (const _e of ref.elements) {
            if (xs.unit === _e.unit) {
                e = _e;
                break;
            }
        }
        if (!e) return;

        if (xs.unit !== e.unit) return;

        return getElementRange(loci.structure.root, e, xs)
    }

    private prevHighlight: StructureElement.Loci | undefined = void 0;

    accumulateInteractiveHighlight(loci: Loci) {
        if (StructureElement.Loci.is(loci)) {
            if (this.prevHighlight) {
                this.prevHighlight = StructureElement.Loci.union(this.prevHighlight, loci);
            } else {
                this.prevHighlight = loci;
            }
        }
        return this.prevHighlight;
    }

    clearInteractiveHighlight() {
        const ret = this.prevHighlight;
        this.prevHighlight = void 0;
        return ret || EmptyLoci;
    }

    /** Count of all selected elements */
    elementCount() {
        let count = 0
        this.entries.forEach(v => {
            count += StructureElement.Loci.size(v.selection)
        })
        return count
    }

    getBoundary() {
        const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE)
        const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE)

        boundaryHelper.reset(0);

        const boundaries: Boundary[] = []
        this.entries.forEach(v => {
            const loci = v.selection
            if (!StructureElement.Loci.isEmpty(loci)) {
                boundaries.push(StructureElement.Loci.getBoundary(loci))
            }
        })

        for (let i = 0, il = boundaries.length; i < il; ++i) {
            const { box, sphere } = boundaries[i];
            Vec3.min(min, min, box.min);
            Vec3.max(max, max, box.max);
            boundaryHelper.boundaryStep(sphere.center, sphere.radius)
        }

        boundaryHelper.finishBoundaryStep();

        for (let i = 0, il = boundaries.length; i < il; ++i) {
            const { sphere } = boundaries[i];
            boundaryHelper.extendStep(sphere.center, sphere.radius);
        }

        return { box: { min, max }, sphere: boundaryHelper.getSphere() };
    }

    getPrincipalAxes(): PrincipalAxes {
        const elementCount = this.elementCount()
        const positions = new Float32Array(3 * elementCount)
        let offset = 0
        this.entries.forEach(v => {
            StructureElement.Loci.toPositionsArray(v.selection, positions, offset)
            offset += StructureElement.Loci.size(v.selection) * 3
        })
        return PrincipalAxes.ofPositions(positions)
    }

    modify(modifier: StructureSelectionModifier, loci: Loci) {
        let changed = false;
        switch (modifier) {
            case 'add': changed = this.add(loci); break;
            case 'remove': changed = this.remove(loci); break;
            case 'set': changed = this.set(loci); break;
        }

        if (changed) {
            this.state.stats = void 0;
            this.events.changed.next();
        }
    }

    private get applicableStructures() {
        // TODO: use "current structures" once implemented
        return this.plugin.state.dataState.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure)).map(s => s.obj!.data)
    }

    private triggerInteraction(modifier: StructureSelectionModifier, loci: Loci, applyGranularity = true) {
        switch (modifier) {
            case 'add':
                this.plugin.interactivity.lociSelects.select({ loci }, applyGranularity)
                break
            case 'remove':
                this.plugin.interactivity.lociSelects.deselect({ loci }, applyGranularity)
                break
            case 'set':
                this.plugin.interactivity.lociSelects.selectOnly({ loci }, applyGranularity)
                break
        }
    }

    fromSelectionQuery(modifier: StructureSelectionModifier, selectionQuery: StructureSelectionQuery, applyGranularity = true) {
        this.plugin.runTask(Task.create('Structure Selection', async runtime => {
            // const loci: Loci[] = [];
            for (const s of this.applicableStructures) {
                const loci = await StructureSelectionQuery.getLoci(this.plugin, runtime, selectionQuery, s);
                this.triggerInteraction(modifier, loci, applyGranularity);
            }
        }))
    }

    constructor(private plugin: PluginContext) {
        super({ entries: new Map(), history: [], stats: SelectionStats() });

        plugin.state.dataState.events.object.removed.subscribe(e => this.onRemove(e.ref));
        plugin.state.dataState.events.object.updated.subscribe(e => this.onUpdate(e.ref, e.oldObj, e.obj));
    }
}

interface SelectionStats {
    structureCount: number,
    elementCount: number,
    label: string
}

function SelectionStats(): SelectionStats { return { structureCount: 0, elementCount: 0, label: 'Nothing Selected' } };

class SelectionEntry {
    private _selection: StructureElement.Loci;
    private _structure?: Structure = void 0;

    get selection() { return this._selection; }
    set selection(value: StructureElement.Loci) {
        this._selection = value;
        this._structure = void 0
    }

    get structure(): Structure | undefined {
        if (this._structure) return this._structure;
        if (Loci.isEmpty(this._selection)) {
            this._structure = void 0;
        } else {
            this._structure = StructureElement.Loci.toStructure(this._selection);
        }
        return this._structure;
    }

    constructor(selection: StructureElement.Loci) {
        this._selection = selection;
    }
}

interface HistoryEntry {
    loci: StructureElement.Loci,
    label: string
}

/** remap `selection-entry` to be related to `structure` if possible */
function remapSelectionEntry(e: SelectionEntry, s: Structure): SelectionEntry {
    return new SelectionEntry(StructureElement.Loci.remap(e.selection, s));
}

/**
 * Assumes `ref` and `ext` belong to the same unit in the same structure
 */
function getElementRange(structure: Structure, ref: StructureElement.Loci['elements'][0], ext: StructureElement.Loci['elements'][0]) {
    const min = Math.min(OrderedSet.min(ref.indices), OrderedSet.min(ext.indices))
    const max = Math.max(OrderedSet.max(ref.indices), OrderedSet.max(ext.indices))

    return StructureElement.Loci(structure, [{
        unit: ref.unit,
        indices: OrderedSet.ofRange(min as StructureElement.UnitIndex, max as StructureElement.UnitIndex)
    }]);
}