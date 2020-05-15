/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../../mol-data/int';
import { BoundaryHelper } from '../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Structure, StructureElement, StructureSelection } from '../../../mol-model/structure';
import { Boundary } from '../../../mol-model/structure/structure/util/boundary';
import { PluginContext } from '../../../mol-plugin/context';
import { StateObjectRef } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { structureElementStatsLabel } from '../../../mol-theme/label';
import { arrayRemoveAtInPlace } from '../../../mol-util/array';
import { StatefulPluginComponent } from '../../component';
import { StructureSelectionQuery } from '../../helpers/structure-selection-query';
import { PluginStateObject as PSO } from '../../objects';
import { UUID } from '../../../mol-util';
import { StructureRef } from './hierarchy-state';

interface StructureSelectionManagerState {
    entries: Map<string, SelectionEntry>,
    additionsHistory: StructureSelectionHistoryEntry[],
    stats?: SelectionStats
}

const boundaryHelper = new BoundaryHelper('98');
const HISTORY_CAPACITY = 24;

export type StructureSelectionModifier = 'add' | 'remove' | 'intersect' | 'set'

export class StructureSelectionManager extends StatefulPluginComponent<StructureSelectionManagerState> {
    readonly events = {
        changed: this.ev<undefined>(),
        additionsHistoryUpdated: this.ev<undefined>(),

        loci: {
            add: this.ev<StructureElement.Loci>(),
            remove: this.ev<StructureElement.Loci>(),
            clear: this.ev<undefined>()
        }
    }

    private referenceLoci: StructureElement.Loci | undefined

    get entries() { return this.state.entries; }
    get additionsHistory() { return this.state.additionsHistory; }
    get stats() {
        if (this.state.stats) return this.state.stats;
        this.state.stats = this.calcStats();
        return this.state.stats;
    }

    private getEntry(s: Structure) {
        // ignore decorators to get stable ref
        const cell = this.plugin.helpers.substructureParent.get(s, true);
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
        let structureCount = 0;
        let elementCount = 0;
        const stats = StructureElement.Stats.create();

        this.entries.forEach(v => {
            const { elements } = v.selection;
            if (elements.length) {
                structureCount += 1;
                for (let i = 0, il = elements.length; i < il; ++i) {
                    elementCount += OrderedSet.size(elements[i].indices);
                }
                StructureElement.Stats.add(stats, stats, StructureElement.Stats.ofLoci(v.selection));
            }
        });

        const label = structureElementStatsLabel(stats, { countsOnly: true });

        return { structureCount, elementCount, label };
    }

    private add(loci: Loci): boolean {
        if (!StructureElement.Loci.is(loci)) return false;

        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = StructureElement.Loci.union(entry.selection, loci);
        this.tryAddHistory(loci);
        this.referenceLoci = loci;
        this.events.loci.add.next(loci);
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    private remove(loci: Loci) {
        if (!StructureElement.Loci.is(loci)) return false;

        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = StructureElement.Loci.subtract(entry.selection, loci);
        // this.addHistory(loci);
        this.referenceLoci = loci;
        this.events.loci.remove.next(loci);
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    private intersect(loci: Loci): boolean {
        if (!StructureElement.Loci.is(loci)) return false;

        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = StructureElement.Loci.intersect(entry.selection, loci);
        // this.addHistory(loci);
        this.referenceLoci = loci;
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    private set(loci: Loci) {
        if (!StructureElement.Loci.is(loci)) return false;

        const entry = this.getEntry(loci.structure);
        if (!entry) return false;

        const sel = entry.selection;
        entry.selection = loci;
        this.tryAddHistory(loci);
        this.referenceLoci = undefined;
        return !StructureElement.Loci.areEqual(sel, entry.selection);
    }

    modifyHistory(entry: StructureSelectionHistoryEntry, action: 'remove' | 'up' | 'down', modulus?: number, groupByStructure = false) {
        const history = this.additionsHistory;
        const idx = history.indexOf(entry);
        if (idx < 0) return;

        let swapWith: number | undefined = void 0;

        switch (action) {
            case 'remove': arrayRemoveAtInPlace(history, idx); break;
            case 'up': swapWith = idx - 1; break;
            case 'down': swapWith = idx + 1; break;
        }

        if (swapWith !== void 0) {
            const mod = modulus ? Math.min(history.length, modulus) : history.length;
            while (true) {
                swapWith = swapWith % mod;
                if (swapWith < 0) swapWith += mod;

                if (!groupByStructure || history[idx].loci.structure === history[swapWith].loci.structure) {
                    const t = history[idx];
                    history[idx] = history[swapWith];
                    history[swapWith] = t;
                    break;
                } else {
                    swapWith += action === 'up' ? -1 : +1;
                }
            }
        }

        this.events.additionsHistoryUpdated.next();
    }

    private tryAddHistory(loci: StructureElement.Loci) {
        if (Loci.isEmpty(loci)) return;

        let idx = 0, entry: StructureSelectionHistoryEntry | undefined = void 0;
        for (const l of this.additionsHistory) {
            if (Loci.areEqual(l.loci, loci)) {
                entry = l;
                break;
            }
            idx++;
        }

        if (entry) {
            // move to top
            arrayRemoveAtInPlace(this.additionsHistory, idx);
            this.additionsHistory.unshift(entry);
            this.events.additionsHistoryUpdated.next();
            return;
        }

        const stats = StructureElement.Stats.ofLoci(loci);
        const label = structureElementStatsLabel(stats, { reverse: true });

        this.additionsHistory.unshift({ id: UUID.create22(), loci, label });
        if (this.additionsHistory.length > HISTORY_CAPACITY) this.additionsHistory.pop();

        this.events.additionsHistoryUpdated.next();
    }

    private clearHistory() {
        if (this.state.additionsHistory.length !== 0) {
            this.state.additionsHistory = [];
            this.events.additionsHistoryUpdated.next();
        }
    }

    private clearHistoryForStructure(structure: Structure) {
        const historyEntryToRemove: StructureSelectionHistoryEntry[] = [];
        for (const e of this.state.additionsHistory) {
            if (e.loci.structure.root === structure.root) {
                historyEntryToRemove.push(e);
            }
        }
        for (const e of historyEntryToRemove) {
            this.modifyHistory(e, 'remove');
        }
        if (historyEntryToRemove.length !== 0) {
            this.events.additionsHistoryUpdated.next();
        }
    }

    private onRemove(ref: string, obj: PSO.Molecule.Structure | undefined) {
        if (this.entries.has(ref)) {
            this.entries.delete(ref);
            if (obj?.data) {
                this.clearHistoryForStructure(obj.data);
            }
            if (this.referenceLoci?.structure === obj?.data) {
                this.referenceLoci = undefined;
            }
            this.state.stats = void 0;
            this.events.changed.next();
        }
    }

    private onUpdate(ref: string, oldObj: PSO.Molecule.Structure | undefined, obj: PSO.Molecule.Structure) {

        // no change to structure
        if (oldObj === obj || oldObj?.data === obj.data) return;

        // ignore decorators to get stable ref
        const cell = this.plugin.helpers.substructureParent.get(obj.data, true);
        if (!cell) return;

        ref = cell.transform.ref;
        if (!this.entries.has(ref)) return;

        // use structure from last decorator as reference
        const structure = this.plugin.helpers.substructureParent.get(obj.data)?.obj?.data;
        if (!structure) return;

        // oldObj is not defined for inserts (e.g. TransformStructureConformation)
        if (!oldObj?.data || Structure.areUnitAndIndicesEqual(oldObj.data, obj.data)) {
            this.entries.set(ref, remapSelectionEntry(this.entries.get(ref)!, structure));

            // remap referenceLoci & prevHighlight if needed and possible
            if (this.referenceLoci?.structure.root === structure.root) {
                this.referenceLoci = StructureElement.Loci.remap(this.referenceLoci, structure);
            }

            // remap history locis if needed and possible
            let changedHistory = false;
            for (const e of this.state.additionsHistory) {
                if (e.loci.structure.root === structure.root) {
                    e.loci = StructureElement.Loci.remap(e.loci, structure);
                    changedHistory = true;
                }
            }
            if (changedHistory) this.events.additionsHistoryUpdated.next();
        } else {
            // clear the selection for ref
            this.entries.set(ref, new SelectionEntry(StructureElement.Loci(structure, [])));

            if (this.referenceLoci?.structure.root === structure.root) {
                this.referenceLoci = undefined;
            }

            this.clearHistoryForStructure(structure);

            this.state.stats = void 0;
            this.events.changed.next();
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
        this.referenceLoci = undefined;
        this.state.stats = void 0;
        this.events.changed.next();
        this.events.loci.clear.next();
        this.clearHistory();
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

    structureHasSelection(structure: StructureRef) {
        const s = structure.cell?.obj?.data;
        if (!s) return false;
        const entry = this.getEntry(s);
        return !!entry && !StructureElement.Loci.isEmpty(entry.selection);
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

        const ref = this.referenceLoci;
        if (!ref || !StructureElement.Loci.is(ref) || ref.structure !== loci.structure) return;

        let e: StructureElement.Loci['elements'][0] | undefined;
        for (const _e of ref.elements) {
            if (xs.unit === _e.unit) {
                e = _e;
                break;
            }
        }
        if (!e) return;

        if (xs.unit !== e.unit) return;

        return getElementRange(loci.structure, e, xs);
    }

    /** Count of all selected elements */
    elementCount() {
        let count = 0;
        this.entries.forEach(v => {
            count += StructureElement.Loci.size(v.selection);
        });
        return count;
    }

    getBoundary() {
        const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE);
        const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE);

        boundaryHelper.reset();

        const boundaries: Boundary[] = [];
        this.entries.forEach(v => {
            const loci = v.selection;
            if (!StructureElement.Loci.isEmpty(loci)) {
                boundaries.push(StructureElement.Loci.getBoundary(loci));
            }
        });

        for (let i = 0, il = boundaries.length; i < il; ++i) {
            const { box, sphere } = boundaries[i];
            Vec3.min(min, min, box.min);
            Vec3.max(max, max, box.max);
            boundaryHelper.includePositionRadius(sphere.center, sphere.radius);
        }
        boundaryHelper.finishedIncludeStep();
        for (let i = 0, il = boundaries.length; i < il; ++i) {
            const { sphere } = boundaries[i];
            boundaryHelper.radiusPositionRadius(sphere.center, sphere.radius);
        }

        return { box: { min, max }, sphere: boundaryHelper.getSphere() };
    }

    getPrincipalAxes(): PrincipalAxes {
        const elementCount = this.elementCount();
        const positions = new Float32Array(3 * elementCount);
        let offset = 0;
        this.entries.forEach(v => {
            StructureElement.Loci.toPositionsArray(v.selection, positions, offset);
            offset += StructureElement.Loci.size(v.selection) * 3;
        });
        return PrincipalAxes.ofPositions(positions);
    }

    modify(modifier: StructureSelectionModifier, loci: Loci) {
        let changed = false;
        switch (modifier) {
            case 'add': changed = this.add(loci); break;
            case 'remove': changed = this.remove(loci); break;
            case 'intersect': changed = this.intersect(loci); break;
            case 'set': changed = this.set(loci); break;
        }

        if (changed) {
            this.state.stats = void 0;
            this.events.changed.next();
        }
    }

    private get applicableStructures() {
        return this.plugin.managers.structure.hierarchy.selection.structures
            .filter(s => !!s.cell.obj)
            .map(s => s.cell.obj!.data);
    }

    private triggerInteraction(modifier: StructureSelectionModifier, loci: Loci, applyGranularity = true) {
        switch (modifier) {
            case 'add':
                this.plugin.managers.interactivity.lociSelects.select({ loci }, applyGranularity);
                break;
            case 'remove':
                this.plugin.managers.interactivity.lociSelects.deselect({ loci }, applyGranularity);
                break;
            case 'intersect':
                this.plugin.managers.interactivity.lociSelects.selectJoin({ loci }, applyGranularity);
                break;
            case 'set':
                this.plugin.managers.interactivity.lociSelects.selectOnly({ loci }, applyGranularity);
                break;
        }
    }

    fromLoci(modifier: StructureSelectionModifier, loci: Loci, applyGranularity = true) {
        this.triggerInteraction(modifier, loci, applyGranularity);
    }

    fromSelectionQuery(modifier: StructureSelectionModifier, query: StructureSelectionQuery, applyGranularity = true) {
        this.plugin.runTask(Task.create('Structure Selection', async runtime => {
            for (const s of this.applicableStructures) {
                const loci = await query.getSelection(this.plugin, runtime, s);
                this.triggerInteraction(modifier, StructureSelection.toLociWithSourceUnits(loci), applyGranularity);
            }
        }));
    }

    fromSelections(ref: StateObjectRef<PSO.Molecule.Structure.Selections>) {
        const cell = StateObjectRef.resolveAndCheck(this.plugin.state.data, ref);
        if (!cell || !cell.obj) return;

        if (!PSO.Molecule.Structure.Selections.is(cell.obj)) {
            console.warn('fromSelections applied to wrong object type.', cell.obj);
            return;
        }

        this.clear();
        for (const s of cell.obj?.data) {
            this.fromLoci('set', s.loci);
        }
    }

    constructor(private plugin: PluginContext) {
        super({ entries: new Map(), additionsHistory: [], stats: SelectionStats() });

        // listen to events from substructureParent helper to ensure it is updated
        plugin.helpers.substructureParent.events.removed.subscribe(e => this.onRemove(e.ref, e.obj));
        plugin.helpers.substructureParent.events.updated.subscribe(e => this.onUpdate(e.ref, e.oldObj, e.obj));
    }
}

interface SelectionStats {
    structureCount: number,
    elementCount: number,
    label: string
}

function SelectionStats(): SelectionStats { return { structureCount: 0, elementCount: 0, label: 'Nothing Selected' }; };

class SelectionEntry {
    private _selection: StructureElement.Loci;
    private _structure?: Structure = void 0;

    get selection() { return this._selection; }
    set selection(value: StructureElement.Loci) {
        this._selection = value;
        this._structure = void 0;
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

export interface StructureSelectionHistoryEntry {
    id: UUID,
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
    const min = Math.min(OrderedSet.min(ref.indices), OrderedSet.min(ext.indices));
    const max = Math.max(OrderedSet.max(ref.indices), OrderedSet.max(ext.indices));

    return StructureElement.Loci(structure, [{
        unit: ref.unit,
        indices: OrderedSet.ofRange(min as StructureElement.UnitIndex, max as StructureElement.UnitIndex)
    }]);
}