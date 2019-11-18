/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../mol-data/int';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Structure, StructureElement } from '../../mol-model/structure';
import { StateObject } from '../../mol-state';
import { PluginContext } from '../context';
import { PluginStateObject } from '../state/objects';
import { structureElementStatsLabel } from '../../mol-theme/label';
import { Vec3 } from '../../mol-math/linear-algebra';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { Boundary } from '../../mol-model/structure/structure/util/boundary';
import { PrincipalAxes } from '../../mol-math/linear-algebra/matrix/principal-axes';
import Matrix from '../../mol-math/linear-algebra/matrix/matrix';

const boundaryHelper = new BoundaryHelper();

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

    /** Count of all selected elements */
    size() {
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
        const elementCount = this.size()
        const positions = new Float32Array(3 * elementCount)
        let offset = 0
        this.entries.forEach(v => {
            StructureElement.Loci.toPositionsArray(v.selection, positions, offset)
            offset += StructureElement.Loci.size(v.selection) * 3
        })
        return PrincipalAxes.ofPoints(Matrix.fromArray(positions, 3, elementCount))
    }

    get stats() {
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

        const label = structureElementStatsLabel(stats, true)

        return { structureCount, elementCount, label }
    }

    add(loci: Loci): Loci {
        if (StructureElement.Loci.is(loci)) {
            const entry = this.getEntry(loci.structure);
            if (entry) {
                entry.selection = StructureElement.Loci.union(entry.selection, loci);
                this.plugin.events.interactivity.selectionUpdated.next()
                return entry.selection;
            }
        }
        return EmptyLoci
    }

    remove(loci: Loci): Loci {
        if (StructureElement.Loci.is(loci)) {
            const entry = this.getEntry(loci.structure);
            if (entry) {
                entry.selection = StructureElement.Loci.subtract(entry.selection, loci);
                this.plugin.events.interactivity.selectionUpdated.next()
                return StructureElement.Loci.isEmpty(entry.selection) ? EmptyLoci : entry.selection;
            }
        }
        return EmptyLoci
    }

    set(loci: Loci): Loci {
        if (StructureElement.Loci.is(loci)) {
            const entry = this.getEntry(loci.structure);
            if (entry) {
                entry.selection = loci;
                this.plugin.events.interactivity.selectionUpdated.next()
                return StructureElement.Loci.isEmpty(entry.selection) ? EmptyLoci : entry.selection;
            }
        }
        return EmptyLoci;
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
        this.plugin.events.interactivity.selectionUpdated.next()
        return selections;
    }

    get(structure: Structure) {
        const entry = this.getEntry(structure);
        if (!entry) return EmptyLoci;
        return entry.selection;
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

        let xs = loci.elements[0];
        if (!xs) return;
        let e: StructureElement.Loci['elements'][0] | undefined;
        for (const _e of entry.selection.elements) {
            if (xs.unit === _e.unit) {
                e = _e;
                break;
            }
        }
        if (!e) return;

        return tryGetElementRange(entry.selection.structure, e, xs)
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

/** remap `selection-entry` to be related to `structure` if possible */
function remapSelectionEntry(e: SelectionEntry, s: Structure): SelectionEntry {
    return {
        selection: StructureElement.Loci.remap(e.selection, s)
    };
}

/**
 * Assumes `ref` and `ext` belong to the same unit in the same structure
 */
function tryGetElementRange(structure: Structure, ref: StructureElement.Loci['elements'][0], ext: StructureElement.Loci['elements'][0]) {

    const refMin = OrderedSet.min(ref.indices)
    const refMax = OrderedSet.max(ref.indices)
    const extMin = OrderedSet.min(ext.indices)
    const extMax = OrderedSet.max(ext.indices)

    let min: number
    let max: number

    if (refMax < extMin) {
        min = refMax + 1
        max = extMax
    } else if (extMax < refMin) {
        min = extMin
        max = refMin - 1
    } else {
        // TODO handle range overlap cases
        return
    }

    return StructureElement.Loci(structure, [{
        unit: ref.unit,
        indices: OrderedSet.ofRange(min as StructureElement.UnitIndex, max as StructureElement.UnitIndex)
    }]);
}