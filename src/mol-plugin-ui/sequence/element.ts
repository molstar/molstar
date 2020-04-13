/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { SequenceWrapper, StructureUnit } from './wrapper';
import { OrderedSet, Interval } from '../../mol-data/int';
import { Loci } from '../../mol-model/loci';
import { ColorNames } from '../../mol-util/color/names';
import { MarkerAction, applyMarkerAction } from '../../mol-util/marker-action';

export class ElementSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private unitIndices: Map<number, Interval<StructureElement.UnitIndex>>

    residueLabel(seqIdx: number) {
        return 'X';
    }
    residueColor(seqIdx: number) {
        return ColorNames.black;
    }
    residueClass(seqIdx: number) {
        return 'msp-sequence-present';
    }

    mark(loci: Loci, action: MarkerAction) {
        let changed = false;
        const { structure, units } = this.data;
        if (StructureElement.Loci.is(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false;
            loci = StructureElement.Loci.remap(loci, structure);

            for (const e of loci.elements) {
                const indices = this.unitIndices.get(e.unit.id);
                if (indices) {
                    if (OrderedSet.isSubset(indices, e.indices)) {
                        if (applyMarkerAction(this.markerArray, e.indices, action)) changed = true;
                    }
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false;

            for (let i = 0, il = units.length; i < il; ++i) {
                const indices = this.unitIndices.get(units[i].id)!;
                if (applyMarkerAction(this.markerArray, indices, action)) changed = true;
            }
        }
        return changed;
    }

    getLoci(seqIdx: number) {
        const { units } = this.data;
        const lociElements: StructureElement.Loci['elements'][0][] = [];
        let offset = 0;
        for (let i = 0, il = units.length; i < il; ++i) {
            const unit = units[i];
            if (seqIdx < offset + unit.elements.length) {
                lociElements.push({ unit, indices: Interval.ofSingleton(seqIdx - offset) });
                break;
            }
            offset += unit.elements.length;
        }
        return StructureElement.Loci(this.data.structure, lociElements);
    }

    constructor(data: StructureUnit) {
        let length = 0;

        const unitIndices = new Map<number, Interval<StructureElement.UnitIndex>>();
        const lociElements: StructureElement.Loci['elements'][0][] = [];

        for (let i = 0, il = data.units.length; i < il; ++i) {
            const unit = data.units[i];
            length += unit.elements.length;

            const indices = Interval.ofBounds(0, unit.elements.length);
            unitIndices.set(unit.id, indices);
            lociElements.push({ unit, indices });
        }
        const markerArray = new Uint8Array(length);

        super(data, markerArray, length);

        this.unitIndices = unitIndices;
    }
}