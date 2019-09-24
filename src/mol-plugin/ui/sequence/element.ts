/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement } from '../../../mol-model/structure';
import { SequenceWrapper, StructureUnit } from './wrapper';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { ColorNames } from '../../../mol-util/color/names';

export class ElementSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private indices: Interval<StructureElement.UnitIndex>

    residueLabel(seqIdx: number) {
        return 'X'
    }
    residueColor(seqIdx: number) {
        return ColorNames.black
    }

    eachResidue(loci: Loci, apply: (set: OrderedSet) => boolean) {
        let changed = false
        const { structure, unit } = this.data
        if (StructureElement.Loci.is(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false
            loci = StructureElement.Loci.remap(loci, structure)

            for (const e of loci.elements) {
                if (e.unit.id === unit.id) {
                    if (OrderedSet.isSubset(this.indices, e.indices)) {
                        if (apply(e.indices)) changed = true
                    }
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false

            if (apply(this.indices)) changed = true
        }
        return changed
    }

    getLoci(seqIdx: number) {
        const elements: StructureElement.Loci['elements'][0][] = [{
            unit: this.data.unit,
            indices: Interval.ofSingleton(seqIdx)
        }]
        return StructureElement.Loci(this.data.structure, elements)
    }

    constructor(data: StructureUnit) {
        const length = data.unit.elements.length
        const markerArray = new Uint8Array(length)

        super(data, markerArray, length)

        this.indices = Interval.ofBounds(0, length)
    }
}