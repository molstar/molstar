/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { SequenceWrapper, StructureUnit } from './wrapper';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { ColorNames } from '../../../mol-util/color/names';

export class ChainSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private label: string
    private indices: Interval<StructureElement.UnitIndex>
    private loci: StructureElement.Loci

    residueLabel(seqIdx: number) {
        return this.label
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
                        if (apply(Interval.ofSingleton(0))) changed = true
                    }
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areRootsEquivalent(loci.structure, structure)) return false

            if (apply(Interval.ofSingleton(0))) changed = true
        }
        return changed
    }

    getLoci(seqIdx: number) {
        return this.loci
    }

    constructor(data: StructureUnit) {
        const counts: string[] = []
        const l = StructureElement.Location.create(data.unit, data.unit.elements[0])
        const entitySeq = data.unit.model.sequence.byEntityKey[StructureProperties.entity.key(l)]
        if (entitySeq) counts.push(`${entitySeq.sequence.length} residues`)
        counts.push(`${data.unit.elements.length} elements`)

        const length = 1
        const markerArray = new Uint8Array(length)

        super(data, markerArray, length)

        this.label = `Whole Unit (${counts.join(', ')})`
        this.indices = Interval.ofBounds(0, data.unit.elements.length)
        this.loci = StructureElement.Loci(this.data.structure, [{
            unit: this.data.unit, indices: this.indices
        }])
    }
}