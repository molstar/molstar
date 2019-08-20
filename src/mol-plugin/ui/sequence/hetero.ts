/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement, ResidueIndex } from '../../../mol-model/structure';
import { SequenceWrapper, StructureUnit } from './wrapper';
import { OrderedSet, Segmentation, Interval, SortedArray } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { ColorNames } from '../../../mol-util/color/names';

export class HeteroSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private readonly sequence: string[]
    private readonly sequenceIndices: Map<ResidueIndex, number>
    private readonly residueIndices: Map<number, ResidueIndex>

    residueLabel(seqIdx: number) {
        return this.sequence[seqIdx]
    }
    residueColor(seqIdx: number) {
        return ColorNames.black
    }

    eachResidue(loci: Loci, apply: (set: OrderedSet) => boolean) {
        let changed = false
        const { structure, unit } = this.data
        if (StructureElement.isLoci(loci)) {
            if (!Structure.areParentsEquivalent(loci.structure, structure)) return false
            loci = StructureElement.Loci.remap(loci, structure)

            for (const e of loci.elements) {
                if (e.unit.id === unit.id) {
                    const { index: residueIndex } = e.unit.model.atomicHierarchy.residueAtomSegments
                    OrderedSet.forEach(e.indices, v => {
                        const seqIdx = this.sequenceIndices.get(residueIndex[unit.elements[v]])
                        if (seqIdx !== undefined && apply(Interval.ofSingleton(seqIdx))) changed = true
                    })
                }
            }
        } else if (Structure.isLoci(loci)) {
            if (!Structure.areParentsEquivalent(loci.structure, structure)) return false

            if (apply(Interval.ofBounds(0, this.length))) changed = true
        }
        return changed
    }

    getLoci(seqIdx: number) {
        const elements: StructureElement.Loci['elements'][0][] = []
        const rI = this.residueIndices.get(seqIdx)
        if (rI !== undefined) {
            const { unit } = this.data
            const { offsets } = unit.model.atomicHierarchy.residueAtomSegments
            const start = SortedArray.findPredecessorIndex(unit.elements, offsets[rI])
            const end = SortedArray.findPredecessorIndex(unit.elements, offsets[rI + 1])
            elements.push({ unit, indices: Interval.ofBounds(start, end) })
        }
        return StructureElement.Loci(this.data.structure, elements)
    }

    constructor(data: StructureUnit) {
        const sequence: string[] = []
        const sequenceIndices = new Map<ResidueIndex, number>()
        const residueIndices = new Map<number, ResidueIndex>()

        const residueIt = Segmentation.transientSegments(data.unit.model.atomicHierarchy.residueAtomSegments, data.unit.elements)
        while (residueIt.hasNext) {
            const { index } = residueIt.move()
            sequenceIndices.set(index, sequence.length)
            residueIndices.set(sequence.length, index)
            sequence.push(data.unit.model.atomicHierarchy.residues.label_comp_id.value(index))
        }

        const length = sequence.length
        const markerArray = new Uint8Array(length)

        super(data, markerArray, length)

        this.sequence = sequence
        this.sequenceIndices = sequenceIndices
        this.residueIndices = residueIndices
    }
}