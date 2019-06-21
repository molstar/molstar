/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureSelection, StructureQuery, Structure, Queries, StructureProperties as SP, StructureElement, Unit } from '../../../mol-model/structure';
import { SequenceWrapper } from './util';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { Sequence } from '../../../mol-model/sequence';
import { MissingResidues } from '../../../mol-model/structure/model/properties/common';
import { ColorNames } from '../../../mol-util/color/tables';

export type StructureUnit = { structure: Structure, unit: Unit }

export class PolymerSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private readonly location: StructureElement
    private readonly sequence: Sequence
    private readonly missing: MissingResidues

    private readonly modelNum: number
    private readonly asymId: string

    seqId(i: number) {
        return this.sequence.offset + i + 1
    }
    residueLabel(i: number) {
        return this.sequence.sequence[i]
    }
    residueColor(i: number) {
        return this.missing.has(this.modelNum, this.asymId, this.seqId(i))
            ? ColorNames.grey
            : ColorNames.black
    }

    eachResidue(loci: Loci, apply: (interval: Interval) => boolean) {
        let changed = false
        const { structure, unit } = this.data
        if (!StructureElement.isLoci(loci)) return false
        if (!Structure.areParentsEquivalent(loci.structure, structure)) return false

        const { location } = this
        for (const e of loci.elements) {
            let rIprev = -1
            location.unit = e.unit

            const { index: residueIndex } = e.unit.model.atomicHierarchy.residueAtomSegments

            OrderedSet.forEach(e.indices, v => {
                location.element = e.unit.elements[v]
                const rI = residueIndex[location.element]
                // avoid checking for the same residue multiple times
                if (rI !== rIprev) {
                    if (SP.unit.id(location) !== unit.id) return

                    if (apply(getSeqIdInterval(location))) changed = true
                    rIprev = rI
                }
            })
        }
        return changed
    }

    getLoci(seqId: number) {
        const query = createResidueQuery(this.data.unit.id, seqId);
        return StructureSelection.toLoci2(StructureQuery.run(query, this.data.structure));
    }

    constructor(data: StructureUnit) {
        const l = StructureElement.create(data.unit, data.unit.elements[0])
        const sequence = data.unit.model.sequence.byEntityKey[SP.entity.key(l)].sequence
        const markerArray = new Uint8Array(sequence.sequence.length)

        super(data, markerArray, sequence.sequence.length)

        this.sequence = sequence
        this.location = StructureElement.create()
        this.missing = data.unit.model.properties.missingResidues

        this.modelNum = data.unit.model.modelNum
        this.asymId = SP.chain.label_asym_id(l)
    }
}

function createResidueQuery(unitId: number, label_seq_id: number) {
    return Queries.generators.atoms({
        unitTest: ctx => {
            return SP.unit.id(ctx.element) === unitId
        },
        residueTest: ctx => {
            if (ctx.element.unit.kind === Unit.Kind.Atomic) {
                return SP.residue.label_seq_id(ctx.element) === label_seq_id
            } else {
                return (
                    SP.coarse.seq_id_begin(ctx.element) <= label_seq_id &&
                    SP.coarse.seq_id_end(ctx.element) >= label_seq_id
                )
            }
        }
    });
}

/** Zero-indexed */
function getSeqIdInterval(location: StructureElement): Interval {
    const { unit, element } = location
    const { model } = unit
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const residueIndex = model.atomicHierarchy.residueAtomSegments.index[element]
            const seqId = model.atomicHierarchy.residues.label_seq_id.value(residueIndex)
            return Interval.ofSingleton(seqId - 1)
        case Unit.Kind.Spheres:
            return Interval.ofRange(
                model.coarseHierarchy.spheres.seq_id_begin.value(element) - 1,
                model.coarseHierarchy.spheres.seq_id_end.value(element) - 1
            )
        case Unit.Kind.Gaussians:
            return Interval.ofRange(
                model.coarseHierarchy.gaussians.seq_id_begin.value(element) - 1,
                model.coarseHierarchy.gaussians.seq_id_end.value(element) - 1
            )
    }
}