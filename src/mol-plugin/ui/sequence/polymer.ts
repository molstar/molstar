/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureSelection, StructureQuery, Structure, Queries, StructureProperties as SP, StructureElement, Unit } from '../../../mol-model/structure';
import { SequenceWrapper } from './util';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';

export type StructureUnit = { structure: Structure, unit: Unit }

export class PolymerSequenceWrapper extends SequenceWrapper<StructureUnit> {
    private readonly location = StructureElement.create()

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

    constructor(readonly data: StructureUnit) {
        super()

        const l = this.location
        l.unit = data.unit
        l.element = data.unit.elements[0]

        this.sequence = data.unit.model.sequence.byEntityKey[SP.entity.key(l)].sequence
        this.markerArray = new Uint8Array(this.sequence.sequence.length)
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