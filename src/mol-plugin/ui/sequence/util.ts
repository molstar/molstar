/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSequence, Queries, StructureProperties as SP, StructureElement, Unit } from '../../../mol-model/structure';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { Loci } from '../../../mol-model/loci';
import { applyMarkerAction, MarkerAction } from '../../../mol-util/marker-action';

export function createResidueQuery(entityId: string, label_seq_id: number) {
    return Queries.generators.atoms({
        entityTest: ctx => {
            return SP.entity.id(ctx.element) === entityId
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
export function getSeqIdInterval(location: StructureElement): Interval {
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

export type StructureSeq = { structure: Structure, seq: StructureSequence.Entity }

export function eachResidue(loci: Loci, structureSeq: StructureSeq, apply: (interval: Interval) => boolean) {
    let changed = false
    const { structure, seq } = structureSeq
    if (!StructureElement.isLoci(loci)) return false
    if (!Structure.areParentsEquivalent(loci.structure, structure)) return false
    const l = StructureElement.create()
    for (const e of loci.elements) {
        l.unit = e.unit
        OrderedSet.forEach(e.indices, v => {
            l.element = e.unit.elements[v]
            const entityId = SP.entity.id(l)
            if (entityId === seq.entityId) {
                if (apply(getSeqIdInterval(l))) changed = true
            }
        })
    }
    return changed
}

export function markResidue(loci: Loci, structureSeq: StructureSeq, array: Uint8Array, action: MarkerAction) {
    const { structure, seq } = structureSeq
    return eachResidue(loci, { structure , seq }, (i: Interval) => {
        return applyMarkerAction(array, i, action)
    })
}