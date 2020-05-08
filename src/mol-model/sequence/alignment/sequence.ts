/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit } from '../../structure/structure';
import { AlignmentOptions, align } from './alignment';
import { OrderedSet } from '../../../mol-data/int';

export { AlignSequences };

namespace AlignSequences {
    export type Input = {
        a: StructureElement.Loci.Element,
        b: StructureElement.Loci.Element
    }
    /** `a` and `b` contain matching pairs, i.e. `a.indices[0]` aligns with `b.indices[0]` */
    export type Result = {
        a: StructureElement.Loci.Element,
        b: StructureElement.Loci.Element,
        score: number
    }

    function createSeqIdIndicesMap(element: StructureElement.Loci.Element) {
        const seqIds = new Map<number, StructureElement.UnitIndex[]>();
        if (Unit.isAtomic(element.unit)) {
            const { label_seq_id } = element.unit.model.atomicHierarchy.residues;
            const { residueIndex } = element.unit;
            for (let i = 0, il = OrderedSet.size(element.indices); i < il; ++i) {
                const uI = OrderedSet.getAt(element.indices, i);
                const eI = element.unit.elements[uI];
                const seqId = label_seq_id.value(residueIndex[eI]);
                if (seqIds.has(seqId)) seqIds.get(seqId)!.push(uI);
                else seqIds.set(seqId, [uI]);
            }
        } else if (Unit.isCoarse(element.unit)) {
            const { seq_id_begin } = Unit.isSpheres(element.unit)
                ? element.unit.model.coarseHierarchy.spheres
                : element.unit.model.coarseHierarchy.gaussians;
            for (let i = 0, il = OrderedSet.size(element.indices); i < il; ++i) {
                const uI = OrderedSet.getAt(element.indices, i);
                const eI = element.unit.elements[uI];
                const seqId = seq_id_begin.value(eI);
                seqIds.set(seqId, [uI]);
            }
        }
        return seqIds;
    }

    export function compute(input: Input, options: Partial<AlignmentOptions> = {}): Result {
        const seqA = getSequence(input.a.unit);
        const seqB = getSequence(input.b.unit);

        const seqIdIndicesA = createSeqIdIndicesMap(input.a);
        const seqIdIndicesB = createSeqIdIndicesMap(input.b);

        const indicesA: StructureElement.UnitIndex[] = [];
        const indicesB: StructureElement.UnitIndex[] = [];
        const { aliA, aliB, score } = align(seqA.code.toArray(), seqB.code.toArray(), options);

        let seqIdxA = 0, seqIdxB = 0;
        for (let i = 0, il = aliA.length; i < il; ++i) {
            if (aliA[i] === '-' || aliB[i] === '-') {
                if (aliA[i] !== '-') seqIdxA += 1;
                if (aliB[i] !== '-') seqIdxB += 1;
                continue;
            }

            const seqIdA = seqA.seqId.value(seqIdxA);
            const seqIdB = seqB.seqId.value(seqIdxB);

            if (seqIdIndicesA.has(seqIdA) && seqIdIndicesB.has(seqIdB)) {
                const iA = seqIdIndicesA.get(seqIdA)!;
                const iB = seqIdIndicesB.get(seqIdB)!;
                // use min length to guard against alternate locations
                for (let j = 0, jl = Math.min(iA.length, iB.length); j < jl; ++j) {
                    indicesA.push(iA[j]);
                    indicesB.push(iB[j]);
                }
            }

            seqIdxA += 1, seqIdxB += 1;
        }

        const outA = OrderedSet.intersect(OrderedSet.ofSortedArray(indicesA), input.a.indices);
        const outB = OrderedSet.intersect(OrderedSet.ofSortedArray(indicesB), input.b.indices);

        return {
            a: { unit: input.a.unit, indices: outA },
            b: { unit: input.b.unit, indices: outB },
            score
        };
    }
}

function entityKey(unit: Unit) {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return unit.model.atomicHierarchy.index.getEntityFromChain(unit.chainIndex[unit.elements[0]]);
        case Unit.Kind.Spheres:
            return unit.model.coarseHierarchy.spheres.entityKey[unit.elements[0]];
        case Unit.Kind.Gaussians:
            return unit.model.coarseHierarchy.gaussians.entityKey[unit.elements[0]];
    }
}

function getSequence(unit: Unit) {
    return unit.model.sequence.byEntityKey[entityKey(unit)].sequence;
}