/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ModelCrossLinkRestraint } from './format';
import { Unit, StructureElement, Structure, CustomPropertyDescriptor} from '../../../mol-model/structure';
import { PairRestraints, PairRestraint } from '../pair-restraints';
import { CustomStructureProperty } from '../../common/custom-structure-property';
import { CustomProperty } from '../../common/custom-property';

export type CrossLinkRestraintValue = PairRestraints<CrossLinkRestraint>

export const CrossLinkRestraintProvider: CustomStructureProperty.Provider<{}, CrossLinkRestraintValue> = CustomStructureProperty.createProvider({
    label: 'Cross Link Restraint',
    descriptor: CustomPropertyDescriptor({
        name: 'integrative-cross-link-restraint',
        // TODO `cifExport` and `symbol`
    }),
    type: 'local',
    defaultParams: {},
    getParams: (data: Structure) => ({}),
    isApplicable: (data: Structure) => data.models.some(m => !!ModelCrossLinkRestraint.Provider.get(m)),
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<{}>) => {
        return extractCrossLinkRestraints(data)
    }
})

export { CrossLinkRestraint }

interface CrossLinkRestraint extends PairRestraint {
    readonly restraintType: 'harmonic' | 'upper bound' | 'lower bound'
    readonly distanceThreshold: number
    readonly psi: number
    readonly sigma1: number
    readonly sigma2: number
}

namespace CrossLinkRestraint {
    export enum Tag {
        CrossLinkRestraint = 'cross-link-restraint'
    }

    export function isApplicable(structure: Structure) {
        return structure.models.some(m => !!ModelCrossLinkRestraint.Provider.get(m))
    }
}


//

function _addRestraints(map: Map<number, number>, unit: Unit, restraints: ModelCrossLinkRestraint) {
    const { elements } = unit;
    const elementCount = elements.length;
    const kind = unit.kind

    for (let i = 0; i < elementCount; i++) {
        const e = elements[i];
        restraints.getIndicesByElement(e, kind).forEach(ri => map.set(ri, i))
    }
}

function extractInter(pairs: CrossLinkRestraint[], unitA: Unit, unitB: Unit) {
    if (unitA.model !== unitB.model) return
    if (unitA.model.sourceData.kind !== 'mmCIF') return

    const restraints = ModelCrossLinkRestraint.Provider.get(unitA.model)
    if (!restraints) return

    const rA = new Map<number, StructureElement.UnitIndex>();
    const rB = new Map<number, StructureElement.UnitIndex>();
    _addRestraints(rA, unitA, restraints)
    _addRestraints(rB, unitB, restraints)

    rA.forEach((indexA, ri) => {
        const indexB = rB.get(ri)
        if (indexB !== undefined) {
            pairs.push(
                createCrossLinkRestraint(unitA, indexA, unitB, indexB, restraints, ri),
                createCrossLinkRestraint(unitB, indexB, unitA, indexA, restraints, ri)
            )
        }
    })
}

function extractIntra(pairs: CrossLinkRestraint[], unit: Unit) {
    if (unit.model.sourceData.kind !== 'mmCIF') return

    const restraints = ModelCrossLinkRestraint.Provider.get(unit.model)
    if (!restraints) return

    const { elements } = unit;
    const elementCount = elements.length;
    const kind = unit.kind

    const r = new Map<number, StructureElement.UnitIndex[]>();

    for (let i = 0; i < elementCount; i++) {
        const e = elements[i];
        restraints.getIndicesByElement(e, kind).forEach(ri => {
            const il = r.get(ri)
            if (il) il.push(i as StructureElement.UnitIndex)
            else r.set(ri, [i as StructureElement.UnitIndex])
        })
    }

    r.forEach((il, ri) => {
        if (il.length < 2) return
        const [ indexA, indexB ] = il
        pairs.push(
            createCrossLinkRestraint(unit, indexA, unit, indexB, restraints, ri),
            createCrossLinkRestraint(unit, indexB, unit, indexA, restraints, ri)
        )
    })
}

function createCrossLinkRestraint(unitA: Unit, indexA: StructureElement.UnitIndex, unitB: Unit, indexB: StructureElement.UnitIndex, restraints: ModelCrossLinkRestraint, row: number): CrossLinkRestraint {
    return {
        unitA, indexA, unitB, indexB,

        restraintType: restraints.data.restraint_type.value(row),
        distanceThreshold: restraints.data.distance_threshold.value(row),
        psi: restraints.data.psi.value(row),
        sigma1: restraints.data.sigma_1.value(row),
        sigma2: restraints.data.sigma_2.value(row),
    }
}

function extractCrossLinkRestraints(structure: Structure): PairRestraints<CrossLinkRestraint> {
    const pairs: CrossLinkRestraint[] = []
    if (!structure.models.some(m => ModelCrossLinkRestraint.Provider.get(m))) {
        return new PairRestraints(pairs)
    }

    const n = structure.units.length
    for (let i = 0; i < n; ++i) {
        const unitA = structure.units[i]
        extractIntra(pairs, unitA)
        for (let j = i + 1; j < n; ++j) {
            const unitB = structure.units[j]
            if (unitA.model === unitB.model) {
                extractInter(pairs, unitA, unitB)
            }
        }
    }

    return new PairRestraints(pairs)
}