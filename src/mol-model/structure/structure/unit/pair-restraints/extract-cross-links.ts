/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../../unit';
import Structure from '../../structure';
import { PairRestraints, CrossLinkRestraint } from './data';
import { StructureElement } from '../../../structure';
import { IHMCrossLinkRestraint } from 'mol-model-formats/structure/mmcif/pair-restraint';

function _addRestraints(map: Map<number, number>, unit: Unit, restraints: IHMCrossLinkRestraint) {
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

    const restraints = IHMCrossLinkRestraint.fromModel(unitA.model)
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

    const restraints = IHMCrossLinkRestraint.fromModel(unit.model)
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

function createCrossLinkRestraint(unitA: Unit, indexA: StructureElement.UnitIndex, unitB: Unit, indexB: StructureElement.UnitIndex, restraints: IHMCrossLinkRestraint, row: number): CrossLinkRestraint {
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

    const n = structure.units.length
    for (let i = 0; i < n; ++i) {
        const unitA = structure.units[i]
        extractIntra(pairs, unitA)
        for (let j = i + 1; j < n; ++j) {
            const unitB = structure.units[j]
            extractInter(pairs, unitA, unitB)
        }
    }

    return new PairRestraints(pairs)
}

export { extractCrossLinkRestraints };
