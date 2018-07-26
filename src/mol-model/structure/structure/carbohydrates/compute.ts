/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../unit';
import { ResidueIndex } from '../../model';
import { Interval, Segmentation } from 'mol-data/int';
import Structure from '../structure';
import { Carbohydrates, CarbohydrateLink, CarbohydrateTerminalLink, CarbohydrateElement } from './data';
import { SaccharideNameMap } from './constants';
import { Vec3 } from 'mol-math/linear-algebra';
import { getCenterAndRadius } from '../../util';

function getResidueIndex(elementIndex: number, unit: Unit.Atomic) {
    return unit.model.atomicHierarchy.residueAtomSegments.index[unit.elements[elementIndex]]
}

function getRingIndices(unit: Unit.Atomic, rI: ResidueIndex) {
    const { offsets } = unit.model.atomicHierarchy.residueAtomSegments
    const { elements } = unit
    const interval = Interval.ofBounds(offsets[rI], offsets[rI + 1])
    const rings = unit.rings.byFingerprint.get('C-C-C-C-C-O') || unit.rings.byFingerprint.get('C-C-C-C-O')
    if (rings) {
        for (let i = 0, il = rings.length; i < il; ++i) {
            let withinIntervalCount = 0
            const ring = unit.rings.all[rings[i]]
            for (let j = 0, jl = ring.length; j < jl; ++j) {
                if (Interval.has(interval, elements[ring[j]])) ++withinIntervalCount
            }
            if (withinIntervalCount === ring.length) return ring
        }
    }
}

export function computeCarbohydrates(structure: Structure): Carbohydrates {
    const links: CarbohydrateLink[] = []
    const terminalLinks: CarbohydrateTerminalLink[] = []
    const elements: CarbohydrateElement[] = []

    const elementsMap = new Map<string, number>()

    function elementKey(residueIndex: number, unitId: number) {
        return `${residueIndex}|${unitId}`
    }

    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        if (!Unit.isAtomic(unit)) continue

        const { model } = unit
        const { chainAtomSegments, residueAtomSegments, residues } = model.atomicHierarchy
        const { label_comp_id } = residues

        const chainIt = Segmentation.transientSegments(chainAtomSegments, unit.elements)
        const residueIt = Segmentation.transientSegments(residueAtomSegments, unit.elements)

        while (chainIt.hasNext) {
            residueIt.setSegment(chainIt.move());

            while (residueIt.hasNext) {
                const { index: residueIndex } = residueIt.move();

                const saccharideComp = SaccharideNameMap.get(label_comp_id.value(residueIndex))
                if (!saccharideComp) continue

                const ringIndices = getRingIndices(unit, residueIndex)
                if (ringIndices) {
                    const center = Vec3.zero()
                    const normal = Vec3.zero()
                    const direction = Vec3.zero()
                    const ringRadius = getCenterAndRadius(center, unit, ringIndices)
                    console.log(ringRadius, center)

                    elementsMap.set(elementKey(residueIndex, unit.id), elements.length)
                    elements.push({ center, normal, direction, unit, residueIndex, component: saccharideComp })
                } else {
                    console.warn('No ring found for carbohydrate')
                }
            }
        }
    }

    elementsMap.forEach((elementIndex, key) => {
        const unit = elements[elementIndex].unit
        structure.links.getLinkedUnits(unit).forEach(pairBonds => {
            pairBonds.linkedElementIndices.forEach(indexA => {
                pairBonds.getBonds(indexA).forEach(bondInfo => {
                    let { unitA, unitB } = pairBonds
                    let indexB = bondInfo.indexB
                    let residueIndexA = getResidueIndex(indexA, unitA)
                    let residueIndexB = getResidueIndex(indexB, unitB)
                    let keyA = elementKey(residueIndexA, unitA.id)
                    let keyB = elementKey(residueIndexB, unitB.id)
                    if (key === keyB) {
                        [keyB, keyA] = [keyA, keyB];
                        [indexB, indexA] = [indexA, indexB];
                        [unitB, unitA] = [unitA, unitB];
                    }
                    const elementIndexB = elementsMap.get(keyB)
                    if (elementIndexB !== undefined) {
                        links.push({
                            carbohydrateIndexA: elementIndex,
                            carbohydrateIndexB: elementIndexB
                        })
                    } else {
                        terminalLinks.push({
                            carbohydrateIndex: elementIndex,
                            elementIndex: indexB,
                            elementUnit: unitB,
                            fromCarbohydrate: true
                        })
                        terminalLinks.push({
                            carbohydrateIndex: elementIndex,
                            elementIndex: indexB,
                            elementUnit: unitB,
                            fromCarbohydrate: false
                        })
                    }
                })
            })
        })
    })

    return { links, terminalLinks, elements }
}