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
import { SaccharideNameMap, UnknownSaccharideComponent } from './constants';
import { Vec3 } from 'mol-math/linear-algebra';
import { getCenterAndRadius, getMoleculeType } from '../../util';
import { MoleculeType } from '../../model/types';
import { areConnected } from 'mol-math/graph';
import { combinations } from 'mol-data/util/combination';
import { fillSerial } from 'mol-util/array';

function getResidueIndex(elementIndex: number, unit: Unit.Atomic) {
    return unit.model.atomicHierarchy.residueAtomSegments.index[unit.elements[elementIndex]]
}

function getRingIndices(unit: Unit.Atomic, rI: ResidueIndex) {
    const { offsets } = unit.model.atomicHierarchy.residueAtomSegments
    const { elements } = unit
    const interval = Interval.ofBounds(offsets[rI], offsets[rI + 1])
    const rings: number[] = []
    rings.push(...unit.rings.byFingerprint.get('C-C-C-C-C-O') || [])
    rings.push(...unit.rings.byFingerprint.get('C-C-C-C-O') || [])
    const sugarRings: ReadonlyArray<number>[] = []
    for (let i = 0, il = rings.length; i < il; ++i) {
        let withinIntervalCount = 0
        const ring = unit.rings.all[rings[i]]
        for (let j = 0, jl = ring.length; j < jl; ++j) {
            if (Interval.has(interval, elements[ring[j]])) ++withinIntervalCount
        }
        if (withinIntervalCount === ring.length) sugarRings.push(ring)
    }
    return sugarRings
}

export function computeCarbohydrates(structure: Structure): Carbohydrates {
    const links: CarbohydrateLink[] = []
    const terminalLinks: CarbohydrateTerminalLink[] = []
    const elements: CarbohydrateElement[] = []

    const elementsMap = new Map<string, number>()

    function elementKey(residueIndex: number, unitId: number) {
        return `${residueIndex}|${unitId}`
    }

    // get carbohydrate elements and carbohydrate links induced by intra-residue bonds
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

                const saccharideComp = SaccharideNameMap.get(label_comp_id.value(residueIndex)) || UnknownSaccharideComponent
                if (saccharideComp === UnknownSaccharideComponent) {
                    if (getMoleculeType(unit.model, residueIndex) !== MoleculeType.saccharide) continue
                }

                const sugarRings = getRingIndices(unit, residueIndex)
                const ringElements: number[] = []
                console.log('sugarRings', sugarRings)

                for (let j = 0, jl = sugarRings.length; j < jl; ++j) {
                    const center = Vec3.zero()
                    const normal = Vec3.zero()
                    const direction = Vec3.zero()
                    const elementIndex = elements.length
                    getCenterAndRadius(center, unit, sugarRings[j])

                    ringElements.push(elementIndex)
                    elementsMap.set(elementKey(residueIndex, unit.id), elementIndex)
                    elements.push({ center, normal, direction, unit, residueIndex, component: saccharideComp })
                }

                // add carbohydrate links induced by intra-residue bonds
                const ringCombinations = combinations(fillSerial(new Array(sugarRings.length)), 2)
                for (let j = 0, jl = ringCombinations.length; j < jl; ++j) {
                    const rc = ringCombinations[j]
                    if (areConnected(sugarRings[rc[0]], sugarRings[rc[1]], unit.links, 2)) {
                        links.push({
                            carbohydrateIndexA: ringElements[rc[0]],
                            carbohydrateIndexB: ringElements[rc[1]]
                        })
                        links.push({
                            carbohydrateIndexA: ringElements[rc[1]],
                            carbohydrateIndexB: ringElements[rc[0]]
                        })
                    }
                }

                if (!sugarRings.length) {
                    console.warn('No ring found for carbohydrate')
                }
            }
        }
    }

    // get carbohydrate links induced by inter-unit bonds
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        if (!Unit.isAtomic(unit)) continue

        structure.links.getLinkedUnits(unit).forEach(pairBonds => {
            pairBonds.linkedElementIndices.forEach(indexA => {
                pairBonds.getBonds(indexA).forEach(bondInfo => {
                    const { unitA, unitB } = pairBonds
                    const indexB = bondInfo.indexB
                    const elementIndexA = elementsMap.get(elementKey(getResidueIndex(indexA, unitA), unitA.id))
                    const elementIndexB = elementsMap.get(elementKey(getResidueIndex(indexB, unitB), unitB.id))

                    if (elementIndexA !== undefined && elementIndexB !== undefined) {
                        links.push({
                            carbohydrateIndexA: elementIndexA,
                            carbohydrateIndexB: elementIndexB
                        })
                    } else if (elementIndexA !== undefined) {
                        terminalLinks.push({
                            carbohydrateIndex: elementIndexA,
                            elementIndex: indexB,
                            elementUnit: unitB,
                            fromCarbohydrate: true
                        })
                    } else if (elementIndexB !== undefined) {
                        terminalLinks.push({
                            carbohydrateIndex: elementIndexB,
                            elementIndex: indexA,
                            elementUnit: unitA,
                            fromCarbohydrate: false
                        })
                    }
                })
            })
        })
    }

    return { links, terminalLinks, elements }
}