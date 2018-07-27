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
import { getMoleculeType, getPositionMatrix } from '../../util';
import { MoleculeType, ElementSymbol } from '../../model/types';
import { areConnected } from 'mol-math/graph';
import { combinations } from 'mol-data/util/combination';
import { fillSerial } from 'mol-util/array';
import PrincipalAxes from 'mol-math/linear-algebra/matrix/principal-axes';

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

const C = ElementSymbol('C')
function getDirection(direction: Vec3, unit: Unit.Atomic, indices: ReadonlyArray<number>, center: Vec3) {
    let indexC1 = -1, indexC1X = -1, indexC = -1
    const { elements } = unit
    const { position } = unit.conformation
    const { label_atom_id, type_symbol } = unit.model.atomicHierarchy.atoms
    for (let i = 0, il = indices.length; i < il; ++i) {
        const ei = elements[indices[i]]
        const atomId = label_atom_id.value(ei)
        if (atomId === 'C1') {
            indexC1 = ei
            break
        } else if (indexC1X === -1 && atomId.startsWith('C1')) {
            indexC1X = ei
        } else if (indexC === -1 && type_symbol.value(ei) === C) {
            indexC = ei
        }
    }
    const index = indexC1 !== -1 ? indexC1
        : indexC1X !== -1 ? indexC1X
        : indexC !== -1 ? indexC
        : elements[indices[0]]
    Vec3.normalize(direction, Vec3.sub(direction, center, position(index, direction)))
    return direction
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

                for (let j = 0, jl = sugarRings.length; j < jl; ++j) {
                    const pa = new PrincipalAxes(getPositionMatrix(unit, sugarRings[j]))
                    const center = Vec3.copy(Vec3.zero(), pa.center)
                    const normal = Vec3.copy(Vec3.zero(), pa.normVecC)
                    const direction = getDirection(Vec3.zero(), unit, sugarRings[j], center)
                    Vec3.orthogonalize(direction, normal, direction)

                    const elementIndex = elements.length
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