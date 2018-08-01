/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation, SortedArray } from 'mol-data/int';
import { combinations } from 'mol-data/util/combination';
import { IntAdjacencyGraph } from 'mol-math/graph';
import { Vec3 } from 'mol-math/linear-algebra';
import PrincipalAxes from 'mol-math/linear-algebra/matrix/principal-axes';
import { fillSerial } from 'mol-util/array';
import { ResidueIndex } from '../../model';
import { ElementSymbol, MoleculeType } from '../../model/types';
import { getMoleculeType, getPositionMatrix } from '../../util';
import StructureElement from '../element';
import Structure from '../structure';
import Unit from '../unit';
import { SaccharideNameMap, UnknownSaccharideComponent } from './constants';
import { CarbohydrateElement, CarbohydrateLink, Carbohydrates, CarbohydrateTerminalLink } from './data';
import { UnitRings, UnitRing } from '../unit/rings';
import { ElementIndex } from '../../model/indexing';

const C = ElementSymbol('C'), O = ElementSymbol('O');
const SugarRingFps = [UnitRing.elementFingerprint([C, C, C, C, C, O]), UnitRing.elementFingerprint([C, C, C, C, O])]

function getAnomericCarbon(unit: Unit.Atomic, ringAtoms: ArrayLike<StructureElement.UnitIndex>): ElementIndex {
    let indexHasTwoOxygen = -1, indexHasOxygenAndCarbon = -1, indexHasC1Name = -1, indexIsCarbon = -1
    const { elements } = unit
    const { type_symbol, label_atom_id } = unit.model.atomicHierarchy.atoms
    const { b: neighbor, offset } = unit.links;
    for (let i = 0, il = ringAtoms.length; i < il; ++i) {
        const ei = elements[ringAtoms[i]]
        if (type_symbol.value(ei) !== C) continue
        let linkedOxygenCount = 0
        let linkedCarbonCount = 0
        for (let j = offset[ringAtoms[i]], jl = offset[ringAtoms[i] + 1]; j < jl; ++j) {
            const ej = elements[neighbor[j]]
            const typeSymbol = type_symbol.value(ej)
            if (typeSymbol === O) ++linkedOxygenCount
            else if (typeSymbol === C) ++linkedCarbonCount
        }
        if (linkedOxygenCount === 2) {
            // found anomeric carbon
            indexHasTwoOxygen = ei
            break
        } else if (linkedOxygenCount === 1 && linkedCarbonCount === 1) {
            // possibly an anomeric carbon if this is a mono-saccharide without a glycosidic bond
            indexHasOxygenAndCarbon = ei
        } else if (label_atom_id.value(ei).startsWith('C1')) {
            // likely the anomeric carbon as it is name C1 by convention
            indexHasC1Name = ei
        } else {
            // use any carbon as a fallback
            indexIsCarbon = ei
        }
    }
    return (indexHasTwoOxygen !== -1 ? indexHasTwoOxygen
                : indexHasOxygenAndCarbon !== -1 ? indexHasOxygenAndCarbon
                    : indexHasC1Name !== -1 ? indexHasC1Name
                        : indexIsCarbon !== -1 ? indexIsCarbon
                            : elements[ringAtoms[0]]) as ElementIndex
}

/** Return first non-empty label_alt_id or an empty string */
function getRingAltId(unit: Unit.Atomic, ringAtoms: SortedArray<StructureElement.UnitIndex>) {
    const { elements } = unit
    const { label_alt_id } = unit.model.atomicHierarchy.atoms
    for (let i = 0, il = ringAtoms.length; i < il; ++i) {
        const ei = elements[ringAtoms[i]]
        const altId = label_alt_id.value(ei)
        if (altId) return altId
    }
    return ''
}

function getAltId(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    const { elements } = unit
    const { label_alt_id } = unit.model.atomicHierarchy.atoms
    return label_alt_id.value(elements[index])
}

function getDirection(direction: Vec3, unit: Unit.Atomic, index: ElementIndex, center: Vec3) {
    const { position } = unit.conformation
    Vec3.normalize(direction, Vec3.sub(direction, center, position(index, direction)))
    return direction
}

function getAtomId(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    const { elements } = unit
    const { label_atom_id } = unit.model.atomicHierarchy.atoms
    return label_atom_id.value(elements[index])
}

function filterFusedRings(unitRings: UnitRings, rings: UnitRings.Index[] | undefined) {
    if (!rings || !rings.length) return

    const fusedRings = new Set<UnitRings.Index>()
    const ringCombinations = combinations(fillSerial(new Array(rings.length) as number[]), 2)
    for (let i = 0, il = ringCombinations.length; i < il; ++i) {
        const rc = ringCombinations[i];
        const r0 = unitRings.all[rings[rc[0]]], r1 = unitRings.all[rings[rc[1]]];
        if (SortedArray.areIntersecting(r0, r1)) {
            fusedRings.add(rings[rc[0]])
            fusedRings.add(rings[rc[1]])
        }
    }

    if (fusedRings.size) {
        const filteredRings: UnitRings.Index[] = []
        for (let i = 0, il = rings.length; i < il; ++i) {
            if (!fusedRings.has(rings[i])) filteredRings.push(rings[i])
        }
        return filteredRings
    } else {
        return rings
    }
}

export function computeCarbohydrates(structure: Structure): Carbohydrates {
    const links: CarbohydrateLink[] = []
    const terminalLinks: CarbohydrateTerminalLink[] = []
    const elements: CarbohydrateElement[] = []

    const elementsWithRingMap = new Map<string, number>()

    function elementKey(residueIndex: number, unitId: number, altId: string) {
        return `${residueIndex}|${unitId}|${altId}`
    }

    function fixLinkDirection(iA: number, iB: number) {
        Vec3.sub(elements[iA].geometry!.direction, elements[iB].geometry!.center, elements[iA].geometry!.center)
        Vec3.normalize(elements[iA].geometry!.direction, elements[iA].geometry!.direction)
    }

    const tmpV = Vec3.zero()
    function fixTerminalLinkDirection(iA: number, indexB: number, unitB: Unit.Atomic) {
        const pos = unitB.conformation.position, geo = elements[iA].geometry!;
        Vec3.sub(geo.direction, pos(unitB.elements[indexB], tmpV), geo.center)
        Vec3.normalize(geo.direction, geo.direction)
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

        let sugarResidueMap: Map<ResidueIndex, UnitRings.Index[]> | undefined = void 0;

        while (chainIt.hasNext) {
            residueIt.setSegment(chainIt.move());

            while (residueIt.hasNext) {
                const { index: residueIndex } = residueIt.move();

                const saccharideComp = SaccharideNameMap.get(label_comp_id.value(residueIndex)) || UnknownSaccharideComponent
                if (saccharideComp === UnknownSaccharideComponent) {
                    if (getMoleculeType(unit.model, residueIndex) !== MoleculeType.saccharide) continue
                }

                if (!sugarResidueMap) {
                    sugarResidueMap = UnitRings.byFingerprintAndResidue(unit.rings, SugarRingFps);
                }

                const sugarRings = filterFusedRings(unit.rings, sugarResidueMap.get(residueIndex));

                if (!sugarRings || !sugarRings.length) {
                    elements.push({
                        hasRing: false,
                        unit, residueIndex, component: saccharideComp
                    })
                    continue;
                }

                const rings = unit.rings;
                const ringElements: number[] = []

                for (let j = 0, jl = sugarRings.length; j < jl; ++j) {
                    const ringAtoms = rings.all[sugarRings[j]];
                    const anomericCarbon = getAnomericCarbon(unit, ringAtoms)

                    const pa = new PrincipalAxes(getPositionMatrix(unit, ringAtoms))
                    const center = Vec3.copy(Vec3.zero(), pa.center)
                    const normal = Vec3.copy(Vec3.zero(), pa.normVecC)
                    const direction = getDirection(Vec3.zero(), unit, anomericCarbon, center)
                    Vec3.orthogonalize(direction, normal, direction)

                    const altId = getRingAltId(unit, ringAtoms)
                    const elementIndex = elements.length
                    ringElements.push(elementIndex)
                    elementsWithRingMap.set(elementKey(residueIndex, unit.id, altId), elementIndex)
                    elements.push({
                        geometry: { center, normal, direction },
                        hasRing: true,
                        unit, residueIndex, component: saccharideComp
                    })
                }

                // add carbohydrate links induced by intra-residue bonds
                const ringCombinations = combinations(fillSerial(new Array(sugarRings.length) as number[]), 2)
                for (let j = 0, jl = ringCombinations.length; j < jl; ++j) {
                    const rc = ringCombinations[j];
                    const r0 = rings.all[sugarRings[rc[0]]], r1 = rings.all[sugarRings[rc[1]]];
                    // 1,6 glycosidic links are distance 3 and 1,4 glycosidic links are distance 2
                    if (IntAdjacencyGraph.areVertexSetsConnected(unit.links, r0, r1, 3)) {
                        // TODO handle better, for now fix both directions as it is unclear where the C1 atom is
                        //  would need to know the path connecting the two rings
                        fixLinkDirection(ringElements[rc[0]], ringElements[rc[1]])
                        fixLinkDirection(ringElements[rc[1]], ringElements[rc[0]])
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
            }
        }
    }

    function getElementIndex(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
        return elementsWithRingMap.get(elementKey(unit.getResidueIndex(index), unit.id, getAltId(unit, index)))
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
                    const elementIndexA = getElementIndex(unitA, indexA)
                    const elementIndexB = getElementIndex(unitB, indexB)

                    if (elementIndexA !== undefined && elementIndexB !== undefined) {
                        const atomIdA = getAtomId(unitA, indexA)
                        if (atomIdA.startsWith('O1') || atomIdA.startsWith('C1')) {
                            fixLinkDirection(elementIndexA, elementIndexB)
                        }
                        links.push({
                            carbohydrateIndexA: elementIndexA,
                            carbohydrateIndexB: elementIndexB
                        })
                    } else if (elementIndexA !== undefined) {
                        const atomIdA = getAtomId(unitA, indexA)
                        if (atomIdA.startsWith('O1') || atomIdA.startsWith('C1')) {
                            fixTerminalLinkDirection(elementIndexA, indexB, unitB)
                        }
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