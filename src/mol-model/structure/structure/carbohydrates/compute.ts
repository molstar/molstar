/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation, SortedArray } from '../../../../mol-data/int';
import { combinations } from '../../../../mol-data/util/combination';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { PrincipalAxes } from '../../../../mol-math/linear-algebra/matrix/principal-axes';
import { fillSerial } from '../../../../mol-util/array';
import { ResidueIndex, Model } from '../../model';
import { ElementSymbol } from '../../model/types';
import { getPositions } from '../../util';
import StructureElement from '../element';
import Structure from '../structure';
import Unit from '../unit';
import { CarbohydrateElement, CarbohydrateLink, Carbohydrates, CarbohydrateTerminalLink, PartialCarbohydrateElement, EmptyCarbohydrates } from './data';
import { UnitRings, UnitRing } from '../unit/rings';
import { ElementIndex } from '../../model/indexing';

const C = ElementSymbol('C'), O = ElementSymbol('O');
const SugarRingFps = [
    UnitRing.elementFingerprint([C, C, C, O]),
    UnitRing.elementFingerprint([C, C, C, C, O]),
    UnitRing.elementFingerprint([C, C, C, C, C, O]),
    UnitRing.elementFingerprint([C, C, C, C, C, C, O]),
]

function getAnomericCarbon(unit: Unit.Atomic, ringAtoms: ArrayLike<StructureElement.UnitIndex>): ElementIndex {
    let indexHasTwoOxygen = -1, indexHasOxygenAndCarbon = -1, indexHasC1Name = -1, indexIsCarbon = -1
    const { elements } = unit
    const { type_symbol, label_atom_id } = unit.model.atomicHierarchy.atoms
    const { b: neighbor, offset } = unit.bonds;
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
            // likely the anomeric carbon as it is named C1 by convention
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

function getSaccharideComp(compId: string, model: Model) {
    return model.properties.saccharideComponentMap.get(compId)
}

export function computeCarbohydrates(structure: Structure): Carbohydrates {
    // skip computation if there are no saccharide components in any model
    if (structure.models.reduce((a, v) => a + v.properties.saccharideComponentMap.size, 0) === 0)
        return EmptyCarbohydrates

    const links: CarbohydrateLink[] = []
    const terminalLinks: CarbohydrateTerminalLink[] = []
    const elements: CarbohydrateElement[] = []
    const partialElements: PartialCarbohydrateElement[] = []

    const elementsWithRingMap = new Map<string, number>()

    function ringElementKey(residueIndex: number, unitId: number, altId: string) {
        return `${residueIndex}|${unitId}|${altId}`
    }

    function fixLinkDirection(iA: number, iB: number) {
        Vec3.sub(elements[iA].geometry.direction, elements[iB].geometry.center, elements[iA].geometry.center)
        Vec3.normalize(elements[iA].geometry.direction, elements[iA].geometry.direction)
    }

    const tmpV = Vec3.zero()
    function fixTerminalLinkDirection(iA: number, indexB: number, unitB: Unit.Atomic) {
        const pos = unitB.conformation.position, geo = elements[iA].geometry;
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

                const saccharideComp = getSaccharideComp(label_comp_id.value(residueIndex), model)
                if (!saccharideComp) continue

                if (!sugarResidueMap) {
                    sugarResidueMap = UnitRings.byFingerprintAndResidue(unit.rings, SugarRingFps);
                }

                const sugarRings = filterFusedRings(unit.rings, sugarResidueMap.get(residueIndex));

                if (!sugarRings || !sugarRings.length) {
                    partialElements.push({ unit, residueIndex, component: saccharideComp })
                    continue;
                }

                const rings = unit.rings;
                const ringElements: number[] = []

                for (let j = 0, jl = sugarRings.length; j < jl; ++j) {
                    const ringAtoms = rings.all[sugarRings[j]];
                    const anomericCarbon = getAnomericCarbon(unit, ringAtoms)

                    const ma = PrincipalAxes.calculateMomentsAxes(getPositions(unit, ringAtoms))
                    const center = Vec3.copy(Vec3.zero(), ma.origin)
                    const normal = Vec3.copy(Vec3.zero(), ma.dirC)
                    const direction = getDirection(Vec3.zero(), unit, anomericCarbon, center)
                    Vec3.orthogonalize(direction, normal, direction)

                    const ringAltId = getRingAltId(unit, ringAtoms)
                    const elementIndex = elements.length
                    ringElements.push(elementIndex)
                    elementsWithRingMap.set(ringElementKey(residueIndex, unit.id, ringAltId), elementIndex)
                    elements.push({
                        geometry: { center, normal, direction },
                        component: saccharideComp,
                        unit, residueIndex, anomericCarbon, ringAltId,
                        ringMemberCount: ringAtoms.length
                    })
                }

                // add carbohydrate links induced by intra-residue bonds
                // (e.g. for structures from the PDB archive __before__ carbohydrate remediation)
                const ringCombinations = combinations(fillSerial(new Array(sugarRings.length) as number[]), 2)
                for (let j = 0, jl = ringCombinations.length; j < jl; ++j) {
                    const rc = ringCombinations[j];
                    const r0 = rings.all[sugarRings[rc[0]]], r1 = rings.all[sugarRings[rc[1]]];
                    // 1,6 glycosidic links are distance 3 and 1,4 glycosidic links are distance 2
                    if (IntAdjacencyGraph.areVertexSetsConnected(unit.bonds, r0, r1, 3)) {
                        const re0 = ringElements[rc[0]]
                        const re1 = ringElements[rc[1]]
                        if (elements[re0].ringAltId === elements[re1].ringAltId) {
                            // TODO handle better, for now fix both directions as it is unclear where the C1 atom is
                            //  would need to know the path connecting the two rings
                            fixLinkDirection(re0, re1)
                            fixLinkDirection(re1, re0)
                            links.push({ carbohydrateIndexA: re0, carbohydrateIndexB: re1 })
                            links.push({ carbohydrateIndexA: re1, carbohydrateIndexB: re0 })
                        }
                    }
                }
            }
        }
    }

    function getRingElementIndex(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
        return elementsWithRingMap.get(ringElementKey(unit.getResidueIndex(index), unit.id, getAltId(unit, index)))
    }

    // add carbohydrate links induced by intra-unit bonds
    // (e.g. for structures from the PDB archive __after__ carbohydrate remediation)
    for (let i = 0, il = elements.length; i < il; ++i) {
        const carbohydrate = elements[i]
        const { unit, residueIndex, anomericCarbon } = carbohydrate
        const { offset, b } = unit.bonds
        const ac = SortedArray.indexOf(unit.elements, anomericCarbon) as StructureElement.UnitIndex

        for (let j = offset[ac], jl = offset[ac + 1]; j < jl; ++j) {
            const bj = b[j] as StructureElement.UnitIndex
            if (residueIndex !== unit.getResidueIndex(bj)) {
                const ringElementIndex = getRingElementIndex(unit, bj)
                if (ringElementIndex !== undefined && ringElementIndex !== i) {
                    fixLinkDirection(i, ringElementIndex)
                    links.push({
                        carbohydrateIndexA: i,
                        carbohydrateIndexB: ringElementIndex
                    })
                    links.push({
                        carbohydrateIndexA: ringElementIndex,
                        carbohydrateIndexB: i
                    })
                }
            }
        }

    }

    // get carbohydrate links induced by inter-unit bonds, that is
    // terminal links plus inter monosaccharide links for structures from the
    // PDB archive __before__ carbohydrate remediation
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        if (!Unit.isAtomic(unit)) continue

        structure.interUnitBonds.getConnectedUnits(unit).forEach(pairBonds => {
            pairBonds.connectedIndices.forEach(indexA => {
                pairBonds.getEdges(indexA).forEach(bondInfo => {
                    const { unitA, unitB } = pairBonds
                    const indexB = bondInfo.indexB
                    const ringElementIndexA = getRingElementIndex(unitA, indexA)
                    const ringElementIndexB = getRingElementIndex(unitB, indexB)

                    if (ringElementIndexA !== undefined && ringElementIndexB !== undefined) {
                        const atomIdA = getAtomId(unitA, indexA)
                        if (atomIdA.startsWith('O1') || atomIdA.startsWith('C1')) {
                            fixLinkDirection(ringElementIndexA, ringElementIndexB)
                        }
                        links.push({
                            carbohydrateIndexA: ringElementIndexA,
                            carbohydrateIndexB: ringElementIndexB
                        })
                    } else if (ringElementIndexA !== undefined) {
                        const atomIdA = getAtomId(unitA, indexA)
                        if (atomIdA.startsWith('O1') || atomIdA.startsWith('C1')) {
                            fixTerminalLinkDirection(ringElementIndexA, indexB, unitB)
                        }
                        terminalLinks.push({
                            carbohydrateIndex: ringElementIndexA,
                            elementIndex: indexB,
                            elementUnit: unitB,
                            fromCarbohydrate: true
                        })
                    } else if (ringElementIndexB !== undefined) {
                        terminalLinks.push({
                            carbohydrateIndex: ringElementIndexB,
                            elementIndex: indexA,
                            elementUnit: unitA,
                            fromCarbohydrate: false
                        })
                    }
                })
            })
        })
    }

    return { links, terminalLinks, elements, partialElements, ...buildLookups(elements, links, terminalLinks) }
}

function buildLookups (elements: CarbohydrateElement[], links: CarbohydrateLink[], terminalLinks: CarbohydrateTerminalLink[]) {
    // element lookup

    function elementKey(unit: Unit, anomericCarbon: ElementIndex) {
        return `${unit.id}|${anomericCarbon}`
    }

    const elementMap = new Map<string, number>()
    for (let i = 0, il = elements.length; i < il; ++i) {
        const { unit, anomericCarbon } = elements[i]
        elementMap.set(elementKey(unit, anomericCarbon), i)
    }

    function getElementIndex(unit: Unit, anomericCarbon: ElementIndex) {
        return elementMap.get(elementKey(unit, anomericCarbon))
    }

    // link lookup

    function linkKey(unitA: Unit, anomericCarbonA: ElementIndex, unitB: Unit, anomericCarbonB: ElementIndex) {
        return `${unitA.id}|${anomericCarbonA}|${unitB.id}|${anomericCarbonB}`
    }

    const linkMap = new Map<string, number>()
    for (let i = 0, il = links.length; i < il; ++i) {
        const l = links[i]
        const { unit: unitA, anomericCarbon: anomericCarbonA } = elements[l.carbohydrateIndexA]
        const { unit: unitB, anomericCarbon: anomericCarbonB } = elements[l.carbohydrateIndexB]
        linkMap.set(linkKey(unitA, anomericCarbonA, unitB, anomericCarbonB), i)
    }

    function getLinkIndex(unitA: Unit, anomericCarbonA: ElementIndex, unitB: Unit, anomericCarbonB: ElementIndex) {
        return linkMap.get(linkKey(unitA, anomericCarbonA, unitB, anomericCarbonB))
    }

    // links lookup

    function linksKey(unit: Unit, anomericCarbon: ElementIndex) {
        return `${unit.id}|${anomericCarbon}`
    }

    const linksMap = new Map<string, number[]>()
    for (let i = 0, il = links.length; i < il; ++i) {
        const l = links[i]
        const { unit, anomericCarbon } = elements[l.carbohydrateIndexA]
        const k = linksKey(unit, anomericCarbon)
        const e = linksMap.get(k)
        if (e === undefined) linksMap.set(k, [i])
        else e.push(i)
    }

    function getLinkIndices(unit: Unit, anomericCarbon: ElementIndex): ReadonlyArray<number> {
        return linksMap.get(linksKey(unit, anomericCarbon)) || []
    }

    // terminal link lookup

    function terminalLinkKey(unitA: Unit, elementA: ElementIndex, unitB: Unit, elementB: ElementIndex) {
        return `${unitA.id}|${elementA}|${unitB.id}|${elementB}`
    }

    const terminalLinkMap = new Map<string, number>()
    for (let i = 0, il = terminalLinks.length; i < il; ++i) {
        const { fromCarbohydrate, carbohydrateIndex, elementUnit, elementIndex } = terminalLinks[i]
        const { unit, anomericCarbon } = elements[carbohydrateIndex]
        if (fromCarbohydrate) {
            terminalLinkMap.set(terminalLinkKey(unit, anomericCarbon, elementUnit, elementUnit.elements[elementIndex]), i)
        } else {
            terminalLinkMap.set(terminalLinkKey(elementUnit, elementUnit.elements[elementIndex], unit, anomericCarbon), i)
        }
    }

    function getTerminalLinkIndex(unitA: Unit, elementA: ElementIndex, unitB: Unit, elementB: ElementIndex) {
        return terminalLinkMap.get(terminalLinkKey(unitA, elementA, unitB, elementB))
    }

    // terminal links lookup

    function terminalLinksKey(unit: Unit, element: ElementIndex) {
        return `${unit.id}|${element}`
    }

    const terminalLinksMap = new Map<string, number[]>()
    for (let i = 0, il = terminalLinks.length; i < il; ++i) {
        const { fromCarbohydrate, carbohydrateIndex, elementUnit, elementIndex } = terminalLinks[i]
        const { unit, anomericCarbon } = elements[carbohydrateIndex]
        let k: string
        if (fromCarbohydrate) {
            k = terminalLinksKey(unit, anomericCarbon)
        } else {
            k = terminalLinksKey(elementUnit, elementUnit.elements[elementIndex])
        }
        const e = terminalLinksMap.get(k)
        if (e === undefined) terminalLinksMap.set(k, [i])
        else e.push(i)
    }

    function getTerminalLinkIndices(unit: Unit, element: ElementIndex): ReadonlyArray<number> {
        return terminalLinksMap.get(terminalLinksKey(unit, element)) || []
    }

    // anomeric carbon lookup

    function anomericCarbonKey(unit: Unit, residueIndex: ResidueIndex) {
        return `${unit.id}|${residueIndex}`
    }

    const anomericCarbonMap = new Map<string, ElementIndex[]>()
    for (let i = 0, il = elements.length; i < il; ++i) {
        const { unit, anomericCarbon } = elements[i]
        const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[anomericCarbon]
        const k = anomericCarbonKey(unit, residueIndex)
        if (anomericCarbonMap.has(k)) {
            anomericCarbonMap.get(k)!.push(anomericCarbon)
        } else {
            anomericCarbonMap.set(k, [anomericCarbon])
        }
    }

    const EmptyArray: ReadonlyArray<any> = []
    function getAnomericCarbons(unit: Unit, residueIndex: ResidueIndex) {
        return anomericCarbonMap.get(anomericCarbonKey(unit, residueIndex)) || EmptyArray
    }

    return { getElementIndex, getLinkIndex, getLinkIndices, getTerminalLinkIndex, getTerminalLinkIndices, getAnomericCarbons }
}