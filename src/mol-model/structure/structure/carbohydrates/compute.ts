/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Segmentation } from 'mol-data/int';
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

const C = ElementSymbol('C'), O = ElementSymbol('O');
const SugarRingFps = [UnitRing.elementFingerprint([C, C, C, C, C, O]), UnitRing.elementFingerprint([C, C, C, C, O])]

function getDirection(direction: Vec3, unit: Unit.Atomic, indices: ArrayLike<StructureElement.UnitIndex>, center: Vec3) {
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

function getAtomId(unit: Unit.Atomic, index: number) {
    const { elements } = unit
    const { label_atom_id } = unit.model.atomicHierarchy.atoms
    return label_atom_id.value(elements[index])
}


export function computeCarbohydrates(structure: Structure): Carbohydrates {
    const links: CarbohydrateLink[] = []
    const terminalLinks: CarbohydrateTerminalLink[] = []
    const elements: CarbohydrateElement[] = []

    const elementsWithRingMap = new Map<string, number>()

    function elementKey(residueIndex: number, unitId: number) {
        return `${residueIndex}|${unitId}`
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

                const sugarRings = sugarResidueMap.get(residueIndex);

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

                    const pa = new PrincipalAxes(getPositionMatrix(unit, ringAtoms))
                    const center = Vec3.copy(Vec3.zero(), pa.center)
                    const normal = Vec3.copy(Vec3.zero(), pa.normVecC)
                    const direction = getDirection(Vec3.zero(), unit, ringAtoms, center)
                    Vec3.orthogonalize(direction, normal, direction)

                    const elementIndex = elements.length
                    ringElements.push(elementIndex)
                    elementsWithRingMap.set(elementKey(residueIndex, unit.id), elementIndex)
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
                    if (IntAdjacencyGraph.areVertexSetsConnected(unit.links, r0, r1, 2)) {
                        // fix both directions as it is unclear where the C1 atom is
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

    // get carbohydrate links induced by inter-unit bonds
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        if (!Unit.isAtomic(unit)) continue

        structure.links.getLinkedUnits(unit).forEach(pairBonds => {
            pairBonds.linkedElementIndices.forEach(indexA => {
                pairBonds.getBonds(indexA).forEach(bondInfo => {
                    const { unitA, unitB } = pairBonds
                    const indexB = bondInfo.indexB
                    const elementIndexA = elementsWithRingMap.get(elementKey(unitA.getResidueIndex(indexA), unitA.id))
                    const elementIndexB = elementsWithRingMap.get(elementKey(unitB.getResidueIndex(indexB), unitB.id))

                    if (elementIndexA !== undefined && elementIndexB !== undefined) {
                        if (getAtomId(unitA, indexA).startsWith('C1')) {
                            fixLinkDirection(elementIndexA, elementIndexB)
                        }
                        links.push({
                            carbohydrateIndexA: elementIndexA,
                            carbohydrateIndexB: elementIndexB
                        })
                    } else if (elementIndexA !== undefined) {
                        if (getAtomId(unitA, indexA).startsWith('C1')) {
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