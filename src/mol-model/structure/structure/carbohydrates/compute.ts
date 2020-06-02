/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation, SortedArray } from '../../../../mol-data/int';
import { combinations } from '../../../../mol-data/util/combination';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { PrincipalAxes } from '../../../../mol-math/linear-algebra/matrix/principal-axes';
import { fillSerial, arraySetAdd } from '../../../../mol-util/array';
import { ResidueIndex, Model } from '../../model';
import { ElementSymbol, BondType } from '../../model/types';
import { getPositions } from '../../util';
import StructureElement from '../element';
import Structure from '../structure';
import Unit from '../unit';
import { CarbohydrateElement, CarbohydrateLink, Carbohydrates, CarbohydrateTerminalLink, PartialCarbohydrateElement, EmptyCarbohydrates } from './data';
import { UnitRings, UnitRing } from '../unit/rings';
import { ElementIndex } from '../../model/indexing';
import { cantorPairing } from '../../../../mol-data/util';

const C = ElementSymbol('C'), O = ElementSymbol('O');
const SugarRingFps = [
    UnitRing.elementFingerprint([C, C, C, O]),
    UnitRing.elementFingerprint([C, C, C, C, O]),
    UnitRing.elementFingerprint([C, C, C, C, C, O]),
    UnitRing.elementFingerprint([C, C, C, C, C, C, O]),
];

function getAnomericCarbon(unit: Unit.Atomic, ringAtoms: ArrayLike<StructureElement.UnitIndex>): ElementIndex {
    let indexHasTwoOxygen = -1, indexHasOxygenAndCarbon = -1, indexHasC1Name = -1, indexIsCarbon = -1;
    const { elements } = unit;
    const { type_symbol, label_atom_id } = unit.model.atomicHierarchy.atoms;
    const { b: neighbor, offset } = unit.bonds;
    for (let i = 0, il = ringAtoms.length; i < il; ++i) {
        const ei = elements[ringAtoms[i]];
        if (type_symbol.value(ei) !== C) continue;
        let linkedOxygenCount = 0;
        let linkedCarbonCount = 0;
        for (let j = offset[ringAtoms[i]], jl = offset[ringAtoms[i] + 1]; j < jl; ++j) {
            const ej = elements[neighbor[j]];
            const typeSymbol = type_symbol.value(ej);
            if (typeSymbol === O) ++linkedOxygenCount;
            else if (typeSymbol === C) ++linkedCarbonCount;
        }
        if (linkedOxygenCount === 2) {
            // found anomeric carbon
            indexHasTwoOxygen = ei;
            break;
        } else if (linkedOxygenCount === 1 && linkedCarbonCount === 1) {
            // possibly an anomeric carbon if this is a mono-saccharide without a glycosidic bond
            indexHasOxygenAndCarbon = ei;
        } else if (label_atom_id.value(ei).startsWith('C1')) {
            // likely the anomeric carbon as it is named C1 by convention
            indexHasC1Name = ei;
        } else {
            // use any carbon as a fallback
            indexIsCarbon = ei;
        }
    }
    return (indexHasTwoOxygen !== -1 ? indexHasTwoOxygen
        : indexHasOxygenAndCarbon !== -1 ? indexHasOxygenAndCarbon
            : indexHasC1Name !== -1 ? indexHasC1Name
                : indexIsCarbon !== -1 ? indexIsCarbon
                    : elements[ringAtoms[0]]) as ElementIndex;
}

function getAltId(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    const { elements } = unit;
    const { label_alt_id } = unit.model.atomicHierarchy.atoms;
    return label_alt_id.value(elements[index]);
}

function getDirection(direction: Vec3, unit: Unit.Atomic, index: ElementIndex, center: Vec3) {
    const { position } = unit.conformation;
    Vec3.normalize(direction, Vec3.sub(direction, center, position(index, direction)));
    return direction;
}

function getAtomId(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    const { elements } = unit;
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    return label_atom_id.value(elements[index]);
}

function filterFusedRings(unitRings: UnitRings, rings: UnitRings.Index[] | undefined) {
    if (!rings || !rings.length) return;

    const { unit, all } = unitRings;
    const fusedRings = new Set<UnitRings.Index>();
    const ringCombinations = combinations(fillSerial(new Array(rings.length) as number[]), 2);
    for (let i = 0, il = ringCombinations.length; i < il; ++i) {
        const rc = ringCombinations[i];
        const r0 = all[rings[rc[0]]], r1 = all[rings[rc[1]]];
        if (SortedArray.areIntersecting(r0, r1) &&
                UnitRing.getAltId(unit, r0) === UnitRing.getAltId(unit, r1)) {
            fusedRings.add(rings[rc[0]]);
            fusedRings.add(rings[rc[1]]);
        }
    }

    if (fusedRings.size) {
        const filteredRings: UnitRings.Index[] = [];
        for (let i = 0, il = rings.length; i < il; ++i) {
            if (!fusedRings.has(rings[i])) filteredRings.push(rings[i]);
        }
        return filteredRings;
    } else {
        return rings;
    }
}

function getSaccharideComp(compId: string, model: Model) {
    return model.properties.saccharideComponentMap.get(compId);
}

export function computeCarbohydrates(structure: Structure): Carbohydrates {
    // skip computation if there are no saccharide components in any model
    if (structure.models.reduce((a, v) => a + v.properties.saccharideComponentMap.size, 0) === 0)
        return EmptyCarbohydrates;

    const links: CarbohydrateLink[] = [];
    const terminalLinks: CarbohydrateTerminalLink[] = [];
    const elements: CarbohydrateElement[] = [];
    const partialElements: PartialCarbohydrateElement[] = [];

    const elementsWithRingMap = new Map<string, number[]>();
    function ringElementKey(residueIndex: number, unitId: number, altId: string) {
        return `${residueIndex}|${unitId}|${altId}`;
    }
    function addRingElement(key: string, elementIndex: number) {
        if (elementsWithRingMap.has(key)) elementsWithRingMap.get(key)!.push(elementIndex);
        else elementsWithRingMap.set(key, [elementIndex]);
    }

    function fixLinkDirection(iA: number, iB: number) {
        Vec3.sub(elements[iA].geometry.direction, elements[iB].geometry.center, elements[iA].geometry.center);
        Vec3.normalize(elements[iA].geometry.direction, elements[iA].geometry.direction);
    }

    const tmpV = Vec3.zero();
    function fixTerminalLinkDirection(iA: number, indexB: number, unitB: Unit.Atomic) {
        const pos = unitB.conformation.position, geo = elements[iA].geometry;
        Vec3.sub(geo.direction, pos(unitB.elements[indexB], tmpV), geo.center);
        Vec3.normalize(geo.direction, geo.direction);
    }

    // get carbohydrate elements and carbohydrate links induced by intra-residue bonds
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i];
        if (!Unit.isAtomic(unit)) continue;

        const { model, rings } = unit;
        const { chainAtomSegments, residueAtomSegments, atoms } = model.atomicHierarchy;
        const { label_comp_id } = atoms;

        const chainIt = Segmentation.transientSegments(chainAtomSegments, unit.elements);
        const residueIt = Segmentation.transientSegments(residueAtomSegments, unit.elements);

        let sugarResidueMap: Map<ResidueIndex, UnitRings.Index[]> | undefined = void 0;

        while (chainIt.hasNext) {
            residueIt.setSegment(chainIt.move());

            while (residueIt.hasNext) {
                const { index: residueIndex } = residueIt.move();

                const saccharideComp = getSaccharideComp(label_comp_id.value(residueAtomSegments.offsets[residueIndex]), model);
                if (!saccharideComp) continue;

                if (!sugarResidueMap) {
                    sugarResidueMap = UnitRings.byFingerprintAndResidue(rings, SugarRingFps);
                }

                const sugarRings = filterFusedRings(rings, sugarResidueMap.get(residueIndex));

                if (!sugarRings || !sugarRings.length) {
                    partialElements.push({ unit, residueIndex, component: saccharideComp });
                    continue;
                }

                const ringElements: number[] = [];

                for (let j = 0, jl = sugarRings.length; j < jl; ++j) {
                    const ringAtoms = rings.all[sugarRings[j]];
                    const anomericCarbon = getAnomericCarbon(unit, ringAtoms);

                    const ma = PrincipalAxes.calculateMomentsAxes(getPositions(unit, ringAtoms));
                    const center = Vec3.copy(Vec3.zero(), ma.origin);
                    const normal = Vec3.copy(Vec3.zero(), ma.dirC);
                    const direction = getDirection(Vec3.zero(), unit, anomericCarbon, center);
                    Vec3.orthogonalize(direction, normal, direction);

                    const ringAltId = UnitRing.getAltId(unit, ringAtoms);
                    const elementIndex = elements.length;
                    ringElements.push(elementIndex);

                    addRingElement(ringElementKey(residueIndex, unit.id, ringAltId), elementIndex);
                    if (ringAltId) addRingElement(ringElementKey(residueIndex, unit.id, ''), elementIndex);

                    elements.push({
                        geometry: { center, normal, direction },
                        component: saccharideComp,
                        ringIndex: sugarRings[j],
                        altId: ringAltId,
                        unit, residueIndex
                    });
                }

                // add carbohydrate links induced by intra-residue bonds
                // (e.g. for structures from the PDB archive __before__ carbohydrate remediation)
                const ringCombinations = combinations(fillSerial(new Array(sugarRings.length) as number[]), 2);
                for (let j = 0, jl = ringCombinations.length; j < jl; ++j) {
                    const rc = ringCombinations[j];
                    const r0 = rings.all[sugarRings[rc[0]]], r1 = rings.all[sugarRings[rc[1]]];
                    // 1,6 glycosidic links are distance 3 and 1,4 glycosidic links are distance 2
                    if (IntAdjacencyGraph.areVertexSetsConnected(unit.bonds, r0, r1, 3)) {
                        const re0 = ringElements[rc[0]];
                        const re1 = ringElements[rc[1]];
                        if (elements[re0].altId === elements[re1].altId) {
                            // TODO handle better, for now fix both directions as it is unclear where the C1 atom is
                            //      would need to know the path connecting the two rings
                            fixLinkDirection(re0, re1);
                            fixLinkDirection(re1, re0);
                            links.push({ carbohydrateIndexA: re0, carbohydrateIndexB: re1 });
                            links.push({ carbohydrateIndexA: re1, carbohydrateIndexB: re0 });
                        }
                    }
                }
            }
        }
    }

    function getRingElementIndices(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
        return elementsWithRingMap.get(ringElementKey(unit.getResidueIndex(index), unit.id, getAltId(unit, index))) || [];
    }

    // add carbohydrate links induced by intra-unit bonds
    // (e.g. for structures from the PDB archive __after__ carbohydrate remediation)
    for (let i = 0, il = elements.length; i < il; ++i) {
        const cA = elements[i];
        const { unit } = cA;

        for (let j = i + 1; j < il; ++j) {
            const cB = elements[j];
            if (unit !== cB.unit || cA.residueIndex === cB.residueIndex) continue;
            const rA = unit.rings.all[cA.ringIndex];
            const rB = unit.rings.all[cB.ringIndex];

            if (IntAdjacencyGraph.areVertexSetsConnected(unit.bonds, rA, rB, 3)) {
                // TODO handle better, for now fix both directions as it is unclear where the C1 atom is
                //      would need to know the path connecting the two rings
                fixLinkDirection(i, j);
                fixLinkDirection(j, i);
                links.push({ carbohydrateIndexA: i, carbohydrateIndexB: j });
                links.push({ carbohydrateIndexA: j, carbohydrateIndexB: i });
            }
        }
    }

    // get carbohydrate links induced by inter-unit bonds, that is
    // inter monosaccharide links for structures from the
    // PDB archive __before__ carbohydrate remediation
    // plus terminal links for __before__ and __after__
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i];
        if (!Unit.isAtomic(unit)) continue;

        structure.interUnitBonds.getConnectedUnits(unit).forEach(pairBonds => {
            pairBonds.connectedIndices.forEach(indexA => {
                pairBonds.getEdges(indexA).forEach(({ props, indexB }) => {
                    if (!BondType.isCovalent(props.flag)) return;

                    const { unitA, unitB } = pairBonds;
                    const ringElementIndicesA = getRingElementIndices(unitA, indexA);
                    const ringElementIndicesB = getRingElementIndices(unitB, indexB);
                    if (ringElementIndicesA.length > 0 && ringElementIndicesB.length > 0) {
                        const lA = ringElementIndicesA.length;
                        const lB = ringElementIndicesB.length;
                        for (let j = 0, jl = Math.max(lA, lB); j < jl; ++j) {
                            const ringElementIndexA = ringElementIndicesA[Math.min(j, lA - 1)];
                            const ringElementIndexB = ringElementIndicesB[Math.min(j, lB - 1)];
                            const atomIdA = getAtomId(unitA, indexA);
                            if (atomIdA.startsWith('O1') || atomIdA.startsWith('C1')) {
                                fixLinkDirection(ringElementIndexA, ringElementIndexB);
                            }
                            links.push({
                                carbohydrateIndexA: ringElementIndexA,
                                carbohydrateIndexB: ringElementIndexB
                            });
                        }
                    } else if (ringElementIndicesB.length === 0) {
                        for (const ringElementIndexA of ringElementIndicesA) {
                            const atomIdA = getAtomId(unitA, indexA);
                            if (atomIdA.startsWith('O1') || atomIdA.startsWith('C1')) {
                                fixTerminalLinkDirection(ringElementIndexA, indexB, unitB);
                            }
                            terminalLinks.push({
                                carbohydrateIndex: ringElementIndexA,
                                elementIndex: indexB,
                                elementUnit: unitB,
                                fromCarbohydrate: true
                            });
                        }
                    } else if (ringElementIndicesA.length === 0) {
                        for (const ringElementIndexB of ringElementIndicesB) {
                            terminalLinks.push({
                                carbohydrateIndex: ringElementIndexB,
                                elementIndex: indexA,
                                elementUnit: unitA,
                                fromCarbohydrate: false
                            });
                        }
                    }
                });
            });
        });
    }

    return { links, terminalLinks, elements, partialElements, ...buildLookups(elements, links, terminalLinks) };
}

function buildLookups (elements: CarbohydrateElement[], links: CarbohydrateLink[], terminalLinks: CarbohydrateTerminalLink[]) {

    function key(unit: Unit, element: ElementIndex) {
        return cantorPairing(unit.id, element);
    }

    function getIndices(map: Map<number, number[]>, unit: Unit.Atomic, index: ElementIndex): ReadonlyArray<number> {
        const indices: number[] = [];
        const il = map.get(key(unit, index));
        if (il !== undefined) {
            for (const i of il) arraySetAdd(indices, i);
        }
        return indices;
    }

    // elements

    const elementsMap = new Map<number, number[]>();
    for (let i = 0, il = elements.length; i < il; ++i) {
        const { unit, ringIndex } = elements[i];
        const ring = unit.rings.all[ringIndex];
        for (let j = 0, jl = ring.length; j < jl; ++j) {
            const k = key(unit, unit.elements[ring[j]]);
            const e = elementsMap.get(k);
            if (e === undefined) elementsMap.set(k, [i]);
            else e.push(i);
        }
    }

    function getElementIndices(unit: Unit.Atomic, index: ElementIndex) {
        return getIndices(elementsMap, unit, index);
    }

    // links

    const linksMap = new Map<number, number[]>();
    for (let i = 0, il = links.length; i < il; ++i) {
        const l = links[i];
        const { unit, ringIndex } = elements[l.carbohydrateIndexA];
        const ring = unit.rings.all[ringIndex];
        for (let j = 0, jl = ring.length; j < jl; ++j) {
            const k = key(unit, unit.elements[ring[j]]);
            const e = linksMap.get(k);
            if (e === undefined) linksMap.set(k, [i]);
            else e.push(i);
        }
    }

    function getLinkIndices(unit: Unit.Atomic, index: ElementIndex) {
        return getIndices(linksMap, unit, index);
    }

    // terminal links

    const terminalLinksMap = new Map<number, number[]>();
    for (let i = 0, il = terminalLinks.length; i < il; ++i) {
        const { fromCarbohydrate, carbohydrateIndex, elementUnit, elementIndex } = terminalLinks[i];
        if (fromCarbohydrate) {
            const { unit, ringIndex } = elements[carbohydrateIndex];
            const ring = unit.rings.all[ringIndex];
            for (let j = 0, jl = ring.length; j < jl; ++j) {
                const k = key(unit, unit.elements[ring[j]]);
                const e = terminalLinksMap.get(k);
                if (e === undefined) terminalLinksMap.set(k, [i]);
                else e.push(i);
            }
        } else {
            const k = key(elementUnit, elementUnit.elements[elementIndex]);
            const e = terminalLinksMap.get(k);
            if (e === undefined) terminalLinksMap.set(k, [i]);
            else e.push(i);
        }
    }

    function getTerminalLinkIndices(unit: Unit.Atomic, index: ElementIndex) {
        return getIndices(terminalLinksMap, unit, index);
    }

    return { getElementIndices, getLinkIndices, getTerminalLinkIndices };
}