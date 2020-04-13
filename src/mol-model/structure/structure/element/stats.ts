/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '../../../../mol-data/int';
import Unit from '../unit';
import { Loci } from './loci';
import { Location } from './location';
import { ChainIndex } from '../../model/indexing';
import Structure from '../structure';

export interface Stats {
    elementCount: number
    conformationCount: number
    residueCount: number
    chainCount: number
    unitCount: number
    structureCount: number

    firstElementLoc: Location
    firstConformationLoc: Location
    firstResidueLoc: Location
    firstChainLoc: Location
    firstUnitLoc: Location
    firstStructureLoc: Location
}

export namespace Stats {
    export function create(): Stats {
        return {
            elementCount: 0,
            conformationCount: 0,
            residueCount: 0,
            chainCount: 0,
            unitCount: 0,
            structureCount: 0,

            firstElementLoc: Location.create(void 0),
            firstConformationLoc: Location.create(void 0),
            firstResidueLoc: Location.create(void 0),
            firstChainLoc: Location.create(void 0),
            firstUnitLoc: Location.create(void 0),
            firstStructureLoc: Location.create(void 0),
        };
    }

    function addCountHelper<K>(map: Map<K, number>, key: K, inc: number) {
        const count = map.get(key) || 0;
        map.set(key, count + inc);
    }

    function handleElement(stats: Stats, structure: Structure, element: Loci['elements'][0]) {
        const { indices, unit } = element;
        const { elements } = unit;
        const size = OrderedSet.size(indices);

        const lociResidueAltIdCounts = new Map<string, number>();
        const residueAltIdCounts = new Map<string, number>();

        if (size > 0) {
            Location.set(stats.firstElementLoc, structure, unit, elements[OrderedSet.start(indices)]);
        }

        // count single element unit as unit not element
        if (size === elements.length) {
            stats.unitCount += 1;
            if (stats.unitCount === 1) {
                Location.set(stats.firstUnitLoc, structure, unit, elements[OrderedSet.start(indices)]);
            }
        } else if (size === 1) {
            if (Unit.Traits.is(unit.traits, Unit.Trait.MultiChain)) {
                return;
            } else {
                stats.elementCount += 1;
                if (stats.elementCount === 1) {
                    Location.set(stats.firstElementLoc, structure, unit, elements[OrderedSet.start(indices)]);
                }
            }
        } else {
            if (Unit.isAtomic(unit)) {
                const { index, offsets } = unit.model.atomicHierarchy.residueAtomSegments;
                const { label_alt_id } = unit.model.atomicHierarchy.atoms;
                let i = 0;
                while (i < size) {
                    lociResidueAltIdCounts.clear();
                    let j = 0;
                    const eI = elements[OrderedSet.getAt(indices, i)];
                    const rI = index[eI];
                    addCountHelper(lociResidueAltIdCounts, label_alt_id.value(eI), 1);
                    ++i;
                    ++j;
                    while (i < size) {
                        const eI = elements[OrderedSet.getAt(indices, i)];
                        if (index[eI] !== rI) break;
                        addCountHelper(lociResidueAltIdCounts, label_alt_id.value(eI), 1);
                        ++i;
                        ++j;
                    }

                    if (offsets[rI + 1] - offsets[rI] === j) {
                        // full residue
                        stats.residueCount += 1;
                        if (stats.residueCount === 1) {
                            Location.set(stats.firstResidueLoc, structure, unit, offsets[rI]);
                        }
                    } else {
                        // partial residue
                        residueAltIdCounts.clear();
                        for (let l = offsets[rI], _l = offsets[rI + 1]; l < _l; ++l) {
                            addCountHelper(residueAltIdCounts, label_alt_id.value(l), 1);
                        }
                        // check if shared atom count match
                        if (residueAltIdCounts.get('') === lociResidueAltIdCounts.get('')) {
                            lociResidueAltIdCounts.forEach((v, k) => {
                                if (residueAltIdCounts.get(k) !== v) return;
                                if (k !== '') {
                                    stats.conformationCount += 1;
                                    if (stats.conformationCount === 1) {
                                        for (let l = offsets[rI], _l = offsets[rI + 1]; l < _l; ++l) {
                                            if (k === label_alt_id.value(l)) {
                                                Location.set(stats.firstConformationLoc, structure, unit, l);
                                                break;
                                            }
                                        }
                                    }
                                }
                                j -= v;
                            });
                        }
                        stats.elementCount += j;
                    }
                }
            } else {
                stats.elementCount += size;
                if (stats.elementCount === 1) {
                    Location.set(stats.firstElementLoc, structure, unit, elements[OrderedSet.start(indices)]);
                }
            }
        }
    }

    function handleUnitChainsSimple(stats: Stats, structure: Structure, element: Loci['elements'][0]) {
        const { indices, unit } = element;
        const size = OrderedSet.size(indices);
        if (size === 0) return;

        const { elements } = unit;

        if (!Unit.Traits.is(unit.traits, Unit.Trait.MultiChain)) {
            if (size === elements.length) {
                stats.chainCount += 1;
                if (stats.chainCount === 1) {
                    Location.set(stats.firstChainLoc, structure, unit, elements[OrderedSet.start(indices)]);
                }
            }
            return;
        }

        const segments = Unit.isAtomic(unit)
            ? unit.model.atomicHierarchy.chainAtomSegments
            : Unit.isSpheres(unit)
                ? unit.model.coarseHierarchy.spheres.chainElementSegments
                : Unit.isGaussians(unit)
                    ? unit.model.coarseHierarchy.gaussians.chainElementSegments
                    : void 0;

        if (!segments) {
            console.warn('StructureElement loci stats: unknown unit type');
            return;
        }

        const { index, offsets } = segments;
        let i = 0;
        while (i < size) {
            let j = 0;
            const eI = elements[OrderedSet.getAt(indices, i)];
            const cI = index[eI];
            ++i;
            ++j;
            while (i < size) {
                const eI = elements[OrderedSet.getAt(indices, i)];
                if (index[eI] !== cI) break;
                ++i;
                ++j;
            }

            if (offsets[cI + 1] - offsets[cI] === j) {
                // full chain
                stats.chainCount += 1;
                if (stats.chainCount === 1) {
                    Location.set(stats.firstChainLoc, structure, unit, offsets[cI]);
                }
            }
        }
    }

    function handleUnitChainsPartitioned(stats: Stats, structure: Structure, lociElements: Loci['elements'], start: number, end: number) {
        let element = lociElements[start];

        // all the elements have the same model since they are part of the same group so this is ok.
        const segments = Unit.isAtomic(element.unit)
            ? element.unit.model.atomicHierarchy.chainAtomSegments
            : Unit.isSpheres(element.unit)
                ? element.unit.model.coarseHierarchy.spheres.chainElementSegments
                : Unit.isGaussians(element.unit)
                    ? element.unit.model.coarseHierarchy.gaussians.chainElementSegments
                    : void 0;

        if (!segments) {
            console.warn('StructureElement loci stats: unknown unit type');
            return;
        }

        const { index, offsets } = segments;

        const chainCounts = new Map<ChainIndex, number>();

        for (let elIndex = start; elIndex < end; elIndex++) {
            element = lociElements[elIndex];

            const { indices, unit } = element;
            const size = OrderedSet.size(indices);
            if (size === 0) continue;

            const { elements } = unit;

            if (!Unit.Traits.is(unit.traits, Unit.Trait.MultiChain)) {
                const eI = elements[OrderedSet.start(indices)];
                addCountHelper(chainCounts, index[eI], elements.length);
                continue;
            }

            let i = 0;
            while (i < size) {
                let j = 0;
                const eI = elements[OrderedSet.getAt(indices, i)];
                const cI = index[eI];
                ++i;
                ++j;
                while (i < size) {
                    const eI = elements[OrderedSet.getAt(indices, i)];
                    if (index[eI] !== cI) break;
                    ++i;
                    ++j;
                }

                addCountHelper(chainCounts, cI, j);
            }
        }

        let firstCI = -1;
        chainCounts.forEach((count, cI) => {
            if (offsets[cI + 1] - offsets[cI] === count) {
                // full chain
                stats.chainCount += 1;
                if (stats.chainCount === 1) {
                    firstCI = cI;
                }
            }
        });

        if (firstCI < 0) return;

        for (let elIndex = start; elIndex < end; elIndex++) {
            element = lociElements[elIndex];

            const { indices, unit } = element;
            const size = OrderedSet.size(indices);
            if (size === 0) continue;

            const { elements } = unit;

            let i = 0;
            while (i < size) {
                const eI = elements[OrderedSet.getAt(indices, i)];
                const cI = index[eI];
                if (cI === firstCI) {
                    Location.set(stats.firstChainLoc, structure, unit, eI);
                    return;
                }
            }
        }
    }

    export function ofLoci(loci: Loci) {
        const stats = create();
        if (Loci.isEmpty(loci)) return stats;

        let hasPartitions = false;
        if (Loci.isWholeStructure(loci)) {
            stats.structureCount += 1;
            if (stats.structureCount === 1) {
                const { unit, indices } = loci.elements[0];
                Location.set(stats.firstStructureLoc, loci.structure, unit, unit.elements[OrderedSet.min(indices)]);
            }
        } else {
            for (const e of loci.elements) {
                handleElement(stats, loci.structure, e);
                if (!Unit.Traits.is(e.unit.traits, Unit.Trait.Partitioned)) {
                    handleUnitChainsSimple(stats, loci.structure, e);
                } else {
                    hasPartitions = true;
                }
            }
        }

        if (hasPartitions) {
            for (let i = 0, len = loci.elements.length; i < len; i++) {
                const e = loci.elements[i];
                if (!Unit.Traits.is(e.unit.traits, Unit.Trait.Partitioned)) continue;

                const start = i;
                while (i < len && Unit.areSameChainOperatorGroup(loci.elements[i].unit, e.unit)) {
                    i++;
                }
                const end = i;
                i--;
                if (end - start === 1) {
                    handleUnitChainsSimple(stats, loci.structure, e);
                } else {
                    handleUnitChainsPartitioned(stats, loci.structure, loci.elements, start, end);
                }
            }
        }

        return stats;
    }

    /** Adds counts of two Stats objects together, assumes they describe different structures */
    export function add(out: Stats, a: Stats, b: Stats) {
        if (a.elementCount === 1 && b.elementCount === 0) {
            Location.copy(out.firstElementLoc, a.firstElementLoc);
        } else if (a.elementCount === 0 && b.elementCount === 1) {
            Location.copy(out.firstElementLoc, b.firstElementLoc);
        }

        if (a.conformationCount === 1 && b.conformationCount === 0) {
            Location.copy(out.firstConformationLoc, a.firstConformationLoc);
        } else if (a.conformationCount === 0 && b.conformationCount === 1) {
            Location.copy(out.firstConformationLoc, b.firstConformationLoc);
        }

        if (a.residueCount === 1 && b.residueCount === 0) {
            Location.copy(out.firstResidueLoc, a.firstResidueLoc);
        } else if (a.residueCount === 0 && b.residueCount === 1) {
            Location.copy(out.firstResidueLoc, b.firstResidueLoc);
        }

        if (a.chainCount === 1 && b.chainCount === 0) {
            Location.copy(out.firstChainLoc, a.firstChainLoc);
        } else if (a.chainCount === 0 && b.chainCount === 1) {
            Location.copy(out.firstChainLoc, b.firstChainLoc);
        }

        if (a.unitCount === 1 && b.unitCount === 0) {
            Location.copy(out.firstUnitLoc, a.firstUnitLoc);
        } else if (a.unitCount === 0 && b.unitCount === 1) {
            Location.copy(out.firstUnitLoc, b.firstUnitLoc);
        }

        if (a.structureCount === 1 && b.structureCount === 0) {
            Location.copy(out.firstStructureLoc, a.firstStructureLoc);
        } else if (a.structureCount === 0 && b.structureCount === 1) {
            Location.copy(out.firstStructureLoc, b.firstStructureLoc);
        }

        out.elementCount = a.elementCount + b.elementCount;
        out.conformationCount = a.conformationCount + b.conformationCount;
        out.residueCount = a.residueCount + b.residueCount;
        out.chainCount = a.chainCount + b.chainCount;
        out.unitCount = a.unitCount + b.unitCount;
        out.structureCount = a.structureCount + b.structureCount;
        return out;
    }
}