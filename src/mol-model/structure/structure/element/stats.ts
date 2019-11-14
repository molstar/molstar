/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../../../mol-data/int';
import Unit from '../unit';
import { Loci } from './loci';
import { Location } from './location';

export interface Stats {
    elementCount: number
    conformationCount: number
    residueCount: number
    unitCount: number

    firstElementLoc: Location
    firstConformationLoc: Location
    firstResidueLoc: Location
    firstUnitLoc: Location
}

export namespace Stats {
    export function create(): Stats {
        return {
            elementCount: 0,
            conformationCount: 0,
            residueCount: 0,
            unitCount: 0,

            firstElementLoc: Location.create(),
            firstConformationLoc: Location.create(),
            firstResidueLoc: Location.create(),
            firstUnitLoc: Location.create(),
        }
    }

    function handleElement(stats: Stats, element: Loci['elements'][0]) {
        const { indices, unit } = element
        const { elements } = unit
        const size = OrderedSet.size(indices)

        const lociResidueAltIdCounts = new Map<string, number>()
        const residueAltIdCounts = new Map<string, number>()
        const addCount = (map: Map<string, number>, altId: string) => {
            const count = map.get(altId) || 0
            map.set(altId, count + 1)
        }

        if (size > 0) {
            Location.set(stats.firstElementLoc, unit, elements[OrderedSet.start(indices)])
        }

        if (size === 1) {
            stats.elementCount += 1
            if (stats.elementCount === 1) {
                Location.set(stats.firstElementLoc, unit, elements[OrderedSet.start(indices)])
            }
        } else if (size === elements.length) {
            stats.unitCount += 1
            if (stats.unitCount === 1) {
                Location.set(stats.firstUnitLoc, unit, elements[OrderedSet.start(indices)])
            }
        } else {
            if (Unit.isAtomic(unit)) {
                const { index, offsets } = unit.model.atomicHierarchy.residueAtomSegments
                const { label_alt_id } = unit.model.atomicHierarchy.atoms;
                let i = 0
                while (i < size) {
                    lociResidueAltIdCounts.clear()
                    let j = 0
                    const eI = elements[OrderedSet.getAt(indices, i)]
                    const rI = index[eI]
                    addCount(lociResidueAltIdCounts, label_alt_id.value(eI))
                    ++i
                    ++j
                    while (i < size) {
                        const eI = elements[OrderedSet.getAt(indices, i)]
                        if (index[eI] !== rI) break
                        addCount(lociResidueAltIdCounts, label_alt_id.value(eI))
                        ++i
                        ++j
                    }

                    if (offsets[rI + 1] - offsets[rI] === j) {
                        // full residue
                        stats.residueCount += 1
                        if (stats.residueCount === 1) {
                            Location.set(stats.firstResidueLoc, unit, elements[offsets[rI]])
                        }
                    } else {
                        // partial residue
                        residueAltIdCounts.clear()
                        for (let l = offsets[rI], _l = offsets[rI + 1]; l < _l; ++l) {
                            addCount(residueAltIdCounts, label_alt_id.value(l))
                        }
                        // check if shared atom count match
                        if (residueAltIdCounts.get('') === lociResidueAltIdCounts.get('')) {
                            lociResidueAltIdCounts.forEach((v, k) => {
                                if (residueAltIdCounts.get(k) !== v) return
                                if (k !== '') {
                                    stats.conformationCount += 1
                                    if (stats.conformationCount === 1) {
                                        for (let l = offsets[rI], _l = offsets[rI + 1]; l < _l; ++l) {
                                            if (k === label_alt_id.value(l)) {
                                                Location.set(stats.firstConformationLoc, unit, l)
                                                break
                                            }
                                        }
                                    }
                                }
                                j -= v
                            })
                        }
                        stats.elementCount += j
                    }
                }
            } else {
                stats.elementCount += size
                if (stats.elementCount === 1) {
                    Location.set(stats.firstElementLoc, unit, elements[OrderedSet.start(indices)])
                }
            }
        }
    }

    export function ofLoci(loci: Loci) {
        const stats = create()
        if (!Loci.isEmpty(loci)) {
            for (const e of loci.elements) handleElement(stats, e)
        }
        return stats
    }

    /** Adds counts of two Stats objects together, assumes they describe different structures */
    export function add(out: Stats, a: Stats, b: Stats) {
        if (a.elementCount === 1 && b.elementCount === 0) {
            Location.copy(out.firstElementLoc, a.firstElementLoc)
        } else if (a.elementCount === 0 && b.elementCount === 1) {
            Location.copy(out.firstElementLoc, b.firstElementLoc)
        }

        if (a.conformationCount === 1 && b.conformationCount === 0) {
            Location.copy(out.firstConformationLoc, a.firstConformationLoc)
        } else if (a.conformationCount === 0 && b.conformationCount === 1) {
            Location.copy(out.firstConformationLoc, b.firstConformationLoc)
        }

        if (a.residueCount === 1 && b.residueCount === 0) {
            Location.copy(out.firstResidueLoc, a.firstResidueLoc)
        } else if (a.residueCount === 0 && b.residueCount === 1) {
            Location.copy(out.firstResidueLoc, b.firstResidueLoc)
        }

        if (a.unitCount === 1 && b.unitCount === 0) {
            Location.copy(out.firstUnitLoc, a.firstUnitLoc)
        } else if (a.unitCount === 0 && b.unitCount === 1) {
            Location.copy(out.firstUnitLoc, b.firstUnitLoc)
        }

        out.elementCount = a.elementCount + b.elementCount
        out.conformationCount = a.conformationCount + b.conformationCount
        out.residueCount = a.residueCount + b.residueCount
        out.unitCount = a.unitCount + b.unitCount
        return out
    }
}