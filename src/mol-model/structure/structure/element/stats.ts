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
    residueCount: number
    unitCount: number

    firstElementLoc: Location
    firstResidueLoc: Location
    firstUnitLoc: Location
}

export namespace Stats {
    export function create(): Stats {
        return {
            elementCount: 0,
            residueCount: 0,
            unitCount: 0,

            firstElementLoc: Location.create(),
            firstResidueLoc: Location.create(),
            firstUnitLoc: Location.create(),
        }
    }

    function handleElement(stats: Stats, element: Loci['elements'][0]) {
        const { indices, unit } = element
        const { elements } = unit
        const size = OrderedSet.size(indices)
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
                let i = 0
                while (i < size) {
                    const eI = elements[OrderedSet.getAt(indices, i)]
                    const rI = index[eI]
                    while (i < size && index[elements[OrderedSet.getAt(indices, i)]] === rI) {
                        ++i
                    }

                    if (offsets[rI + 1] - offsets[rI] === i) {
                        // full residue
                        stats.residueCount += 1
                        if (stats.residueCount === 1) {
                            Location.set(stats.firstResidueLoc, unit, elements[OrderedSet.start(indices)])
                        }
                    } else {
                        // partial residue
                        stats.elementCount += i
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

    export function add(out: Stats, a: Stats, b: Stats) {
        if (a.elementCount === 1 && b.elementCount === 0) {
            Location.copy(out.firstElementLoc, a.firstElementLoc)
        } else if (a.elementCount === 0 && b.elementCount === 1) {
            Location.copy(out.firstElementLoc, b.firstElementLoc)
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
        out.residueCount = a.residueCount + b.residueCount
        out.unitCount = a.unitCount + b.unitCount
        return out
    }
}