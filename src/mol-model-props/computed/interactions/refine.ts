/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * based in part on NGL (https://github.com/arose/ngl)
 */

import { Interactions } from './interactions';
import { InteractionType, InteractionFlag, InteractionsIntraContacts, FeatureType } from './common';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Unit, Structure } from '../../../mol-model/structure';
import { Features } from './features';

interface ContactRefiner {
    isApplicable: (type: InteractionType) => boolean
    handleInterContact: (index: number, infoA: Features.Info, infoB: Features.Info) => void
    startUnit: (contacts: InteractionsIntraContacts) => void
    handleIntraContact: (index: number, infoA: Features.Info, infoB: Features.Info) => void
}

export function refineInteractions(structure: Structure, interactions: Interactions) {
    const { contacts, unitsContacts, unitsFeatures } = interactions

    const contactRefiners: ContactRefiner[] = [
        hydrophobicRefiner(structure, interactions),
        weakHydrogenBondsRefiner(structure, interactions),
    ]

    for (let i = 0, il = contacts.edgeCount; i < il; ++i) {
        const e = contacts.edges[i]

        const infoA = Features.Info(structure, e.unitA as Unit.Atomic, unitsFeatures.get(e.unitA.id))
        infoA.feature = e.indexA
        const infoB = Features.Info(structure, e.unitB as Unit.Atomic, unitsFeatures.get(e.unitB.id))
        infoB.feature = e.indexB

        for (const refiner of contactRefiners) {
            if (refiner.isApplicable(e.props.type)) refiner.handleInterContact(i, infoA, infoB)
        }
    }

    //

    const ucKeys = unitsContacts.keys()

    while (true) {
        const { done, value } = ucKeys.next();
        if (done) break;

        const contacts = unitsContacts.get(value)
        const features = unitsFeatures.get(value)
        const unit = structure.unitMap.get(value)
        if (!Unit.isAtomic(unit)) continue

        const infoA = Features.Info(structure, unit, features)
        const infoB = Features.Info(structure, unit, features)

        for (const refiner of contactRefiners) refiner.startUnit(contacts)

        for (let i = 0, il = contacts.edgeCount * 2; i < il; ++i) {
            infoA.feature = contacts.a[i]
            infoB.feature = contacts.b[i]

            for (const refiner of contactRefiners) {
                if (refiner.isApplicable(contacts.edgeProps.type[i])) refiner.handleIntraContact(i, infoA, infoB)
            }
        }
    }
}

/**
 * For atoms interacting with several atoms in the same residue
 * only the one with the closest distance is kept.
 */
function hydrophobicRefiner(structure: Structure, interactions: Interactions): ContactRefiner {
    const { contacts } = interactions

    /* keep only closest contact between residues */
    const handleResidueContact = function (dist: number, edge: number, key: string, map: Map<string, [number, number]>, set: (i: number) => void) {
        const [minDist, minIndex] = map.get(key) || [Infinity, -1]
        if (dist < minDist) {
            if (minIndex !== -1) set(minIndex)
            map.set(key, [dist, edge])
        } else {
            set(edge)
        }
    }

    const pA = Vec3()
    const pB = Vec3()

    function handleEdge(edge: number, infoA: Features.Info, infoB: Features.Info, map: Map<string, [number, number]>, set: (i: number) => void) {
        const elementA = infoA.members[infoA.offsets[infoA.feature]]
        const elementB = infoB.members[infoB.offsets[infoB.feature]]
        const residueA = infoA.unit.getResidueIndex(elementA)
        const residueB = infoB.unit.getResidueIndex(elementB)

        infoA.unit.conformation.position(infoA.unit.elements[elementA], pA)
        infoB.unit.conformation.position(infoB.unit.elements[elementB], pB)
        const dist = Vec3.distance(pA, pB)

        const keyA = `${elementA}|${infoA.unit.id}|${residueB}|${infoB.unit.id}|A`
        const keyB = `${elementB}|${infoB.unit.id}|${residueA}|${infoA.unit.id}|B`

        handleResidueContact(dist, edge, keyA, map, set)
        handleResidueContact(dist, edge, keyB, map, set)
    }

    const residueInterMap = new Map<string, [number, number]>()
    const setInterFiltered = (i: number) => contacts.edges[i].props.flag = InteractionFlag.Filtered

    let residueIntraMap: Map<string, [number, number]>
    let setIntraFiltered: (i: number) => void

    return {
        isApplicable: (type: InteractionType) => type === InteractionType.Hydrophobic,
        handleInterContact: (index: number, infoA: Features.Info, infoB: Features.Info) => {
            handleEdge(index, infoA, infoB, residueInterMap, setInterFiltered)
        },
        startUnit: (contacts: InteractionsIntraContacts) => {
            residueIntraMap = new Map<string, [number, number]>()
            setIntraFiltered = (i: number) => contacts.edgeProps.flag[i] = InteractionFlag.Filtered
        },
        handleIntraContact: (index: number, infoA: Features.Info, infoB: Features.Info) => {
            handleEdge(index, infoA, infoB, residueIntraMap, setIntraFiltered)
        }
    }
}

/**
 * Remove weak hydrogen bonds when the acceptor is involved in
 * a normal/strong hydrogen bond
 */
function weakHydrogenBondsRefiner(structure: Structure, interactions: Interactions): ContactRefiner {
    const { contacts } = interactions

    const hasHydrogenBond = (infoA: Features.Info, infoB: Features.Info) => {
        const acc = infoA.types[infoA.feature] === FeatureType.WeakHydrogenDonor ? infoB : infoA

        // check intra
        const eI = acc.members[acc.offsets[acc.feature]]
        const { offsets, indices } = interactions.unitsFeatures.get(acc.unit.id).elementsIndex
        const { type } = interactions.unitsContacts.get(acc.unit.id).edgeProps
        for (let i = offsets[eI], il = offsets[eI + 1]; i < il; ++i) {
            if (type[indices[i]] === InteractionType.HydrogenBond) return true
        }

        // check inter
        const interIndices = contacts.getEdgeIndices(acc.feature, acc.unit)
        for (let i = 0, il = interIndices.length; i < il; ++i) {
            if (contacts.edges[i].props.type === InteractionType.HydrogenBond) return true
        }

        return false
    }

    return {
        isApplicable: (type: InteractionType) => type === InteractionType.WeakHydrogenBond,
        handleInterContact: (index: number, infoA: Features.Info, infoB: Features.Info) => {
            if (hasHydrogenBond(infoA, infoB)) {
                contacts.edges[index].props.flag = InteractionFlag.Filtered
            }
        },
        startUnit: (contacts: InteractionsIntraContacts) => { },
        handleIntraContact: (index: number, infoA: Features.Info, infoB: Features.Info) => {
            if (hasHydrogenBond(infoA, infoB)) {
                const { flag } = interactions.unitsContacts.get(infoA.unit.id).edgeProps
                flag[index] = InteractionFlag.Filtered
            }
        }
    }
}