/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * based in part on NGL (https://github.com/arose/ngl)
 */

import { Interactions } from './interactions';
import { InteractionType, InteractionFlag } from './common';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Unit, Structure } from '../../../mol-model/structure';

export function refineInteractions(structure: Structure, interactions: Interactions) {
    refineHydrophobicContacts(structure, interactions)
}

/**
 * For atoms interacting with several atoms in the same residue
 * only the one with the closest distance is kept.
 */
function refineHydrophobicContacts(structure: Structure, interactions: Interactions) {
    const { contacts, unitsContacts, unitsFeatures } = interactions

    /* keep only closest contact between residues */
    const handleResidueContact = function (dist: number, i: number, key: string, map: Map<string, [number, number]>, set: (i: number) => void) {
        const [minDist, minIndex] = map.get(key) || [Infinity, -1]
        if (dist < minDist) {
            if (minIndex !== -1) set(i)
            map.set(key, [dist, i])
        } else {
            set(i)
        }
    }

    const pA = Vec3()
    const pB = Vec3()

    //

    const residueInterMap = new Map<string, [number, number]>()
    const setInterFiltered = (i: number) => contacts.edges[i].props.flag = InteractionFlag.Filtered

    for (let i = 0, il = contacts.edgeCount; i < il; ++i) {
        const e = contacts.edges[i]
        if (e.props.type !== InteractionType.Hydrophobic) continue

        const featureA = e.indexA
        const featureB = e.indexB

        const { offsets: offsetsA, members: membersA } = unitsFeatures.get(e.unitA.id)
        const { offsets: offsetsB, members: membersB } = unitsFeatures.get(e.unitB.id)

        const elementA = membersA[offsetsA[featureA]]
        const elementB = membersB[offsetsB[featureB]]
        const residueA = (e.unitA as Unit.Atomic).getResidueIndex(elementA)
        const residueB = (e.unitB as Unit.Atomic).getResidueIndex(elementB)

        e.unitA.conformation.position(e.unitA.elements[elementA], pA)
        e.unitB.conformation.position(e.unitB.elements[elementB], pB)
        const dist = Vec3.distance(pA, pB)

        const keyA = `${elementA}|${e.unitA.id}|${residueB}|${e.unitB.id}`
        const keyB = `${elementB}|${e.unitB.id}|${residueA}|${e.unitA.id}`

        handleResidueContact(dist, i, keyA, residueInterMap, setInterFiltered)
        handleResidueContact(dist, i, keyB, residueInterMap, setInterFiltered)
    }

    //

    const uc = unitsContacts.keys()

    while (true) {
        const { done, value } = uc.next();
        if (done) break;

        const contacts = unitsContacts.get(value)
        const { offsets, members } = unitsFeatures.get(value)
        const unit = structure.unitMap.get(value)

        const residueIntraMap = new Map<string, [number, number]>()
        const setIntraFiltered = (i: number) => contacts.edgeProps.flag[i] = InteractionFlag.Filtered

        for (let i = 0, il = contacts.edgeCount * 2; i < il; ++i) {
            if (contacts.edgeProps.type[i] !== InteractionType.Hydrophobic) continue

            const featureA = contacts.a[i]
            const featureB = contacts.b[i]
            const elementA = members[offsets[featureA]]
            const elementB = members[offsets[featureB]]
            const residueA = (unit as Unit.Atomic).getResidueIndex(elementA)
            const residueB = (unit as Unit.Atomic).getResidueIndex(elementB)

            unit.conformation.position(unit.elements[elementA], pA)
            unit.conformation.position(unit.elements[elementB], pB)
            const dist = Vec3.distance(pA, pB)

            const keyA = `${elementA}|${residueB}`
            const keyB = `${elementB}|${residueA}}`

            handleResidueContact(dist, i, keyA, residueIntraMap, setIntraFiltered)
            handleResidueContact(dist, i, keyB, residueIntraMap, setIntraFiltered)
        }
    }
}