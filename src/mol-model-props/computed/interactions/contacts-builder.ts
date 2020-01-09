/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Features } from './features';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { InteractionType, InteractionsIntraContacts, InteractionsInterContacts, InteractionFlag } from './common';
import { Unit } from '../../../mol-model/structure/structure';
import { UniqueArray } from '../../../mol-data/generic';
import { NumberArray } from '../../../mol-util/type-helpers';
import { IntMap } from '../../../mol-data/int';

export { IntraContactsBuilder }

interface IntraContactsBuilder {
    add: (indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) => void
    getContacts: () => InteractionsIntraContacts
}

namespace IntraContactsBuilder {
    export function create(features: Features, elementsCount: number): IntraContactsBuilder {
        const aIndices: Features.FeatureIndex[] = []
        const bIndices: Features.FeatureIndex[] = []
        const types: number[] = []

        return {
            add(indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) {
                aIndices[aIndices.length] = indexA
                bIndices[bIndices.length] = indexB
                types[types.length] = type
            },
            getContacts() {
                const builder = new IntAdjacencyGraph.EdgeBuilder(features.count, aIndices, bIndices)
                const type = new Int8Array(builder.slotCount) as ArrayLike<InteractionType>
                const flag = new Int8Array(builder.slotCount) as NumberArray
                for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
                    builder.addNextEdge()
                    builder.assignProperty(type, types[i])
                }
                const graph = builder.createGraph({ type, flag })

                let elementsIndex: InteractionsIntraContacts.ElementsIndex
                const contacts: InteractionsIntraContacts = Object.defineProperty(graph, 'elementsIndex', {
                    get: () => {
                        return elementsIndex || (elementsIndex = InteractionsIntraContacts.createElementsIndex(graph, features, elementsCount))
                    }
                })
                return contacts
            }
        }
    }
}

export { InterContactsBuilder }

interface InterContactsBuilder {
    startUnitPair: (unitA: Unit, unitB: Unit) => void
    finishUnitPair: () => void
    add: (indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) => void
    getContacts: (unitsFeatures: IntMap<Features>) => InteractionsInterContacts
}

namespace InterContactsBuilder {
    function addMapEntry<A, B>(map: Map<A, B[]>, a: A, b: B) {
        if (map.has(a)) map.get(a)!.push(b);
        else map.set(a, [b]);
    }

    export function create(): InterContactsBuilder {
        let uA: Unit
        let uB: Unit
        let mapAB: Map<number, InteractionsInterContacts.Info[]>
        let mapBA: Map<number, InteractionsInterContacts.Info[]>
        let linkedA: UniqueArray<Features.FeatureIndex, Features.FeatureIndex>
        let linkedB: UniqueArray<Features.FeatureIndex, Features.FeatureIndex>
        let linkCount: number

        const map = new Map<number, InteractionsInterContacts.Pair[]>();

        return {
            startUnitPair(unitA: Unit, unitB: Unit) {
                uA = unitA
                uB = unitB
                mapAB = new Map()
                mapBA = new Map()
                linkedA = UniqueArray.create()
                linkedB = UniqueArray.create()
                linkCount = 0
            },
            finishUnitPair() {
                if (linkCount === 0) return
                addMapEntry(map, uA.id, new InteractionsInterContacts.Pair(uA, uB, linkCount, linkedA.array, mapAB))
                addMapEntry(map, uB.id, new InteractionsInterContacts.Pair(uB, uA, linkCount, linkedB.array, mapBA))
            },
            add(indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) {
                addMapEntry(mapAB, indexA, { indexB, props: { type, flag: InteractionFlag.None } })
                addMapEntry(mapBA, indexB, { indexB: indexA, props: { type, flag: InteractionFlag.None } })
                UniqueArray.add(linkedA, indexA, indexA)
                UniqueArray.add(linkedB, indexB, indexB)
                linkCount += 1
            },
            getContacts(unitsFeatures: IntMap<Features>) {
                return new InteractionsInterContacts(map, unitsFeatures);
            }
        }
    }
}