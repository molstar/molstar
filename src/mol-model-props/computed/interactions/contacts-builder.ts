/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Features } from './features';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { InteractionType, InteractionsIntraContacts, InteractionsInterContacts, InteractionFlag } from './common';
import { Unit } from '../../../mol-model/structure/structure';
import { NumberArray } from '../../../mol-util/type-helpers';
import { IntMap } from '../../../mol-data/int';
import { InterUnitGraph } from '../../../mol-math/graph/inter-unit-graph';

export { IntraContactsBuilder };

interface IntraContactsBuilder {
    add: (indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) => void
    getContacts: () => InteractionsIntraContacts
}

namespace IntraContactsBuilder {
    export function create(features: Features, elementsCount: number): IntraContactsBuilder {
        const aIndices: Features.FeatureIndex[] = [];
        const bIndices: Features.FeatureIndex[] = [];
        const types: number[] = [];

        return {
            add(indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) {
                aIndices[aIndices.length] = indexA;
                bIndices[bIndices.length] = indexB;
                types[types.length] = type;
            },
            getContacts() {
                const builder = new IntAdjacencyGraph.EdgeBuilder(features.count, aIndices, bIndices);
                const type = new Int8Array(builder.slotCount) as ArrayLike<InteractionType>;
                const flag = new Int8Array(builder.slotCount) as NumberArray;
                for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
                    builder.addNextEdge();
                    builder.assignProperty(type, types[i]);
                }
                const graph = builder.createGraph({ type, flag });

                let elementsIndex: InteractionsIntraContacts.ElementsIndex;
                const contacts: InteractionsIntraContacts = Object.defineProperty(graph, 'elementsIndex', {
                    get: () => {
                        return elementsIndex || (elementsIndex = InteractionsIntraContacts.createElementsIndex(graph, features, elementsCount));
                    }
                });
                return contacts;
            }
        };
    }
}

export { InterContactsBuilder };

interface InterContactsBuilder {
    startUnitPair: (unitA: Unit, unitB: Unit) => void
    finishUnitPair: () => void
    add: (indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) => void
    getContacts: (unitsFeatures: IntMap<Features>) => InteractionsInterContacts
}

namespace InterContactsBuilder {
    export function create(): InterContactsBuilder {
        const builder = new InterUnitGraph.Builder<Unit, Features.FeatureIndex, InteractionsInterContacts.Props>();

        return {
            startUnitPair(unitA: Unit, unitB: Unit) {
                builder.startUnitPair(unitA, unitB);
            },
            finishUnitPair() {
                builder.finishUnitPair();
            },
            add(indexA: Features.FeatureIndex, indexB: Features.FeatureIndex, type: InteractionType) {
                builder.add(indexA, indexB, { type, flag: InteractionFlag.None });
            },
            getContacts(unitsFeatures: IntMap<Features>) {
                return new InteractionsInterContacts(builder.getMap(), unitsFeatures);
            }
        };
    }
}