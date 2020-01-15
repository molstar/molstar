/**
 * Copyright (c) 2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../../mol-model/structure/model/model'
import { CustomPropertyDescriptor } from '../../../../mol-model/structure';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Column } from '../../../../mol-data/db';

export type IndexPairBonds = IntAdjacencyGraph<number, { readonly order: ArrayLike<number> }>

function getGraph(indexA: ArrayLike<number>, indexB: ArrayLike<number>, _order: ArrayLike<number>, count: number): IndexPairBonds {
    const builder = new IntAdjacencyGraph.EdgeBuilder(count, indexA, indexB);
    const order = new Int8Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(order, _order[i]);
    }

    return builder.createGraph({ order });
}

export namespace IndexPairBonds {
    export const Descriptor: CustomPropertyDescriptor = {
        isStatic: true,
        name: 'index_pair_bonds',
    }

    export type Data = {
        pairs: {
            indexA: Column<number>,
            indexB: Column<number>
            order: Column<number>
        },
        count: number
    }

    export function attachFromData(model: Model, data: Data): boolean {
        if (model.customProperties.has(Descriptor)) return true;

        model.customProperties.add(Descriptor);
        model._staticPropertyData.__IndexPairBondsData__ = data;
        return true;
    }

    function getIndexPairBonds(model: Model) {
        return model._staticPropertyData.__IndexPairBondsData__ as Data;
    }

    export const PropName = '__IndexPairBonds__';
    export function get(model: Model): IndexPairBonds | undefined {
        if (model._staticPropertyData[PropName]) return model._staticPropertyData[PropName];
        if (!model.customProperties.has(Descriptor)) return void 0;

        const data = getIndexPairBonds(model);
        if (!data) return void 0;
        const { pairs, count } = data

        const indexA = pairs.indexA.toArray()
        const indexB = pairs.indexB.toArray()
        const order = pairs.order.toArray()

        const indexPairBonds = getGraph(indexA, indexB, order, count);
        model._staticPropertyData[PropName] = indexPairBonds;
        return indexPairBonds;
    }
}