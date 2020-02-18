/**
 * Copyright (c) 2019-2020 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor } from '../../../../mol-model/structure';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Column } from '../../../../mol-data/db';
import { FormatPropertyProvider } from '../../common/property';

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
        name: 'index_pair_bonds',
    }

    export const Provider = FormatPropertyProvider.create<IndexPairBonds>(Descriptor)

    export type Data = {
        pairs: {
            indexA: Column<number>,
            indexB: Column<number>
            order: Column<number>
        },
        count: number
    }

    export function fromData(data: Data) {
        const { pairs, count } = data
        const indexA = pairs.indexA.toArray()
        const indexB = pairs.indexB.toArray()
        const order = pairs.order.toArray()
        return getGraph(indexA, indexB, order, count);
    }
}