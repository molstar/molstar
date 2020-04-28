/**
 * Copyright (c) 2019-2020 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor } from '../../../../mol-model/custom-property';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Column } from '../../../../mol-data/db';
import { FormatPropertyProvider } from '../../common/property';

export type IndexPairBondsProps = {
    readonly order: ArrayLike<number>
    readonly symmetryA: ArrayLike<string>
    readonly symmetryB: ArrayLike<string>
}
export type IndexPairBonds = IntAdjacencyGraph<number, IndexPairBondsProps>

function getGraph(indexA: ArrayLike<number>, indexB: ArrayLike<number>, props: Partial<IndexPairBondsProps>, count: number): IndexPairBonds {
    const builder = new IntAdjacencyGraph.EdgeBuilder(count, indexA, indexB);
    const order = new Int8Array(builder.slotCount);
    const symmetryA = new Array(builder.slotCount);
    const symmetryB = new Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(order, props.order ? props.order[i] : 1);
        builder.assignProperty(symmetryA, props.symmetryA ? props.symmetryA[i] : '');
        builder.assignProperty(symmetryB, props.symmetryB ? props.symmetryB[i] : '');
    }

    return builder.createGraph({ order, symmetryA, symmetryB });
}

export namespace IndexPairBonds {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'index_pair_bonds',
    };

    export const Provider = FormatPropertyProvider.create<IndexPairBonds>(Descriptor);

    export type Data = {
        pairs: {
            indexA: Column<number>,
            indexB: Column<number>
            order?: Column<number>,
            symmetryA?: Column<string>,
            symmetryB?: Column<string>,
        },
        count: number
    }

    export function fromData(data: Data) {
        const { pairs, count } = data;
        const indexA = pairs.indexA.toArray();
        const indexB = pairs.indexB.toArray();
        const order = pairs.order && pairs.order.toArray();
        const symmetryA = pairs.symmetryA && pairs.symmetryA.toArray();
        const symmetryB = pairs.symmetryB && pairs.symmetryB.toArray();
        return getGraph(indexA, indexB, { order, symmetryA, symmetryB }, count);
    }
}