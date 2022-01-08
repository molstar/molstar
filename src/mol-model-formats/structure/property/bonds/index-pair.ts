/**
 * Copyright (c) 2019-2022 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor } from '../../../../mol-model/custom-property';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Column } from '../../../../mol-data/db';
import { FormatPropertyProvider } from '../../common/property';
import { BondType } from '../../../../mol-model/structure/model/types';
import { ElementIndex } from '../../../../mol-model/structure';

export type IndexPairsProps = {
    readonly id: ArrayLike<number>
    readonly order: ArrayLike<number>
    readonly distance: ArrayLike<number>
    readonly flag: ArrayLike<BondType.Flag>
}
export type IndexPairs = IntAdjacencyGraph<ElementIndex, IndexPairsProps>
export type IndexPairBonds = { bonds: IndexPairs, maxDistance: number }

function getGraph(indexA: ArrayLike<ElementIndex>, indexB: ArrayLike<ElementIndex>, props: Partial<IndexPairsProps>, count: number): IndexPairs {
    const builder = new IntAdjacencyGraph.EdgeBuilder(count, indexA, indexB);
    const id = new Int32Array(builder.slotCount);
    const order = new Int8Array(builder.slotCount);
    const distance = new Array(builder.slotCount);
    const flag = new Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(id, props.id ? props.id[i] : -1);
        builder.assignProperty(order, props.order ? props.order[i] : 1);
        builder.assignProperty(distance, props.distance ? props.distance[i] : -1);
        builder.assignProperty(flag, props.flag ? props.flag[i] : BondType.Flag.Covalent);
    }

    return builder.createGraph({ id, order, distance, flag });
}

export namespace IndexPairBonds {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'index_pair_bonds',
    };

    export const Provider = FormatPropertyProvider.create<IndexPairBonds>(Descriptor);

    export type Data = {
        pairs: {
            indexA: Column<number>,
            indexB: Column<number>,
            id?: Column<number>,
            order?: Column<number>,
            /**
             * Useful for bonds in periodic cells. That is, only bonds within the given
             * distance are added. This allows for bond between periodic image but
             * avoids unwanted bonds with wrong distances.
             */
            distance?: Column<number>,
            flag?: Column<BondType.Flag>,
        },
        count: number
    }

    export const DefaultProps = {
        /**
         * If -1, test using element-based threshold, otherwise distance in Angstrom.
         *
         * This option exists to handle bonds in periodic cells. For systems that are
         * made from beads (as opposed to atomic elements), set to a spicific distance.
         */
        maxDistance: -1
    };
    export type Props = typeof DefaultProps

    export function fromData(data: Data, props: Partial<Props> = {}): IndexPairBonds {
        const p = { ...DefaultProps, ...props };
        const { pairs, count } = data;
        const indexA = pairs.indexA.toArray() as ArrayLike<ElementIndex>;
        const indexB = pairs.indexB.toArray() as ArrayLike<ElementIndex>;
        const id = pairs.id && pairs.id.toArray();
        const order = pairs.order && pairs.order.toArray();
        const distance = pairs.distance && pairs.distance.toArray();
        const flag = pairs.flag && pairs.flag.toArray();
        return {
            bonds: getGraph(indexA, indexB, { id, order, distance, flag }, count),
            maxDistance: p.maxDistance
        };
    }
}