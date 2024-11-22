/**
 * Copyright (c) 2019-2024 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CustomPropertyDescriptor } from '../../../../mol-model/custom-property';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Column } from '../../../../mol-data/db';
import { FormatPropertyProvider } from '../../common/property';
import { BondType } from '../../../../mol-model/structure/model/types';
import { ElementIndex } from '../../../../mol-model/structure';

export type IndexPairsProps = {
    readonly key: ArrayLike<number>
    readonly operatorA: ArrayLike<number>
    readonly operatorB: ArrayLike<number>
    readonly order: ArrayLike<number>
    readonly distance: ArrayLike<number>
    readonly flag: ArrayLike<BondType.Flag>
}
export type IndexPairs = IntAdjacencyGraph<ElementIndex, IndexPairsProps>
export type IndexPairBonds = {
    bonds: IndexPairs
    /** Distance in Angstrom. If negative, element-based threshold is used. */
    maxDistance: number
    /** Can be cached in `ElementSetIntraBondCache` */
    cacheable: boolean
    /** Has operatorA & operatorB set for each bond */
    hasOperators: boolean
    /** Bonds with same operator grouped together */
    bySameOperator: Map<number, ArrayLike<number>>
}

function getGraph(indexA: ArrayLike<ElementIndex>, indexB: ArrayLike<ElementIndex>, props: Partial<IndexPairsProps>, count: number): IndexPairs {
    const builder = new IntAdjacencyGraph.EdgeBuilder(count, indexA, indexB);
    const key = new Int32Array(builder.slotCount);
    const operatorA = new Array(builder.slotCount);
    const operatorB = new Array(builder.slotCount);
    const order = new Int8Array(builder.slotCount);
    const distance = new Array(builder.slotCount);
    const flag = new Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(key, props.key ? props.key[i] : -1);
        builder.assignDirectedProperty(
            operatorA, props.operatorA ? props.operatorA[i] : -1,
            operatorB, props.operatorB ? props.operatorB[i] : -1
        );
        builder.assignProperty(order, props.order ? props.order[i] : 1);
        builder.assignProperty(distance, props.distance ? props.distance[i] : -1);
        builder.assignProperty(flag, props.flag ? props.flag[i] : BondType.Flag.Covalent);
    }

    return builder.createGraph({ key, operatorA, operatorB, order, distance, flag });
}

export namespace IndexPairBonds {
    export const Descriptor: CustomPropertyDescriptor = {
        name: 'index_pair_bonds',
    };

    export const Provider = FormatPropertyProvider.create<IndexPairBonds>(Descriptor, { asDynamic: true });

    export type Data = {
        pairs: {
            indexA: Column<number>,
            indexB: Column<number>,
            key?: Column<number>,
            /** Operator key for indexA. Used in bond computation. */
            operatorA?: Column<number>,
            /** Operator key for indexB. Used in bond computation. */
            operatorB?: Column<number>,
            order?: Column<number>,
            /**
             * Useful for bonds in periodic cells. That is, only bonds within the given
             * distance are added. This allows for bond between periodic image but
             * avoids unwanted bonds with wrong distances. If negative, test using the
             * `maxDistance` option from `Props`.
             */
            distance?: Column<number>,
            flag?: Column<BondType.Flag>,
        },
        count: number
    }

    export const DefaultProps = {
        /**
         * If negative, test using element-based threshold, otherwise distance in Angstrom.
         *
         * This option exists to handle bonds in periodic cells. For systems that are
         * made from beads (as opposed to atomic elements), set to a specific distance.
         *
         * Note that `Data` has a `distance` field which allows specifying a distance
         * for each bond individually which takes precedence over this option.
         */
        maxDistance: -1,
        /** Can be cached in `ElementSetIntraBondCache` */
        cacheable: true,
    };
    export type Props = typeof DefaultProps

    export function fromData(data: Data, props: Partial<Props> = {}): IndexPairBonds {
        const p = { ...DefaultProps, ...props };
        const { pairs, count } = data;
        const indexA = pairs.indexA.toArray() as ArrayLike<ElementIndex>;
        const indexB = pairs.indexB.toArray() as ArrayLike<ElementIndex>;
        const key = pairs.key && pairs.key.toArray();
        const operatorA = pairs.operatorA && pairs.operatorA.toArray();
        const operatorB = pairs.operatorB && pairs.operatorB.toArray();
        const order = pairs.order && pairs.order.toArray();
        const distance = pairs.distance && pairs.distance.toArray();
        const flag = pairs.flag && pairs.flag.toArray();

        let hasOperators = false;
        if (operatorA && operatorB) {
            hasOperators = true;
            for (let i = 0, il = operatorA.length; i < il; ++i) {
                if (operatorA[i] === -1 || operatorB[i] === -1) {
                    hasOperators = false;
                    break;
                }
            }
        }

        const bonds = getGraph(indexA, indexB, { key, operatorA, operatorB, order, distance, flag }, count);

        const bySameOperator = new Map<number, number[]>();
        if (hasOperators) {
            const { operatorA, operatorB } = bonds.edgeProps;
            for (let i = 0, il = operatorA.length; i < il; ++i) {
                if (operatorA[i] === operatorB[i]) {
                    const op = operatorA[i];
                    if (bySameOperator.has(op)) bySameOperator.get(op)!.push(i);
                    else bySameOperator.set(op, [i]);
                }
            }
        }

        return {
            bonds,
            maxDistance: p.maxDistance,
            cacheable: p.cacheable,
            hasOperators,
            bySameOperator,
        };
    }

    /** Like `getEdgeIndex` but taking `edgeProps.operatorA` and `edgeProps.operatorB` into account */
    export function getEdgeIndexForOperators(bonds: IndexPairs, i: ElementIndex, j: ElementIndex, opI: number, opJ: number): number {
        let a, b, opA, opB;
        if (i < j) {
            a = i; b = j;
            opA = opI; opB = opJ;
        } else {
            a = j; b = i;
            opA = opJ; opB = opI;
        }
        for (let t = bonds.offset[a], _t = bonds.offset[a + 1]; t < _t; t++) {
            if (bonds.b[t] === b && bonds.edgeProps.operatorA[t] === opA && bonds.edgeProps.operatorB[t] === opB) return t;
        }
        return -1;
    }
}