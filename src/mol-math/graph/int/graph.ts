/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Represent a graph using vertex adjacency list.
 *
 * Edges of the i-th vertex are stored in the arrays a and b
 * for indices in the range [offset[i], offset[i+1]).
 *
 * Edge properties are indexed same as in the arrays a and b.
 */
type IntGraph<EdgeProperties extends object = { }> = {
    readonly offset: ArrayLike<number>,
    readonly a: ArrayLike<number>,
    readonly b: ArrayLike<number>,
    readonly vertexCount: number,
    readonly edgeCount: number,

    /**
     * Get the edge index between i-th and j-th vertex.
     * -1 if the edge does not exist.
     *
     * Because the a and b arrays contains each edge twice,
     * this always returns the smaller of the indices.
     */
    getEdgeIndex(i: number, j: number): number,
    getVertexEdgeCount(i: number): number
} & EdgeProperties

namespace IntGraph {
    class Impl implements IntGraph<any> {
        readonly vertexCount: number;

        getEdgeIndex(i: number, j: number): number {
            let a, b;
            if (i < j) { a = i; b = j; }
            else { a = j; b = i; }
            for (let t = this.offset[a], _t = this.offset[a + 1]; t < _t; t++) {
                if (this.b[t] === b) return t;
            }
            return -1;
        }

        getVertexEdgeCount(i: number): number {
            return this.offset[i + 1] - this.offset[i];
        }

        constructor(public offset: ArrayLike<number>, public a: ArrayLike<number>, public b: ArrayLike<number>, public edgeCount: number, props?: any) {
            this.vertexCount = offset.length - 1;
            if (props) {
                for (const p of Object.keys(props)) {
                    (this as any)[p] = props[p];
                }
            }
        }
    }

    export function create<EdgeProps extends object = { }>(offset: ArrayLike<number>, a: ArrayLike<number>, b: ArrayLike<number>, edgeCount: number, edgeProps?: EdgeProps): IntGraph<EdgeProps> {
        return new Impl(offset, a, b, edgeCount, edgeProps) as IntGraph<EdgeProps>;
    }

    export class EdgeBuilder {
        private bucketFill: Int32Array;
        private current = 0;
        private curA: number = 0;
        private curB: number = 0;

        offsets: Int32Array;
        edgeCount: number;
        /** the size of the A and B arrays */
        slotCount: number;
        a: Int32Array;
        b: Int32Array;

        createGraph<EdgeProps extends object = { }>(edgeProps?: EdgeProps) {
            return create(this.offsets, this.a, this.b, this.edgeCount, edgeProps);
        }

        /**
         * @example
         *   const property = new Int32Array(builder.slotCount);
         *   for (let i = 0; i < builder.edgeCount; i++) {
         *     builder.addNextEdge();
         *     builder.assignProperty(property, srcProp[i]);
         *   }
         * return builder.createGraph({ property });
         */
        addNextEdge() {
            const a = this.xs[this.current], b = this.ys[this.current];

            const oa = this.offsets[a] + this.bucketFill[a];
            const ob = this.offsets[b] + this.bucketFill[b];

            this.a[oa] = a;
            this.b[oa] = b;
            this.bucketFill[a]++;

            this.a[ob] = b;
            this.b[ob] = a;
            this.bucketFill[b]++;

            this.current++;
            this.curA = oa;
            this.curB = ob;
        }

        assignProperty<T>(prop: { [i: number]: T }, value: T) {
            prop[this.curA] = value;
            prop[this.curB] = value;
        }

        constructor(public vertexCount: number, public xs: ArrayLike<number>, public ys: ArrayLike<number>) {
            this.edgeCount = xs.length;
            this.offsets = new Int32Array(this.vertexCount + 1);
            this.bucketFill = new Int32Array(this.vertexCount);

            const bucketSizes = new Int32Array(this.vertexCount);
            for (let i = 0, _i = this.xs.length; i < _i; i++) bucketSizes[this.xs[i]]++;
            for (let i = 0, _i = this.ys.length; i < _i; i++) bucketSizes[this.ys[i]]++;

            let offset = 0;
            for (let i = 0; i < this.vertexCount; i++) {
                this.offsets[i] = offset;
                offset += bucketSizes[i];
            }
            this.offsets[this.vertexCount] = offset;
            this.slotCount = offset;
            this.a = new Int32Array(offset);
            this.b = new Int32Array(offset);
        }
    }
}

export { IntGraph }