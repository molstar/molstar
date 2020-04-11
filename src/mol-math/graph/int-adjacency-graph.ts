/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { arrayPickIndices, cantorPairing } from '../../mol-data/util';
import { LinkedIndex, SortedArray } from '../../mol-data/int';
import { AssignableArrayLike } from '../../mol-util/type-helpers';

/**
 * Represent a graph using vertex adjacency list.
 *
 * Edges of the i-th vertex are stored in the arrays a and b
 * for indices in the range [offset[i], offset[i+1]).
 *
 * Edge properties are indexed same as in the arrays a and b.
 */
export interface IntAdjacencyGraph<VertexIndex extends number, EdgeProps extends IntAdjacencyGraph.EdgePropsBase> {
    readonly offset: ArrayLike<number>,
    readonly a: ArrayLike<VertexIndex>,
    readonly b: ArrayLike<VertexIndex>,
    readonly vertexCount: number,
    readonly edgeCount: number,
    readonly edgeProps: Readonly<EdgeProps>

    /**
     * Get the edge index between i-th and j-th vertex.
     * -1 if the edge does not exist.
     *
     * Because the a and b arrays contains each edge twice,
     * this always returns the smaller of the indices.
     *
     * `getEdgeIndex(i, j) === getEdgeIndex(j, i)`
     */
    getEdgeIndex(i: VertexIndex, j: VertexIndex): number,
    /**
     * Get the edge index between i-th and j-th vertex.
     * -1 if the edge does not exist.
     *
     * `getEdgeIndex(i, j) !== getEdgeIndex(j, i)`
     */
    getDirectedEdgeIndex(i: VertexIndex, j: VertexIndex): number,
    getVertexEdgeCount(i: VertexIndex): number
}

export namespace IntAdjacencyGraph {
    export type EdgePropsBase = { [name: string]: ArrayLike<any> }

    class IntGraphImpl<VertexIndex extends number, EdgeProps extends IntAdjacencyGraph.EdgePropsBase> implements IntAdjacencyGraph<VertexIndex, EdgeProps> {
        readonly vertexCount: number;
        readonly edgeProps: EdgeProps;

        getEdgeIndex(i: VertexIndex, j: VertexIndex): number {
            let a, b;
            if (i < j) {
                a = i; b = j;
            } else {
                a = j; b = i;
            }
            for (let t = this.offset[a], _t = this.offset[a + 1]; t < _t; t++) {
                if (this.b[t] === b) return t;
            }
            return -1;
        }

        getDirectedEdgeIndex(i: VertexIndex, j: VertexIndex): number {
            for (let t = this.offset[i], _t = this.offset[i + 1]; t < _t; t++) {
                if (this.b[t] === j) return t;
            }
            return -1;
        }

        getVertexEdgeCount(i: VertexIndex): number {
            return this.offset[i + 1] - this.offset[i];
        }

        constructor(public offset: ArrayLike<number>, public a: ArrayLike<VertexIndex>, public b: ArrayLike<VertexIndex>, public edgeCount: number, edgeProps?: EdgeProps) {
            this.vertexCount = offset.length - 1;
            this.edgeProps = (edgeProps || {}) as EdgeProps;
        }
    }

    export function create<VertexIndex extends number, EdgeProps extends IntAdjacencyGraph.EdgePropsBase>(offset: ArrayLike<number>, a: ArrayLike<VertexIndex>, b: ArrayLike<VertexIndex>, edgeCount: number, edgeProps?: EdgeProps): IntAdjacencyGraph<VertexIndex, EdgeProps> {
        return new IntGraphImpl(offset, a, b, edgeCount, edgeProps) as IntAdjacencyGraph<VertexIndex, EdgeProps>;
    }

    export class EdgeBuilder<VertexIndex extends number> {
        private bucketFill: Int32Array;
        private current = 0;
        private curA: number = 0;
        private curB: number = 0;

        offsets: Int32Array;
        edgeCount: number;
        /** the size of the A and B arrays */
        slotCount: number;
        a: AssignableArrayLike<VertexIndex>;
        b: AssignableArrayLike<VertexIndex>;

        createGraph<EdgeProps extends IntAdjacencyGraph.EdgePropsBase>(edgeProps: EdgeProps) {
            return create<VertexIndex, EdgeProps>(this.offsets, this.a, this.b, this.edgeCount, edgeProps);
        }

        /**
         * @example
         *   const property = new Int32Array(builder.slotCount);
         *   for (let i = 0; i < builder.edgeCount; i++) {
         *     builder.addNextEdge();
         *     builder.assignProperty(property, srcProp[i]);
         *   }
         *   return builder.createGraph({ property });
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

        /** Builds property-less graph */
        addAllEdges() {
            for (let i = 0; i < this.edgeCount; i++) {
                this.addNextEdge();
            }
        }

        assignProperty<T>(prop: { [i: number]: T }, value: T) {
            prop[this.curA] = value;
            prop[this.curB] = value;
        }

        constructor(public vertexCount: number, public xs: ArrayLike<VertexIndex>, public ys: ArrayLike<VertexIndex>) {
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
            this.a = new Int32Array(offset) as unknown as AssignableArrayLike<VertexIndex>;
            this.b = new Int32Array(offset) as unknown as AssignableArrayLike<VertexIndex>;
        }
    }

    export class DirectedEdgeBuilder<VertexIndex extends number> {
        private bucketFill: Int32Array;
        private current = 0;
        private curA: number = 0;

        offsets: Int32Array;
        edgeCount: number;
        /** the size of the A and B arrays */
        slotCount: number;
        a: Int32Array;
        b: Int32Array;

        createGraph<EdgeProps extends IntAdjacencyGraph.EdgePropsBase>(edgeProps: EdgeProps) {
            return create(this.offsets, this.a, this.b, this.edgeCount, edgeProps);
        }

        /**
         * @example
         *   const property = new Int32Array(builder.slotCount);
         *   for (let i = 0; i < builder.edgeCount; i++) {
         *     builder.addNextEdge();
         *     builder.assignProperty(property, srcProp[i]);
         *   }
         *   return builder.createGraph({ property });
         */
        addNextEdge() {
            const a = this.xs[this.current], b = this.ys[this.current];

            const oa = this.offsets[a] + this.bucketFill[a];

            this.a[oa] = a;
            this.b[oa] = b;
            this.bucketFill[a]++;

            this.current++;
            this.curA = oa;
        }

        /** Builds property-less graph */
        addAllEdges() {
            for (let i = 0; i < this.edgeCount; i++) {
                this.addNextEdge();
            }
        }

        assignProperty<T>(prop: { [i: number]: T }, value: T) {
            prop[this.curA] = value;
        }

        constructor(public vertexCount: number, public xs: ArrayLike<VertexIndex>, public ys: ArrayLike<VertexIndex>) {
            this.edgeCount = xs.length;
            this.offsets = new Int32Array(this.vertexCount + 1);
            this.bucketFill = new Int32Array(this.vertexCount);

            const bucketSizes = new Int32Array(this.vertexCount);
            for (let i = 0, _i = this.xs.length; i < _i; i++) bucketSizes[this.xs[i]]++;

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

    export class UniqueEdgeBuilder<VertexIndex extends number> {
        private xs: VertexIndex[] = [];
        private ys: VertexIndex[] = [];
        private included = new Set<number>();

        addEdge(i: VertexIndex, j: VertexIndex) {
            let u = i, v = j;
            if (i > j) { u = j; v = i; }
            const id = cantorPairing(u, v);
            if (this.included.has(id)) return false;
            this.included.add(id);
            this.xs[this.xs.length] = u;
            this.ys[this.ys.length] = v;
            return true;
        }

        getGraph(): IntAdjacencyGraph<VertexIndex, {}> {
            return fromVertexPairs(this.vertexCount, this.xs, this.ys);
        }

        // if we cant to add custom props as well
        getEdgeBuiler() {
            return new EdgeBuilder(this.vertexCount, this.xs, this.ys);
        }

        constructor(public vertexCount: number) {
        }
    }

    export function fromVertexPairs<V extends number>(vertexCount: number, xs: V[], ys: V[]) {
        const graphBuilder = new IntAdjacencyGraph.EdgeBuilder(vertexCount, xs, ys);
        graphBuilder.addAllEdges();
        return graphBuilder.createGraph({});
    }

    export function induceByVertices<V extends number, P extends IntAdjacencyGraph.EdgePropsBase>(graph: IntAdjacencyGraph<V, P>, vertexIndices: ArrayLike<number>): IntAdjacencyGraph<V, P> {
        const { b, offset, vertexCount, edgeProps } = graph;
        const vertexMap = new Int32Array(vertexCount);
        for (let i = 0, _i = vertexIndices.length; i < _i; i++) vertexMap[vertexIndices[i]] = i + 1;

        let newEdgeCount = 0;
        for (let i = 0; i < vertexCount; i++) {
            if (vertexMap[i] === 0) continue;
            for (let j = offset[i], _j = offset[i + 1]; j < _j; j++) {
                if (b[j] > i && vertexMap[b[j]] !== 0) newEdgeCount++;
            }
        }

        const newOffsets = new Int32Array(vertexIndices.length + 1);
        const edgeIndices = new Int32Array(2 * newEdgeCount);
        const newA = new Int32Array(2 * newEdgeCount) as unknown as AssignableArrayLike<V>;
        const newB = new Int32Array(2 * newEdgeCount) as unknown as AssignableArrayLike<V>;
        let eo = 0, vo = 0;
        for (let i = 0; i < vertexCount; i++) {
            if (vertexMap[i] === 0) continue;
            const aa = vertexMap[i] - 1;
            for (let j = offset[i], _j = offset[i + 1]; j < _j; j++) {
                const bb = vertexMap[b[j]];
                if (bb === 0) continue;

                newA[eo] = aa as V;
                newB[eo] = bb - 1 as V;
                edgeIndices[eo] = j;
                eo++;
            }
            newOffsets[++vo] = eo;
        }

        const newEdgeProps = {} as P;
        for (const key of Object.keys(edgeProps) as (keyof P)[]) {
            newEdgeProps[key] = arrayPickIndices(edgeProps[key], edgeIndices) as P[keyof P];
        }

        return create(newOffsets, newA, newB, newEdgeCount, newEdgeProps);
    }

    export function connectedComponents(graph: IntAdjacencyGraph<any, any>): { componentCount: number, componentIndex: Int32Array } {
        const vCount = graph.vertexCount;

        if (vCount === 0) return { componentCount: 0, componentIndex: new Int32Array(0) };
        if (graph.edgeCount === 0) {
            const componentIndex = new Int32Array(vCount);
            for (let i = 0, _i = vCount; i < _i; i++) {
                componentIndex[i] = i;
            }
            return { componentCount: vCount, componentIndex };
        }

        const componentIndex = new Int32Array(vCount);
        for (let i = 0, _i = vCount; i < _i; i++) componentIndex[i] = -1;
        let currentComponent = 0;
        componentIndex[0] = currentComponent;

        const { offset, b: neighbor } = graph;
        const stack = [0];
        const list = LinkedIndex(vCount);
        list.remove(0);

        while (stack.length > 0) {
            const v = stack.pop()!;
            const cIdx = componentIndex[v];

            for (let eI = offset[v], _eI = offset[v + 1]; eI < _eI; eI++) {
                const n = neighbor[eI];
                if (!list.has(n)) continue;
                list.remove(n);
                stack.push(n);
                componentIndex[n] = cIdx;
            }

            // check if we visited all vertices.
            // If not, create a new component and continue.
            if (stack.length === 0 && list.head >= 0) {
                stack.push(list.head);
                componentIndex[list.head] = ++currentComponent;
                list.remove(list.head);
            }
        }

        return { componentCount: vCount, componentIndex };
    }

    /**
     * Check if any vertex in `verticesA` is connected to any vertex in `verticesB`
     * via at most `maxDistance` edges.
     *
     * Returns true if verticesA and verticesB are intersecting.
     */
    export function areVertexSetsConnected(graph: IntAdjacencyGraph<any, any>, verticesA: SortedArray<number>, verticesB: SortedArray<number>, maxDistance: number): boolean {
        // check if A and B are intersecting, this handles maxDistance = 0
        if (SortedArray.areIntersecting(verticesA, verticesB)) return true;
        if (maxDistance < 1) return false;

        const visited = new Set<number>();
        for (let i = 0, il = verticesA.length; i < il; ++i) {
            visited.add(verticesA[i]);
        }

        return areVertexSetsConnectedImpl(graph, verticesA, verticesB, maxDistance, visited);
    }
}

function areVertexSetsConnectedImpl(graph: IntAdjacencyGraph<any, any>, frontier: ArrayLike<number>, target: SortedArray<number>, distance: number, visited: Set<number>): boolean {
    const { b: neighbor, offset } = graph;
    const newFrontier: number[] = [];

    for (let i = 0, il = frontier.length; i < il; ++i) {
        const src = frontier[i];

        for (let j = offset[src], jl = offset[src + 1]; j < jl; ++j) {
            const other = neighbor[j];
            if (visited.has(other)) continue;
            if (SortedArray.has(target, other)) return true;

            visited.add(other);
            newFrontier[newFrontier.length] = other;
        }
    }

    return distance > 1 ? areVertexSetsConnectedImpl(graph, newFrontier, target, distance - 1, visited) : false;
}