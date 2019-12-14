/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export { InterUnitGraph }

class InterUnitGraph<Unit extends InterUnitGraph.UnitBase, VertexIndex extends number, EdgeProps extends InterUnitGraph.EdgePropsBase = {}> {
    /** Number of inter-unit edges */
    readonly edgeCount: number
    /** Array of inter-unit edges */
    readonly edges: ReadonlyArray<InterUnitGraph.Edge<Unit, VertexIndex, EdgeProps>>
    private readonly edgeKeyIndex: Map<string, number>
    private readonly vertexKeyIndex: Map<string, number[]>

    /** Get an array of unit-pair-edges that are linked to the given unit */
    getLinkedUnits(unit: Unit): ReadonlyArray<InterUnitGraph.UnitPairEdges<Unit, VertexIndex, EdgeProps>> {
        if (!this.map.has(unit.id)) return emptyArray;
        return this.map.get(unit.id)!;
    }

    /** Index into this.edges */
    getEdgeIndex(indexA: VertexIndex, unitA: Unit, indexB: VertexIndex, unitB: Unit): number {
        const edgeKey = InterUnitGraph.getEdgeKey<Unit, VertexIndex>(indexA, unitA, indexB, unitB)
        const index = this.edgeKeyIndex.get(edgeKey)
        return index !== undefined ? index : -1
    }

    /** Check if edge exists */
    hasEdge(indexA: VertexIndex, unitA: Unit, indexB: VertexIndex, unitB: Unit): boolean {
        return this.getEdgeIndex(indexA, unitA, indexB, unitB) !== -1
    }

    /** Get inter-unit edge given a pair of indices and units */
    getEdge(indexA: VertexIndex, unitA: Unit, indexB: VertexIndex, unitB: Unit): InterUnitGraph.Edge<Unit, VertexIndex, EdgeProps> | undefined {
        const index = this.getEdgeIndex(indexA, unitA, indexB, unitB)
        return index !== -1 ? this.edges[index] : undefined
    }

    /** Indices into this.edges */
    getEdgeIndices(index: VertexIndex, unit: Unit): ReadonlyArray<number> {
        return this.vertexKeyIndex.get(InterUnitGraph.getVertexKey(index, unit)) || []
    }

    constructor(private map: Map<number, InterUnitGraph.UnitPairEdges<Unit, VertexIndex, EdgeProps>[]>) {
        let count = 0
        const edges: (InterUnitGraph.Edge<Unit, VertexIndex, EdgeProps>)[] = []
        const edgeKeyIndex = new Map<string, number>()
        const elementKeyIndex = new Map<string, number[]>()

        this.map.forEach(pairEdgesArray => {
            pairEdgesArray.forEach(pairEdges => {
                count += pairEdges.edgeCount
                pairEdges.connectedIndices.forEach(indexA => {
                    pairEdges.getEdges(indexA).forEach(edgeInfo => {
                        const { unitA, unitB } = pairEdges

                        const edgeKey = InterUnitGraph.getEdgeKey<Unit, VertexIndex>(indexA, unitA, edgeInfo.indexB, unitB)
                        edgeKeyIndex.set(edgeKey, edges.length)

                        const elementKey = InterUnitGraph.getVertexKey(indexA, unitA)
                        const e = elementKeyIndex.get(elementKey)
                        if (e === undefined) elementKeyIndex.set(elementKey, [edges.length])
                        else e.push(edges.length)

                        edges.push({ ...edgeInfo, indexA, unitA, unitB })
                    })
                })
            })
        })

        this.edgeCount = count
        this.edges = edges
        this.edgeKeyIndex = edgeKeyIndex
        this.vertexKeyIndex = elementKeyIndex
    }
}

namespace InterUnitGraph {
    export class UnitPairEdges<Unit extends UnitBase, VertexIndex extends number, EdgeProps extends EdgePropsBase = {}> {
        hasEdges(indexA: VertexIndex) {
            return this.edgeMap.has(indexA);
        }

        getEdges(indexA: VertexIndex): ReadonlyArray<EdgeInfo<VertexIndex, EdgeProps>> {
            if (!this.edgeMap.has(indexA)) return emptyArray;
            return this.edgeMap.get(indexA)!;
        }

        get areUnitsOrdered() {
            return this.unitA.id < this.unitB.id;
        }

        constructor(public unitA: Unit, public unitB: Unit,
            public edgeCount: number, public connectedIndices: ReadonlyArray<VertexIndex>,
            private edgeMap: Map<number, EdgeInfo<VertexIndex, EdgeProps>[]>) {
        }
    }

    export type UnitBase = { id: number }
    export type EdgePropsBase = { [name: string]: any }

    export interface EdgeInfo<VertexIndex extends number, EdgeProps extends EdgePropsBase = {}> {
        /** indexInto */
        readonly indexB: VertexIndex,
        readonly props: EdgeProps
    }

    export interface Edge<Unit extends UnitBase, VertexIndex extends number, EdgeProps extends EdgePropsBase = {}> {
        readonly unitA: Unit,
        readonly unitB: Unit,
        readonly indexA: VertexIndex,
        readonly indexB: VertexIndex,
        readonly props: EdgeProps
    }

    export function getEdgeKey<Unit extends UnitBase, VertexIndex extends number>(indexA: VertexIndex, unitA: Unit, indexB: VertexIndex, unitB: Unit) {
        return `${indexA}|${unitA.id}|${indexB}|${unitB.id}`
    }

    export function getVertexKey<Unit extends UnitBase, VertexIndex extends number>(index: VertexIndex, unit: Unit) {
        return `${index}|${unit.id}`
    }
}

const emptyArray: any[] = [];