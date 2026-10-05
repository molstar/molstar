/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { IntAdjacencyGraph } from '../int-adjacency-graph';

describe('IntGraph', () => {
    const vc = 3;
    const xs = [0, 1, 2];
    const ys = [1, 2, 0];
    const _prop = [10, 11, 12];

    const builder = new IntAdjacencyGraph.EdgeBuilder(vc, xs, ys);
    const prop: number[] = new Array(builder.slotCount);
    for (let i = 0; i < builder.edgeCount; i++) {
        builder.addNextEdge();
        builder.assignProperty(prop, _prop[i]);
    }
    const graph = builder.createGraph({ prop });

    it('triangle-edgeCount', () => expect(graph.edgeCount).toBe(3));
    it('triangle-vertexEdgeCounts', () => {
        expect(graph.getVertexEdgeCount(0)).toBe(2);
        expect(graph.getVertexEdgeCount(1)).toBe(2);
        expect(graph.getVertexEdgeCount(2)).toBe(2);
    });

    it('triangle-propAndEdgeIndex', () => {
        const prop = graph.edgeProps.prop;
        expect(prop[graph.getEdgeIndex(0, 1)]).toBe(10);
        expect(prop[graph.getEdgeIndex(1, 2)]).toBe(11);
        expect(prop[graph.getEdgeIndex(2, 0)]).toBe(12);
    });

    it('induce', () => {
        const induced = IntAdjacencyGraph.induceByVertices(graph, [1, 2]);
        expect(induced.vertexCount).toBe(2);
        expect(induced.edgeCount).toBe(1);
        expect(induced.edgeProps.prop[induced.getEdgeIndex(0, 1)]).toBe(11);
    });

    describe('connectedComponents', () => {
        function fromEdges(vertexCount: number, edges: number[][]) {
            const builder = new IntAdjacencyGraph.UniqueEdgeBuilder(vertexCount);
            for (const [i, j] of edges) builder.addEdge(i, j);
            return builder.getGraph();
        }

        /** Component ids must cover 0..componentCount-1 with no gaps, which is what callers rely on to size buffers. */
        function expectConsistentIndex(vertexCount: number, componentCount: number, componentIndex: Int32Array) {
            expect(componentIndex.length).toBe(vertexCount);

            const used = new Array(componentCount);
            let maxId = -1;
            for (let i = 0; i < vertexCount; i++) {
                const c = componentIndex[i];
                expect(c).toBeGreaterThanOrEqual(0);
                expect(c).toBeLessThan(componentCount);
                used[c] = true;
                if (c > maxId) maxId = c;
            }
            expect(maxId + 1).toBe(componentCount);
            for (let i = 0; i < componentCount; i++) expect(used[i]).toBe(true);
        }

        it('empty graph', () => {
            const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(fromEdges(0, []));
            expect(componentCount).toBe(0);
            expect(componentIndex.length).toBe(0);
        });

        it('no edges, each vertex is its own component', () => {
            const vertexCount = 5;
            const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(fromEdges(vertexCount, []));
            expect(componentCount).toBe(vertexCount);
            expectConsistentIndex(vertexCount, componentCount, componentIndex);
        });

        it('single fully connected component', () => {
            const vertexCount = 3;
            const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(fromEdges(vertexCount, [[0, 1], [1, 2], [0, 2]]));
            expect(componentCount).toBe(1);
            expectConsistentIndex(vertexCount, componentCount, componentIndex);
        });

        it('two disjoint triangles', () => {
            const vertexCount = 6;
            const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(fromEdges(vertexCount, [[0, 1], [1, 2], [0, 2], [3, 4], [4, 5], [3, 5]]));
            expect(componentCount).toBe(2);
            expectConsistentIndex(vertexCount, componentCount, componentIndex);

            expect(componentIndex[0]).toBe(componentIndex[2]);
            expect(componentIndex[3]).toBe(componentIndex[5]);
            expect(componentIndex[0]).not.toBe(componentIndex[3]);
        });

        it('two chains plus an isolated vertex', () => {
            const vertexCount = 6;
            const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(fromEdges(vertexCount, [[0, 1], [1, 2], [3, 4]]));
            expect(componentCount).toBe(3);
            expectConsistentIndex(vertexCount, componentCount, componentIndex);
        });

        it('many vertices, few components', () => {
            const chainLength = 33;
            const chainCount = 3;
            const edges: number[][] = [];
            for (let c = 0; c < chainCount; c++) {
                for (let i = 0; i < chainLength - 1; i++) {
                    edges.push([c * chainLength + i, c * chainLength + i + 1]);
                }
            }
            const vertexCount = chainCount * chainLength;

            const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(fromEdges(vertexCount, edges));
            expect(componentCount).toBe(chainCount);
            expectConsistentIndex(vertexCount, componentCount, componentIndex);
        });
    });
});