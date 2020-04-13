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
});