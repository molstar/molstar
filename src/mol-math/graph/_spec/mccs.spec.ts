/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { IntAdjacencyGraph } from '../int-adjacency-graph';
import { Mccs } from '../mccs';

type Edge = [number, number, number?]; // [a, b, order = 1]

/** Build an IntAdjacencyGraph with a per-edge `order` property from an undirected edge list. */
function graph(vertexCount: number, edges: Edge[]) {
    const xs = edges.map(e => e[0]);
    const ys = edges.map(e => e[1]);
    const builder = new IntAdjacencyGraph.EdgeBuilder(vertexCount, xs, ys);
    const order = new Int32Array(builder.slotCount);
    for (let i = 0; i < builder.edgeCount; i++) {
        builder.addNextEdge();
        builder.assignProperty(order, edges[i][2] ?? 1);
    }
    return builder.createGraph({ order });
}

function ring(n: number, order = 1): Edge[] {
    const edges: Edge[] = [];
    for (let i = 0; i < n; i++) edges.push([i, (i + 1) % n, order]);
    return edges;
}

function path(n: number): Edge[] {
    const edges: Edge[] = [];
    for (let i = 0; i + 1 < n; i++) edges.push([i, i + 1]);
    return edges;
}

describe('Mccs', () => {
    it('finds the maximum common connected subgraph (6-ring in ring vs ring+pendant)', () => {
        const a = graph(6, ring(6));
        const b = graph(7, [...ring(6), [0, 6]]); // 6-ring plus a pendant on vertex 0
        const cliques = Mccs.find(a, b);
        expect(cliques[0].length).toBe(6);
    });

    it('returns symmetry-degenerate maxima for an identical ring', () => {
        const a = graph(6, ring(6));
        const b = graph(6, ring(6));
        const cliques = Mccs.find(a, b);
        expect(cliques.length).toBeGreaterThan(1); // ring symmetry -> several equally-sized mappings
        expect(cliques[0].length).toBe(6);
    });

    it('vertexTest restricts which vertices may pair', () => {
        const labelsA = ['C', 'C', 'C', 'N'];
        const labelsB = ['C', 'C', 'C', 'O'];
        const a = graph(4, path(4));
        const b = graph(4, path(4));
        const cliques = Mccs.find(a, b, { vertexTest: (i, j) => labelsA[i] === labelsB[j] });
        expect(cliques[0].length).toBe(3); // shared C-C-C prefix; terminal N/O incompatible
    });

    it('edgeTest gates c-edges (mismatched edge labels prevent a connected match)', () => {
        const a = graph(6, ring(6, 2)); // "order 2"
        const b = graph(6, ring(6, 1)); // "order 1"
        const aOrder = a.edgeProps.order, bOrder = b.edgeProps.order;
        const cliques = Mccs.find(a, b, {
            edgeTest: (ea, eb) => aOrder[ea] === bOrder[eb],
            minMatchedVertices: 3
        });
        // no c-edges survive the edge test, so the only connected cliques are singletons
        expect(cliques).toHaveLength(0);
    });

    it('pathCutoff 0 disables d-edges, collapsing a ring to a single edge', () => {
        const a = graph(6, ring(6));
        const b = graph(6, ring(6));
        const cliques = Mccs.find(a, b, { pathCutoff: 0, minMatchedVertices: 1 });
        // without d-edges a clique must be a complete subgraph; a 6-cycle has no triangle
        expect(cliques[0].length).toBe(2);
    });

    it('does not bridge disconnected components (connected-only)', () => {
        const a = graph(4, [[0, 1], [2, 3]]);
        const b = graph(4, [[0, 1], [2, 3]]);
        const cliques = Mccs.find(a, b, { minMatchedVertices: 1 });
        // cross-component pairs are unreachable -> no d-edges across the gap -> max connected match is one edge
        expect(cliques[0].length).toBe(2);
    });

    it('is deterministic across runs', () => {
        const a = graph(6, ring(6));
        const b = graph(7, [...ring(6), [0, 6]]);
        expect(JSON.stringify(Mccs.find(a, b))).toEqual(JSON.stringify(Mccs.find(a, b)));
    });

    it('honors a tiny search budget and reports truncation via the status out-parameter', () => {
        const a = graph(6, ring(6));
        const b = graph(6, ring(6));

        const clipped: Mccs.Status = { truncated: false };
        expect(() => Mccs.find(a, b, { maxIterations: 1, minMatchedVertices: 1 }, clipped)).not.toThrow();
        expect(clipped.truncated).toBe(true); // budget hit -> best-so-far

        const full: Mccs.Status = { truncated: false };
        Mccs.find(a, b, { minMatchedVertices: 1 }, full);
        expect(full.truncated).toBe(false); // completed within budget
    });

    it('returns nothing below minMatchedVertices', () => {
        const a = graph(2, [[0, 1]]);
        const b = graph(2, [[0, 1]]);
        expect(Mccs.find(a, b, { minMatchedVertices: 3 })).toHaveLength(0);
    });
});