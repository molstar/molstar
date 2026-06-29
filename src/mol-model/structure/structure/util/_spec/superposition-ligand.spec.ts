/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { ElementIndex } from '../../../model';
import { BondType } from '../../../model/types';
import {
    DefaultLigandMccsOptions, findMccsCliques, findMccsPairs,
    LigandGraph, LigandGraphVertex
} from '../superposition-ligand';

/**
 * The MCCS matcher operates purely on a LigandGraph's `vertices` and `bonds`, so it can be exercised
 * with hand-built graphs (no Structure/Unit needed). The subsequent superpose/RMSD step reuses the
 * already-tested `MinimizeRmsd` machinery and is best covered by an in-app integration test against
 * a real structure.
 */

type Bond = [number, number, number?, number?]; // [a, b, order = 1, flags = 0]

function makeGraph(elements: string[], bonds: Bond[]): LigandGraph {
    const vertices: LigandGraphVertex[] = elements.map((elementSymbol, index) => ({
        index,
        element: index as unknown as ElementIndex,
        atomName: `${elementSymbol}${index}`,
        elementSymbol
    }));

    const xs = bonds.map(e => e[0]);
    const ys = bonds.map(e => e[1]);
    const builder = new IntAdjacencyGraph.EdgeBuilder(elements.length, xs, ys);
    const order = new Int32Array(builder.slotCount);
    const flags = new Int32Array(builder.slotCount);
    for (let i = 0; i < builder.edgeCount; i++) {
        builder.addNextEdge();
        builder.assignProperty(order, bonds[i][2] ?? 1);
        builder.assignProperty(flags, bonds[i][3] ?? 0);
    }
    const graph = builder.createGraph({ order, flags });

    return { structure: undefined, unit: undefined, compId: 'LIG', residueKey: 0, vertices, bonds: graph } as unknown as LigandGraph;
}

const AR = BondType.Flag.Aromatic;

/** A monocycle of `n` aromatic carbons with the given per-bond Kekulé orders (length n). */
function aromaticRing(kekuleOrders: number[]): LigandGraph {
    const n = kekuleOrders.length;
    const elements = new Array<string>(n).fill('C');
    const bonds: Bond[] = [];
    for (let i = 0; i < n; i++) bonds.push([i, (i + 1) % n, kekuleOrders[i], AR]);
    return makeGraph(elements, bonds);
}

/** A monocycle of `n` carbons with the given bond order. */
function carbonRing(n: number, order: number): LigandGraph {
    const elements = new Array<string>(n).fill('C');
    const bonds: Bond[] = [];
    for (let i = 0; i < n; i++) bonds.push([i, (i + 1) % n, order]);
    return makeGraph(elements, bonds);
}

/** A linear chain of the given element symbols (single bonds). */
function chain(elements: string[]): LigandGraph {
    const bonds: Bond[] = [];
    for (let i = 0; i + 1 < elements.length; i++) bonds.push([i, i + 1]);
    return makeGraph(elements, bonds);
}

const strictBondOrder = (orderA: number, flagsA: number, orderB: number, flagsB: number) => (orderA | 0) === (orderB | 0);

describe('Ligand MCCS', () => {
    it('benzene vs toluene: ring MCCS of size 6', () => {
        const benzene = carbonRing(6, 1);
        const toluene = makeGraph(
            ['C', 'C', 'C', 'C', 'C', 'C', 'C'],
            [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6]]
        );

        const pairs = findMccsPairs(benzene, toluene);
        expect(pairs).toBeDefined();
        expect(pairs!.length).toBe(6);
    });

    it('benzene vs cyclohexane: no MCCS when bond orders must match (aromatic vs aliphatic)', () => {
        const benzene = carbonRing(6, 2); // stand-in for aromatic
        const cyclohexane = carbonRing(6, 1);

        const cliques = findMccsCliques(benzene, cyclohexane, {
            vertexTest: (a, b) => a.elementSymbol === b.elementSymbol,
            edgeTest: strictBondOrder
        });
        // no c-edges survive the bond-order test, so the only c-connected cliques are singletons
        expect(cliques).toHaveLength(0);
    });

    it('identical benzene: full ring match with symmetry-related candidates', () => {
        const a = carbonRing(6, 1);
        const b = carbonRing(6, 1);

        const cliques = findMccsCliques(a, b);
        expect(cliques.length).toBeGreaterThan(1); // ring symmetry -> several equally-sized mappings
        expect(cliques[0]).toHaveLength(6);
    });

    it('linear chain C-C-C-O matches itself fully', () => {
        const a = chain(['C', 'C', 'C', 'O']);
        const b = chain(['C', 'C', 'C', 'O']);

        const pairs = findMccsPairs(a, b);
        expect(pairs!.length).toBe(4);
    });

    it('element mismatch trims the correspondence (C-C-C-N vs C-C-C-O)', () => {
        const a = chain(['C', 'C', 'C', 'N']);
        const b = chain(['C', 'C', 'C', 'O']);

        const pairs = findMccsPairs(a, b);
        expect(pairs!.length).toBe(3); // shared C-C-C backbone; terminal N/O are incompatible
    });

    it('pathCutoff 0 disables d-edges, collapsing a ring to a single bond', () => {
        const a = carbonRing(6, 1);
        const b = carbonRing(6, 1);

        const cliques = findMccsCliques(a, b, { pathCutoff: 0, minMatchedAtoms: 1 });
        // without d-edges a clique must be a complete subgraph; a 6-cycle has no triangle
        expect(cliques[0]).toHaveLength(2);
    });

    it('does not bridge disconnected fragments (connected-only)', () => {
        // two separate C-C dimers in each graph
        const a = makeGraph(['C', 'C', 'C', 'C'], [[0, 1], [2, 3]]);
        const b = makeGraph(['C', 'C', 'C', 'C'], [[0, 1], [2, 3]]);

        const cliques = findMccsCliques(a, b, { minMatchedAtoms: 1 });
        // cross-fragment pairs are unreachable -> no d-edges across the gap -> max connected match is one dimer
        expect(cliques[0]).toHaveLength(2);
    });

    it('is deterministic across runs', () => {
        const a = carbonRing(6, 1);
        const toluene = makeGraph(
            ['C', 'C', 'C', 'C', 'C', 'C', 'C'],
            [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6]]
        );

        const first = JSON.stringify(findMccsPairs(a, toluene));
        const second = JSON.stringify(findMccsPairs(a, toluene));
        expect(first).toEqual(second);
    });

    it('honors a tiny search budget without throwing', () => {
        const a = carbonRing(6, 1);
        const b = carbonRing(6, 1);

        expect(() => findMccsCliques(a, b, { maxTimeMs: 0, maxIterations: 1 })).not.toThrow();
    });

    it('matches aromatic rings across differing Kekulé bond-order assignments', () => {
        // two benzene rings flagged aromatic but with opposite single/double placements
        const a = aromaticRing([2, 1, 2, 1, 2, 1]);
        const b = aromaticRing([1, 2, 1, 2, 1, 2]);

        const pairs = findMccsPairs(a, b); // default edgeTest is aromatic-aware
        expect(pairs!.length).toBe(6);
    });

    it('does not match an aromatic ring to an aliphatic ring of equal bond orders', () => {
        // identical integer orders, but only one ring is flagged aromatic
        const aromatic = aromaticRing([1, 1, 1, 1, 1, 1]);
        const aliphatic = carbonRing(6, 1); // no aromatic flag

        const cliques = findMccsCliques(aromatic, aliphatic);
        expect(cliques).toHaveLength(0);
    });

    it('returns nothing when the overlap is below minMatchedAtoms', () => {
        const a = chain(['C', 'O']); // only 2 atoms
        const b = chain(['C', 'O']);

        expect(findMccsPairs(a, b, { minMatchedAtoms: 3 })).toBeUndefined();
        expect(DefaultLigandMccsOptions.minMatchedAtoms).toBe(3);
    });
});