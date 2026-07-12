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
 * These cover the ligand-domain layer: the default vertex/edge tests (element identity, aromatic-aware
 * bond matching), hydrogen skipping, and the atom-count gate. The generic graph algorithm underneath
 * (connectivity, symmetry, d-edges, budget/truncation) is tested in mol-math/graph/_spec/mccs.spec.ts.
 * Graphs are hand-built LigandGraphs (no Structure/Unit needed); the RMSD/superpose step reuses the
 * already-tested MinimizeRmsd machinery and is best covered by an in-app integration test.
 */

type Bond = [a: number, b: number, order?: number, flags?: number];

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

/** A monocycle of `n` carbons with the given bond order (not aromatic-flagged). */
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

describe('Ligand MCCS', () => {
    it('matches a shared ring (benzene vs toluene) via the default options', () => {
        const benzene = carbonRing(6, 1);
        const toluene = makeGraph(
            ['C', 'C', 'C', 'C', 'C', 'C', 'C'],
            [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [0, 6]]
        );

        const pairs = findMccsPairs(benzene, toluene);
        expect(pairs).toBeDefined();
        expect(pairs!.length).toBe(6);
    });

    it('default vertex test is element identity (C-C-C-N vs C-C-C-O trims to the shared backbone)', () => {
        const a = chain(['C', 'C', 'C', 'N']);
        const b = chain(['C', 'C', 'C', 'O']);

        const pairs = findMccsPairs(a, b);
        expect(pairs!.length).toBe(3); // shared C-C-C; terminal N/O are different elements
    });

    it('default edge test matches aromatic rings across differing Kekulé assignments', () => {
        // two benzene rings flagged aromatic but with opposite single/double placements
        const a = aromaticRing([2, 1, 2, 1, 2, 1]);
        const b = aromaticRing([1, 2, 1, 2, 1, 2]);

        expect(findMccsPairs(a, b)!.length).toBe(6);
    });

    it('default edge test does not match an aromatic ring to an aliphatic ring of equal bond orders', () => {
        const aromatic = aromaticRing([1, 1, 1, 1, 1, 1]); // flagged aromatic
        const aliphatic = carbonRing(6, 1); // same integer orders, not aromatic

        expect(findMccsCliques(aromatic, aliphatic)).toHaveLength(0);
    });

    it('ignoreHydrogens (default) excludes hydrogens from the correspondence', () => {
        const a = makeGraph(['C', 'C', 'H'], [[0, 1], [1, 2]]);
        const b = makeGraph(['C', 'C', 'H'], [[0, 1], [1, 2]]);

        expect(findMccsPairs(a, b, { minMatchedAtoms: 1 })!.length).toBe(2); // H skipped
        expect(findMccsPairs(a, b, { minMatchedAtoms: 1, ignoreHydrogens: false })!.length).toBe(3); // H included
    });

    it('returns nothing when the overlap is below minMatchedAtoms', () => {
        const a = chain(['C', 'O']); // only 2 atoms
        const b = chain(['C', 'O']);

        expect(findMccsPairs(a, b, { minMatchedAtoms: 3 })).toBeUndefined();
        expect(DefaultLigandMccsOptions.minMatchedAtoms).toBe(3);
    });
});
