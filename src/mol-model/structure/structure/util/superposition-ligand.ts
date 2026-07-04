/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { SortedArray } from '../../../../mol-data/int';
import { ElementIndex } from '../../model';
import { BondType } from '../../model/types';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Mat4 } from '../../../../mol-math/linear-algebra';
import { MinimizeRmsd, RmsdTransformState } from '../../../../mol-math/linear-algebra/3d/minimize-rmsd';
import { RuntimeContext } from '../../../../mol-task';
import { BitSet } from '../../../../mol-util/bit-set';
import { Structure, StructureElement, StructureProperties, Unit } from '../../structure';

/**
 * Ligand MCCS utilities:
 * - fast path: same compId -> match by atom name
 * - MCCS search: connected common subgraph search with configurable vertex/edge matching lambdas
 */
export interface LigandGraphVertex {
    readonly index: number;
    readonly element: ElementIndex;
    readonly atomName: string;
    readonly elementSymbol: string;
}

/** Bond order and BondType flags, stored as parallel per-edge arrays on the adjacency graph. */
export type LigandBondProps = { readonly order: ArrayLike<number>, readonly flags: ArrayLike<number> };

/** Local-vertex-indexed adjacency of one ligand residue (an induced subgraph of `unit.bonds`). */
export type LigandBondGraph = IntAdjacencyGraph<number, LigandBondProps>;

export interface LigandGraph {
    readonly structure: Structure;
    readonly unit: Unit.Atomic;
    readonly compId: string;
    readonly residueKey: number;
    readonly vertices: readonly LigandGraphVertex[];
    readonly bonds: LigandBondGraph;
}

export type LigandMccsVertexTest = (a: LigandGraphVertex, b: LigandGraphVertex) => boolean;
export type LigandMccsEdgeTest = (orderA: number, flagsA: number, orderB: number, flagsB: number) => boolean;

export interface LigandMccsOptions {
    /** Whether atom i of A may be paired with atom j of B (default: same element symbol). */
    readonly vertexTest: LigandMccsVertexTest;
    /** Whether a bond in A may be paired with a bond in B (default: equal order, unknown = wildcard). */
    readonly edgeTest: LigandMccsEdgeTest;
    /** skip hydrogen/deuterium atoms */
    readonly ignoreHydrogens: boolean;
    /**
     * Topological-distance tolerance for "d-edges" (pairs of atoms that are non-bonded in both
     * ligands). These are required for the clique formulation to represent non-complete-graph
     * substructures (e.g. rings), so the default is 1: two non-bonded atom pairs may be matched
     * iff their shortest-path distances differ by at most this much. Higher values are more
     * permissive about topological stretch; 0 disables d-edges entirely (only complete-subgraph
     * matches survive, rarely useful).
     */
    readonly pathCutoff: number;
    readonly maxTimeMs: number;
    readonly minMatchedAtoms: number;
    /** hard backstop on clique-search recursion steps */
    readonly maxIterations: number;
    /** number of equally-sized maximum cliques to retain and RMSD-rank when picking the pose */
    readonly maxCliquesForPose: number;
}

export const DefaultLigandMccsOptions: LigandMccsOptions = {
    vertexTest: (a, b) => a.elementSymbol === b.elementSymbol,
    edgeTest: (orderA, flagsA, orderB, flagsB) => {
        // aromatic bonds match each other regardless of the (arbitrary) Kekulé order the CCD assigns,
        // but an aromatic bond is not compatible with an aliphatic one
        const aArom = (flagsA & BondType.Flag.Aromatic) !== 0;
        const bArom = (flagsB & BondType.Flag.Aromatic) !== 0;
        if (aArom || bArom) return aArom === bArom;
        // unknown/0 bond orders are treated as compatible
        const ao = orderA | 0, bo = orderB | 0;
        if (ao > 0 && bo > 0) return ao === bo;
        return true;
    },
    ignoreHydrogens: true,
    pathCutoff: 1,
    maxTimeMs: 400,
    minMatchedAtoms: 3,
    maxIterations: 10000000,
    maxCliquesForPose: 64
};

function isAtomicUnit(u: Unit): u is Unit.Atomic {
    return u.kind === Unit.Kind.Atomic;
}

/**
 * Returns residueKey/compId iff:
 * - loci is non-empty
 * - all atoms are in the same residue (one ligand residue)
 * - entity type is not polymer
 */
export function isSingleLigandLoci(loci: StructureElement.Loci): { residueKey: number, compId: string } | undefined {
    if (StructureElement.Loci.isEmpty(loci)) return;

    const loc = StructureElement.Location.create(loci.structure);
    const first = StructureElement.Loci.getFirstLocation(loci, loc);
    if (!first) return;

    if (StructureProperties.entity.type(first) === 'polymer') return;

    const firstResidueKey = StructureProperties.residue.key(first);
    const firstCompId = StructureProperties.atom.label_comp_id(first);

    let ok = true;
    StructureElement.Loci.forEachLocation(loci, l => {
        if (!ok) return;
        if (StructureProperties.entity.type(l) === 'polymer') ok = false;
        if (StructureProperties.residue.key(l) !== firstResidueKey) ok = false;
    });

    if (!ok) return;
    return { residueKey: firstResidueKey, compId: firstCompId };
}

export function buildLigandGraphFromLoci(loci: StructureElement.Loci): LigandGraph | undefined {
    const info = isSingleLigandLoci(loci);
    if (!info) return;

    if (loci.elements.length !== 1) return;
    const e0 = loci.elements[0];
    const unit = e0.unit;
    if (!isAtomicUnit(unit)) return;

    const { residueKey, compId } = info;

    const indices: ElementIndex[] = [];
    StructureElement.Loci.forEachLocation(loci, l => {
        indices.push(l.element);
    });
    if (indices.length === 0) return;

    const loc = StructureElement.Location.create(loci.structure, unit);
    for (let i = 0; i < indices.length; i++) {
        loc.element = indices[i];
        if (StructureProperties.residue.key(loc) !== residueKey) return;
    }

    indices.sort((a, b) => a - b);

    // intra-unit bonds are keyed by unit-local UnitIndex, not the global ElementIndex, so track both
    const unitIndices: number[] = new Array(indices.length);
    const unitIndexToVertex = new Map<number, number>();
    const vertices: LigandGraphVertex[] = indices.map((el, i) => {
        loc.element = el;
        const uIndex = SortedArray.indexOf(unit.elements, el);
        unitIndices[i] = uIndex;
        unitIndexToVertex.set(uIndex, i);
        return {
            index: i,
            element: el,
            atomName: StructureProperties.atom.label_atom_id(loc),
            elementSymbol: StructureProperties.atom.type_symbol(loc)
        };
    });

    // collect the ligand's intra-residue bonds as local-vertex edge pairs (each undirected edge once)
    const { offset, b, edgeProps } = unit.bonds;
    const srcOrder = edgeProps.order, srcFlags = edgeProps.flags;

    const xs: number[] = [], ys: number[] = [], orders: number[] = [], flags: number[] = [];
    for (let vi = 0; vi < vertices.length; vi++) {
        const uIndex = unitIndices[vi];
        for (let ei = offset[uIndex], end = offset[uIndex + 1]; ei < end; ei++) {
            const vj = unitIndexToVertex.get(b[ei]); // b[ei] is a UnitIndex
            if (vj === void 0 || vj <= vi) continue; // neighbor outside the ligand, or already added (vj < vi)
            xs.push(vi); ys.push(vj);
            orders.push(srcOrder[ei]); flags.push(srcFlags[ei]);
        }
    }

    const builder = new IntAdjacencyGraph.EdgeBuilder(vertices.length, xs, ys);
    const order = new Int32Array(builder.slotCount);
    const flag = new Int32Array(builder.slotCount);
    for (let i = 0; i < builder.edgeCount; i++) {
        builder.addNextEdge();
        builder.assignProperty(order, orders[i]);
        builder.assignProperty(flag, flags[i]);
    }
    const bonds = builder.createGraph({ order, flags: flag });

    return { structure: loci.structure, unit, compId, residueKey, vertices, bonds };
}

/** Fast path: same compId -> match by atom name (preferring same element symbol). */
export function matchLigandsByAtomName(a: LigandGraph, b: LigandGraph): Array<[number, number]> | undefined {
    if (!a.compId || !b.compId) return;
    if (a.compId.toUpperCase() !== b.compId.toUpperCase()) return;

    const byName = new Map<string, number[]>();
    for (const v of b.vertices) {
        const name = v.atomName || '';
        if (!name) continue;
        const arr = byName.get(name);
        if (arr) arr.push(v.index);
        else byName.set(name, [v.index]);
    }

    const usedB = new Uint8Array(b.vertices.length);
    const pairs: Array<[number, number]> = [];

    for (const va of a.vertices) {
        const name = va.atomName || '';
        if (!name) continue;

        const cands = byName.get(name);
        if (!cands) continue;

        let chosen = -1;
        for (let i = 0; i < cands.length; i++) {
            const bi = cands[i];
            if (usedB[bi]) continue;
            if (b.vertices[bi].elementSymbol === va.elementSymbol) { chosen = bi; break; }
        }
        if (chosen < 0) {
            for (let i = 0; i < cands.length; i++) {
                const bi = cands[i];
                if (!usedB[bi]) { chosen = bi; break; }
            }
        }
        if (chosen < 0) continue;

        usedB[chosen] = 1;
        pairs.push([va.index, chosen]);
    }

    if (!pairs.length) return;
    pairs.sort((p, q) => p[0] - q[0]);
    return pairs;
}

/** LIFO pool of fixed-size BitSets so the (heavily recursive) clique search does not allocate per call. */
class BitSetPool {
    private readonly free: BitSet[] = [];
    constructor(private readonly nbits: number) {}
    acquire(): BitSet {
        const b = this.free.pop();
        return b !== void 0 ? b : new BitSet(this.nbits);
    }
    acquireCopy(src: BitSet): BitSet {
        const b = this.acquire();
        b.words.set(src.words);
        return b;
    }
    release(b: BitSet) { this.free.push(b); }
}

/**
 * Compatibility (modular product) graph. Vertices are atom pairs (pairA[k] in A, pairB[k] in B)
 * that pass the vertex test. Two vertices are joined by a c-edge when their atoms are bonded in
 * both ligands with compatible bonds, or by a d-edge when they are non-bonded in both with similar
 * topological distance. A c-connected clique of this graph is a connected common substructure.
 */
interface CompatGraph {
    readonly n: number;
    readonly pairA: Int32Array;
    readonly pairB: Int32Array;
    readonly adj: BitSet[]; // c-edges and d-edges
    readonly cadj: BitSet[]; // c-edges only
}

function isHydrogenVertex(v: LigandGraphVertex): boolean {
    const s = v.elementSymbol.toUpperCase();
    return s === 'H' || s === 'D';
}

/** BFS all-pairs shortest paths (edge count) within a ligand graph; -1 marks unreachable. */
function ligandShortestPaths(g: LigandGraph): Int32Array[] {
    const n = g.vertices.length;
    const { offset, b } = g.bonds;
    const dist: Int32Array[] = new Array(n);
    const queue = new Int32Array(n);
    for (let s = 0; s < n; s++) {
        const d = new Int32Array(n).fill(-1);
        let head = 0, tail = 0;
        queue[tail++] = s;
        d[s] = 0;
        while (head < tail) {
            const u = queue[head++];
            const du = d[u];
            for (let t = offset[u], end = offset[u + 1]; t < end; t++) {
                const w = b[t];
                if (d[w] < 0) { d[w] = du + 1; queue[tail++] = w; }
            }
        }
        dist[s] = d;
    }
    return dist;
}

function buildCompatGraph(a: LigandGraph, b: LigandGraph, opts: LigandMccsOptions): CompatGraph {
    const nA = a.vertices.length, nB = b.vertices.length;
    const ignoreH = opts.ignoreHydrogens;

    const pa: number[] = [], pb: number[] = [];
    for (let i = 0; i < nA; i++) {
        if (ignoreH && isHydrogenVertex(a.vertices[i])) continue;
        for (let j = 0; j < nB; j++) {
            if (ignoreH && isHydrogenVertex(b.vertices[j])) continue;
            if (opts.vertexTest(a.vertices[i], b.vertices[j])) { pa.push(i); pb.push(j); }
        }
    }

    const n = pa.length;
    const pairA = Int32Array.from(pa);
    const pairB = Int32Array.from(pb);

    const adj: BitSet[] = new Array(n);
    const cadj: BitSet[] = new Array(n);
    for (let i = 0; i < n; i++) { adj[i] = new BitSet(n); cadj[i] = new BitSet(n); }

    const useD = opts.pathCutoff > 0;
    const distA = useD ? ligandShortestPaths(a) : undefined;
    const distB = useD ? ligandShortestPaths(b) : undefined;

    const aBonds = a.bonds, bBonds = b.bonds;
    const aOrder = aBonds.edgeProps.order, aFlags = aBonds.edgeProps.flags;
    const bOrder = bBonds.edgeProps.order, bFlags = bBonds.edgeProps.flags;

    for (let i = 0; i < n; i++) {
        const ai = pairA[i], bi = pairB[i];
        for (let j = i + 1; j < n; j++) {
            const aj = pairA[j], bj = pairB[j];

            // 1:1 constraint: a given atom may not be matched twice within one clique
            if (ai === aj || bi === bj) continue;

            const eaIdx = aBonds.getDirectedEdgeIndex(ai, aj);
            const ebIdx = bBonds.getDirectedEdgeIndex(bi, bj);
            const adjA = eaIdx >= 0;
            const adjB = ebIdx >= 0;

            let cEdge = false;
            if (adjA && adjB) cEdge = opts.edgeTest(aOrder[eaIdx], aFlags[eaIdx], bOrder[ebIdx], bFlags[ebIdx]);

            let dEdge = false;
            if (!cEdge && useD && !adjA && !adjB) {
                const dA = distA![ai][aj];
                const dB = distB![bi][bj];
                if (dA >= 0 && dB >= 0 && Math.abs(dA - dB) <= opts.pathCutoff) dEdge = true;
            }

            if (cEdge || dEdge) {
                adj[i].set(j); adj[j].set(i);
                if (cEdge) { cadj[i].set(j); cadj[j].set(i); }
            }
        }
    }

    return { n, pairA, pairB, adj, cadj };
}

/** Retains every maximum-size clique encountered (up to a cap) for downstream RMSD pose selection. */
class CliqueCollector {
    bestSize = 0;
    cliques: number[][] = [];
    private readonly cap: number;
    constructor(cap: number) { this.cap = Math.max(1, cap); }
    offer(clique: number[]) {
        const size = clique.length;
        if (size < this.bestSize) return;
        if (size > this.bestSize) { this.bestSize = size; this.cliques = [clique]; return; }
        if (this.cliques.length < this.cap) this.cliques.push(clique);
    }
}

/**
 * Enumerate maximal c-connected cliques of the compatibility graph (Koch's connected-clique variant
 * of Bron-Kerbosch). c-edges keep the match a single connected fragment; d-edges may appear inside a
 * clique but never seed connectivity. Bounded by a wall-clock and a recursion-step budget.
 */
function enumerateConnectedCliques(g: CompatGraph, opts: LigandMccsOptions, collector: CliqueCollector) {
    const n = g.n;
    if (n === 0) return;

    const minK = Math.max(1, opts.minMatchedAtoms);
    const adj = g.adj, cadj = g.cadj;

    const dAdj: BitSet[] = new Array(n);
    for (let i = 0; i < n; i++) { const d = adj[i].clone(); d.andNot(cadj[i]); dAdj[i] = d; }

    const all = BitSet.full(n);
    const t = new BitSet(n); // already-processed seeds, so each maximal clique is found once
    const pool = new BitSetPool(n);

    const deadline = Date.now() + opts.maxTimeMs;
    const maxIter = opts.maxIterations;
    let iters = 0;
    let aborted = false;

    const checkAbort = () => {
        if (aborted) return true;
        if ((maxIter > 0 && iters > maxIter) || Date.now() > deadline) { aborted = true; return true; }
        return false;
    };

    const cClique = (R: BitSet, P: BitSet, Q: BitSet, X: BitSet, Y: BitSet) => {
        if (checkAbort()) return;
        iters++;

        if (P.isEmpty() && X.isEmpty()) {
            const size = R.cardinality();
            if (size >= minK) {
                const clique: number[] = new Array(size);
                let idx = 0;
                for (let v = R.nextSetBit(0); v >= 0; v = R.nextSetBit(v + 1)) clique[idx++] = v;
                collector.offer(clique);
            }
            return;
        }

        // scratch sets are drawn from a pool and returned at the end of each iteration/frame
        const pIter = pool.acquireCopy(P);
        for (let ui = pIter.nextSetBit(0); ui >= 0; ui = pIter.nextSetBit(ui + 1)) {
            P.clear(ui);

            const dN = dAdj[ui], cN = cadj[ui], neigh = adj[ui];

            const rNew = pool.acquireCopy(R); rNew.set(ui);
            const qNew = pool.acquireCopy(Q); qNew.and(dN);
            const yNew = pool.acquireCopy(Y); yNew.and(dN);
            const pNew = pool.acquireCopy(P); pNew.and(neigh);
            const xNew = pool.acquireCopy(X); xNew.and(neigh);
            const tmp = pool.acquireCopy(Q); tmp.and(cN); pNew.or(tmp); // P' = (P ∩ neigh) ∪ (Q ∩ cN)
            tmp.words.set(Y.words); tmp.and(cN); xNew.or(tmp); // X' = (X ∩ neigh) ∪ (Y ∩ cN)
            pool.release(tmp);

            cClique(rNew, pNew, qNew, xNew, yNew);

            pool.release(xNew); pool.release(pNew); pool.release(yNew); pool.release(qNew); pool.release(rNew);

            X.set(ui);
            if (checkAbort()) break;
        }
        pool.release(pIter);
    };

    for (let ui = 0; ui < n; ui++) {
        const dN = dAdj[ui], cN = cadj[ui];

        const r = pool.acquire(); r.zero(); r.set(ui);
        const q = pool.acquireCopy(all); q.andNot(t); q.and(dN);
        const p = pool.acquireCopy(all); p.andNot(t); p.and(cN);
        const y = pool.acquireCopy(dN); y.and(t);
        const x = pool.acquireCopy(cN); x.and(t);

        cClique(r, p, q, x, y);

        pool.release(x); pool.release(y); pool.release(p); pool.release(q); pool.release(r);

        t.set(ui);
        if (checkAbort()) break;
    }
}

function cliqueToPairs(g: CompatGraph, clique: number[]): Array<[number, number]> {
    const pairs: Array<[number, number]> = clique.map(ci => [g.pairA[ci], g.pairB[ci]] as [number, number]);
    pairs.sort((p, q) => p[0] - q[0]);
    return pairs;
}

/**
 * Enumerate the maximum-size MCCS atom correspondences between two ligand graphs. Several equally
 * sized candidates may be returned (symmetry-related mappings, e.g. ring flips) so the caller can
 * pick the best pose by RMSD. Each correspondence is a list of [aVertexIndex, bVertexIndex] pairs.
 */
export function findMccsCliques(a: LigandGraph, b: LigandGraph, options: Partial<LigandMccsOptions> = {}): Array<Array<[number, number]>> {
    const opts: LigandMccsOptions = { ...DefaultLigandMccsOptions, ...options };
    const g = buildCompatGraph(a, b, opts);
    if (g.n === 0) return [];

    const collector = new CliqueCollector(opts.maxCliquesForPose);
    enumerateConnectedCliques(g, opts, collector);
    if (collector.bestSize < opts.minMatchedAtoms) return [];

    return collector.cliques.map(clique => cliqueToPairs(g, clique));
}

/** Convenience wrapper returning a single (the first) maximum MCCS correspondence. */
export function findMccsPairs(a: LigandGraph, b: LigandGraph, options: Partial<LigandMccsOptions> = {}): Array<[number, number]> | undefined {
    const cliques = findMccsCliques(a, b, options);
    return cliques.length ? cliques[0] : undefined;
}

export type LigandSuperposeMethod = 'atom-name' | 'mccs';

export interface LigandSuperposeResult {
    method: LigandSuperposeMethod;
    atomCount: number;
    rmsd: number;
    bTransform: Mat4;
    /**
     * matched atoms as model ElementIndex, paired 1:1 (by array position) with `targetElements`
     */
    referenceElements: ElementIndex[];
    /** matched atoms in the target ligand as model ElementIndex, paired 1:1 with `referenceElements` */
    targetElements: ElementIndex[];
}

type MutablePositions = { x: Float64Array, y: Float64Array, z: Float64Array };

/** Write the paired atom coordinates into reused position buffers, in correspondence order. */
function fillPairPositions(gA: LigandGraph, gB: LigandGraph, pairs: Array<[number, number]>, posA: MutablePositions, posB: MutablePositions) {
    const ca = gA.unit.conformation, cb = gB.unit.conformation;
    for (let i = 0; i < pairs.length; i++) {
        const ea = gA.vertices[pairs[i][0]].element;
        const eb = gB.vertices[pairs[i][1]].element;
        posA.x[i] = ca.x(ea); posA.y[i] = ca.y(ea); posA.z[i] = ca.z(ea);
        posB.x[i] = cb.x(eb); posB.y[i] = cb.y(eb); posB.z[i] = cb.z(eb);
    }
}

function pairsToElements(gA: LigandGraph, gB: LigandGraph, pairs: Array<[number, number]>) {
    const referenceElements: ElementIndex[] = new Array(pairs.length);
    const targetElements: ElementIndex[] = new Array(pairs.length);
    for (let i = 0; i < pairs.length; i++) {
        referenceElements[i] = gA.vertices[pairs[i][0]].element;
        targetElements[i] = gB.vertices[pairs[i][1]].element;
    }
    return { referenceElements, targetElements };
}

/**
 * Given two ligand loci (each one ligand residue), compute the bTransform that superposes the
 * target onto the reference. Tries the atom-name fast path first (identical compId), then falls
 * back to an MCCS search whose equally-sized maximum correspondences are each superposed so the
 * lowest-RMSD pose wins (this resolves symmetry-related mappings the largest-clique alone cannot).
 */
export function superposeLigandsByMccs(
    reference: StructureElement.Loci,
    target: StructureElement.Loci,
    options: Partial<LigandMccsOptions> = {},
    ctx?: RuntimeContext
): LigandSuperposeResult | undefined {
    const opts: LigandMccsOptions = { ...DefaultLigandMccsOptions, ...options };

    const gA = buildLigandGraphFromLoci(reference);
    const gB = buildLigandGraphFromLoci(target);
    if (!gA || !gB) return;

    // fast path: identical compId -> match by atom name
    const namePairs = matchLigandsByAtomName(gA, gB);
    if (namePairs && namePairs.length >= opts.minMatchedAtoms) {
        const n = namePairs.length;
        const posA = MinimizeRmsd.Positions.empty(n);
        const posB = MinimizeRmsd.Positions.empty(n);
        fillPairPositions(gA, gB, namePairs, posA, posB);
        const { bTransform, rmsd } = MinimizeRmsd.compute({ a: posA, b: posB, length: n });
        const { referenceElements, targetElements } = pairsToElements(gA, gB, namePairs);
        return { method: 'atom-name', atomCount: n, rmsd, bTransform, referenceElements, targetElements };
    }

    // MCCS: superpose each maximum-size correspondence and keep the lowest-RMSD pose. Every clique
    // shares the maximum size, so one set of reused buffers (and a reused RMSD state) fits all of them.
    const cliques = findMccsCliques(gA, gB, opts);
    if (cliques.length === 0) return;

    const n = cliques[0].length;
    const posA = MinimizeRmsd.Positions.empty(n);
    const posB = MinimizeRmsd.Positions.empty(n);
    const result: MinimizeRmsd.Result = { bTransform: Mat4.zero(), rmsd: 0, nAlignedElements: 0 };
    const state = new RmsdTransformState({ a: posA, b: posB, length: n }, result);

    let bestRmsd = Number.POSITIVE_INFINITY;
    let bestPairs: Array<[number, number]> | undefined;
    const bestTransform = Mat4.zero();
    for (const pairs of cliques) {
        fillPairPositions(gA, gB, pairs, posA, posB);
        MinimizeRmsd.compute({ a: posA, b: posB, length: n }, result, state);
        if (result.rmsd < bestRmsd) {
            bestRmsd = result.rmsd;
            Mat4.copy(bestTransform, result.bTransform);
            bestPairs = pairs;
        }
    }
    if (!bestPairs) return;

    const { referenceElements, targetElements } = pairsToElements(gA, gB, bestPairs);
    return { method: 'mccs', atomCount: bestPairs.length, rmsd: bestRmsd, bTransform: bestTransform, referenceElements, targetElements };
}