/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { BitSet } from '../../mol-util/bit-set';
import { IntAdjacencyGraph } from './int-adjacency-graph';

export { Mccs };

/**
 * Maximum common connected subgraph (MCCS) matching between two graphs.
 *
 * Builds the modular-product (compatibility) graph and enumerates its connected cliques (Koch's
 * connected-clique variant of Bron-Kerbosch): product vertices are compatible vertex pairs
 * (`a` in graph A, `b` in graph B); two product vertices are joined by a "c-edge" when the
 * underlying vertices are adjacent in both graphs with compatible edges, or by a "d-edge" when they
 * are non-adjacent in both with similar topological distance. A clique connected via c-edges is a
 * common *connected* subgraph, and the maximum such clique is the MCCS. The vertex/edge
 * compatibility predicates are supplied by the caller, so the routine is domain-agnostic.
 */
namespace Mccs {
    export interface Options {
        /** Whether vertex `a` of graph A may be paired with vertex `b` of graph B. */
        vertexTest: (a: number, b: number) => boolean;
        /** Whether edge `edgeA` of A may be paired with edge `edgeB` of B (indices into the graphs' edge arrays). */
        edgeTest: (edgeA: number, edgeB: number) => boolean;
        /**
         * Topological-distance tolerance for "d-edges" (vertex pairs non-adjacent in both graphs).
         * These are required for the clique formulation to represent non-complete subgraphs (e.g.
         * rings), so the default is 1: two non-adjacent pairs may be matched iff their shortest-path
         * distances differ by at most this much. Higher values allow more topological stretch; 0
         * disables d-edges (only complete-subgraph matches survive, rarely useful).
         */
        pathCutoff: number;
        /** Wall-clock budget; the search returns the best cliques found so far once exceeded. */
        maxTimeMs: number;
        /** Hard backstop on recursion steps. */
        maxIterations: number;
        /** Discard cliques smaller than this. */
        minMatchedVertices: number;
        /** Number of equally-sized maximum cliques to retain. */
        maxCliques: number;
    }

    export const DefaultOptions: Options = {
        vertexTest: () => true,
        edgeTest: () => true,
        pathCutoff: 1,
        maxTimeMs: 400,
        maxIterations: 10000000,
        minMatchedVertices: 1,
        maxCliques: 64
    };

    /** Optional out-parameter reporting whether the search hit its time/iteration budget. */
    export interface Status {
        truncated: boolean;
    }

    /**
     * Enumerate the maximum-size common-connected-subgraph correspondences between `a` and `b`.
     * Several equally sized candidates may be returned when the maximum is symmetry-degenerate; each
     * is a list of `[vertexA, vertexB]` pairs sorted by `vertexA`. Empty if nothing meets
     * `minMatchedVertices`.
     *
     * Best-effort: the search stops at `maxTimeMs`/`maxIterations` and returns the best cliques found
     * so far. Pass `status` to learn whether that budget was hit (`status.truncated`).
     */
    export function find(a: IntAdjacencyGraph<any, any>, b: IntAdjacencyGraph<any, any>, options: Partial<Options> = {}, status?: Status): Array<Array<[number, number]>> {
        const opts: Options = { ...DefaultOptions, ...options };
        const g = buildCompatGraph(a, b, opts);
        if (g.n === 0) return [];

        const collector = new CliqueCollector(opts.maxCliques);
        const truncated = enumerateConnectedCliques(g, opts, collector);
        if (status) status.truncated = truncated;
        if (collector.bestSize < opts.minMatchedVertices) return [];

        return collector.cliques.map(clique => cliqueToPairs(g, clique));
    }
}

type Graph = IntAdjacencyGraph<any, any>;

/**
 * Compatibility (modular product) graph. Vertices are vertex pairs (pairA[k] in A, pairB[k] in B)
 * that pass the vertex test; c-edges join pairs adjacent in both graphs with compatible edges,
 * d-edges join pairs non-adjacent in both with similar topological distance.
 */
interface CompatGraph {
    readonly n: number;
    readonly pairA: Int32Array;
    readonly pairB: Int32Array;
    /** c-edges and d-edges */
    readonly adj: BitSet[];
    /** c-edges only */
    readonly cadj: BitSet[];
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

/** BFS all-pairs shortest paths (edge count) within a graph; -1 marks unreachable. */
function graphShortestPaths(g: Graph): Int32Array[] {
    const n = g.vertexCount;
    const { offset, b } = g;
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

function buildCompatGraph(a: Graph, b: Graph, opts: Mccs.Options): CompatGraph {
    const nA = a.vertexCount, nB = b.vertexCount;

    const pa: number[] = [], pb: number[] = [];
    for (let i = 0; i < nA; i++) {
        for (let j = 0; j < nB; j++) {
            if (opts.vertexTest(i, j)) { pa.push(i); pb.push(j); }
        }
    }

    const n = pa.length;
    const pairA = Int32Array.from(pa);
    const pairB = Int32Array.from(pb);

    const adj: BitSet[] = new Array(n);
    const cadj: BitSet[] = new Array(n);
    for (let i = 0; i < n; i++) { adj[i] = new BitSet(n); cadj[i] = new BitSet(n); }

    const useD = opts.pathCutoff > 0;
    const distA = useD ? graphShortestPaths(a) : undefined;
    const distB = useD ? graphShortestPaths(b) : undefined;

    for (let i = 0; i < n; i++) {
        const ai = pairA[i], bi = pairB[i];
        for (let j = i + 1; j < n; j++) {
            const aj = pairA[j], bj = pairB[j];

            // 1:1 constraint: a given vertex may not be matched twice within one clique
            if (ai === aj || bi === bj) continue;

            const eaIdx = a.getDirectedEdgeIndex(ai, aj);
            const ebIdx = b.getDirectedEdgeIndex(bi, bj);
            const adjA = eaIdx >= 0;
            const adjB = ebIdx >= 0;

            let cEdge = false;
            if (adjA && adjB) cEdge = opts.edgeTest(eaIdx, ebIdx);

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

/** Retains every maximum-size clique encountered (up to a cap). */
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
 * clique but never seed connectivity. Bounded by a wall-clock and a recursion-step budget; returns
 * true if that budget was hit (the result is a best-so-far, not proven maximal).
 */
function enumerateConnectedCliques(g: CompatGraph, opts: Mccs.Options, collector: CliqueCollector): boolean {
    const n = g.n;
    if (n === 0) return false;

    const minK = Math.max(1, opts.minMatchedVertices);
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
    return aborted;
}

function cliqueToPairs(g: CompatGraph, clique: number[]): Array<[number, number]> {
    const pairs: Array<[number, number]> = clique.map(ci => [g.pairA[ci], g.pairB[ci]] as [number, number]);
    pairs.sort((p, q) => p[0] - q[0]);
    return pairs;
}
