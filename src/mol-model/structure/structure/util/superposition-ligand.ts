/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { OrderedSet, SortedArray } from '../../../../mol-data/int';
import { ElementIndex } from '../../model';
import { BondType } from '../../model/types';
import { Mat4 } from '../../../../mol-math/linear-algebra';
import { RuntimeContext } from '../../../../mol-task';
import { Structure, StructureElement, StructureProperties, Unit } from '../../structure';
import { superpose } from './superposition';
import { UnitIndex } from '../element/util';

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

export interface LigandGraphEdge {
    readonly a: number;
    readonly b: number;
    readonly order: number;
    readonly flags: number;
    readonly key: number;
}

export interface LigandGraph {
    readonly structure: Structure;
    readonly unit: Unit.Atomic;
    readonly compId: string;
    readonly residueKey: number;
    readonly vertices: readonly LigandGraphVertex[];
    readonly adj: readonly Map<number, LigandGraphEdge>[];
}

export type LigandMccsVertexTest = (a: LigandGraphVertex, b: LigandGraphVertex) => boolean;
export type LigandMccsEdgeTest = (a: LigandGraphEdge, b: LigandGraphEdge) => boolean;

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
    edgeTest: (a, b) => {
        // aromatic bonds match each other regardless of the (arbitrary) Kekulé order the CCD assigns,
        // but an aromatic bond is not compatible with an aliphatic one
        const aArom = (a.flags & BondType.Flag.Aromatic) !== 0;
        const bArom = (b.flags & BondType.Flag.Aromatic) !== 0;
        if (aArom || bArom) return aArom === bArom;
        // unknown/0 bond orders are treated as compatible
        const ao = a.order | 0, bo = b.order | 0;
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

    const adj: Map<number, LigandGraphEdge>[] = [];
    for (let i = 0; i < vertices.length; i++) adj.push(new Map());

    const { offset, b, edgeProps } = unit.bonds;
    const { order, flags, key } = edgeProps;

    for (let vi = 0; vi < vertices.length; vi++) {
        const uIndex = unitIndices[vi];
        const start = offset[uIndex];
        const end = offset[uIndex + 1];

        for (let ei = start; ei < end; ei++) {
            const vj = unitIndexToVertex.get(b[ei]); // b[ei] is a UnitIndex
            if (vj === void 0) continue; // neighbor outside the selected ligand

            const edge: LigandGraphEdge = {
                a: vi,
                b: vj,
                order: order[ei],
                flags: flags[ei],
                key: key ? key[ei] : -1
            };

            if (!adj[vi].has(vj)) adj[vi].set(vj, edge);
            if (!adj[vj].has(vi)) adj[vj].set(vi, { ...edge, a: vj, b: vi });
        }
    }

    return { structure: loci.structure, unit, compId, residueKey, vertices, adj };
}

function lociFromVertexOrder(graph: LigandGraph, vertexOrder: readonly number[]): StructureElement.Loci {
    // one element per atom so the A/B atom correspondence (and thus the RMSD pairing) is preserved
    // in order; grouping all indices into a single OrderedSet would re-sort them and break it
    const elements: StructureElement.Loci['elements'][0][] = [];
    for (let i = 0; i < vertexOrder.length; i++) {
        const v = graph.vertices[vertexOrder[i]];
        const uIdx = SortedArray.indexOf(graph.unit.elements, v.element) as UnitIndex;
        elements.push({ unit: graph.unit, indices: OrderedSet.ofSingleton(uIdx) });
    }
    return StructureElement.Loci(graph.structure, elements);
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

/** Minimal fixed-size bit set backed by a Uint32Array, used by the clique search. */
class BitSet {
    readonly nbits: number;
    readonly words: Uint32Array;

    constructor(nbits: number) {
        this.nbits = nbits;
        this.words = new Uint32Array((nbits + 31) >>> 5);
    }

    static full(nbits: number): BitSet {
        const b = new BitSet(nbits);
        b.words.fill(0xffffffff);
        const rem = nbits & 31;
        if (rem !== 0) b.words[b.words.length - 1] = ~(0xffffffff << rem);
        return b;
    }

    clone(): BitSet {
        const b = new BitSet(this.nbits);
        b.words.set(this.words);
        return b;
    }

    set(i: number) { this.words[i >>> 5] |= (1 << (i & 31)); }
    clear(i: number) { this.words[i >>> 5] &= ~(1 << (i & 31)); }
    get(i: number): boolean { return (this.words[i >>> 5] & (1 << (i & 31))) !== 0; }

    and(o: BitSet) { const w = this.words, x = o.words; for (let k = 0, n = w.length; k < n; k++) w[k] &= x[k]; }
    andNot(o: BitSet) { const w = this.words, x = o.words; for (let k = 0, n = w.length; k < n; k++) w[k] &= ~x[k]; }
    or(o: BitSet) { const w = this.words, x = o.words; for (let k = 0, n = w.length; k < n; k++) w[k] |= x[k]; }

    isEmpty(): boolean {
        const w = this.words;
        for (let k = 0, n = w.length; k < n; k++) if (w[k] !== 0) return false;
        return true;
    }

    cardinality(): number {
        const w = this.words;
        let c = 0;
        for (let k = 0, n = w.length; k < n; k++) c += bitCount(w[k]);
        return c;
    }

    nextSetBit(from: number): number {
        if (from < 0) from = 0;
        if (from >= this.nbits) return -1;
        const w = this.words;
        let wi = from >>> 5;
        let word = w[wi] & (0xffffffff << (from & 31));
        while (true) {
            if (word !== 0) {
                const bit = (wi << 5) + (31 - Math.clz32(word & -word));
                return bit < this.nbits ? bit : -1;
            }
            if (++wi >= w.length) return -1;
            word = w[wi];
        }
    }
}

function bitCount(n: number): number {
    n = n - ((n >>> 1) & 0x55555555);
    n = (n & 0x33333333) + ((n >>> 2) & 0x33333333);
    return (((n + (n >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
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
            for (const w of g.adj[u].keys()) {
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

    for (let i = 0; i < n; i++) {
        const ai = pairA[i], bi = pairB[i];
        for (let j = i + 1; j < n; j++) {
            const aj = pairA[j], bj = pairB[j];

            // 1:1 constraint: a given atom may not be matched twice within one clique
            if (ai === aj || bi === bj) continue;

            const ea = a.adj[ai].get(aj);
            const eb = b.adj[bi].get(bj);
            const adjA = ea !== void 0;
            const adjB = eb !== void 0;

            let cEdge = false;
            if (adjA && adjB) cEdge = opts.edgeTest(ea!, eb!);

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

        const pIter = P.clone();
        for (let ui = pIter.nextSetBit(0); ui >= 0; ui = pIter.nextSetBit(ui + 1)) {
            P.clear(ui);

            const dN = dAdj[ui], cN = cadj[ui], neigh = adj[ui];

            const rNew = R.clone(); rNew.set(ui);
            const qNew = Q.clone(); qNew.and(dN);
            const pNew = P.clone(); pNew.and(neigh);
            const qCapC = Q.clone(); qCapC.and(cN); pNew.or(qCapC);
            const yNew = Y.clone(); yNew.and(dN);
            const xNew = X.clone(); xNew.and(neigh);
            const yCapC = Y.clone(); yCapC.and(cN); xNew.or(yCapC);

            cClique(rNew, pNew, qNew, xNew, yNew);

            X.set(ui);
            if (checkAbort()) return;
        }
    };

    for (let ui = 0; ui < n; ui++) {
        const r = new BitSet(n); r.set(ui);
        const dN = dAdj[ui], cN = cadj[ui];

        const q = all.clone(); q.andNot(t); q.and(dN);
        const p = all.clone(); p.andNot(t); p.and(cN);
        const y = dN.clone(); y.and(t);
        const x = cN.clone(); x.and(t);

        cClique(r, p, q, x, y);

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
}

/** Superpose target onto reference using a fixed atom correspondence; returns the RMSD-minimizing transform. */
function superposeFromPairs(gA: LigandGraph, gB: LigandGraph, pairs: Array<[number, number]>): { rmsd: number, bTransform: Mat4 } | undefined {
    if (pairs.length === 0) return;

    const aOrder: number[] = [];
    const bOrder: number[] = [];
    for (let i = 0; i < pairs.length; i++) {
        aOrder.push(pairs[i][0]);
        bOrder.push(pairs[i][1]);
    }

    const lociA = lociFromVertexOrder(gA, aOrder);
    const lociB = lociFromVertexOrder(gB, bOrder);

    const results = superpose([lociA, lociB]);
    if (!results.length) return;
    const { bTransform, rmsd } = results[0];
    return { bTransform, rmsd };
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
        const res = superposeFromPairs(gA, gB, namePairs);
        if (res) return { method: 'atom-name', atomCount: namePairs.length, rmsd: res.rmsd, bTransform: res.bTransform };
    }

    // MCCS: superpose each maximum-size correspondence and keep the lowest-RMSD pose
    let best: LigandSuperposeResult | undefined;
    for (const pairs of findMccsCliques(gA, gB, opts)) {
        if (pairs.length < opts.minMatchedAtoms) continue;
        const res = superposeFromPairs(gA, gB, pairs);
        if (!res) continue;
        if (!best || res.rmsd < best.rmsd) {
            best = { method: 'mccs', atomCount: pairs.length, rmsd: res.rmsd, bTransform: res.bTransform };
        }
    }
    return best;
}