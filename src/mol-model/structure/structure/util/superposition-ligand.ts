/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { SortedArray } from '../../../../mol-data/int';
import { ElementIndex } from '../../model';
import { BondType } from '../../model/types';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Mccs } from '../../../../mol-math/graph/mccs';
import { Mat4 } from '../../../../mol-math/linear-algebra';
import { MinimizeRmsd, RmsdTransformState } from '../../../../mol-math/linear-algebra/3d/minimize-rmsd';
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

function isHydrogenVertex(v: LigandGraphVertex): boolean {
    const s = v.elementSymbol.toUpperCase();
    return s === 'H' || s === 'D';
}

/** Optional out-parameter reporting whether the MCCS search hit its time/iteration budget. */
export interface LigandMccsStatus {
    truncated: boolean;
}

/**
 * Enumerate the maximum-size MCCS atom correspondences between two ligand graphs. Several equally
 * sized candidates may be returned (symmetry-related mappings, e.g. ring flips) so the caller can
 * pick the best pose by RMSD. Each correspondence is a list of [aVertexIndex, bVertexIndex] pairs.
 *
 * The search is best-effort: it stops at `maxTimeMs`/`maxIterations` and returns the best cliques
 * found so far. Pass `status` to learn whether that budget was hit (`status.truncated`).
 */
export function findMccsCliques(a: LigandGraph, b: LigandGraph, options: Partial<LigandMccsOptions> = {}, status?: LigandMccsStatus): Array<Array<[number, number]>> {
    const opts: LigandMccsOptions = { ...DefaultLigandMccsOptions, ...options };

    // adapt the ligand-domain options to the generic MCCS predicates: fold hydrogen skipping + the
    // element-based vertex test into vertexTest(i, j), and read bond order/flags for edgeTest(ea, eb)
    const ignoreH = opts.ignoreHydrogens;
    const aOrder = a.bonds.edgeProps.order, aFlags = a.bonds.edgeProps.flags;
    const bOrder = b.bonds.edgeProps.order, bFlags = b.bonds.edgeProps.flags;

    return Mccs.find(a.bonds, b.bonds, {
        vertexTest: (ia, ib) => {
            if (ignoreH && (isHydrogenVertex(a.vertices[ia]) || isHydrogenVertex(b.vertices[ib]))) return false;
            return opts.vertexTest(a.vertices[ia], b.vertices[ib]);
        },
        edgeTest: (ea, eb) => opts.edgeTest(aOrder[ea], aFlags[ea], bOrder[eb], bFlags[eb]),
        pathCutoff: opts.pathCutoff,
        maxTimeMs: opts.maxTimeMs,
        maxIterations: opts.maxIterations,
        minMatchedVertices: opts.minMatchedAtoms,
        maxCliques: opts.maxCliquesForPose
    }, status);
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
    /** true if the MCCS search hit its time/iteration budget, so the result may not be the true maximum (always false for the atom-name fast path) */
    truncated: boolean;
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
    options: Partial<LigandMccsOptions> = {}
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
        return { method: 'atom-name', atomCount: n, rmsd, bTransform, referenceElements, targetElements, truncated: false };
    }

    // MCCS: superpose each maximum-size correspondence and keep the lowest-RMSD pose. Every clique
    // shares the maximum size, so one set of reused buffers (and a reused RMSD state) fits all of them.
    const status: LigandMccsStatus = { truncated: false };
    const cliques = findMccsCliques(gA, gB, opts, status);
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
    return { method: 'mccs', atomCount: bestPairs.length, rmsd: bestRmsd, bTransform: bestTransform, referenceElements, targetElements, truncated: status.truncated };
}