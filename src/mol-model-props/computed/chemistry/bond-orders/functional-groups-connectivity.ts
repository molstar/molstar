/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Based on Roger Sayle (2001) method. Functional groups patterns are used to assign bond orders
 * after a geometry assignment step. The patterns are supposed to be unambiguous.
 * Most of the patterns duplicate the ones already in chemistry/functional-group.ts, but the properties added
 * in the declaration allows assignment of bond orders based on information collected during matching.
 * Note the following quirks:
 * - It's possible that 2 functional groups overlap. e.g. HBI C[C@H](O)[C@H](O)C1=NC2=C(NC1)N=C(N)NC2=O has a N in common between an amide and a guanidine.
 *   No backtracking is done to find the best match: for now, order gives priority.
 * - Amide pattern is treated differently from Sayle's depiction: no single bond is enforced to the N partners. This addresses a limitation discussed in Labute (2005).
 */

import { Elements } from '../../../../mol-model/structure/model/properties/atomic/types';
import { BondType, ElementSymbol } from '../../../../mol-model/structure/model/types';
import { isDebugMode } from '../../../../mol-util/debug';
import { AtomGeometry } from '../geometry';
import { atomId } from '../util';
import { computeOpenValence, distSq, getOrder, hasMultipleBond, isAssignedBond, isHydrogenElement, setAssignedBond, setBond, State } from './common';

/** Element spec for functional-group matching: an element symbol or '*' = any heavy atom. */
type ElSpec = Elements | '*';

type NeighborConnectivity = 'terminal' | 'connected' | 'either';
interface NeighbourSpec { el: ElSpec; connectivity: NeighborConnectivity; order: number; }
interface FunctionalGroup { name: string; center: { el: ElSpec; geometry: AtomGeometry }; neighbours: NeighbourSpec[]; fp: number;}

const FunctionalGroups: FunctionalGroup[] = [
    // *N(*)C=O
    { name: 'formamide', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'either', order: 1 },
        // H?
    ], fp: 0 },
    // CC(=O)O
    { name: 'carboxylic acid', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: Elements.C, connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 1 },
    ], fp: 0 },
    // *C(=O)O*
    { name: 'ester', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'connected', order: 1 },
    ], fp: 0 },
    // *C(=O)S*
    { name: 'thioester', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.S, connectivity: 'connected', order: 1 },
    ], fp: 0 },
    // *C(=O)N* — covers amide (C-C(=O)-N, incl. acetamide with a terminal methyl),
    // urea (N-C(=O)-N, e.g. the C2 carbonyl of uracil/thymine) and carbamate (O-C(=O)-N).
    // Note: this differs from the depiction in Sayle (2001) which hints at a carbonyl substituent
    { name: 'amide', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'either', order: 1 },
    ], fp: 0 },
    // *C(=S)S*
    { name: 'dithioester', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.S, connectivity: 'terminal', order: 2 },
        { el: Elements.S, connectivity: 'connected', order: 1 },
    ], fp: 0 },
    // *C(=S)N*
    { name: 'thioamide', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.S, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'connected', order: 1 },
    ], fp: 0 },
    // *NC(=N*)N*
    { name: 'guanidinium', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: Elements.N, connectivity: 'either', order: 2 },
        { el: Elements.N, connectivity: 'either', order: 1 },
        { el: Elements.N, connectivity: 'connected', order: 1 },
    ], fp: 0 },
    // *C(=N)N
    { name: 'amidine', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.N, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'terminal', order: 1 },
    ], fp: 0 },
    // *N=[N+]=[N-]
    { name: 'azide', center: { el: Elements.N, geometry: AtomGeometry.Linear }, neighbours: [
        { el: Elements.N, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'connected', order: 2 },
    ], fp: 0 },
    // *P(=O)(*)*
    { name: 'phosphoryl', center: { el: Elements.P, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'either', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
    ], fp: 0 },
    // *S(=O)(=O)*
    { name: 'sulfonyl', center: { el: Elements.S, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
    ], fp: 0 },
    // *S(=N)(=N)*
    { name: 'thiodiimine', center: { el: Elements.S, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.N, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'terminal', order: 2 },
    ], fp: 0 },
    // *[N+](=O)[O-]
    { name: 'nitro', center: { el: Elements.N, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 1 },
    ], fp: 0 },
    // *C#N
    { name: 'nitrile', center: { el: Elements.C, geometry: AtomGeometry.Linear }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.N, connectivity: 'terminal', order: 3 },
    ], fp: 0 },
    // *[Se](=O)O
    { name: 'seleninic acid', center: { el: Elements.SE, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 1 },
    ], fp: 0 },
];

const FingerprintBits = {
    C: 1 << 0,
    N1: 1 << 1,
    N2: 1 << 2,
    N3: 1 << 3,
    O1: 1 << 4,
    O2: 1 << 5,
    S1: 1 << 6,
    S2: 1 << 7,
};

FunctionalGroups.forEach(g => g.fp = computeFingerprint(g.neighbours.map(n => n.el)));

function computeFingerprint(elements: (ElementSymbol|ElSpec)[]) {
    let fp = 0x0;
    for (const el of elements) {
        switch (el) {
            case Elements.C: fp |= FingerprintBits.C; break;
            case Elements.N:
                if (fp & FingerprintBits.N1) {
                    if (fp & FingerprintBits.N2) fp |= FingerprintBits.N3;
                    else fp |= FingerprintBits.N2;
                } else fp |= FingerprintBits.N1;
                break;
            case Elements.O:
                if (fp & FingerprintBits.O1) {
                    fp |= FingerprintBits.O2;
                } else fp |= FingerprintBits.O1;
                break;
            case Elements.S:
                if (fp & FingerprintBits.S1) {
                    fp |= FingerprintBits.S2;
                } else fp |= FingerprintBits.S1;
                break;
        }
    }
    return fp;
}

function matchEl(spec: ElSpec, el: ElementSymbol) {
    return spec === '*' ? !isHydrogenElement(el) : spec === el;
}

function connectivityOk(c: NeighborConnectivity, isTerminal: boolean) {
    return c === 'either' || (c === 'terminal' ? isTerminal : !isTerminal);
}

export function applyFunctionalGroups(state: State) {
    const { unitIndices, n, geometry, heavyNeighbours } = state;
    for (let i = 0; i < n; i++) {
         // If the center already has a multiple bond (e.g. a backbone C=O or phosphate
        // P=O from the order table), its pi system is already placed - don't add another.
        if (hasMultipleBond(state, i)) continue;
        const el = state.el[i];
        const geom = geometry[i];

        // Skip common cases that are not in patterns
        if (el === Elements.C && (geom !== AtomGeometry.Trigonal && geom !== AtomGeometry.Linear)) continue;
        if (el === Elements.O) continue;
        if (el === Elements.H) continue;
        if (geom === AtomGeometry.Terminal) continue;

        const nb = heavyNeighbours[i];
        const u = unitIndices[i];
        const fp = computeFingerprint(nb.map(j => state.el[j]));
        isDebugMode && console.log('Functional group search for ', atomId(state.unit, u), ' el=', el, ' geom=', geom, ' heavyDeg=', nb.length);
        for (const g of FunctionalGroups) {
            if (g.center.el !== el) continue;
            if (geom !== g.center.geometry) continue;
            if (nb.length !== g.neighbours.length) continue;
            if ((fp & g.fp) !== g.fp) continue; // fingerprint mismatch: some required element is missing
            const slots = assignSlots(state, i, g.neighbours);
            if (!slots) continue;

            isDebugMode && console.log('Functional group ', g.name, ' matched at ', atomId(state.unit, u));
            for (let k = 0; k < slots.length; k++) {
                const v = unitIndices[nb[slots[k]]];
                setAssignedBond(state, u, v);
                const order = g.neighbours[k].order;
                isDebugMode && console.log('Functional group ', g.name, ' setting bond ', atomId(state.unit, u), '-', atomId(state.unit, v), ' order ', order);
                if (order < 2) continue;
                setBond(state.bonds, u, v, order, BondType.Flag.Computed);
            }
            break;
        }
    }
    computeOpenValence(state);
}

/**
 * Match the center's neighbours with the neigbours specs.
 * Returns an array of local indices ordered by the slots in `specs`, or undefined if no complete assignment exists.
 *
 * Perfect bipartite matching of slots to neighbours. A per-pair cost ranks the assignments.
 * Heuristics:
 * - higher-order prefer short bonds (Labute 2005)
 * - higher-order prefer in-ring bonds (Sayle 2001)
 */
function assignSlots(state: State, i: number, specs: NeighbourSpec[]): number[] | undefined {
    const nb = state.heavyNeighbours[i];
    const m = nb.length;
    const s = specs.length;
    const u = state.unitIndices[i];

    // Precompute the slot×neighbour cost matrix once (Infinity = element/connectivity mismatch, forbidden).
    // Heuristics: higher-order prefer short (Labute 2005) and in-ring (Sayle 2001) bonds.
    const cost = new Array<number>(s * m);
    for (let k = 0; k < m; k++) {
        const j = nb[k];
        const elJ = state.el[j];
        const isTerminalHSuppressed = state.heavyNeighbours[j].length === 1;
        const d2 = distSq(state, i, j);
        const inRingJ = state.inRing[j];
        const v = state.unitIndices[j];
        const assigned = isAssignedBond(state, u, v);
        const curOrder = assigned ? getOrder(state.bonds, u, v) : 0;
        for (let si = 0; si < s; si++) {
            const spec = specs[si];
            let c: number;
            if (!matchEl(spec.el, elJ) || !connectivityOk(spec.connectivity, isTerminalHSuppressed)) {
                c = Infinity;
            } else {
                const ringBonus = (spec.order > 1 && inRingJ) ? -2 : 0;
                // If already assigned bond and conflicts with the spec, penalize to allow picking a better slot
                // (e.g. guanidine double bond assignment after amide recognition)
                const preAssignmentPenalty = (assigned && spec.order !== curOrder) ? 1 : 0;
                c = spec.order > 1 ? d2 + ringBonus - preAssignmentPenalty : -d2 + preAssignmentPenalty;
            }
            cost[si * m + k] = c;
        }
    }

    const used = new Array<boolean>(m).fill(false);
    const curAssign = new Array<number>(s).fill(-1);
    let bestCost = Infinity;
    let bestAssign: number[] | undefined;

    const search = (si: number, acc: number) => {
        if (si === s) {
            if (acc < bestCost) { bestCost = acc; bestAssign = curAssign.slice(); }
            return;
        }
        const base = si * m;
        for (let k = 0; k < m; k++) {
            if (used[k]) continue;
            const c = cost[base + k];
            if (!Number.isFinite(c)) continue;
            used[k] = true; curAssign[si] = k;
            search(si + 1, acc + c);
            used[k] = false;
        }
    };
    search(0, 0);
    return bestAssign;
}