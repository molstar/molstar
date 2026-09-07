/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Double/Triple bond assignment methods. At this stage, some double/triple bonds have been assigned after
 * functional group detection, Hückel aromatic rings have been assigned and some unambiguous single
 * bonds have been enforced.
 * The algorithm follows steps from Sayle (2001) "PDB: Cruft to Content", but the details of the kekulization
 * algorithm are not described in the paper ("The actual kekule form assignment algorithm is complex").
 * Labute (2005) "On the Perception of Molecules from 3D Atomic Coordinates" describes a maximum weighted matching
 * algorithm, which is loosely implemented here (used only at certain steps, not for the entire bipartite graph matching).
 * - kekulizaton of aromatic rings
 * - assignment of carbonyls (C=O) by distance
 * - kekulization of remaining sp2 centers in all rings
 * - assignment of localized multiple bonds (sp/sp, sp2/sp2, sp/terminal, sp2/terminal)
 * - final pass over ambiguous cases to recover planar double bonds (sp2/sp2, sp2/ambiguous)
 */

import { Elements } from '../../../../mol-model/structure/model/properties/atomic/types';
import { BondType } from '../../../../mol-model/structure/model/types';
import { getElementIdx } from '../../../../mol-model/structure/structure/unit/bonds/common';
import { AtomGeometry } from '../geometry';
import { State, computeOpenValence, distSq, getFlags, hasMultipleBond, isAssignedBond, setBond } from './common';
import { getUnambiguousSingleBondThreshold, hasPlanarDihedral, carbonylBondMaxLengthSq, multipleBondMaxSq } from './thresholds';

const ShortLengthBonus = 0.1;

/**
 * Score of a candidate double i-j (residue-local indices) for the weighted matching.
 * Composite of 2 parts:
 * - number of interior endpoints (vs terminal endpoints): a ring bond scores 2, a ring-terminal bond scores 1.
 * - bonus for bond's length being shorter than the single bond threshold reference (formula from Labute 2005).
 */
function edgeWeight(state: State, i: number, j: number): number {
    const base = (isSp2OrAmbiguousDoubleDemand(state, i) ? 1 : 0) + (isSp2OrAmbiguousDoubleDemand(state, j) ? 1 : 0);
    const t = getUnambiguousSingleBondThreshold(getElementIdx(state.el[i]), getElementIdx(state.el[j]));
    let bonus = 0;
    if (t >= 0) {
        const diff = t - Math.sqrt(distSq(state, i, j)); // how far below the single-bond length
        bonus = diff > 0.25 ? 3 * ShortLengthBonus : diff > 0.1 ? 2 * ShortLengthBonus : 0;
    }
    return base + bonus;
}

/** Ensure the final result of the kekulization is chemically correct by not leaving
 * any sp2 carbon with only single bonds.
 * TODO: add also atoms from aromatic rings typed as contributing 1e to the ring
 */
function mustCoverDouble(state: State, i: number): boolean {
    return state.geometry[i] === AtomGeometry.Trigonal && state.el[i] === Elements.C && state.degree[i] === 3;
}

/**
 * Filter graph nodes/edges to satisfy the predicates and then run a maximum matching algorithm.
 */
function matchDoubleBondsBy(state: State, isDemand: (state: State, i: number) => boolean, canPartner: (state: State, i: number) => boolean, edgeOk: (state: State, i: number, j: number) => boolean) {
    computeOpenValence(state);
    const { bonds, n, unitIndices, heavyNeighbours, assignedBonds } = state;

    const vertices: number[] = [];
    const adj = new Map<number, number[]>();
    for (let i = 0; i < n; i++) {
        if (!isDemand(state, i)) continue;
        vertices.push(i);
        const list: number[] = [];
        for (const j of heavyNeighbours[i]) {
            if (!canPartner(state, j)) continue;
            if (!edgeOk(state, i, j)) continue;
            list.push(j);
        }
        adj.set(i, list);
    }

    const matched = matchDemand(vertices, adj, (i, j) => edgeWeight(state, i, j), (i) => mustCoverDouble(state, i));
    for (const [i, j] of matched) {
        if (i < j) setBond(bonds, unitIndices[i], unitIndices[j], 2, BondType.Flag.Computed, assignedBonds);
    }
}

/**
 * Maximum matching of the demand vertices over `adj`.
 * First step is to split the graph into connected components and then solve each component separately.
 */
function matchDemand(demand: number[], adj: Map<number, number[]>, weightOf: (i: number, j: number) => number, mustCover: (i: number) => boolean): Map<number, number> {
    const matched = new Map<number, number>();
    const inDemand = new Set(demand);
    const seen = new Set<number>();
    for (const s of demand) {
        if (seen.has(s)) continue;
        // collect the connected component of `s` (over demand-internal edges)
        const comp: number[] = [];
        const stack = [s];
        seen.add(s);
        while (stack.length) {
            const i = stack.pop()!;
            comp.push(i);
            for (const j of adj.get(i)!) {
                if (inDemand.has(j) && !seen.has(j)) { seen.add(j); stack.push(j); }
            }
        }
        matchComponent(comp, adj, matched, weightOf, mustCover);
    }
    return matched;
}

/**
 * Exact matching of one component that maximises the number of **demand** atoms matched, by constraint
 * propagation + backtracking. `comp` holds the demand vertices; `adj` may point a demand atom at a
 * partner outside `comp` (a terminal — see `matchDoubleBondsBy`). The objective counts covered demand
 * atoms, not pairs, so a demand–demand double (covers 2) is preferred over a demand–terminal one
 * (covers 1) when they compete; terminal can be left unmatched.
 *
 * Objective is a single scalar `score = Σ edgeWeight(pair)`: its integer part is the covered
 * demand-atom count (so coverage stays primary), and a fractional distance bonus breaks
 * coverage ties toward shorter (more double-like) bonds. The best matching maximises `score`.
 */
function matchComponent(comp: number[], adj: Map<number, number[]>, matched: Map<number, number>, weightOf: (i: number, j: number) => number, mustCover: (i: number) => boolean) {
    if (comp.length > 60) { // safety valve for pathological components
        for (const i of comp) {
            if (matched.has(i)) continue;
            for (const j of adj.get(i)!) { if (!matched.has(j)) { matched.set(i, j); matched.set(j, i); break; } }
        }
        return;
    }

    let best = new Map<number, number>();
    let bestScore = -Infinity; // max score found: (covered demand atoms) + Σ bonus − invalid penalties
    const maxCovered = comp.length; // every demand atom matched
    // a matching that strands a `mustCover` atom (sp2 carbon → chemically invalid) is penalised beyond
    // any achievable coverage+bonus, so it never wins over one that covers all coverable carbons.
    const invalidPenalty = comp.length + 1;
    const cur = new Map<number, number>(); // matched pairs (both directions), may include terminal partners
    const skipped = new Set<number>(); // demand vertices deliberately left unmatched in this branch

    // available partners of `u`: candidate neighbours neither matched nor skipped
    const avail = (i: number) => adj.get(i)!.filter(j => !cur.has(j) && !skipped.has(j));
    // demand atoms (comp members) currently matched
    const covered = () => { let c = 0; for (const i of comp) if (cur.has(i)) c++; return c; };

    // is `i` wanted by an undecided demand other than `j`? Such a contested sole-partner must be branched.
    const contested = (i: number, j: number) => {
        for (const k of comp) {
            if (k === i || cur.has(k) || skipped.has(k)) continue;
            if (avail(k).includes(j)) return true;
        }
        return false;
    };

    function rec() {
        if (bestScore >= maxCovered) return; // all demand covered — coverage can't be beaten

        // unit propagation: force every undecided demand vertex whose sole available partner is
        // uncontested (contested sole-partners are left to weighted branching — see the docstring).
        const forced: number[] = [];
        let changed = true;
        while (changed) {
            changed = false;
            for (const i of comp) {
                if (cur.has(i) || skipped.has(i)) continue;
                const av = avail(i);
                if (av.length === 1 && !contested(i, av[0])) {
                    const j = av[0];
                    cur.set(i, j); cur.set(j, i);
                    forced.push(i, j);
                    changed = true;
                }
            }
        }

        // pick the lowest-available-degree undecided demand vertex with at least one partner left
        let pick = -1, pickDeg = Infinity, pending = 0;
        for (const i of comp) {
            if (cur.has(i) || skipped.has(i)) continue;
            const d = avail(i).length;
            if (d === 0) continue; // already unmatchable; just left unmatched
            pending++;
            if (d < pickDeg) { pickDeg = d; pick = i; }
        }

        if (pick === -1) {
            let score = 0;
            for (const [i, j] of cur) if (i < j) score += weightOf(i, j);
            // reject chemically-invalid results: an sp2 carbon left without a double is invalid
            for (const i of comp) if (!cur.has(i) && mustCover(i)) score -= invalidPenalty;
            if (score > bestScore) { bestScore = score; best = new Map(cur); }
            for (const f of forced) cur.delete(f);
            return;
        }

        // upper bound: even covering every still-matchable demand (`pending`) plus a full bonus (< 1)
        // cannot reach `bestScore` — the +1 slack keeps equal-coverage, higher-bonus branches alive.
        if (covered() + pending + 1 <= bestScore) {
            for (const f of forced) cur.delete(f);
            return;
        }

        for (const i of avail(pick)) {
            cur.set(pick, i); cur.set(i, pick);
            rec();
            cur.delete(pick); cur.delete(i);
            if (bestScore >= maxCovered) { for (const f of forced) cur.delete(f); return; }
        }
        // leave `pick` unmatched
        skipped.add(pick);
        rec();
        skipped.delete(pick);

        for (const f of forced) cur.delete(f);
    }
    rec();
    for (const [i, j] of best) matched.set(i, j);
}

/**
 * Sayle (2001) Step 9c: multiple-bond assignment between sp/sp/terminal or
 * sp2/sp2/terminal that have unfilled valence and no assigned multiple bond.
 * The distance threshold from Sayle (2001) must be satisfied.
 */
function assignLocalizedMultipleBonds(state: State) {
    const { n, bonds, geometry, heavyNeighbours, open, unitIndices } = state;
    computeOpenValence(state);
    type Kind = 'sp' | 'sp2' | 'term' | 'none';
    const kind = (i: number): Kind => {
        if (geometry[i] === AtomGeometry.Linear) return 'sp';
        if (geometry[i] === AtomGeometry.Trigonal || state.ambiguous[i] !== 0) return 'sp2';
        if (geometry[i] === AtomGeometry.Terminal) return 'term';
        return 'none';
    };

    for (let i = 0; i < n; i++) {
        const u = unitIndices[i];
        if (open[i] === 0 || hasMultipleBond(state, i)) continue;
        const ki = kind(i);
        if (ki === 'none') continue;

        for (const j of heavyNeighbours[i]) {
            if (i >= j) continue;
            if (open[j] === 0 || hasMultipleBond(state, j)) continue;
            const v = unitIndices[j];
            const kj = kind(j);
            if (kj === 'none') continue;
            if (isAssignedBond(state, u, v)) continue;

            let orders: number[];
            if (ki === 'sp' && kj === 'sp') orders = [3];
            else if ((ki === 'sp' && kj === 'term') || (ki === 'term' && kj === 'sp')) orders = [3];
            else if (ki === 'sp2' && kj === 'sp2') orders = [2];
            else if ((ki === 'sp2' && kj === 'term') || (ki === 'term' && kj === 'sp2')) orders = [2];
            else if (ki === 'term' && kj === 'term') orders = [3, 2];
            else continue;

            let found = false;
            for (const order of orders) {
                if (open[i] < order - 1 || open[j] < order - 1) continue;
                const maxSq = multipleBondMaxSq(state.el[i], state.el[j], order);
                if (maxSq === undefined || distSq(state, i, j) > maxSq) continue;
                setBond(bonds, u, v, order, BondType.Flag.Computed, state.assignedBonds);
                computeOpenValence(state);
                found = true;
                break;
            }
            if (found) break;
        }
    }
}

/**
 * Sayle (2001) step 9d. Pass that recovers bonds between 2 sp2 that are less certain:
 * - Both atoms are classified as sp2 and the bond was not rejected previously (i.e. not unambiguous single bond).
 * - One atom is classified as sp2, the other is ambiguous (sp3 with 2 neighbours) and the bond length is below the unambiguous double bond threshold.
 */
function recoverPlanarDoubleBonds(state: State) {
    matchDoubleBondsBy(state, isSp2DoubleDemand, isSp2OrAmbiguousDoubleDemand, bondHasSingleGeometricViolation);
}

/**
 * Sayle (2001) Step 8b: place unambiguous double bonds with terminal oxygen using distance threshold.
 * This favors keto vs enol.
 */
function assignCarbonylsByDistance(state: State) {
    computeOpenValence(state);
    const { n, unitIndices, heavyNeighbours, assignedBonds, bonds, geometry, open, el } = state;
    for (let i = 0; i < n; i++) {
        if (el[i] !== Elements.C || !isSp2DoubleDemand(state, i)) continue;
        for (const j of heavyNeighbours[i]) {
            if (geometry[j] !== AtomGeometry.Terminal || open[j] <= 0 || el[j] !== Elements.O) continue;
            if (distSq(state, i, j) > carbonylBondMaxLengthSq) continue;

            setBond(bonds, unitIndices[i], unitIndices[j], 2, BondType.Flag.Computed, assignedBonds);
            computeOpenValence(state);
            break;
        }
    }
}

/**
 * This pass must be run after the aromatic rings kekulization.
 * A pyramidal nitrogen hints at a non-conjugated system. sp2 atoms next to it are probably misassigned.
 * A noticeable exception is a Nitrogen next to an aromatic ring: Nitrogen is not planar due to
 * weaker conjugation.
 */
function preventConjugationNextToPyramidalNitrogen(state: State) {
    const { n, heavyNeighbours, geometry, ambiguous, degree, open, pyramidalNitrogen } = state;
    let found = false;
    for (let i = 0; i < n; i++) {
        if (!pyramidalNitrogen[i]) continue;
        for (const j of heavyNeighbours[i]) {
            const isSp2 = geometry[j] === AtomGeometry.Trigonal || ambiguous[j] > 0;
            if (!isSp2 || open[j] === 0) continue;
            if (degree[j] > 2) continue; // sp2 perception with 3 substituents is probably correct. Degree 2 is ambiguous.
            geometry[j] = AtomGeometry.Tetrahedral;
            state.ambiguous[j] = 0;
            found = true;
        }
    }
    if (found) computeOpenValence(state);
}

// --- atom / edge predicates for the kekulization passes ------------------

/**
 * A *partner*: an atom that can receive one (more) double bond — open valence, no existing multiple
 *  bond, and a suitable geometry: a genuine Trigonal centre, a geometrically `ambiguous` one, or a
 *  terminal atom.
 */
function canAcceptDouble(state: State, i: number) {
    return state.open[i] > 0 && !hasMultipleBond(state, i) &&
        (state.geometry[i] === AtomGeometry.Trigonal || state.ambiguous[i] !== 0 ||
            state.geometry[i] === AtomGeometry.Terminal);
}

/**
 * A *demand*: an atom the matching must try to satisfy with a double — a **non-terminal** sp2 centre
 *  (genuine Trigonal or geometrically `ambiguous`). This is `canAcceptDouble` minus the terminal case
 */
function isSp2OrAmbiguousDoubleDemand(state: State, i: number) {
    return state.open[i] > 0 && !hasMultipleBond(state, i) &&
        (state.geometry[i] === AtomGeometry.Trigonal || state.ambiguous[i] !== 0);
}

/**
 * A *demand*: an atom the matching must try to satisfy with a double — a **non-terminal** sp2 centre
 *  (genuine Trigonal). This is `isSp2OrAmbiguousDoubleDemand` minus the ambiguous case
 */
function isSp2DoubleDemand(state: State, i: number) {
    return state.open[i] > 0 && !hasMultipleBond(state, i) && state.geometry[i] === AtomGeometry.Trigonal;
}

/** The bond was flagged Hückel-aromatic by `perceiveAromaticRings`. */
function isAromaticBond(state: State, i: number, j: number) {
    const u = state.unitIndices[i], v = state.unitIndices[j];
    return (getFlags(state.bonds, u, v) & BondType.Flag.AromaticHuckel) !== 0 && !isAssignedBond(state, u, v);
}

/** Both endpoints lie in a ring within the residue. */
function isRingBond(state: State, i: number, j: number) {
    // TODO: only allow the same ring! Not inter-rings bonds
    const u = state.unitIndices[i], v = state.unitIndices[j];
    return state.inRing[i] !== 0 && state.inRing[j] !== 0 && !isAssignedBond(state, u, v);
}

/**
 * Reevaluates bonds that may have a single geometric violation on the condition
 * that other characteristics are relatively unambiguous.
 * - A bond length violation (too long) between sp2 atoms that have a degree of 3 (no ambiguous dimension) and which is planar.
 * - A bond between one sp2 atom and one ambiguous (checked by the atom predicate), with no geometric violation,
 *   and with a length shorter than the double bond threshold including a 0.05 Å tolerance.
 */
function bondHasSingleGeometricViolation(state: State, i: number, j: number) {
    const { unitIndices, el, ambiguous, degree } = state;
    const u = unitIndices[i], v = unitIndices[j];
    const noAmbiguous = ambiguous[i] === 0 && ambiguous[j] === 0 && degree[i] === 3 && degree[j] === 3;
    if (isAssignedBond(state, u, v)) {
        return noAmbiguous && hasPlanarDihedral(state, u, v);
    }
    return noAmbiguous
        ? true
        : distSq(state, i, j) <= (multipleBondMaxSq(el[i], el[j], 2, 0.05) ?? 0);
}

export function assignBondOrders(state: State) {
    // step 8: kekulize the genuine (Hückel-flagged) aromatic rings first
    matchDoubleBondsBy(state, isSp2OrAmbiguousDoubleDemand, canAcceptDouble, isAromaticBond);

    // step 8b: place distance-confirmed terminal carbonyl-type doubles before the "any ring" kekulé.
    assignCarbonylsByDistance(state);
    preventConjugationNextToPyramidalNitrogen(state);

    // step 9a: kekulize the remaining in-ring sp2 atoms
    matchDoubleBondsBy(state, isSp2DoubleDemand, canAcceptDouble, isRingBond);

    // step 9b/9c: localized terminal/triple bonds (carbonyls, nitriles, etc.).
    assignLocalizedMultipleBonds(state);

    // step 9d: planar fallback for any non-ring / stretched conjugated doubles.
    recoverPlanarDoubleBonds(state);
}