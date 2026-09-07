/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Aromaticity perception based on Huckel theory, following the algorithm from Sayle (2001).
 * A ring system with sp2 or sp3-ambiguous atoms is candidate for aromaticity perception.
 * The algorithm tries to find whether a combination of atomic π-electron contributions
 * can satisfy the 4n+2 rule, and if so:
 * - assigns the bonds as aromatic.
 * - assigns the atoms as sp2 (trigonal).
 * - kekulizes the ring system by taking into account the degree of freedom of certain configurations:
 *   e.g. the contribution of a nitrogen can be 1 or 2 electrons, and only one may be compatible with
 *   the 4n+2 rule.
 */

import { Elements } from '../../../../mol-model/structure/model/properties/atomic/types';
import { BondType, ElementSymbol } from '../../../../mol-model/structure/model/types';
import { isDebugMode } from '../../../../mol-util/debug';
import { AtomGeometry } from '../geometry';
import { atomId } from '../util';
import { Ambiguity, State, computeOpenValence, getOrder, hasMultipleBond, isAssignedBond, setBond } from './common';
import { exceedsMaxDoubleLength, hasPlanarDihedral } from './thresholds';

type RingAtomType = { electrons: number, ambiguous: 'C' | 'N' | null };

type ExoSpec =
    | { order: 0 } // no exocyclic heavy bond (implicit H / 2-connected)
    | { order: 1 | 2; el: 'O' | '*' }; // exocyclic single/double; partner O-specific or any heavy

interface AromaticAtomConfig {
    name: string;
    el: Elements; // C | N | O | S | SE
    charge?: 1; // the two N+ motifs
    ringDoubles: 0 | 1; // # in-ring double bonds (each ring atom has exactly 2 ring bonds)
    exo: ExoSpec;
    electrons: number | [number, number]; // π contribution; [lo,hi] = ambiguous
}

// Order matters: most restrictive motifs come before equivalent wildcards
const AromaticAtomConfigs: AromaticAtomConfig[] = [
    // carbon, ringDoubles 0
    { name: '*C*', el: Elements.C, ringDoubles: 0, exo: { order: 0 }, electrons: 1 },
    { name: '*C(O)*', el: Elements.C, ringDoubles: 0, exo: { order: 1, el: 'O' }, electrons: [0, 1] }, // ambiguous if the bond order to exocyclic O is unknown
    { name: '*C(*)*', el: Elements.C, ringDoubles: 0, exo: { order: 1, el: '*' }, electrons: 1 },
    { name: '*C(=O)*', el: Elements.C, ringDoubles: 0, exo: { order: 2, el: 'O' }, electrons: 0 },
    { name: '*C(=*)*', el: Elements.C, ringDoubles: 0, exo: { order: 2, el: '*' }, electrons: 1 },
    // carbon, ringDoubles 1
    { name: '*C=*', el: Elements.C, ringDoubles: 1, exo: { order: 0 }, electrons: 1 },
    { name: '*C(O)=*', el: Elements.C, ringDoubles: 1, exo: { order: 1, el: 'O' }, electrons: 1 },
    { name: '*C(*)=*', el: Elements.C, ringDoubles: 1, exo: { order: 1, el: '*' }, electrons: 1 },
    // nitrogen, ringDoubles 0  (neutral *N(*)* before charged *[N+](*)*)
    { name: '*N*', el: Elements.N, ringDoubles: 0, exo: { order: 0 }, electrons: [1, 2] },
    { name: '*N(=O)*', el: Elements.N, ringDoubles: 0, exo: { order: 2, el: 'O' }, electrons: 1 }, // To be clarified: there are consequences on the N charge
    { name: '*N(*)*', el: Elements.N, ringDoubles: 0, exo: { order: 1, el: '*' }, electrons: 2 },
    // The following seems like a standardization issue, i.e. this config should be post-processed to '*[N+](*)=*' (with a ring double)
    // { name: '*[N+](*)*', el: Elements.N, ringDoubles: 0, exo: { order: 1, el: '*' }, charge: 1, electrons: 1 },
    // nitrogen, ringDoubles 1
    { name: '*N=*', el: Elements.N, ringDoubles: 1, exo: { order: 0 }, electrons: 1 },
    // The following is mentionned in Sayle (2001), but it probably refers to a standardization issue
    // and cannot happen in this code path.
    // See: https://doi.org/10.1186/s13321-018-0293-8
    // { name: '*N(=*)=*', el: Elements.N, ringDoubles: 1, exo: { order: 2, el: '*' }, electrons: 1 },
    { name: '*[N+](*)=*', el: Elements.N, ringDoubles: 1, exo: { order: 1, el: '*' }, charge: 1, electrons: 1 },
    // chalcogens
    { name: '*O*', el: Elements.O, ringDoubles: 0, exo: { order: 0 }, electrons: 2 },
    { name: '*S*', el: Elements.S, ringDoubles: 0, exo: { order: 0 }, electrons: 2 },
    { name: '*[Se]*', el: Elements.SE, ringDoubles: 0, exo: { order: 0 }, electrons: 2 },
];

/**
 * Type a ring atom against `AromaticAtomConfigs` and return its π-electron contribution (Sayle 2001, step 8).
 * There are 2 ambiguous cases:
 * - C with exo-cyclic O. The ketone (C=O) contributes 0, the enol (C-OH) contributes 1.
 * - N with no exo-cyclic heavy atom. With no H, it contributes 2 (lone pair), with H it contributes 1 (pyrrole-like).
 * Returns null when no configuration matches.
 */
function typeRingAtom(state: State, ringSet: Set<number>, i: number): RingAtomType | null {
    const el = state.el[i];

    // honor an explicit formal charge: a pyridinium N⁺ donates 1 electron (one ring double)
    if (el === Elements.N && state.charge[i] === 1) return { electrons: 1, ambiguous: null };

    // in-ring double bonds already on this atom (e.g. placed by a step-7 functional group)
    let ringDoubles = 0;
    // exocyclic bond (at most one for an sp2 ring atom); prefer a double if present
    let exoOrder = 0;
    let exoPartner: ElementSymbol | undefined;
    const u = state.unitIndices[i];

    for (const j of state.allNeighbours[i]) {
        const v = state.unitIndices[j];
        const bo = getOrder(state.bonds, u, v);
        if (ringSet.has(j)) {
            if (bo > 1) ringDoubles++;
            continue;
        }
        // exo-cyclic bond
        exoOrder = bo;
        exoPartner = state.el[j];
    }
    const partnerIsO = exoPartner === Elements.O;

    for (const cfg of AromaticAtomConfigs) {
        if (cfg.el !== el || cfg.ringDoubles !== ringDoubles || cfg.exo.order !== exoOrder) continue;
        if (cfg.exo.order !== 0 && cfg.exo.el !== '*' && !(cfg.exo.el === 'O' && partnerIsO)) continue;
        const e = cfg.electrons;
        if (Array.isArray(e)) return { electrons: 1, ambiguous: el === Elements.C ? 'C' : 'N' };
        return { electrons: e, ambiguous: null };
    }
    return null;
}

/**
 * Filter out ring systems that cannot be aromatic candidates:
 *  - every atom of the set is sp2 (planar) or geometrically ambiguous.
 *  - every ring bond is not an unambiguous single bond.
 */
function isAromaticCandidate(state: State, atoms: number[]): boolean {
    const { geometry, unitIndices, ambiguous, bonds } = state;
    const ringSet = new Set(atoms);
    const seen = new Set<number>();
    const n = atoms.length;
    for (const i of atoms) {
        const u = unitIndices[i];
        if (geometry[i] !== AtomGeometry.Trigonal && ambiguous[i] === Ambiguity.None) return false;

        for (const j of state.heavyNeighbours[i]) {
            if (!ringSet.has(j) || j <= i) continue; // each ring bond once
            const localBondId = n * i + j;
            if (seen.has(localBondId)) continue;
            seen.add(localBondId);
            const v = unitIndices[j];
            // Note: this is more lenient than the hybridization code: both exclusion criteria must be satisfied
            if (isAssignedBond(state,u, v) && getOrder(bonds, u, v) === 1
                && exceedsMaxDoubleLength(state, u, v) && !hasPlanarDihedral(state, u, v)
            ) {
                return false;
            }
        }
    }
    return true;
}

/**
 * Original Sayle (2001) algorithm considers only single aromatic rings (vs aromatic systems). The
 * resolution of the 4n+2 rule depends on the total degree of freedom coming from ambiguous atom types.
 * This cannot be resolved by a greedy single-atom choice, so the algorithm is extended to consider the
 * whole ring system as a unit.
 */
function resolveAromaticDOF(state: State, atoms: number[], ringSet: Set<number>, ringAtomTypes: RingAtomType[]): RingAtomType[] | null {
    const ambIdx: number[] = [];
    let sumFixed = 0, positiveDOF = 0, negativeDOF = 0;
    for (let i = 0; i < atoms.length; i++) {
        if (ringAtomTypes[i].ambiguous) {
            ambIdx.push(i);
            if (ringAtomTypes[i].ambiguous === 'N') positiveDOF ++;
            else negativeDOF ++;
        } else {
            sumFixed += ringAtomTypes[i].electrons;
        }
    }
    const k = ambIdx.length;
    if (k === 0 || k > 16) return null; // nothing to search or too large
    const totalDOF = positiveDOF + negativeDOF;
    if (totalDOF < 4) {
        const targetElectronsDiff = (sumFixed + k) % 4 - 2;
        if (targetElectronsDiff < -negativeDOF || targetElectronsDiff > positiveDOF) return null;
    }

    const altElectrons = (i: number) => ringAtomTypes[i].ambiguous === 'N' ? 2 : 0; // N contributes 1e or 2e, C contributes 1e or 0e

    for (let mask = 0; mask < (1 << k); mask++) {
        const electrons = atoms.map((_, i) => ringAtomTypes[i].ambiguous ? 1 : ringAtomTypes[i].electrons);
        let total = sumFixed;
        for (let j = 0; j < k; j++) {
            const e = (mask & (1 << j)) ? altElectrons(ambIdx[j]) : 1;
            electrons[ambIdx[j]] = e; total += e;
        }
        if (total % 4 !== 2) continue;
        if (!feasibleRingKekule(state, atoms, ringSet, electrons)) continue;
        return atoms.map((_, i) => ({ electrons: electrons[i], ambiguous: null }));
    }
    return null;
}

/**
 * Each atom that contributes 1 electron to the ring must be perfectly matched with another such atom over ring bonds.
 */
function feasibleRingKekule(state: State, atoms: number[], ringSet: Set<number>, electrons: number[]): boolean {
    const idxOf = new Map<number, number>();
    atoms.forEach((a, idx) => idxOf.set(a, idx));
    const needsDouble = (k: number) => electrons[k] === 1 && state.open[atoms[k]] > 0 && !hasMultipleBond(state, atoms[k]);
    const demand: number[] = [];
    for (let k = 0; k < atoms.length; k++) if (needsDouble(k)) demand.push(k);
    if (demand.length % 2 !== 0) return false;
    if (demand.length > 60) return true; // safety valve: assume feasible, let the matcher decide

    const adj = new Map<number, number[]>(); // acceptors for each demand atom
    for (const i of demand) {
        const nbrs: number[] = [];
        for (const v of state.heavyNeighbours[atoms[i]]) {
            if (!ringSet.has(v)) continue;
            const j = idxOf.get(v);
            if (j !== undefined && needsDouble(j)) nbrs.push(j);
        }
        if (nbrs.length === 0) return false; // stranded demand atom
        adj.set(i, nbrs);
    }

    const matched = new Set<number>();
    const rec = (): boolean => {
        let pick = -1, pickDeg = Infinity;
        for (const i of demand) {
            if (matched.has(i)) continue;
            let d = 0; for (const j of adj.get(i)!) if (!matched.has(j)) d++;
            if (d < pickDeg) { pickDeg = d; pick = i; }
        }
        if (pick === -1) return true; // all matched
        if (pickDeg === 0) return false; // stranded demand atom
        for (const j of adj.get(pick)!) {
            if (matched.has(j)) continue;
            matched.add(pick); matched.add(j);
            if (rec()) return true;
            matched.delete(pick); matched.delete(j);
        }
        return false;
    };
    return rec();
}

function perceiveAromaticSet(state: State, atoms: number[]): boolean {
    const { open, unitIndices } = state;
    const ringSet = new Set(atoms);
    // ring neighbours of each atom = its heavy neighbours that are also in the set
    const ringNeighbours = (i: number) => state.heavyNeighbours[i].filter(j => ringSet.has(j));

    const typed: RingAtomType[] = [];
    for (const i of atoms) {
        const t = typeRingAtom(state, ringSet, i);
        if (t === null) return false;
        typed.push(t);
    }

    // A ring atom that cannot accept a double bond with the current source atom
    const cannotAcceptDoubleBond = (acceptor: number, src: number) => open[acceptor] === 0 || hasMultipleBond(state, acceptor) || isAssignedBond(state, unitIndices[acceptor], unitIndices[src]);

    // Sayle (2001) step 8a: Try to resolve ambiguous case. If no neighbour can accept a double bond, then:
    // - ambiguous C-O resolves to the keto form C=O
    // - ambiguous N resolves the pyrole form N-H
    for (let k = 0; k < atoms.length; k++) {
        const t = typed[k];
        if (!t.ambiguous) continue;
        const i = atoms[k];
        const nbs = ringNeighbours(i);
        if (nbs.every(j => cannotAcceptDoubleBond(j, i))) {
            t.electrons = t.ambiguous === 'N' ? 2 : 0;
            t.ambiguous = null;
        }
    }
    let sum = typed.reduce((s, t) => s + t.electrons, 0);

    // Sayle (2001) step 8b: count ≡ 1 (mod 4) and an ambiguous N → set it to pyrrole-like (1→2)
    if (sum % 4 === 1) {
        const t = typed.find(t => t.ambiguous === 'N');
        if (t) { t.electrons = 2; t.ambiguous = null; sum++; }
    }

    // Sayle (2001) step 8c: count ≡ 3 (mod 4) and an ambiguous C → set to keto (1→0)
    if (sum % 4 === 3) {
        const t = typed.find(t => t.ambiguous === 'C');
        if (t) { t.electrons = 0; t.ambiguous = null; sum--; }
    }

    // Sayle (2001) step 8d: count ≡ 3 (mod 4) after 8c and the ring holds an uncharged pyrrole-like N → charge it to a
    // pyridinium (2→1 electron).
    if (sum % 4 === 3) {
        const i = typed.findIndex((t, k) => t.ambiguous === null && t.electrons === 2
            && state.el[atoms[k]] === Elements.N && state.charge[atoms[k]] === 0);
        if (i > -1) { typed[i].electrons = 1; state.charge[atoms[i]] = 1; sum--; }
    }

    // In ring systems, multiple combinations of ambiguous atom assignments may exist to reach 4n+2.
    if (sum % 4 !== 2) {
        const resolved = resolveAromaticDOF(state, atoms, ringSet, typed);
        if (!resolved) return false; // no assignment reaches 4n+2 with a valid Kekulé → not aromatic
        sum = 0;
        for (let i = 0; i < atoms.length; i++) { typed[i] = resolved[i]; sum += typed[i].electrons; }
    }

    // Sayle (2001) step 8e: Hückel 4n+2 → flag every intra-set ring bond as perceived-aromatic. A plain perceivable single
    // bond is (re)set to single so the aromatic Kekulé pass (`matchDoubleBondsBy(isAromaticBond)`) can
    // place the doubles freely; a double already placed by a functional group (e.g. guanidine's in-ring
    // C=N) or an authoritative order is kept, so the matcher fills only the rest around it. Exocyclic
    // doubles (e.g. carbonyls) are never touched - this loop only visits ring bonds.
    isDebugMode && console.log('Aromatic set ', atoms.map(a => atomId(state.unit, state.unitIndices[a])).join('-'), ' π=', sum);
    for (const i of atoms) {
        const u = state.unitIndices[i];
        for (const j of state.heavyNeighbours[i]) {
            if (!ringSet.has(j) || j < i) continue;
            const v = state.unitIndices[j];
            const keptOrder = getOrder(state.bonds, u, v);
            setBond(state.bonds, u, v, keptOrder, BondType.Flag.AromaticHuckel);
        }
    }

    // Enforce the typings with 0 Degree Of Freedom, before yielding to the kekulization procedure.
    // Each atom can contribute 0, 1 or 2 π-electrons to the ring system.
    // The cases of 0 (Carbonyl) or 2 (lone-pair donor: pyrrole N / furan O / thiophene S) are unambiguous,
    // and constrain the ring bonds to be single.
    // Remaining atoms accept a ring double and are resolved by the kekulization procedure.
    for (let k = 0; k < atoms.length; k++) {
        const t = typed[k];
        if (t.electrons === 1) continue; // double bond to be placed by kekulization
        const i = atoms[k];
        const u = state.unitIndices[i];
        for (const j of ringNeighbours(i)) {
            setBond(state.bonds, u, state.unitIndices[j], 1, BondType.Flag.AromaticHuckel, state.assignedBonds);
        }
        if (t.electrons === 2) {
            // pyrrole-like N: enforce implicit H
            if (state.el[i] === Elements.N && state.heavyNeighbours[i].length === 2) state.hCount[i] = 1;
        } else if (!hasMultipleBond(state, i)) {
            // carbonyl (electrons 0): enforce C=O double
            const o = state.heavyNeighbours[i].find(j => !ringSet.has(j) && state.el[j] === Elements.O);
            if (o !== undefined) setBond(state.bonds, u, state.unitIndices[o], 2, BondType.Flag.Computed, state.assignedBonds);
        }
    }
    return true;
}


/**
 * Aromatic perception is tried first on the whole ring system, and if not satisfied, individual rings are considered.
 * isAromaticCandidate ensures that planarity and sp2 requirements are met, before testing that the
 * Hückel 4n+2 rule can be satisfied.
 */
export function perceiveAromaticRings(state: State) {
    computeOpenValence(state);
    for (const system of state.localRingSystems) {
        let ringSystemIsAromatic = false;
        if (isAromaticCandidate(state, system.atoms)) {
            ringSystemIsAromatic = perceiveAromaticSet(state, system.atoms);
        }

        if (!ringSystemIsAromatic && system.rings.length > 1) {
            for (const ring of system.rings) {
                if (isAromaticCandidate(state, ring)) perceiveAromaticSet(state, ring);
            }
        }
    }
}