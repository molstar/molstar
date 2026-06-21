/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Perception of bond orders (double / triple / aromatic) from 3D coordinates,
 * following Roger Sayle's "PDB: Cruft to Content" algorithm, steps 6-9
 * (hybridization, functional groups, aromatic perception, bond-order assignment).
 *
 * Connectivity (which atoms are bonded) and ring detection are *not* recomputed here
 * - the existing single-bond graph and `computeRings` are reused. This routine only
 * upgrades the orders of bonds whose order is not authoritative, for residues whose
 * orders are not otherwise known. "Order not authoritative" bonds are flagged
 * `BondType.Flag.Computed`; this covers both distance-computed bonds and PDB `CONECT`
 * / `struct_conn` records that provide basic connectivity with no explicit bond order
 * (these are flagged `Computed` by the struct_conn parser - see `struct_conn.ts`).
 *
 * Only bond orders / aromatic flags are written. Formal charge and implicit-H stay
 * with the `ValenceModel`, which derives them from the perceived orders.
 */

import { Segmentation } from '../../../../../mol-data/int';
import { Vec3 } from '../../../../../mol-math/linear-algebra';
import { degToRad, radToDeg } from '../../../../../mol-math/misc';
import { AtomGeometry } from '../../../../../mol-model-props/computed/chemistry/geometry';
import { type Model } from '../../../model/model';
import { hasIntraBondOrderFromTable } from '../../../model/properties/atomic/bonds';
import { Elements } from '../../../model/properties/atomic/types';
import { BondType, ElementSymbol } from '../../../model/types';
import { type Unit } from '../../unit';
import { UnitRings } from '../rings';
import { type IntraUnitBonds } from './data';
import { type NumberArray } from '../../../../../mol-util/type-helpers';

/** Element spec for functional-group matching: an element symbol or '*' = any heavy atom. */
type ElSpec = Elements | '*';

function isHydrogenElement(el: ElementSymbol) {
    return el === Elements.H || el === Elements.D || el === Elements.T;
}

/** Standard (Daylight-like) heavy valence used for open-valence bookkeeping. -1 = unknown. */
function defaultValence(el: ElementSymbol): number {
    switch (el) {
        case Elements.C: case Elements.SI: return 4;
        case Elements.N: case Elements.P: case Elements.B: return 3;
        case Elements.O: case Elements.S: case Elements.SE: return 2;
        case Elements.F: case Elements.CL: case Elements.BR: case Elements.I:
        case Elements.H: case Elements.D: case Elements.T: return 1;
        default: return -1;
    }
}

// Atom geometry (from bond angles) is represented with the shared `AtomGeometry` enum:
//   Terminal = terminal, Linear = sp, Trigonal = sp2, Tetrahedral = sp3.

// Loose double-bond length cap (1.42^2) for the kekulization passes (step 8 aromatic, step 9a
// conjugated): those bonds run up to ~1.40 A, so the per-element step-9c table below would
// wrongly reject them.
const DoubleMaxSq = 2.0164;

// Step 9c per-bond-type reference lengths (Angstrom), from Sayle "PDB: Cruft to Content". A bond
// of the given order is accepted when its length is at or below the reference. Pairs absent here
// (e.g. C=N) are never assigned a localized multiple bond by 9c.
//   triples: C#C 1.25, C#N 1.22 ; doubles: C=C 1.38, C=O 1.28, C=S 1.70, N=N 1.32
function multipleBondMaxSq(elA: ElementSymbol, elB: ElementSymbol, order: number): number | undefined {
    const a = elA < elB ? elA : elB;
    const b = elA < elB ? elB : elA;
    if (order === 3) {
        if (a === Elements.C && b === Elements.C) return 1.5625; // 1.25^2
        if (a === Elements.C && b === Elements.N) return 1.4884; // 1.22^2
        return undefined;
    }
    if (order === 2) {
        if (a === Elements.C && b === Elements.C) return 1.9044; // 1.38^2
        if (a === Elements.C && b === Elements.O) return 1.6384; // 1.28^2
        if (a === Elements.C && b === Elements.S) return 2.8900; // 1.70^2
        if (a === Elements.N && b === Elements.N) return 1.7424; // 1.32^2
        return undefined;
    }
    return undefined;
}

// Sayle thresholds for geometries (minimal angles rather than ideal angles)
const SP_ANGLE = degToRad(155);
const SP2_ANGLE = degToRad(115);

// Max average |torsion| (degrees) for a ring to count as planar (sp2).
const Torsion5MaxDeg = 7.5;
const Torsion6MaxDeg = 12;
// Max deviation from 0°/180° of the substituent dihedral across a candidate double bond for it
// to count as planar (`isBondPlanar`). Looser than the ring threshold: a real double inside a
// distorted ring still twists locally (e.g. the chlorin C1=C12 dihedral is ~10.6°), while a
// pyramidal sp3 centre sits near ±60°.
const BondPlanarMaxDeg = 30;

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpVecC = Vec3();
const tmpVecD = Vec3();

type Pattern = { a: string, b: string, order: number, flags: number }[];

const PerceivedCache = new WeakMap<Model, Map<string, Pattern>>();

function getModelCache(model: Model) {
    let c = PerceivedCache.get(model);
    if (!c) { c = new Map(); PerceivedCache.set(model, c); }
    return c;
}

/** Set the order (and optionally OR-in flags) of the undirected bond u-v in both directed slots. */
function setBond(bonds: IntraUnitBonds, u: number, v: number, order: number, addFlags: number) {
    const { offset, b, edgeProps } = bonds;
    const ord = edgeProps.order as NumberArray;
    const flg = edgeProps.flags as NumberArray;
    for (let t = offset[u], _t = offset[u + 1]; t < _t; t++) {
        if (b[t] === v) { ord[t] = order; if (addFlags) flg[t] |= addFlags; break; }
    }
    for (let t = offset[v], _t = offset[v + 1]; t < _t; t++) {
        if (b[t] === u) { ord[t] = order; if (addFlags) flg[t] |= addFlags; break; }
    }
}

function getOrder(bonds: IntraUnitBonds, u: number, v: number) {
    const { offset, b, edgeProps: { order } } = bonds;
    for (let t = offset[u], _t = offset[u + 1]; t < _t; t++) {
        if (b[t] === v) return order[t];
    }
    return 0;
}

function getFlags(bonds: IntraUnitBonds, u: number, v: number) {
    const { offset, b, edgeProps: { flags } } = bonds;
    for (let t = offset[u], _t = offset[u + 1]; t < _t; t++) {
        if (b[t] === v) return flags[t];
    }
    return BondType.Flag.None;
}

/**
 * A bond is eligible for order perception if it is a single covalent bond whose order
 * is not authoritative. Such bonds are flagged `BondType.Flag.Computed`, which is set
 * both for distance-computed bonds and for PDB `CONECT` / `struct_conn` records that
 * give only basic connectivity without an explicit order (see `struct_conn.ts`). Bonds
 * with an explicit order (`chem_comp_bond`, struct_conn `doub`/`sing`, mol/mol2/sdf)
 * are not flagged `Computed` and are therefore left untouched.
 */
function isPerceivable(flags: number, order: number) {
    return order === 1 && (flags & BondType.Flag.Computed) !== 0 && BondType.isCovalent(flags);
}

interface State {
    unit: Unit.Atomic
    bonds: IntraUnitBonds
    start: number
    end: number
    x: ArrayLike<number>
    y: ArrayLike<number>
    z: ArrayLike<number>
    /** element symbol per local atom (offset by `start`) */
    el: ElementSymbol[]
    /** attached hydrogen count per local atom */
    hCount: Int8Array
    /** total covalent heavy degree (incl. cross-residue) per local atom */
    heavyDegree: Int8Array
    /** intra-residue heavy neighbours (local UnitIndex) per local atom */
    neighbours: number[][]
    /** geometry (from bond angles) per local atom, as `AtomGeometry` */
    geometry: Int8Array
    /** average heavy-neighbour bond angle (radians) per local atom; 0 for degree < 2 */
    angle: Float32Array
    /** remaining open valence per local atom */
    open: Int8Array
    /** 1 if the atom belongs to a ring within the residue */
    inRing: Uint8Array
    /** `Ambiguity` bit-flags per local atom (geometry that can't be decided) */
    ambiguous: Uint8Array
}

/**
 * Cases where an atom's sp2/sp3 hybridization can't be decided from coordinates alone, so it
 * may take a double if a genuine sp2 neighbour and the bond geometry support it. Bit-flags so
 * more cases can be added later.
 */
const enum Ambiguity {
    None = 0,
    /** degree-2 atom in a non-planar 5-ring: the ~108° angle can't tell sp2 from sp3 */
    Sp2OrSp3InRing5 = 1,
    /** degree-2 atom whose single angle sits just below the sp2 threshold (e.g. an aromatic-ring
     *  aldehyde C at 114.7°): too close to call sp2 vs sp3 from one angle */
    Sp2OrSp3Borderline = 2,
}

function pos(state: State, local: number, out: Vec3) {
    const eI = state.unit.elements[state.start + local];
    return Vec3.set(out, state.x[eI], state.y[eI], state.z[eI]);
}

function distSq(state: State, aIdx: number, bIdx: number) {
    pos(state, aIdx, tmpVecA);
    pos(state, bIdx, tmpVecB);
    return Vec3.squaredDistance(tmpVecA, tmpVecB);
}

/** Returns undefined if there is nothing to perceive in this residue. */
function State(unit: Unit.Atomic, bonds: IntraUnitBonds, start: number, end: number): State | undefined {
    const { offset, b, edgeProps } = bonds;
    const { flags, order } = edgeProps;
    const { type_symbol } = unit.model.atomicHierarchy.atoms;
    const { x, y, z } = unit.model.atomicConformation;
    const n = end - start;

    const el = new Array<ElementSymbol>(n);
    const hCount = new Int8Array(n);
    const heavyDegree = new Int8Array(n);
    const neighbours: number[][] = [];
    let perceivable = false;

    for (let i = 0; i < n; i++) {
        const u = start + i;
        el[i] = type_symbol.value(unit.elements[u]);
        const list: number[] = [];
        for (let t = offset[u], _t = offset[u + 1]; t < _t; t++) {
            if (!BondType.isCovalent(flags[t])) continue;
            const v = b[t];
            const ev = type_symbol.value(unit.elements[v]);
            if (isHydrogenElement(ev)) { hCount[i]++; continue; }
            heavyDegree[i]++;
            if (v >= start && v < end) {
                list.push(v - start);
                if (isPerceivable(flags[t], order[t])) perceivable = true;
            }
        }
        neighbours.push(list);
    }

    if (!perceivable) return undefined;

    return {
        unit, bonds, start, end, x, y, z,
        el, hCount, heavyDegree, neighbours,
        geometry: new Int8Array(n),
        angle: new Float32Array(n),
        open: new Int8Array(n),
        inRing: new Uint8Array(n),
        ambiguous: new Uint8Array(n),
    };
}

// --- Step 6: hybridization -------------------------------------------------

function assignHybridization(state: State) {
    const n = state.end - state.start;
    for (let i = 0; i < n; i++) {
        const nb = state.neighbours[i];
        const deg = state.heavyDegree[i];
        if (deg <= 1) {
            state.geometry[i] = AtomGeometry.Terminal;
            continue;
        }

        // 6a: average / max bond angle over heavy neighbours
        pos(state, i, tmpVecA);
        let sum = 0, max = 0, count = 0;
        for (let a = 0; a < nb.length; a++) {
            pos(state, nb[a], tmpVecB);
            Vec3.sub(tmpVecB, tmpVecB, tmpVecA);
            for (let bIdx = a + 1; bIdx < nb.length; bIdx++) {
                pos(state, nb[bIdx], tmpVecC);
                Vec3.sub(tmpVecC, tmpVecC, tmpVecA);
                const ang = Vec3.angle(tmpVecB, tmpVecC);
                sum += ang; count++;
                if (ang > max) max = ang;
            }
        }
        const avg = count > 0 ? sum / count : 0;
        state.angle[i] = avg;
        if (max > SP_ANGLE) state.geometry[i] = AtomGeometry.Linear;
        else if (avg > SP2_ANGLE) state.geometry[i] = AtomGeometry.Trigonal;
        else state.geometry[i] = AtomGeometry.Tetrahedral;
    }
}

/** 6b: planarity override for 5/6-membered rings entirely within the residue. */
function applyRingPlanarity(state: State, rings: number[][]) {
    for (const seq of rings) {
        const len = seq.length;
        if (len !== 5 && len !== 6) continue;

        // average absolute in-ring torsion
        let sum = 0;
        for (let i = 0; i < len; i++) {
            pos(state, seq[i], tmpVecA);
            pos(state, seq[(i + 1) % len], tmpVecB);
            pos(state, seq[(i + 2) % len], tmpVecC);
            pos(state, seq[(i + 3) % len], tmpVecD);
            sum += Math.abs(radToDeg(Vec3.dihedralAngle(tmpVecA, tmpVecB, tmpVecC, tmpVecD)));
        }
        const avg = sum / len;
        const threshold = len === 5 ? Torsion5MaxDeg : Torsion6MaxDeg;
        if (avg < threshold) {
            for (const a of seq) state.geometry[a] = AtomGeometry.Trigonal;
        }
    }
}

// A degree-2 atom is ambiguous when its single angle falls in this band just below the sp2
// threshold: too close to call sp2 vs sp3 from one angle (an aromatic-ring aldehyde C is ~114.7°,
// a hydroxymethyl ~115.2°). 5° wide; the step-9c distance then decides.
const BorderlineAngleMin = SP2_ANGLE - degToRad(2);

/**
 * 6c: flag degree-2 atoms whose sp2/sp3 hybridization can't be decided from one bond angle, so a
 * missing H hides the real coordination. Two cases:
 *  - in a non-planar 5-ring: the ~108° angle is dictated by ring geometry (e.g. the chlorin
 *    methine C12 in CHR, a real sp2 =CH- that looks Tetrahedral);
 *  - a degree-2 atom whose angle sits just below the sp2 threshold (e.g. an aromatic-ring aldehyde
 *    C at 114.7°, mis-read Tetrahedral).
 * A *planar* 5-ring was already promoted to Trigonal by `applyRingPlanarity`, so anything still
 * Tetrahedral is genuinely undecidable. Run after `applyRingPlanarity`.
 */
function markAmbiguous(state: State, rings: number[][]) {
    for (const seq of rings) {
        if (seq.length !== 5) continue;
        for (const a of seq) {
            if (state.heavyDegree[a] === 2 && state.geometry[a] === AtomGeometry.Tetrahedral) {
                state.ambiguous[a] |= Ambiguity.Sp2OrSp3InRing5;
            }
        }
    }
    const n = state.end - state.start;
    for (let i = 0; i < n; i++) {
        if (state.heavyDegree[i] === 2 && state.geometry[i] === AtomGeometry.Tetrahedral &&
            state.angle[i] >= BorderlineAngleMin) {
            state.ambiguous[i] |= Ambiguity.Sp2OrSp3Borderline;
        }
    }
}

/**
 * Whether the bond u-v looks like a planar (sp2-sp2) double: the substituents on either end are
 * roughly coplanar across the bond. Picks one heavy neighbour on each side (other than the bond
 * partner) and tests the n_u-u-v-n_v torsion against ~0°/180° within `BondPlanarMaxDeg`. Returns
 * true when a side has no such neighbour (planarity can't be disproven, e.g. a terminal partner).
 */
function isBondPlanar(state: State, u: number, v: number) {
    const nu = state.neighbours[u].find(w => w !== v);
    const nv = state.neighbours[v].find(w => w !== u);
    if (nu === undefined || nv === undefined) return true;
    pos(state, nu, tmpVecA);
    pos(state, u, tmpVecB);
    pos(state, v, tmpVecC);
    pos(state, nv, tmpVecD);
    const deg = Math.abs(radToDeg(Vec3.dihedralAngle(tmpVecA, tmpVecB, tmpVecC, tmpVecD)));
    return deg <= BondPlanarMaxDeg || deg >= 180 - BondPlanarMaxDeg;
}

/**
 * order ring atoms indices by cyclic walk
 */
function ringTraversalSort(state: State, ring: number[]): number[] | undefined {
    const len = ring.length;

    const seq: number[] = [ring[0]];
    let current = ring[0];
    for (let step = 1; step < len; step++) {
        let found = false;
        for (const nb of state.neighbours[current]) {
            if (!ring.includes(nb) || seq.includes(nb)) continue;
            found = true;
            seq.push(nb);
            current = nb;
            break;
        }
        if (!found) return undefined; // shoud never happen
    }
    return seq;
}

// --- open valence ----------------------------------------------------------

function computeOpenValence(state: State) {
    const n = state.end - state.start;
    const { bonds, start } = state;
    for (let i = 0; i < n; i++) {
        const val = defaultValence(state.el[i]);
        if (val < 0) { state.open[i] = 0; continue; }
        // sum current bond orders to heavy neighbours (intra + cross-residue counted as >=1)
        let used = state.hCount[i];
        // cross-residue heavy bonds: count as single each
        used += (state.heavyDegree[i] - state.neighbours[i].length);
        for (const nb of state.neighbours[i]) {
            used += getOrder(bonds, start + i, start + nb);
        }
        state.open[i] = Math.max(0, val - used);
    }
}

// --- Step 7: functional groups --------------------------------------------

type NeighborConnectivity = 'terminal' | 'connected' | 'either';
interface neighbourspec { el: ElSpec; connectivity: NeighborConnectivity; order: number; }
// `geometry` is the center's hybridization required by the depiction (Sayle): a carbonyl/
// nitro/amidine center is sp2 (Trigonal), a nitrile / azide-middle is sp (Linear), a
// phosphate / sulfonyl is sp3 (Tetrahedral). It is enforced before slot assignment so e.g. an
// aminomethyl R-CH2-NH2 (Trigonal C) cannot match the nitrile R-C#N (Linear C).
interface FunctionalGroup { name: string; center: { el: ElSpec; geometry: AtomGeometry }; neighbours: neighbourspec[]; }

const FunctionalGroups: FunctionalGroup[] = [
    // *N(*)C=O
    { name: 'formamide', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'either', order: 1 },
        // H?
    ] },
    // CC(=O)O
    { name: 'carboxylic acid', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: Elements.C, connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 1 },
    ] },
    // *C(=O)O*
    { name: 'ester', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'connected', order: 1 },
    ] },
    // *C(=O)S*
    { name: 'thioester', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.S, connectivity: 'connected', order: 1 },
    ] },
    // *C(=O)N* — covers amide (C-C(=O)-N, incl. acetamide with a terminal methyl),
    // urea (N-C(=O)-N, e.g. the C2 carbonyl of uracil/thymine) and carbamate (O-C(=O)-N).
    // Note: this differs from the depiction in the original article which hinted at a carbonyl substituent
    { name: 'amide', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'either', order: 1 },
    ] },
    // *C(=S)S*
    { name: 'dithioester', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.S, connectivity: 'terminal', order: 2 },
        { el: Elements.S, connectivity: 'connected', order: 1 },
    ] },
    // *C(=S)N*
    { name: 'thioamide', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.S, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'connected', order: 1 },
    ] },
    // *NC(=N*)N*
    { name: 'guanidinium', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: Elements.N, connectivity: 'either', order: 2 },
        { el: Elements.N, connectivity: 'either', order: 1 },
        { el: Elements.N, connectivity: 'connected', order: 1 },
    ] },
    // *C(=N)N
    { name: 'amidine', center: { el: Elements.C, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.N, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'terminal', order: 1 },
    ] },
    // *N=[N+]=[N-]
    { name: 'azide', center: { el: Elements.N, geometry: AtomGeometry.Linear }, neighbours: [
        { el: Elements.N, connectivity: 'terminal', order: 2 },
        { el: Elements.N, connectivity: 'connected', order: 2 },
    ] },
    // *P(=O)(*)*
    { name: 'phosphoryl', center: { el: Elements.P, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'either', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
    ] },
    // *S(=O)(*)=O
    { name: 'sulfonyl', center: { el: Elements.S, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: '*', connectivity: 'either', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
    ] },
    // *[N+](=O)[O-]
    { name: 'nitro', center: { el: Elements.N, geometry: AtomGeometry.Trigonal }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 1 },
    ] },
    // *C#N
    { name: 'nitrile', center: { el: Elements.C, geometry: AtomGeometry.Linear }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.N, connectivity: 'terminal', order: 3 },
    ] },
    // *[Se](=O)O
    { name: 'seleninic acid', center: { el: Elements.SE, geometry: AtomGeometry.Tetrahedral }, neighbours: [
        { el: '*', connectivity: 'connected', order: 1 },
        { el: Elements.O, connectivity: 'terminal', order: 2 },
        { el: Elements.O, connectivity: 'terminal', order: 1 },
    ] },
];

type ExoSpec =
    | { order: 0 } // no exocyclic heavy bond (implicit H / 2-connected)
    | { order: 1 | 2; el: 'O' | '*' }; // exocyclic single/double; partner O-specific or any heavy

interface AromaticAtomConfig {
    name: string;
    el: Elements; // C | N | O | S | SE
    charge?: 1; // the two N+ motifs (recorded, not acted on)
    ringDoubles: 0 | 1; // # in-ring double bonds (each ring atom has exactly 2 ring bonds)
    exo: ExoSpec;
    electrons: number | [number, number]; // π contribution; [lo,hi] = ambiguous
}

const AromaticAtomConfigs: AromaticAtomConfig[] = [
    { name: '*C(*)*', el: Elements.C, ringDoubles: 0, exo: { order: 1, el: '*' }, electrons: 1 },
    { name: '*C(*)=*', el: Elements.C, ringDoubles: 1, exo: { order: 1, el: '*' }, electrons: 1 },
    { name: '*C(=*)*', el: Elements.C, ringDoubles: 0, exo: { order: 2, el: '*' }, electrons: 1 },
    { name: '*C*', el: Elements.C, ringDoubles: 0, exo: { order: 0 }, electrons: 1 },
    { name: '*C=*', el: Elements.C, ringDoubles: 1, exo: { order: 0 }, electrons: 1 },
    { name: '*C(O)*', el: Elements.C, ringDoubles: 0, exo: { order: 1, el: 'O' }, electrons: [0,1] },
    { name: '*C(=O)*', el: Elements.C, ringDoubles: 0, exo: { order: 2, el: 'O' }, electrons: 0 },
    { name: '*C(O)=*', el: Elements.C, ringDoubles: 1, exo: { order: 1, el: 'O' }, electrons: 1 },
    { name: '*N=*', el: Elements.N, ringDoubles: 1, exo: { order: 0 }, electrons: 1 },
    { name: '*N(=O)N', el: Elements.N, ringDoubles: 0, exo: { order: 2, el: 'O' }, electrons: 1 },
    { name: '*N(=*)=*', el: Elements.N, ringDoubles: 1, exo: { order: 2, el: '*' }, electrons: 1 },
    { name: '*[N+](*)*', el: Elements.N, ringDoubles: 0, exo: { order: 1, el: '*' }, charge: 1, electrons: 1 },
    { name: '*[N+](*)=*', el: Elements.N, ringDoubles: 1, exo: { order: 1, el: '*' }, charge: 1,electrons: 1 },
    { name: '*N*', el: Elements.N, ringDoubles: 0, exo: { order: 0 }, electrons: [1,2] },
    { name: '*N(*)*', el: Elements.N, ringDoubles: 0, exo: { order: 1, el: '*' }, electrons: 2 },
    { name: '*O*', el: Elements.O, ringDoubles: 0, exo: { order: 0 }, electrons: 2 },
    { name: '*S*', el: Elements.S, ringDoubles: 0, exo: { order: 0 }, electrons: 2 },
    { name: '*[Se]*', el: Elements.SE, ringDoubles: 0, exo: { order: 0 }, electrons: 2 },
];

function matchEl(spec: ElSpec, el: ElementSymbol) { return spec === '*' ? !isHydrogenElement(el) : spec === el; }

function connectivityOk(c: NeighborConnectivity, isTerminal: boolean) {
    return c === 'either' || (c === 'terminal' ? isTerminal : !isTerminal);
}

function applyFunctionalGroups(state: State) {
    const n = state.end - state.start;
    for (let i = 0; i < n; i++) {
         // If the center already has a multiple bond (e.g. a backbone C=O or phosphate
        // P=O from the order table), its pi system is already placed - don't add another.
        if (hasMultipleBond(state, i)) continue;
        const el = state.el[i];
        const nb = state.neighbours[i];

        // Skip common cases that are not in patterns
        if (el === Elements.C && state.geometry[i] !== AtomGeometry.Trigonal) continue;
        if (el === Elements.O) continue;

        for (const g of FunctionalGroups) {
            if (g.center.el !== el) continue;
            // The depiction fixes the center hybridization; enforce it so e.g. an aminomethyl
            // R-CH2-NH2 (Trigonal C) cannot match the nitrile R-C#N (Linear C).
            if (state.geometry[i] !== g.center.geometry) continue;
            if (nb.length !== g.neighbours.length) continue; // exact heavy degree
            const slots = assignSlots(state, i, g.neighbours);
            if (!slots) continue;
            for (let k = 0; k < slots.length; k++) {
                const v = nb[slots[k]];
                const order = g.neighbours[k].order;
                // only upgrade bonds whose order is not authoritative
                if (order <= 1 || !isPerceivableBond(state, i, v)) continue;
                setBond(state.bonds, state.start + i, state.start + v, order, BondType.Flag.Computed);
            }
            break;
        }
    }
    computeOpenValence(state);
}

/**
 * Assign the center's real neighbours to ordered slots by (element, connectivity).
 * Returns, for each slot, the index into `state.neighbours[i]`; interchangeable slots are
 * disambiguated by giving higher orders to shorter bonds. undefined if no full match.
 */
function assignSlots(state: State, i: number, specs: neighbourspec[]): number[] | undefined {
    const nb = state.neighbours[i];
    const m = nb.length;
    const used = new Array<boolean>(m).fill(false);
    const result = new Array<number>(specs.length).fill(-1);

    // order slots so that higher-order, more-constrained slots are filled first
    const slotOrder = specs.map((_, idx) => idx).sort((a, bIdx) => specs[bIdx].order - specs[a].order);

    for (const si of slotOrder) {
        const spec = specs[si];
        // candidate neighbours matching element + connectivity, not yet used
        let best = -1, bestDistSq = Infinity;
        for (let k = 0; k < m; k++) {
            if (used[k]) continue;
            const v = nb[k];
            if (!matchEl(spec.el, state.el[v])) continue;
            const isTerminal = state.heavyDegree[v] <= 1;
            if (!connectivityOk(spec.connectivity, isTerminal)) continue;
            const d2 = distSq(state, i, v);
            // higher orders prefer shorter bonds; lower orders prefer longer.
            // Ring neighbours beat non-ring neighbours for double-bond slots (Sayle step 7
            // heuristic): a ring atom takes its pi bond in-ring. -2 undercuts any d² > 0.
            const ringBonus = (spec.order > 1 && state.inRing[v]) ? -2 : 0;
            const key = spec.order > 1 ? d2 + ringBonus : -d2;
            if (key < bestDistSq) { bestDistSq = key; best = k; }
        }
        if (best < 0) return undefined;
        used[best] = true;
        result[si] = best;
    }
    return result;
}

// --- Step 8 / 9a: aromatic perception + Kekule ----------------------------

function perceiveAromaticRings(state: State, rings: number[][]) {
    for (const seq of rings) {
        const len = seq.length;
        if (len !== 5 && len !== 6) continue;
        // all ring atoms must be sp2 (trigonal) and ring bonds still perceivable
        let ok = true;
        for (const a of seq) {
            if (state.geometry[a] !== AtomGeometry.Trigonal) { ok = false; break; }
        }
        if (!ok) continue;
        // mark ring bonds aromatic (only the perceivable ones)
        for (let i = 0; i < len; i++) {
            const u = seq[i], v = seq[(i + 1) % len];
            const flags = getFlags(state.bonds, state.start + u, state.start + v);
            if (isPerceivable(flags, getOrder(state.bonds, state.start + u, state.start + v))) {
                setBond(state.bonds, state.start + u, state.start + v, 1, BondType.Flag.Aromatic);
            }
        }
    }
}

// --- Kekule + remaining double/triple via maximum matching ----------------

/** Whether atom `i` already has an incident multiple (double/triple) intra bond. */
function hasMultipleBond(state: State, i: number) {
    for (const v of state.neighbours[i]) {
        if (getOrder(state.bonds, state.start + i, state.start + v) > 1) return true;
    }
    return false;
}

function isPerceivableBond(state: State, u: number, v: number) {
    const flags = getFlags(state.bonds, state.start + u, state.start + v);
    return isPerceivable(flags, getOrder(state.bonds, state.start + u, state.start + v));
}

/**
 * One maximum-matching pass that assigns double bonds. `bondEligible` restricts which
 * single bonds may become double, letting us run aromatic ring bonds first (step 8)
 * before the remaining bonds (step 9a).
 */
function matchDoubleBonds(state: State, bondEligible: (u: number, v: number) => boolean, requireSp2 = true) {
    computeOpenValence(state);
    const n = state.end - state.start;

    const isSp2 = (i: number) =>
        state.geometry[i] === AtomGeometry.Trigonal || state.geometry[i] === AtomGeometry.Terminal;

    // An atom can accept a double bond if it has open valence and no existing multiple bond.
    // For step 8 (aromatic Kekule) both endpoints must be sp2/terminal. For step 9a (conjugated /
    // ring kekulization) an endpoint must be a genuine sp2 (Trigonal) or a geometrically
    // undecidable ambiguous atom — Terminal is excluded here (carbonyls and other terminal-atom
    // doubles are handled by the step-9c distance table), as are genuine sp3 (Tetrahedral) and sp
    // (Linear) atoms.
    const canAccept = (i: number) =>
        state.open[i] > 0 && !hasMultipleBond(state, i) &&
        (requireSp2
            ? isSp2(i)
            : state.geometry[i] === AtomGeometry.Trigonal || state.ambiguous[i] !== 0);

    // vertices = all atoms that can accept one double bond; carbons first so
    // heteroatoms tend to be left unmatched (lone-pair donors, e.g. pyrrole N / furan O)
    const vertices: number[] = [];
    for (let i = 0; i < n; i++) if (canAccept(i)) vertices.push(i);
    vertices.sort((a, bIdx) => (state.el[a] === Elements.C ? 0 : 1) - (state.el[bIdx] === Elements.C ? 0 : 1));

    // candidate adjacency: perceivable bonds of double-bond length between two acceptors
    const adj = new Map<number, number[]>();
    for (const u of vertices) {
        const list: number[] = [];
        for (const v of state.neighbours[u]) {
            if (!canAccept(v)) continue;
            if (!isPerceivableBond(state, u, v)) continue;
            const d2 = distSq(state, u, v);
            if (d2 > DoubleMaxSq) continue;
            if (!bondEligible(u, v)) continue;
            if (!requireSp2) {
                // A π bond needs a genuine sp2 (Trigonal) anchor; an ambiguous atom is never an
                // anchor, so two undecidable atoms can't pair up on their own.
                const uTri = state.geometry[u] === AtomGeometry.Trigonal;
                const vTri = state.geometry[v] === AtomGeometry.Trigonal;
                if (!uTri && !vTri) continue;
                // An ambiguous endpoint must additionally look like a planar (sp2-sp2) double:
                // its substituents are coplanar across the bond (e.g. the chlorin C1=C12).
                if ((state.ambiguous[u] || state.ambiguous[v]) && !isBondPlanar(state, u, v)) continue;
            }
            list.push(v);
        }
        adj.set(u, list);
    }

    const matched = matchDemand(vertices, adj);
    for (const [u, v] of matched) {
        if (u < v) setBond(state.bonds, state.start + u, state.start + v, 2, BondType.Flag.Computed);
    }
}

function kekulize(state: State) {
    const { start, end, bonds, geometry, neighbours, open } = state;
    const n = end - start;

    // step 8: assign aromatic ring bonds first so a ring atom consumes its pi bond
    // in-ring (and cannot then donate a spurious double to an exocyclic neighbor)
    const isAromaticBond = (u: number, v: number) =>
        (getFlags(bonds, start + u, start + v) & BondType.Flag.Aromatic) !== 0;
    matchDoubleBonds(state, isAromaticBond);

    // step 9a: kekulize remaining conjugated / non-aromatic-ring double bonds (non-terminal sp2
    // or ambiguous atoms) by maximum matching.
    matchDoubleBonds(state, () => true, false);

    // step 9c: final distance pass for bonds still with unfilled valence and no incident multiple
    // bond. Which reference length to test is chosen by hybridization (Sayle): sp–sp / sp–terminal
    // → triple; sp2–sp2 / sp2–terminal → double; terminal–terminal → both. Element-specific
    // thresholds come from `multipleBondMaxSq`; pairs absent from the table are never assigned.
    // (This approximates Sayle's full refinement of remaining valences, not an exhaustive search.)
    computeOpenValence(state);
    type Kind = 'sp' | 'sp2' | 'term' | 'none';
    const kind = (i: number): Kind => {
        if (geometry[i] === AtomGeometry.Linear) return 'sp';
        if (geometry[i] === AtomGeometry.Trigonal || state.ambiguous[i] !== 0) return 'sp2';
        if (geometry[i] === AtomGeometry.Terminal) return 'term';
        return 'none';
    };
    for (let u = 0; u < n; u++) {
        for (const v of neighbours[u]) {
            if (u >= v) continue;
            if (hasMultipleBond(state, u) || hasMultipleBond(state, v)) continue;
            if (!isPerceivableBond(state, u, v)) continue;
            const ku = kind(u), kv = kind(v);
            if (ku === 'none' || kv === 'none') continue;
            // candidate orders by hybridization pair
            let orders: number[];
            if (ku === 'sp' && kv === 'sp') orders = [3];
            else if ((ku === 'sp' && kv === 'term') || (ku === 'term' && kv === 'sp')) orders = [3];
            else if (ku === 'sp2' && kv === 'sp2') orders = [2];
            else if ((ku === 'sp2' && kv === 'term') || (ku === 'term' && kv === 'sp2')) orders = [2];
            else if (ku === 'term' && kv === 'term') orders = [3, 2];
            else continue; // mixed sp/sp2 — not covered by the table rules
            for (const order of orders) {
                if (open[u] < order - 1 || open[v] < order - 1) continue;
                const maxSq = multipleBondMaxSq(state.el[u], state.el[v], order);
                if (maxSq === undefined || distSq(state, u, v) > maxSq) continue;
                // an ambiguous endpoint forming a double must look like a planar sp2-sp2 bond
                // (no-op when the partner is terminal — there is no second substituent to test)
                if (order === 2 && (state.ambiguous[u] || state.ambiguous[v]) && !isBondPlanar(state, u, v)) continue;
                setBond(bonds, start + u, start + v, order, BondType.Flag.Computed);
                computeOpenValence(state);
                break;
            }
        }
    }
}

/**
 * Maximum matching of the demand vertices over `adj`. Solved per connected component: a maximum
 * matching of a disjoint union is the union of the components' maximum matchings, so this keeps
 * the exact result while reducing the exponential backtracking to small independent subproblems
 * (e.g. a ligand's several separate aromatic rings instead of all their atoms at once).
 */
function matchDemand(demand: number[], adj: Map<number, number[]>): Map<number, number> {
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
            const u = stack.pop()!;
            comp.push(u);
            for (const v of adj.get(u)!) {
                if (inDemand.has(v) && !seen.has(v)) { seen.add(v); stack.push(v); }
            }
        }
        matchComponent(comp, adj, matched);
    }
    return matched;
}

/** Exact maximum matching of one component via depth-first backtracking with a (loose) bound. */
function matchComponent(comp: number[], adj: Map<number, number[]>, matched: Map<number, number>) {
    if (comp.length > 60) { // safety valve for pathological components
        for (const u of comp) {
            if (matched.has(u)) continue;
            for (const v of adj.get(u)!) { if (!matched.has(v)) { matched.set(u, v); matched.set(v, u); break; } }
        }
        return;
    }

    let best = new Map<number, number>();
    let bestSize = -1;
    const maxSize = Math.floor(comp.length / 2);
    const cur = new Map<number, number>();

    function rec(idx: number, size: number) {
        if (bestSize >= maxSize) return; // can't do better than a perfect matching
        if (size + (comp.length - idx) <= bestSize) return; // prune
        if (idx >= comp.length) {
            if (size > bestSize) { bestSize = size; best = new Map(cur); }
            return;
        }
        const u = comp[idx];
        if (cur.has(u)) { rec(idx + 1, size); return; }
        for (const v of adj.get(u)!) {
            if (cur.has(v)) continue;
            cur.set(u, v); cur.set(v, u);
            rec(idx + 1, size + 1);
            cur.delete(u); cur.delete(v);
            if (bestSize >= maxSize) break;
        }
        rec(idx + 1, size); // leave u unmatched
    }
    rec(0, 0);
    for (const [u, v] of best) matched.set(u, v);
}

// --- pattern cache (compute once per residue type) ------------------------

function getCompHash(unit: Unit.Atomic, start: number, end: number, compId: string) {
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const names: string[] = [];
    for (let i = start; i < end; i++) names.push(label_atom_id.value(unit.elements[i]));
    names.sort();
    return `${compId}|${names.join(',')}`;
}

function extractPattern(state: State): Pattern {
    const { label_atom_id } = state.unit.model.atomicHierarchy.atoms;
    const pattern: Pattern = [];
    const n = state.end - state.start;
    for (let i = 0; i < n; i++) {
        for (const v of state.neighbours[i]) {
            if (i >= v) continue;
            const order = getOrder(state.bonds, state.start + i, state.start + v);
            const flags = getFlags(state.bonds, state.start + i, state.start + v);
            if (order > 1 || (flags & BondType.Flag.Aromatic)) {
                pattern.push({
                    a: label_atom_id.value(state.unit.elements[state.start + i]),
                    b: label_atom_id.value(state.unit.elements[state.start + v]),
                    order, flags: flags & (BondType.Flag.Aromatic),
                });
            }
        }
    }
    return pattern;
}

function applyPattern(unit: Unit.Atomic, bonds: IntraUnitBonds, start: number, end: number, pattern: Pattern) {
    if (pattern.length === 0) return;
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const nameToLocal = new Map<string, number>();
    for (let i = start; i < end; i++) nameToLocal.set(label_atom_id.value(unit.elements[i]), i);
    for (const p of pattern) {
        const u = nameToLocal.get(p.a);
        const v = nameToLocal.get(p.b);
        if (u === undefined || v === undefined) continue;
        const flags = getFlags(bonds, u, v);
        if (!isPerceivable(flags, getOrder(bonds, u, v))) continue;
        setBond(bonds, u, v, p.order, BondType.Flag.Computed | p.flags);
    }
}

// --- entry point -----------------------------------------------------------

/**
 * Perceive and assign bond orders in place on `bonds` for residues of `unit` whose
 * orders are not otherwise known. Mutates the `order`/`flags` edge properties.
 */
export function perceiveBondOrders(unit: Unit.Atomic, bonds: IntraUnitBonds) {
    if (unit.elements.length <= 1) return;

    const model = unit.model;
    const { label_comp_id } = model.atomicHierarchy.atoms;
    const cache = getModelCache(model);

    // Rings of the whole unit, computed from the freshly-built graph (not unit.bonds,
    // which isn't assigned yet). Cache the result in `unit.props.rings` so the later
    // `unit.rings` getter reuses it instead of recomputing the same rings.
    const unitRings = UnitRings.create(unit, bonds);
    unit.props.rings = unitRings;
    const allRings = unitRings.all;
    // Build atom → ring-index map directly from allRings to avoid triggering
    // the lazy UnitRings.index getter (which needs aromaticRings → unit.bonds, circular).
    const atomToRings = new Map<number, number[]>();
    for (let ri = 0; ri < allRings.length; ri++) {
        const ring = allRings[ri];
        for (let i = 0; i < ring.length; i++) {
            const u = ring[i];
            let list = atomToRings.get(u);
            if (!list) atomToRings.set(u, list = []);
            list.push(ri);
        }
    }
    const seenRings = new Set<number>();

    const residuesIt = Segmentation.transientSegments(model.atomicHierarchy.residueAtomSegments, unit.elements);
    while (residuesIt.hasNext) {
        const { start, end } = residuesIt.move();
        if (end - start < 2) continue;

        const compId = label_comp_id.value(unit.elements[start]);
        if (hasIntraBondOrderFromTable(compId)) continue;

        const compHash = getCompHash(unit, start, end, compId);
        const cached = cache.get(compHash);
        if (cached) {
            applyPattern(unit, bonds, start, end, cached);
            continue;
        }

        const state = State(unit, bonds, start, end);
        if (!state) {
            cache.set(compHash, []); continue;
        }

        // rings restricted to this residue, as cyclic-walk-ordered local indices;
        const localRings5or6: number[][] = [];
        for (let u = start; u < end; u++) {
            const ringsForAtom = atomToRings.get(u);
            if (!ringsForAtom) continue;
            for (const ri of ringsForAtom) {
                if (seenRings.has(ri)) continue;
                seenRings.add(ri);
                const ring = allRings[ri];
                // residue containment: ring is sorted, [start, end] is a contiguous range
                if (ring[0] < start || ring[ring.length - 1] >= end) continue;
                const local = new Array<number>(ring.length);
                for (let i = 0; i < ring.length; i++) {
                    const resIdx = ring[i] - start;
                    local[i] = resIdx;
                    state.inRing[resIdx] = 1;
                }
                if (ring.length < 5 || ring.length > 6) continue;
                const seq = ringTraversalSort(state, local);
                if (seq) localRings5or6.push(seq);
            }
        }

        assignHybridization(state);
        applyRingPlanarity(state, localRings5or6);
        markAmbiguous(state, localRings5or6);
        applyFunctionalGroups(state);
        perceiveAromaticRings(state, localRings5or6);
        kekulize(state);

        cache.set(compHash, extractPattern(state));
    }
}
