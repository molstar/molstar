/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Assigns hybridization to atoms based on local geometry. Follows Sayle (2001) method based on angles
 * and Labute (2005) method for measuring atoms' local point cloud dimensionality (PCA of the atom + its heavy and hydrogen neighbours).
 * This identifies straightforward cases (linear, trigonal and tetrahedral with 3+ neighbours).
 * The remaining cases (e.g. a C/N atom with 2 neighbours) are marked as ambiguous instead of
 * assigning hybridization based on angle only. This diverges from Sayle's method.
 * Further assessments are made to reject unambiguous single bonds:
 * - bond length above threshold (Labute 2005)
 * - non planar dihedral across the bond (Labute 2005)
 * A special case is made for planar rings, following Sayle (2001) method to promote ambiguous sp3 atoms
 * to sp2 if the average ring torsion is below a threshold (7.5°). Note that in the original work
 * 2 different thresholds were used for 5- and 6-membered rings (7.5° and 12° respectively). It
 * seems that the 12º threshold was overfitting to accomodate the disputed 2TRM benzamidine structure,
 * which is actually non-planar (see discussion in Labute 2005).
 */

import { isDebugMode } from '../../../../mol-util/debug';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { PrincipalAxes } from '../../../../mol-math/linear-algebra/matrix/principal-axes';
import { degToRad, radToDeg } from '../../../../mol-math/misc';
import { UnitIndex } from '../../../../mol-model/structure/structure/element/element';
import { AtomGeometry, geometryLabel } from '../geometry';
import { atomId, eachBondedAtom, eachIntraBondedAtom } from '../util';
import { Ambiguity, State, isBondFromMainAltLoc, isAssignedBond, isHydrogenElement, pos, setAssignedBond } from './common';
import { exceedsMaxDoubleLength, hasPlanarDihedral } from './thresholds';
import { Elements } from '../../../../mol-model/structure/model/properties/atomic/types';

// Per-atom dimensionality: threshold from which the principal component of the coordinates of
// connected atoms is deemed significant: 1D, 2D or 3D.
const PcaDimThresholdSquared = 0.2 * 0.2;

// Max average |torsion| (degrees) for a ring to count as planar (sp2).
const Torsion5MaxDeg = 7.5;

// Sayle (2001) thresholds for geometries (minimal angles rather than ideal angles)
const SP_ANGLE = degToRad(155);
const SP2_ANGLE = degToRad(115);

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpVecC = Vec3();
const tmpVecD = Vec3();

/**
 * Labute per-atom dimensionality: PCA of the local point cloud (atom `i` + its heavy and hydrogen
 * neighbours). Returns the number of significant principal components (1 = linear, 2 = planar, 3 = 3-D)
 * using the *scatter*-matrix √-eigenvalues (RMS spread × √nPoints) against `PcaDimThreshold`.
 */
function atomDimensionality(state: State, i: number): number {
    const { structure, unit, unitIndices } = state;
    const u = unitIndices[i];
    const coords: number[] = [];
    pos(state, i, tmpVecA);
    Vec3.toArray(tmpVecA, coords, 0);
    let nPoints = 1;
    eachBondedAtom(structure, unit, u, (otherUnit, v) => {
        if (!isBondFromMainAltLoc(state, unit, u, otherUnit, v)) return;
        const vElementIndex = otherUnit.elements[v];
        otherUnit.conformation.position(vElementIndex, tmpVecA);
        Vec3.toArray(tmpVecA, coords, coords.length);
        nPoints++;
    });
    const { dirA, dirB, dirC } = PrincipalAxes.calculateMomentsAxes(coords);
    // `dir*` magnitudes are sqrt(W / (nPoints/3)); scale back
    const n3 = nPoints / 3;
    let dims = 0;
    if (Vec3.squaredMagnitude(dirA) * n3 > PcaDimThresholdSquared) dims++;
    if (Vec3.squaredMagnitude(dirB) * n3 > PcaDimThresholdSquared) dims++;
    if (Vec3.squaredMagnitude(dirC) * n3 > PcaDimThresholdSquared) dims++;
    return dims;
}

function assignGeometry(state: State) {
    const { structure, unit, n, degree, geometry, el, unitIndices } = state;
    for (let i = 0; i < n; i++) {
        const u = unitIndices[i];
        const deg = degree[i];
        if (deg === 0) {
            geometry[i] = AtomGeometry.Spherical;
            continue;
        }

        if (deg === 1) {
            geometry[i] = AtomGeometry.Terminal;
            if (isHydrogenElement(el[i])) {
                setAssignedBond(state, u);
            }
            continue;
        }

        // With ≥3 total neighbours (heavy + present H) the local point cloud can be 3-D, so PCA
        // dimensionality (Labute) reliably separates pyramidal sp3 from planar sp2.
        if (deg >= 3) {
            const dims = atomDimensionality(state, i);
            if (dims > 2) {
                geometry[i] = AtomGeometry.Tetrahedral;
                if (el[i] === Elements.N) markPyramidalNitrogen(state, i);
            } else if (dims > 1) {
                geometry[i] = deg > 3 ? AtomGeometry.SquarePlanar : AtomGeometry.Trigonal;
            } else {
                geometry[i] = AtomGeometry.Linear;
            }
            continue;
        }

        // With 2 neighbours dimensionality is pointless,
        // so fall back to the angle test with ambiguities.
        // Sayle's 6a step is a measure of the average bond angle over neighbours, but here
        // we are down to a single possible angle.
        pos(state, i, tmpVecA);
        let p = 0;
        const tmpVecs: Vec3[] = [tmpVecB, tmpVecC];
        eachBondedAtom(structure, unit, u, (otherUnit, v) => {
            if (!isBondFromMainAltLoc(state, unit, u, otherUnit, v)) return;
            const tmpVec = tmpVecs[p++];
            if (tmpVec === void 0) return;
            otherUnit.conformation.position(otherUnit.elements[v], tmpVec);
            Vec3.sub(tmpVec, tmpVec, tmpVecA);
        });
        const avg = Vec3.angle(tmpVecB, tmpVecC);

        if (avg > SP_ANGLE) state.geometry[i] = AtomGeometry.Linear;
        else if (avg > SP2_ANGLE) state.geometry[i] = AtomGeometry.Trigonal;
        else {
            state.geometry[i] = AtomGeometry.Tetrahedral;
            state.ambiguous[i] |= Ambiguity.Sp2OrSp3Borderline;
        }
    }
}

/**
 * A Nitrogen located in a conjugated system tends to have a planar geometry.
 * If a Nitrogen is tetrahedral, this could be used as an argument to prevent
 * over-perception of nearby atoms' hybridization.
 * Note that a nearby aromatic ring forms a weaker conjugation effect with an
 * adjacent Nitrogen: it is not planar in that case.
 * Note that a Nitrogen that is strained (e.g. in bridged rings or fused saturated rings),
 * may be pyramidal even if it is conjugated.
 */
function markPyramidalNitrogen(state: State, i: number) {
    const { geometry, degree, el, inRing, unit, unitIndices } = state;
    if (geometry[i] !== AtomGeometry.Tetrahedral && degree[i] < 3 && el[i] !== Elements.N) return;
    if (inRing[i]) {
        const u = unitIndices[i];
        const ringsCount = unit.rings.elementRingIndices.get(u)?.length ?? 0;
        if (ringsCount > 1) return; // N may be strained in fused or bridged rings: pyramidalization does not reflect hybridization
    }
    state.pyramidalNitrogen[i] = 1;
}

/**
 * Sayle (2001) step 6b: planarity override for 5/6-membered rings.
 * Hybridization based on angles is unreliable in rings. E.g. in a 5 membered all-sp2 ring,
 * the average bond angle is ~108° which falls in the tetrahedral range.
 * Sayle recovers those by checking the average torsion across the rings bonds. If it is
 * below 7.5º, the ring is deemed planar and all atoms are promoted to sp2.
 * In this implementation, if any of the atoms is unambiguous sp3, the ring is ignored.
 */
function applyRingPlanarity(state: State) {
    const { localRings5or6, unit, unitIndices, geometry, degree } = state;
    for (const seq of localRings5or6) {
        const len = seq.length;
        if (len !== 5 && len !== 6) continue;

        // Check that all rings atoms are either trigonal or tetrahedral with degree 2 (hence ambiguous)
        let foundNonTrigonal = false;
        let nbTrigonal = 0;
        isDebugMode && console.log('Ring ', seq.map(i => atomId(unit, unitIndices[i])).join('-'), ' planarity check');
        for (const a of seq) {
            const geom = geometry[a];
            if (geom === AtomGeometry.Trigonal) {
                nbTrigonal++;
                continue;
            }
            if (geom === AtomGeometry.Tetrahedral && degree[a] === 2) continue;
            isDebugMode && console.log('Ring has non-trigonal atom ', atomId(unit, unitIndices[a]), ' geometry ', geometryLabel(geom), ' degree ', degree[a]);
            foundNonTrigonal = true;
            break;
        }
        if (foundNonTrigonal) continue;
        if (nbTrigonal === len) continue; // all trigonal, nothing to do

        // Average torsion across the ring's bonds.
        let sum = 0;
        pos(state, seq[0], tmpVecA);
        pos(state, seq[1], tmpVecB);
        pos(state, seq[2], tmpVecC);

        for (let i = 0; i < len; i++) {
            pos(state, seq[(i + 3) % len], tmpVecD);
            const deg = Math.abs(radToDeg(Vec3.dihedralAngle(tmpVecA, tmpVecB, tmpVecC, tmpVecD)));
            sum += Math.min(deg, 180 - deg);
            // Rotate
            Vec3.copy(tmpVecA, tmpVecB);
            Vec3.copy(tmpVecB, tmpVecC);
            Vec3.copy(tmpVecC, tmpVecD);
        }
        const avg = sum / len;
        // Note: Sayle's 6-ring threshold is 12° which allows 2TRM disputed benzamidine to be perceived as planar, when it has an actual cyclo-hexene chair like conformation.
        // Here we use the more conservative 5-ring threshold for both 5-and 6-rings.
        isDebugMode && console.log('Ring average torsion ', avg.toFixed(2), '°');
        if (avg < Torsion5MaxDeg) {
            isDebugMode && console.log('Ring ', seq.map(i => atomId(unit, unitIndices[i])).join('-'), ' is planar, setting trigonal geometry');
            for (const a of seq) geometry[a] = AtomGeometry.Trigonal;
        }
    }
}


/**
 * Apply marking over atoms and bonds that are ambiguous or unabiguous:
 * - Rings of 5 atoms not passing the planarity test: angle is ~108° which is assigned to tetrahedral at the geometry stage. If the atom
 * has only 2 neighbours, it is marked as ambiguous.
 * - Bonds too long or for which no planar dihedral can be found are marked as unambiguous single bonds.
 */
function markAmbiguous(state: State) {
    for (const seq of state.localRings5or6) {
        if (seq.length !== 5) continue;
        for (const a of seq) {
            if (state.degree[a] === 2 && state.geometry[a] === AtomGeometry.Tetrahedral) {
                state.ambiguous[a] |= Ambiguity.Sp2OrSp3InRing5;
            }
        }
    }

    // Labute's conservative tests: planarity OR single-bond-length test to all perceivable bonds, and mark those that fail as "assigned" (excluded from double/aromatic candidacy).
    const seen = new Set<UnitIndex>();
    const { unit, n, unitIndices } = state;
    for (let i = 0; i < n; i++) {
        const u = unitIndices[i];
        if (state.geometry[i] !== AtomGeometry.Trigonal && state.ambiguous[i] === Ambiguity.None) {
            seen.add(u);
            isDebugMode && console.log('Atom ', atomId(state.unit, u), ' is not trigonal, skipping');
            continue;
        }
        eachIntraBondedAtom(unit, u, (unitB, v) => {
            if (unitB.id !== unit.id) return; // TODO: inter-unit bonds
            if (seen.has(v)) return;
            if (isAssignedBond(state, u, v)) return;
            if (!isBondFromMainAltLoc(state, unit, u, unitB, v)) return;

            if (exceedsMaxDoubleLength(state, u, v)) {
                setAssignedBond(state, u, v);
                isDebugMode && console.log('Assigned bond ', atomId(state.unit, u), '-', atomId(state.unit, v), ' as single (too long)');
                return;
            }

            const j = unitIndices.indexOf(v);
            if (state.degree[j] > 1 // Can't compute dihedral to terminal atom.
                && !hasPlanarDihedral(state, u, v)
            ) {
                setAssignedBond(state, u, v);
                isDebugMode && console.log('Assigned bond ', atomId(state.unit, u), '-', atomId(state.unit, v), ' as single (non-planar)');
                return;
            }
        });
        seen.add(u);
    }

    // "Antialiasing": check for stranded sp or sp2 atoms. Demote to sp2 or sp3.
    for (let i = 0; i < n; i++) {
        if (state.geometry[i] !== AtomGeometry.Linear) continue;
        let found = false;
        for (const j of state.heavyNeighbours[i]) {
            if (state.geometry[j] === AtomGeometry.Terminal || state.geometry[j] === AtomGeometry.Linear) {
                found = true;
                break;
            }
        }
        if (!found) {
            state.geometry[i] = AtomGeometry.Trigonal;
        }
    }

    for (let i = 0; i < n; i++) {
        if (state.geometry[i] !== AtomGeometry.Trigonal && (state.geometry[i] !== AtomGeometry.Tetrahedral && state.ambiguous[i] === 0)) continue;
        let found = false;
        for (const j of state.heavyNeighbours[i]) {
            if (state.geometry[j] === AtomGeometry.Terminal || state.geometry[j] === AtomGeometry.Trigonal || (state.geometry[j] === AtomGeometry.Tetrahedral && state.ambiguous[j] > 0)) {
                found = true;
                break;
            }
        }
        if (!found) {
            // No trigonal or suitable tetrahedral neighbor found, demote to sp3
            state.geometry[i] = AtomGeometry.Tetrahedral;
            state.ambiguous[i] = 0;
        }
    }
}

export function assignHybridization(state: State) {
    assignGeometry(state);
    applyRingPlanarity(state);
    markAmbiguous(state);
}