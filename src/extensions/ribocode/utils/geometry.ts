/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
import { Vec3 } from '../../../mol-math/linear-algebra';
import { QCProt } from './qcprot';
import Big from 'big.js';

/**
 * Compute the centroid of a set of coordinates.
 * @param coordsOrX Array of Vec3 coordinates or X coordinates.
 * @param y Array of Y coordinates.
 * @param z Array of Z coordinates.
 * @returns Centroid as Vec3.
 */
export function computeCentroid(
    coordsOrX: Vec3[] | number[],
    y?: number[],
    z?: number[]): Vec3 {
    const centroid: Vec3 = Vec3.zero();
    let n: number;
    if (Array.isArray(coordsOrX[0])) {
        // coords: Vec3[]
        const coords = coordsOrX as Vec3[];
        n = coords.length;
        for (let i = 0; i < n; i++) {
            centroid[0] += coords[i][0];
            centroid[1] += coords[i][1];
            centroid[2] += coords[i][2];
        }
    } else {
        const x = coordsOrX as number[];
        n = x.length;
        for (let i = 0; i < n; i++) {
            centroid[0] += x[i];
            centroid[1] += y![i];
            centroid[2] += z![i];
        }
    }
    centroid[0] /= n;
    centroid[1] /= n;
    centroid[2] /= n;
    return centroid;
}

/**
 * Compute the centroid of a set of coordinates using Big.js for precision.
 * Accepts either an array of Vec3 or three arrays (x, y, z).
 * @param coordsOrX Array of Vec3 coordinates or X coordinates.
 * @param y Array of Y coordinates (optional).
 * @param z Array of Z coordinates (optional).
 * @returns Centroid as Vec3
 */
export function computeBigCentroid(
    coordsOrX: Vec3[] | number[],
    y?: number[],
    z?: number[]
): Vec3 {
    let sumX = new Big(0), sumY = new Big(0), sumZ = new Big(0);
    let n: number;
    if (Array.isArray(coordsOrX[0])) {
        // coords: Vec3[]
        const coords = coordsOrX as Vec3[];
        n = coords.length;
        for (let i = 0; i < n; i++) {
            sumX = sumX.plus(coords[i][0]);
            sumY = sumY.plus(coords[i][1]);
            sumZ = sumZ.plus(coords[i][2]);
        }
    } else {
        const x = coordsOrX as number[];
        n = x.length;
        for (let i = 0; i < n; i++) {
            sumX = sumX.plus(x[i]);
            sumY = sumY.plus(y![i]);
            sumZ = sumZ.plus(z![i]);
        }
    }
    const len = new Big(n);
    const centroid = Vec3.zero();
    centroid[0] = sumX.div(len).toNumber();
    centroid[1] = sumY.div(len).toNumber();
    centroid[2] = sumZ.div(len).toNumber();
    return centroid;
}

/**
 * Compute the mean of an array of numbers using Big.js for precision.
 * @param x Array of numbers
 * @returns Mean as number
 */
export function computeBigMean(x: number[]): number {
    let sumX = new Big(0);
    for (const xi of x) {
        sumX = sumX.plus(xi);
    }
    const len = new Big(x.length);
    return sumX.div(len).toNumber();
}

/**
 * Calculate rotated coordinates given a rotation matrix.
 * @param rotmat The rotation matrix as a flat array of 9 numbers.
 * @param x The x coordinates.
 * @param y The y coordinates.
 * @param z The z coordinates.
 * @returns An object containing the rotated coordinates {xR, yR, zR}.
 */
export function getRotatedCoordinates(rotmat: number[], x: number[], y: number[], z: number[]): { xR: number[]; yR: number[]; zR: number[] } {
    const xR: number[] = [];
    const yR: number[] = [];
    const zR: number[] = [];
    for (let i = 0; i < x.length; i++) {
        xR[i] = rotmat[0] * x[i] + rotmat[1] * y[i] + rotmat[2] * z[i];
        yR[i] = rotmat[3] * x[i] + rotmat[4] * y[i] + rotmat[5] * z[i];
        zR[i] = rotmat[6] * x[i] + rotmat[7] * y[i] + rotmat[8] * z[i];
    }
    return { xR, yR, zR };
}

/**
 * Calculate rotated coordinates given a rotation matrix.
 * @param rotmat The rotation matrix as a flat array of 9 numbers.
 * @param x The x coordinates.
 * @param y The y coordinates.
 * @param z The z coordinates.
 * @returns An object containing the rotated coordinates {xR, yR, zR}.
 */
export function getBigRotatedCoordinates(rotmat: number[], x: number[], y: number[], z: number[]): { xR: number[]; yR: number[]; zR: number[] } {
    const xR: number[] = [];
    const yR: number[] = [];
    const zR: number[] = [];
    for (let i = 0; i < x.length; i++) {
        xR.push(new Big(rotmat[0]).times(x[i])
            .plus(new Big(rotmat[1]).times(y[i]))
            .plus(new Big(rotmat[2]).times(z[i])).toNumber());
        yR.push(new Big(rotmat[3]).times(x[i])
            .plus(new Big(rotmat[4]).times(y[i]))
            .plus(new Big(rotmat[5]).times(z[i])).toNumber());
        zR.push(new Big(rotmat[6]).times(x[i])
            .plus(new Big(rotmat[7]).times(y[i]))
            .plus(new Big(rotmat[8]).times(z[i])).toNumber());
    }
    return { xR, yR, zR };
}

/**
 * Generalized alignment function for datasets.
 * @param moving Array of atom objects for the moving dataset
 * @param reference Array of atom objects for the reference dataset
 * @param selector Function or object to select atoms for alignment
 * @param coordKeys Object specifying coordinate keys (default: {x:'x',y:'y',z:'z'})
 * @returns Aligned coordinates as {alignedX, alignedY, alignedZ, rotmat, centroid}
 */
export function alignDatasetsGeneral<T extends Record<string, any>>(
    moving: T[],
    reference: T[],
    selector: ((atom: T) => boolean) | { [key: string]: any },
    coordKeys: { x: string; y: string; z: string } = { x: 'x', y: 'y', z: 'z' }
): { alignedX: number[]; alignedY: number[]; alignedZ: number[]; rotmat: number[]; centroidReference: Vec3; centroid: Vec3 } {
    // Build selector function
    let select: (atom: T) => boolean;
    if (typeof selector === 'function') {
        select = selector as (atom: T) => boolean;
    } else {
        select = (atom: T) => Object.entries(selector).every(([k, v]) => atom[k] === v);
    }
    // Indices of selected atoms
    const movingSel = moving.map((a, i) => select(a) ? i : -1).filter(i => i !== -1);
    const refSel = reference.map((a, i) => select(a) ? i : -1).filter(i => i !== -1);
    const newCount = movingSel.length;
    const count = refSel.length;
    // Extract coordinates
    const getCoords = (arr: T[], idxs: number[]) => [
        idxs.map(i => arr[i][coordKeys.x]),
        idxs.map(i => arr[i][coordKeys.y]),
        idxs.map(i => arr[i][coordKeys.z])
    ];
    const [newXS, newYS, newZS] = getCoords(moving, movingSel);
    const [xS, yS, zS] = getCoords(reference, refSel);
    if (newCount === 0 || count === 0) {
        throw new Error('No atoms found for alignment with the specified selector.');
    }
    const nAtomsNew = moving.length;
    // Prepare full coordinate arrays
    const fullX = moving.map(a => a[coordKeys.x]);
    const fullY = moving.map(a => a[coordKeys.y]);
    const fullZ = moving.map(a => a[coordKeys.z]);
    // Alignment logic (same as before, generalized)
    let rotmat: number[] = [];
    if (newCount !== count) {
        const filteredXS: number[] = [];
        const filteredYS: number[] = [];
        const filteredZS: number[] = [];
        if (newCount > count) {
            for (let i = 0; i < count; i++) {
                filteredXS.push(newXS[i]);
                filteredYS.push(newYS[i]);
                filteredZS.push(newZS[i]);
            }
            const filteredCentroid = computeBigCentroid(filteredXS, filteredYS, filteredZS);
            for (let i = 0; i < count; i++) {
                filteredXS[i] -= filteredCentroid[0];
                filteredYS[i] -= filteredCentroid[1];
                filteredZS[i] -= filteredCentroid[2];
            }
            const centroid = computeBigCentroid(xS, yS, zS);
            for (let i = 0; i < count; i++) {
                xS[i] -= centroid[0];
                yS[i] -= centroid[1];
                zS[i] -= centroid[2];
            }
            const qcprot = new QCProt(filteredXS, filteredYS, filteredZS, xS, yS, zS);
            rotmat = qcprot.rotmat;
            for (let i = 0; i < nAtomsNew; i++) {
                fullX[i] -= filteredCentroid[0];
                fullY[i] -= filteredCentroid[1];
                fullZ[i] -= filteredCentroid[2];
            }
            const aligned = getRotatedCoordinates(qcprot.rotmat, fullX, fullY, fullZ);
            for (let i = 0; i < nAtomsNew; i++) {
                aligned.xR[i] += filteredCentroid[0];
                aligned.yR[i] += filteredCentroid[1];
                aligned.zR[i] += filteredCentroid[2];
            }
            return { alignedX: aligned.xR, alignedY: aligned.yR, alignedZ: aligned.zR, rotmat, centroid: filteredCentroid, centroidReference: centroid };
        } else {
            for (let i = 0; i < newCount; i++) {
                filteredXS.push(xS[i]);
                filteredYS.push(yS[i]);
                filteredZS.push(zS[i]);
            }
            const filteredCentroid = computeBigCentroid(filteredXS, filteredYS, filteredZS);
            for (let i = 0; i < newCount; i++) {
                filteredXS[i] -= filteredCentroid[0];
                filteredYS[i] -= filteredCentroid[1];
                filteredZS[i] -= filteredCentroid[2];
            }
            const newCentroid = computeBigCentroid(newXS, newYS, newZS);
            for (let i = 0; i < newCount; i++) {
                newXS[i] -= newCentroid[0];
                newYS[i] -= newCentroid[1];
                newZS[i] -= newCentroid[2];
            }
            const qcprot = new QCProt(newXS, newYS, newZS, filteredXS, filteredYS, filteredZS);
            rotmat = qcprot.rotmat;
            for (let i = 0; i < nAtomsNew; i++) {
                fullX[i] += newCentroid[0];
                fullY[i] += newCentroid[1];
                fullZ[i] += newCentroid[2];
            }
            const aligned = getRotatedCoordinates(qcprot.rotmat, fullX, fullY, fullZ);
            for (let i = 0; i < nAtomsNew; i++) {
                aligned.xR[i] -= newCentroid[0];
                aligned.yR[i] -= newCentroid[1];
                aligned.zR[i] -= newCentroid[2];
            }
            return { alignedX: aligned.xR, alignedY: aligned.yR, alignedZ: aligned.zR, rotmat, centroid: newCentroid, centroidReference: filteredCentroid };
        }
    } else {
        const newCentroid = computeBigCentroid(newXS, newYS, newZS);
        for (let i = 0; i < newCount; i++) {
            newXS[i] -= newCentroid[0];
            newYS[i] -= newCentroid[1];
            newZS[i] -= newCentroid[2];
        }
        const centroid = computeBigCentroid(xS, yS, zS);
        for (let i = 0; i < count; i++) {
            xS[i] -= centroid[0];
            yS[i] -= centroid[1];
            zS[i] -= centroid[2];
        }
        const qcprot = new QCProt(newXS, newYS, newZS, xS, yS, zS);
        rotmat = qcprot.rotmat;
        const aligned = getRotatedCoordinates(qcprot.rotmat, fullX, fullY, fullZ);
        return { alignedX: aligned.xR, alignedY: aligned.yR, alignedZ: aligned.zR, rotmat, centroidReference: newCentroid, centroid: centroid };
    }
}

/**
 * Align incoming dataset to existing dataset using centroid translation and atom type subsetting.
 * @param symbol_type Array of atom types of atoms to align (moving)
 * @param newX X coordinates to align (moving)
 * @param newY Y coordinates to align (moving)
 * @param newZ Z coordinates to align (moving)
 * @param type Array of atom types of data to align with (reference)
 * @param x X coordinates to align with (reference)
 * @param y Y coordinates to align with (reference)
 * @param z Z coordinates to align with (reference)
 * @param selectedAtomTypes Object of atom types to use for alignment (e.g., { 'P': true })
 * @return Aligned coordinates as {alignedX, alignedY, alignedZ, rotmat}
 */
export function alignDataset(
    symbol_type: string[],
    newX: number[],
    newY: number[],
    newZ: number[],
    type: string[],
    x: number[],
    y: number[],
    z: number[],
    selectedAtomTypes: { [atomType: string]: boolean } = { 'P': true }
): { alignedX: number[]; alignedY: number[]; alignedZ: number[]; rotmat: number[]; centroid: Vec3 } {
    // Build atom objects for both datasets
    const moving = symbol_type.map((atomType, i) => ({ type: atomType, x: newX[i], y: newY[i], z: newZ[i] }));
    const reference = type.map((atomType, i) => ({ type: atomType, x: x[i], y: y[i], z: z[i] }));
    return alignDatasetsGeneral(moving, reference, (a) => selectedAtomTypes[a.type]);
}

/**
 * Align incoming dataset to existing dataset using centroid translation and chain/atom type subsetting.
 * @param symbol_type Array of atom types of atoms to align (moving)
 * @param chain_ids_new Array of chain IDs for atoms to align (moving)
 * @param newX X coordinates to align (moving)
 * @param newY Y coordinates to align (moving)
 * @param newZ Z coordinates to align (moving)
 * @param type Array of atom types of data to align with (reference)
 * @param chain_ids_ref Array of chain IDs for atoms to align with (reference)
 * @param x X coordinates to align with (reference)
 * @param y Y coordinates to align with (reference)
 * @param z Z coordinates to align with (reference)
 * @param selectedChainId The chain ID to use for alignment
 * @param selectedAtomTypes Object of atom types to use for alignment (e.g., { 'P': true })
 * @return Aligned coordinates as {alignedX, alignedY, alignedZ, rotmat, centroid}
 */
export function alignDatasetUsingChains(
    selectedAtomTypes: { [atomType: string]: boolean } = { 'P': true },
    toAlignSelectedChainID: string,
    toAlignTypes: string[],
    toAlignChainIDs: string[],
    toAlignXs: number[],
    toAlignYs: number[],
    toAlignZs: number[],
    selectedChainId: string,
    symbolTypes: string[],
    chainIds: string[],
    xs: number[],
    ys: number[],
    zs: number[]
): { alignedX: number[]; alignedY: number[]; alignedZ: number[]; rotmat: number[]; centroid: Vec3; centroidReference: Vec3 } {
    const indexesToAlign: number[] = [];
    // Collect all atom indices from all atomic units
    // Get indexes of atoms to align
    for (let i = 0; i < toAlignChainIDs.length; i++) {
        if (toAlignChainIDs[i] === toAlignSelectedChainID) {
            indexesToAlign.push(i);
        }
    }
    // Get indexes of atoms aligned with
    const indexesRef: number[] = [];
    for (let i = 0; i < chainIds.length; i++) {
        if (chainIds[i] === selectedChainId) {
            indexesRef.push(i);
        }
    }
    return alignDatasetsGeneral(
        indexesToAlign.map(i => ({ type: toAlignTypes[i], x: toAlignXs[i], y: toAlignYs[i], z: toAlignZs[i] })),
        indexesRef.map(i => ({ type: symbolTypes[i], x: xs[i], y: ys[i], z: zs[i] })),
        (a) => selectedAtomTypes[a.type]
    );
}