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
 * @param coords Array of Vec3 coordinates
 * @returns Centroid as Vec3
 */
export function computeCentroid(coords: Vec3[]): Vec3 {
    const centroid = Vec3.zero();
    for (let i = 0; i < coords.length; i++) {
        centroid[0] += coords[i][0];
        centroid[1] += coords[i][1];
        centroid[2] += coords[i][2];
    }
    centroid[0] /= coords.length;
    centroid[1] /= coords.length;
    centroid[2] /= coords.length;
    return centroid;
}

/**
 * Compute the centroid of a set of coordinates.
 * @param coords Array of BigVec3 coordinates.
 * @returns Centroid as Vec3
 */
export function computeBigCentroid(coords: Vec3[]): Vec3 {
    let sumX = new Big(0), sumY = new Big(0), sumZ = new Big(0);
    for (const v of coords) {
        sumX = sumX.plus(v[0]);
        sumY = sumY.plus(v[1]);
        sumZ = sumZ.plus(v[2]);
    }
    const len = new Big(coords.length);
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
export function computeBigMean(x: number[]) : number {
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
    let xR: number[] = [];
    let yR: number[] = [];
    let zR: number[] = [];
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
    let xR: number[] = [];
    let yR: number[] = [];
    let zR: number[] = [];
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
 * Align incoming dataset to existing dataset using centroid translation.
 * @param symbol_type Array of atom types of atoms to align
 * @param newX X coordinates to align
 * @param newY Y coordinates to align
 * @param newZ Z coordinates to align
 * @param type Array of atom types of data to align with
 * @param x X coordinates to align with
 * @param y Y coordinates to align with
 * @param z Z coordinates to align with
 * @return Aligned coordinates as {alignedX, alignedY, alignedZ}
 */
export function alignDataset(
    symbol_type: string[],
    newX: number[],
    newY: number[],
    newZ: number[],
    type: string[],
    x: number[],
    y: number[],
    z: number[]
): { alignedX: number[]; alignedY: number[]; alignedZ: number[] } {
    const nAtomsNew = newX.length;
    //const nAtoms = x.length;
    // Get indices of selected atoms
    const newSIndices: number[] = [];
    let newCount: number = 0;
    // Change this to the desired atom types for alignment
    const selectedAtomTypes = { 'P': true };
    //const selectedAtomTypes = { 'P': true, 'S': true };
    for (let i = 0; i < symbol_type.length; i++) {
        if (selectedAtomTypes.hasOwnProperty(symbol_type[i])) {
            newSIndices.push(i);
            //console.log(`Selected atom found at index ${i}`);
            newCount++;
        }
    }
    console.log(`Total selected atoms found in data to align: ${newCount}`);
    const sIndices: number[] = [];
    let count: number = 0;
    for (let i = 0; i < type.length; i++) {
        if (selectedAtomTypes.hasOwnProperty(type[i])) {
            sIndices.push(i);
            //console.log(`Selected atom found at index ${i}`);
            count++;
        }
    }
    console.log(`Total selected atoms found in data to align from: ${count}`);
    // Create arrays of selected coordinates
    const newXS: number[] = [];
    const newYS: number[] = [];
    const newZS: number[] = [];
    for (const idx of newSIndices) {
        newXS.push(newX[idx]);
        newYS.push(newY[idx]);
        newZS.push(newZ[idx]);
    }
    const xS: number[] = [];
    const yS: number[] = [];
    const zS: number[] = [];
    for (const idx of sIndices) {
        xS.push(x[idx]);
        yS.push(y[idx]);
        zS.push(z[idx]);
    }
    if (newCount !== count) {
        const filteredXS: number[] = [];
        const filteredYS: number[] = [];
        const filteredZS: number[] = [];
        if (newCount > count) {
            // // Select the count closest new selected types of atoms to the origin.
            // // 1. Calculate distances to origin of new selected atoms.
            // const distances: { index: number; distance: number }[] = [];
            // for (let i = 0; i < newCount; i++) {
            //     const distance = newXS[i] * newXS[i] + newYS[i] * newYS[i] + newZS[i] * newZS[i];
            //     distances.push({ index: i, distance: distance });
            // }
            // // 2. Sort distances and select closest/furthest count.
            // distances.sort((a, b) => a.distance - b.distance); //closest
            // //distances.sort((a, b) => b.distance - a.distance); // furthest
            // const selectedIndices = distances.slice(0, count).map(d => d.index);
            // //console.log(`Selected types of atom indices for alignment: ${selectedIndices}`);
            // // 3. Add to filtered.
            // for (const idx of selectedIndices) {
            //     filteredXS.push(newXS[idx]);
            //     filteredYS.push(newYS[idx]);
            //     filteredZS.push(newZS[idx]);
            // }
            // console.log(`Filtered selected types of atoms to ${filteredXS.length} for alignment.`);
            
            // Simply select the first count atoms.
            for (let i = 0; i < count; i++) {
                filteredXS.push(newXS[i]);
                filteredYS.push(newYS[i]);
                filteredZS.push(newZS[i]);
            }

            // 4. Recentralise filtered.
            const filteredXMean = computeBigMean(filteredXS);
            const filteredYMean = computeBigMean(filteredYS);
            const filteredZMean = computeBigMean(filteredZS);
            for (let i = 0; i < count; i++) {
                filteredXS[i] -= filteredXMean;
                filteredYS[i] -= filteredYMean;
                filteredZS[i] -= filteredZMean;
            }
            // 5. Recentralise.
            const xSMean = computeBigMean(xS);
            const ySMean = computeBigMean(yS);
            const zSMean = computeBigMean(zS);
            for (let i = 0; i < count; i++) {
                xS[i] -= xSMean;
                yS[i] -= ySMean;
                zS[i] -= zSMean;
            }
            console.log(`Recentralised ${xS.length} atoms for alignment.`); 
            // 6. Create new aligned coordinates arrays.
            const qcprot = new QCProt(filteredXS, filteredYS, filteredZS, xS, yS, zS);
            // Recentre new coordinates.
            for (let i = 0; i < nAtomsNew; i++) {
                newX[i] -= filteredXMean;
                newY[i] -= filteredYMean;
                newZ[i] -= filteredZMean;
            }
            const aligned = getRotatedCoordinates(qcprot.rotmat, newX, newY, newZ);
            // Recentre new coordinates.
            for (let i = 0; i < nAtomsNew; i++) {
                aligned.xR[i] += filteredXMean;
                aligned.yR[i] += filteredYMean;
                aligned.zR[i] += filteredZMean;
            }
            return { alignedX: aligned.xR, alignedY: aligned.yR, alignedZ: aligned.zR };
        } else {
            // // Select the newCount closest selected types of atoms to the centroid.
            // // 1. Calculate distances to origin of selected types of atoms.
            // const distances: { index: number; distance: number }[] = [];
            // for (let i = 0; i < count; i++) {
            //     const distance = xS[i] * xS[i] + yS[i] * yS[i] + zS[i] * zS[i];
            //     distances.push({ index: i, distance: distance });
            // }
            // // 2. Sort distances and select closest/furthest newCount.
            // //distances.sort((a, b) => a.distance - b.distance); // closest
            // distances.sort((a, b) => b.distance - a.distance); // furthest
            // const selectedIndices = distances.slice(0, newCount).map(d => d.index);
            // //console.log(`Selected selected atom indices for alignment: ${selectedIndices}`);
            // // 3. Add to filtered.
            // for (const idx of selectedIndices) {
            //     filteredXS.push(xS[idx]);
            //     filteredYS.push(yS[idx]);
            //     filteredZS.push(zS[idx]);
            // }
            // console.log(`Filtered selected types of atoms to ${filteredXS.length} for alignment.`);
            
            // Simply select the first newCount atoms.
            for (let i = 0; i < newCount; i++) {
                filteredXS.push(xS[i]);
                filteredYS.push(yS[i]);
                filteredZS.push(zS[i]);
            }

            // 4. Recentralise filtered.
            const filteredXMean = computeBigMean(filteredXS);
            const filteredYMean = computeBigMean(filteredYS);
            const filteredZMean = computeBigMean(filteredZS);
            for (let i = 0; i < newCount; i++) {
                filteredXS[i] -= filteredXMean;
                filteredYS[i] -= filteredYMean;
                filteredZS[i] -= filteredZMean;
            }
            // 5. Recentralise.
            const newXPMean = computeBigMean(newXS);
            const newYPMean = computeBigMean(newYS);
            const newZPMean = computeBigMean(newZS);
            for (let i = 0; i < newCount; i++) {
                newXS[i] -= newXPMean;
                newYS[i] -= newYPMean;
                newZS[i] -= newZPMean;
            }
            console.log(`Recentralised ${newXS.length} atoms for alignment.`);
            // 5. Create new aligned coordinates arrays.
            const qcprot = new QCProt(newXS, newYS, newZS, filteredXS, filteredYS, filteredZS);
            // Recentre new coordinates.
            for (let i = 0; i < nAtomsNew; i++) {
                newX[i] += newXPMean;
                newY[i] += newYPMean;
                newZ[i] += newZPMean;
            }
            //const aligned = getBigRotatedCoordinates(qcprot.rotmat, newX, newY, newZ);
            const aligned = getRotatedCoordinates(qcprot.rotmat, newX, newY, newZ);
            // Recentre new coordinates.
            for (let i = 0; i < nAtomsNew; i++) {
                aligned.xR[i] -= newXPMean;
                aligned.yR[i] -= newYPMean;
                aligned.zR[i] -= newZPMean;
            }
            return { alignedX: aligned.xR, alignedY: aligned.yR, alignedZ: aligned.zR };
        }
    } else {
        // Recentralise.
        const newXMean = computeBigMean(newXS);
        const newYMean = computeBigMean(newYS);
        const newZMean = computeBigMean(newZS);
        for (let i = 0; i < newCount; i++) {
            newXS[i] -= newXMean;
            newYS[i] -= newYMean;
            newZS[i] -= newZMean;
        }
        const xSMean = computeBigMean(xS);
        const ySMean = computeBigMean(yS);
        const zSMean = computeBigMean(zS);
        for (let i = 0; i < count; i++) {
            xS[i] -= xSMean;
            yS[i] -= ySMean;
            zS[i] -= zSMean;
        }
        const qcprot = new QCProt(newXS, newYS, newZS, xS, yS, zS);
        const aligned = getRotatedCoordinates(qcprot.rotmat, newX, newY, newZ);
        return { alignedX: aligned.xR, alignedY: aligned.yR, alignedZ: aligned.zR };
    }
}