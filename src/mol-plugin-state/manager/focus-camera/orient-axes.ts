/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { Mat3, Vec3 } from '../../../mol-math/linear-algebra';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { Structure } from '../../../mol-model/structure';


/** Minimum number of atoms necessary for running PCA.
 * If enough atoms cannot be selected, XYZ axes will be used instead of PCA axes. */
const MIN_ATOMS_FOR_PCA = 3;

/** Rotation matrices for the basic rotations by 90 degrees */
export const ROTATION_MATRICES = {
    // The order of elements in the matrices in column-wise (F-style)
    eye: Mat3.create(1, 0, 0, 0, 1, 0, 0, 0, 1),
    rotX90: Mat3.create(1, 0, 0, 0, 0, 1, 0, -1, 0),
    rotY90: Mat3.create(0, 0, -1, 0, 1, 0, 1, 0, 0),
    rotZ90: Mat3.create(0, 1, 0, -1, 0, 0, 0, 0, 1),
    rotX270: Mat3.create(1, 0, 0, 0, 0, -1, 0, 1, 0),
    rotY270: Mat3.create(0, 0, 1, 0, 1, 0, -1, 0, 0),
    rotZ270: Mat3.create(0, -1, 0, 1, 0, 0, 0, 0, 1),
    rotX180: Mat3.create(1, 0, 0, 0, -1, 0, 0, 0, -1),
    rotY180: Mat3.create(-1, 0, 0, 0, 1, 0, 0, 0, -1),
    rotZ180: Mat3.create(-1, 0, 0, 0, -1, 0, 0, 0, 1),
};


/** Return transformation which will align the PCA axes of an atomic structure
 * (or multiple structures) to the Cartesian axes x, y, z
 * (transformed = rotation * (coords - origin)).
 *
 * There are always 4 equally good rotations to do this (4 flips).
 * If `referenceRotation` is provided, select the one nearest to `referenceRotation`.
 * Otherwise use arbitrary rules to ensure the orientation after transform does not depend on the original orientation.
 */
export function structureLayingTransform(structures: Structure[], referenceRotation?: Mat3): { rotation: Mat3, origin: Vec3 } {
    const coords = selectCoords(structures, MIN_ATOMS_FOR_PCA);
    return layingTransform(coords, referenceRotation);
}

/** Return transformation which will align the PCA axes of a sequence
 * of points to the Cartesian axes x, y, z
 * (transformed = rotation * (coords - origin)).
 *
 * `coords` is a flattened array of 3D coordinates (i.e. the first 3 values are x, y, and z of the first point etc.).
 *
 * There are always 4 equally good rotations to do this (4 flips).
 * If `referenceRotation` is provided, select the one nearest to `referenceRotation`.
 * Otherwise use arbitrary rules to ensure the orientation after transform does not depend on the original orientation.
 */
export function layingTransform(coords: number[], referenceRotation?: Mat3): { rotation: Mat3, origin: Vec3 } {
    if (coords.length === 0) {
        console.warn('Skipping PCA, no atoms');
        return { rotation: Mat3.identity(), origin: Vec3.zero() };
    }
    const axes = PrincipalAxes.calculateMomentsAxes(coords);
    const normAxes = PrincipalAxes.calculateNormalizedAxes(axes);
    const R = mat3FromRows(normAxes.dirA, normAxes.dirB, normAxes.dirC);
    avoidMirrorRotation(R); // The SVD implementation seems to always provide proper rotation, but just to be sure
    const flip = referenceRotation ? minimalFlip(R, referenceRotation) : canonicalFlip(coords, R, axes.origin);
    Mat3.mul(R, flip, R);
    return { rotation: R, origin: normAxes.origin };
}

/** Try these selection strategies until having at least `minAtoms` atoms:
 * 1. only "polymer" atoms (e.g. C-alpha and O3')
 * 2. all non-hydrogen atoms with exception of water (HOH)
 * 3. all atoms
 * Return the coordinates in a flattened array (in triples).
 * If the total number of atoms is less than `minAtoms`, return only those. */
function selectCoords(structs: Structure[], minAtoms: number): number[] {
    let coords = selectCACoords(structs);
    if (coords.length >= 3 * minAtoms) return coords;

    coords = selectHeavyCoords(structs);
    if (coords.length >= 3 * minAtoms) return coords;

    coords = selectAllCoords(structs);
    return coords;
}

/** Select coordinates of C-alpha and O3' atoms */
function selectCACoords(structs: Structure[]): number[] {
    const coords: number[] = [];
    for (const struct of structs) {
        for (const unit of struct.units) {
            const { x, y, z } = unit.conformation;
            for (let i = 0; i < unit.polymerElements.length; i++) {
                const index = unit.polymerElements[i];
                coords.push(x(index), y(index), z(index));
            }
        }
    }
    return coords;
}

/** Select coordinates of non-hydrogen atoms, excluding water */
function selectHeavyCoords(structs: Structure[]): number[] {
    const coords: number[] = [];
    for (const struct of structs) {
        for (const unit of struct.units) {
            const { x, y, z } = unit.model.atomicConformation;
            for (let i = 0; i < unit.elements.length; i++) {
                const index = unit.elements[i];
                const compound = unit.model.atomicHierarchy.atoms.label_comp_id.value(index);
                const element = unit.model.atomicHierarchy.atoms.type_symbol.value(index);
                if (element !== 'H' && compound !== 'HOH') {
                    coords.push(x[index], y[index], z[index]);
                }
            }
        }
    }
    return coords;
}

/** Select coordinates of all atoms */
function selectAllCoords(structs: Structure[]): number[] {
    const coords: number[] = [];
    for (const struct of structs) {
        for (const unit of struct.units) {
            const { x, y, z } = unit.model.atomicConformation;
            for (let i = 0; i < unit.elements.length; i++) {
                const index = unit.elements[i];
                coords.push(x[index], y[index], z[index]);
            }
        }
    }
    return coords;
}

/** Return a flip around XYZ axes which minimizes the difference between flip*rotation and referenceRotation. */
function minimalFlip(rotation: Mat3, referenceRotation: Mat3): Mat3 {
    let bestFlip = ROTATION_MATRICES.eye;
    let bestScore = 0; // there will always be at least one positive score
    const aux = Mat3();
    for (const flip of [ROTATION_MATRICES.eye, ROTATION_MATRICES.rotX180, ROTATION_MATRICES.rotY180, ROTATION_MATRICES.rotZ180]) {
        const score = mat3Dot(Mat3.mul(aux, flip, rotation), referenceRotation);
        if (score > bestScore) {
            bestFlip = flip;
            bestScore = score;
        }
    }
    return bestFlip;
}

/** Measure of how similar two rotation matrices are (dot product of flattened matrices) */
function mat3Dot(a: Mat3, b: Mat3) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4] + a[5] * b[5] + a[6] * b[6] + a[7] * b[7] + a[8] * b[8];
}

/** Return a rotation matrix (flip) that should be applied to `coords` (after being rotated by `rotation`)
 * to ensure a deterministic "canonical" rotation.
 * There are 4 flips to choose from (one identity and three 180-degree rotations around the X, Y, and Z axes).
 * One of these 4 possible results is selected so that:
 *   1) starting and ending coordinates tend to be more in front (z > 0), middle more behind (z < 0).
 *   2) starting coordinates tend to be more left-top (x < y), ending more right-bottom (x > y).
 * These rules are arbitrary, but try to avoid ties for at least some basic symmetries.
 * Provided `origin` parameter MUST be the mean of the coordinates, otherwise it will not work!
 */
function canonicalFlip(coords: number[], rotation: Mat3, origin: Vec3): Mat3 {
    const pcaX = Vec3.create(Mat3.getValue(rotation, 0, 0), Mat3.getValue(rotation, 0, 1), Mat3.getValue(rotation, 0, 2));
    const pcaY = Vec3.create(Mat3.getValue(rotation, 1, 0), Mat3.getValue(rotation, 1, 1), Mat3.getValue(rotation, 1, 2));
    const pcaZ = Vec3.create(Mat3.getValue(rotation, 2, 0), Mat3.getValue(rotation, 2, 1), Mat3.getValue(rotation, 2, 2));
    const n = Math.floor(coords.length / 3);
    const v = Vec3();
    let xCum = 0;
    let yCum = 0;
    let zCum = 0;
    for (let i = 0; i < n; i++) {
        Vec3.fromArray(v, coords, 3 * i);
        Vec3.sub(v, v, origin);
        xCum += i * Vec3.dot(v, pcaX);
        yCum += i * Vec3.dot(v, pcaY);
        zCum += veeSlope(i, n) * Vec3.dot(v, pcaZ);
        // Thanks to subtracting `origin` from `coords` the slope functions `i` and `veeSlope(i, n)`
        // don't have to have zero sum (can be shifted up or down):
        //     sum{(slope[i]+shift)*(coords[i]-origin).PCA} =
        //     = sum{slope[i]*coords[i].PCA - slope[i]*origin.PCA + shift*coords[i].PCA - shift*origin.PCA} =
        //     = sum{slope[i]*(coords[i]-origin).PCA} + shift*sum{coords[i]-origin}.PCA =
        //     = sum{slope[i]*(coords[i]-origin).PCA}
    }
    const wrongFrontBack = zCum < 0;
    const wrongLeftTopRightBottom = wrongFrontBack ? xCum + yCum < 0 : xCum - yCum < 0;
    if (wrongLeftTopRightBottom && wrongFrontBack) {
        return ROTATION_MATRICES.rotY180; // flip around Y = around X then Z
    } else if (wrongFrontBack) {
        return ROTATION_MATRICES.rotX180; // flip around X
    } else if (wrongLeftTopRightBottom) {
        return ROTATION_MATRICES.rotZ180; // flip around Z
    } else {
        return Mat3.identity(); // do not flip
    }
}

/** Auxiliary function defined for i in [0, n), linearly decreasing from 0 to n/2
 * and then increasing back from n/2 to n, resembling letter V. */
function veeSlope(i: number, n: number) {
    const mid = Math.floor(n / 2);
    if (i < mid) {
        if (n % 2) return mid - i;
        else return mid - i - 1;
    } else {
        return i - mid;
    }
}

function mat3FromRows(row0: Vec3, row1: Vec3, row2: Vec3): Mat3 {
    const m = Mat3();
    Mat3.setValue(m, 0, 0, row0[0]);
    Mat3.setValue(m, 0, 1, row0[1]);
    Mat3.setValue(m, 0, 2, row0[2]);
    Mat3.setValue(m, 1, 0, row1[0]);
    Mat3.setValue(m, 1, 1, row1[1]);
    Mat3.setValue(m, 1, 2, row1[2]);
    Mat3.setValue(m, 2, 0, row2[0]);
    Mat3.setValue(m, 2, 1, row2[1]);
    Mat3.setValue(m, 2, 2, row2[2]);
    return m;
}

/** Check if a rotation matrix includes mirroring and invert Z axis in such case, to ensure a proper rotation (in-place). */
function avoidMirrorRotation(rot: Mat3) {
    if (Mat3.determinant(rot) < 0) {
        Mat3.setValue(rot, 2, 0, -Mat3.getValue(rot, 2, 0));
        Mat3.setValue(rot, 2, 1, -Mat3.getValue(rot, 2, 1));
        Mat3.setValue(rot, 2, 2, -Mat3.getValue(rot, 2, 2));
    }
}

/** Return a new camera snapshot with the same target and camera distance from the target as `old`
 * but with diferent orientation.
 * The actual rotation applied to the camera is the inverse of `rotation`,
 * which creates the same effect as if `rotation` were applied to the whole scene without moving the camera.
 * The rotation is relative to the default camera orientation (not to the current orientation). */
export function changeCameraRotation(old: Camera.Snapshot, rotation: Mat3): Camera.Snapshot {
    const cameraRotation = Mat3.invert(Mat3(), rotation);
    const dist = Vec3.distance(old.position, old.target);
    const relPosition = Vec3.transformMat3(Vec3(), Vec3.create(0, 0, dist), cameraRotation);
    const newUp = Vec3.transformMat3(Vec3(), Vec3.create(0, 1, 0), cameraRotation);
    const newPosition = Vec3.add(Vec3(), old.target, relPosition);
    return { ...old, position: newPosition, up: newUp };
}
