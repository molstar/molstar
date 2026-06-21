/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { ParticleList } from '../../mol-model/particles/particle-list';
import { CustomProperties } from '../../mol-model/custom-property';
import { ArtiatomiEmFile } from '../../mol-io/reader/artiatomi/em';
import { degToRad } from '../../mol-math/misc';
import { ModelFormat } from '../format';

// Artiatomi motivelist column layout (20 rows, 0-indexed):
// https://github.com/uermel/Artiatomi/wiki/Motivelist
//
// NamingConvention TomoParticle (preferred):
//   Row 0:  cc        – cross-correlation coefficient
//   Row 1:  x_legacy  – legacy x-coordinate (unused)
//   Row 2:  y_legacy  – legacy y-coordinate (unused)
//   Row 3:  custom    – tomogram-specific particle number (customisable)
//   Row 4:  tomo      – tomogram number
//   Row 5:  particle  – global particle number
//   Row 6:  wedge     – running number of wedge volume
//   Row 7:  x         – x-coordinate in reconstruction volume
//   Row 8:  y         – y-coordinate in reconstruction volume
//   Row 9:  z         – z-coordinate in reconstruction volume
//   Row 10: dx        – x-shift after rotation of template
//   Row 11: dy        – y-shift after rotation of template
//   Row 12: dz        – z-shift after rotation of template
//   Row 13: dx_pre    – x-shift before rotation (legacy, unused)
//   Row 14: dy_pre    – y-shift before rotation (legacy, unused)
//   Row 15: dz_pre    – z-shift before rotation (legacy, unused)
//   Row 16: phi       – Euler angle Phi [deg]
//   Row 17: psi       – Euler angle Psi [deg]
//   Row 18: theta     – Euler angle Theta [deg]
//   Row 19: class     – class number

const ROW_CC = 0;
const ROW_TOMO = 4;
const ROW_X = 7;
const ROW_Y = 8;
const ROW_Z = 9;
const ROW_DX = 10;
const ROW_DY = 11;
const ROW_DZ = 12;
const ROW_PHI = 16;
const ROW_PSI = 17;
const ROW_THETA = 18;
const ROW_CLASS = 19;

export const ArtiatomiMotivelistRowCount = 20;

const tmpRotZphi = Mat4();
const tmpRotXtheta = Mat4();
const tmpRotZpsi = Mat4();

/**
 * Artiatomi/TOM ZXZ Euler rotation: R = Rz(psi) * Rx(theta) * Rz(phi).
 *
 * Derived from Artiatomi's CudaRot::computeRotMat (src/SubTomogramAverageMPI/CudaRot.cpp):
 *   rotMat[i][j] == M[j][i]  (column-major storage)
 *
 * The resulting 3×3 matrix is:
 *   [ cospsi*cosphi - costheta*sinpsi*sinphi ,  sinpsi*cosphi + costheta*cospsi*sinphi ,  sintheta*sinphi ]
 *   [ -cospsi*sinphi - costheta*sinpsi*cosphi, -sinpsi*sinphi + costheta*cospsi*cosphi,  sintheta*cosphi ]
 *   [ sintheta*sinpsi                        , -sintheta*cospsi                        ,  costheta        ]
 */
function artiatomiEulerToRotation(out: Mat4, phi: number, psi: number, theta: number): Mat4 {
    const rotZphi = Mat4.fromRotation(tmpRotZphi, phi, Vec3.unitZ);
    const rotXtheta = Mat4.fromRotation(tmpRotXtheta, theta, Vec3.unitX);
    const rotZpsi = Mat4.fromRotation(tmpRotZpsi, psi, Vec3.unitZ);
    // R = Rz(psi) * Rx(theta) * Rz(phi)
    Mat4.mul(out, rotXtheta, rotZphi);
    Mat4.mul(out, rotZpsi, out);
    return out;
}

export interface ArtiatomiMotivelistOptions {
    /**
     * Pixel size (Å/pixel) used to convert voxel-space coordinates to angstrom.
     * Artiatomi EM motivelists do not encode distance units, so this must be supplied.
     */
    readonly pixelSize: number
    /** If given, only particles whose tomogram number (row 5 in TomoParticle convention) matches are included. */
    readonly tomos?: ReadonlyArray<number>
    readonly label?: string
}

/** Return the sorted set of unique tomogram IDs from the motivelist (row 5, TomoParticle convention). */
export function getArtiatomiMotivelistTomogramIds(data: ArtiatomiEmFile): number[] {
    if (data.header.dimX !== ArtiatomiMotivelistRowCount) return [];
    const n = data.header.dimY * data.header.dimZ;
    const tomos = new Set<number>();
    for (let p = 0; p < n; p++) {
        const tomo = data.data[p * ArtiatomiMotivelistRowCount + ROW_TOMO];
        tomos.add(tomo);
    }
    return Array.from(tomos).sort((a, b) => a - b);
}

export function createParticleListFromArtiatomiEm(data: ArtiatomiEmFile, options: ArtiatomiMotivelistOptions): ParticleList {
    if (data.header.dimX !== ArtiatomiMotivelistRowCount) {
        throw new Error(
            `EM file does not appear to be a motivelist: expected DimX = ${ArtiatomiMotivelistRowCount}, got ${data.header.dimX}.`
        );
    }

    const pixelSize = options.pixelSize > 0 ? options.pixelSize : 1;
    const particleCount = data.header.dimY * data.header.dimZ;
    const tomoFilter = options.tomos?.length ? new Set<number>(options.tomos) : undefined;
    const vals = data.data;

    // Pre-count accepted particles so we can allocate exact-sized typed arrays.
    let acceptedCount = 0;
    for (let p = 0; p < particleCount; p++) {
        if (tomoFilter !== undefined) {
            const tomo = vals[p * ArtiatomiMotivelistRowCount + ROW_TOMO];
            if (!tomoFilter.has(tomo)) continue;
        }
        acceptedCount++;
    }

    if (acceptedCount === 0) {
        throw new Error(tomoFilter !== undefined
            ? `No Artiatomi motivelist particles matched tomos '${options.tomos!.join(', ')}'.`
            : 'No readable Artiatomi motivelist particles were found.');
    }

    const keys = new Int32Array(acceptedCount);
    const coordinates = new Float32Array(acceptedCount * 3);
    const rotations = new Float32Array(acceptedCount * 4);

    const rotation = Mat4();
    const quaternion = Quat();

    let count = 0;
    for (let p = 0; p < particleCount; p++) {
        const base = p * ArtiatomiMotivelistRowCount;

        if (tomoFilter !== undefined) {
            if (!tomoFilter.has(vals[base + ROW_TOMO])) continue;
        }

        // Coordinates (voxels) with post-rotation shifts applied in the lab frame.
        // Shifts (rows 11–13) are the offsets by which the template was displaced after
        // rotation to align with the particle, so the particle centre = coord − shift.
        const x = vals[base + ROW_X];
        const y = vals[base + ROW_Y];
        const z = vals[base + ROW_Z];
        const dx = vals[base + ROW_DX];
        const dy = vals[base + ROW_DY];
        const dz = vals[base + ROW_DZ];

        const cOffset = count * 3;
        coordinates[cOffset + 0] = (x - dx) * pixelSize;
        coordinates[cOffset + 1] = (y - dy) * pixelSize;
        coordinates[cOffset + 2] = (z - dz) * pixelSize;

        // Orientation: R = Rz(psi) * Rx(theta) * Rz(phi)
        const phi = degToRad(vals[base + ROW_PHI]);
        const psi = degToRad(vals[base + ROW_PSI]);
        const theta = degToRad(vals[base + ROW_THETA]);

        artiatomiEulerToRotation(rotation, phi, psi, theta);
        Quat.normalize(quaternion, Quat.fromMat4(quaternion, rotation));

        const qOffset = count * 4;
        rotations[qOffset + 0] = quaternion[0];
        rotations[qOffset + 1] = quaternion[1];
        rotations[qOffset + 2] = quaternion[2];
        rotations[qOffset + 3] = quaternion[3];

        keys[count] = p;
        count++;
    }

    const finalKeys = count === particleCount ? keys : keys.slice(0, count);
    const finalCoords = count === particleCount ? coordinates : coordinates.slice(0, count * 3);
    const finalRotations = count === particleCount ? rotations : rotations.slice(0, count * 4);

    return {
        label: buildArtiatomiLabel(options.label, options.tomos),
        count,
        keys: finalKeys,
        targets: new Int32Array(count),
        coordinates: finalCoords,
        rotations: finalRotations,
        getParticleLabel: (index: number) => {
            const p = finalKeys[index];
            const base = p * ArtiatomiMotivelistRowCount;
            const parts: string[] = [`#${p + 1}`];
            const tomo = vals[base + ROW_TOMO];
            if (tomo > 0) parts.push(`tomo ${tomo}`);
            const classNo = vals[base + ROW_CLASS];
            if (classNo > 0) parts.push(`class ${classNo}`);
            const cc = vals[base + ROW_CC];
            parts.push(`cc ${cc.toFixed(4)}`);
            return parts.join(' | ');
        },
        sourceData: ArtiatomiEmFormat.create(data),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

function buildArtiatomiLabel(label?: string, tomos?: ReadonlyArray<number>): string {
    if (!label) label = 'Particles';
    if (tomos?.length) {
        return tomos.length === 1
            ? `${label} (tomo ${tomos[0]})`
            : `${label} (tomos ${tomos.join(', ')})`;
    }
    return label;
}

//

export { ArtiatomiEmFormat };

type ArtiatomiEmFormat = ModelFormat<ArtiatomiEmFile>

namespace ArtiatomiEmFormat {
    export function is(x?: ModelFormat): x is ArtiatomiEmFormat {
        return x?.kind === 'artiatomi-em';
    }

    export function create(em: ArtiatomiEmFile): ArtiatomiEmFormat {
        return { kind: 'artiatomi-em', name: 'artiatomi-em', data: em };
    }
}
