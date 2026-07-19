/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat } from '../../mol-math/linear-algebra';
import { ParticleList } from '../../mol-model/particles/particle-list';
import { CustomProperties } from '../../mol-model/custom-property';
import { CryoEtDataPortalNdjsonFile } from '../../mol-io/reader/cryoet/ndjson';
import { ModelFormat } from '../format';

function setRotationFromRowMajor3x3(out: Mat4, values: ReadonlyArray<ReadonlyArray<number>>) {
    Mat4.setIdentity(out);
    Mat4.setValue(out, 0, 0, values[0][0]);
    Mat4.setValue(out, 0, 1, values[0][1]);
    Mat4.setValue(out, 0, 2, values[0][2]);
    Mat4.setValue(out, 1, 0, values[1][0]);
    Mat4.setValue(out, 1, 1, values[1][1]);
    Mat4.setValue(out, 1, 2, values[1][2]);
    Mat4.setValue(out, 2, 0, values[2][0]);
    Mat4.setValue(out, 2, 1, values[2][1]);
    Mat4.setValue(out, 2, 2, values[2][2]);
    return out;
}

export interface CryoEtDataPortalParticleListOptions {
    /**
     * Pixel size (Å/pixel) used to convert pixel-space NDJSON coordinates to angstrom.
     * CryoET Data Portal NDJSON does not encode distance units, so this must be supplied.
     */
    readonly pixelSize: number
    readonly type?: string
    /** Uniform particle radius in angstrom assigned to every particle. Leave 0 or undefined to omit radii. */
    readonly particleRadius?: number
}

function buildCryoEtLabel(type?: string) {
    if (type) return `CryoET Data Portal particles (${type})`;
    return 'CryoET Data Portal particles';
}

export function createParticleListFromCryoEtDataPortalNdjson(data: CryoEtDataPortalNdjsonFile, options: CryoEtDataPortalParticleListOptions): ParticleList {
    const pixelSize = options.pixelSize || 1;

    const recordCount = data.records.length;
    const typeFilter = options.type;

    let acceptedCount = 0;
    let hasAnyRotation = false;
    for (let i = 0; i < recordCount; ++i) {
        const record = data.records[i];
        if (typeFilter && record.type !== typeFilter) continue;
        ++acceptedCount;
        if (record.type === 'orientedPoint') hasAnyRotation = true;
    }

    if (acceptedCount === 0) {
        throw new Error(typeFilter
            ? `No CryoET Data Portal ndjson records matched type '${typeFilter}'.`
            : 'No readable CryoET Data Portal ndjson particle records were found.');
    }

    const keys = new Int32Array(acceptedCount);
    const coordinates = new Float32Array(acceptedCount * 3);
    const rotations = hasAnyRotation ? new Float32Array(acceptedCount * 4) : undefined;

    const rotation = Mat4();
    const quaternion = Quat();

    let count = 0;
    for (let i = 0; i < recordCount; ++i) {
        const record = data.records[i];
        if (typeFilter && record.type !== typeFilter) continue;

        const cOffset = count * 3;
        coordinates[cOffset + 0] = record.location.x * pixelSize;
        coordinates[cOffset + 1] = record.location.y * pixelSize;
        coordinates[cOffset + 2] = record.location.z * pixelSize;

        if (rotations) {
            if (record.type === 'orientedPoint') {
                setRotationFromRowMajor3x3(rotation, record.xyz_rotation_matrix);
                Quat.normalize(quaternion, Quat.fromMat4(quaternion, rotation));
            } else {
                Quat.setIdentity(quaternion);
            }
            const qOffset = count * 4;
            rotations[qOffset + 0] = quaternion[0];
            rotations[qOffset + 1] = quaternion[1];
            rotations[qOffset + 2] = quaternion[2];
            rotations[qOffset + 3] = quaternion[3];
        }

        keys[count] = i;
        ++count;
    }

    const radii = options.particleRadius && options.particleRadius > 0
        ? new Float32Array(count).fill(options.particleRadius)
        : undefined;

    return {
        label: buildCryoEtLabel(options.type),
        count,
        keys,
        targets: new Int32Array(count),
        coordinates,
        rotations,
        radii,
        getParticleLabel: (index: number) => {
            const i = keys[index];
            const record = data.records[i];
            return `#${i + 1} | ${record.type}`;
        },
        sourceData: CryoEtDataPortalNdjsonFormat.create(data),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

//

export { CryoEtDataPortalNdjsonFormat };

type CryoEtDataPortalNdjsonFormat = ModelFormat<CryoEtDataPortalNdjsonFile>

namespace CryoEtDataPortalNdjsonFormat {
    export function is(x?: ModelFormat): x is CryoEtDataPortalNdjsonFormat {
        return x?.kind === 'cryoet-data-portal-ndjson';
    }

    export function create(ndjson: CryoEtDataPortalNdjsonFile): CryoEtDataPortalNdjsonFormat {
        return { kind: 'cryoet-data-portal-ndjson', name: 'cryoet-data-portal-ndjson', data: ndjson };
    }
}
