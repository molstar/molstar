/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { ParticleList } from '../../mol-model/particles/particle-list';
import { CustomProperties } from '../../mol-model/custom-property';
import { DynamoTblFile } from '../../mol-io/reader/dynamo/tbl';
import { Column } from '../../mol-data/db';
import { degToRad } from '../../mol-math/misc';
import { ModelFormat } from '../format';

function dynamoEulerToRotation(out: Mat4, tdrot: number, tilt: number, narot: number) {
    const tdrotZ = Mat4.fromRotation(Mat4(), degToRad(-tdrot), Vec3.unitZ);
    const tiltX = Mat4.fromRotation(Mat4(), degToRad(tilt), Vec3.unitX);
    const narotZ = Mat4.fromRotation(Mat4(), degToRad(-narot), Vec3.unitZ);

    Mat4.mul(out, tiltX, tdrotZ);
    Mat4.mul(out, narotZ, out);
    return out;
}

export interface DynamoParticleListOptions {
    readonly label?: string
    readonly tomos?: ReadonlyArray<number>
    /** Override pixel size (Å/pixel) used to convert pixel-space coordinates to angstrom. */
    readonly pixelSize?: number
    /** Uniform particle radius in angstrom assigned to every particle. Leave 0 or undefined to omit radii. */
    readonly particleRadius?: number
}

function buildDynamoLabel(tomos?: ReadonlyArray<number>) {
    if (tomos !== void 0 && tomos.length > 0) {
        return tomos.length === 1
            ? `Dynamo particles (tomo ${tomos[0]})`
            : `Dynamo particles (tomos ${tomos.join(', ')})`;
    }
    return 'Dynamo particles';
}

function buildTblAttrMaps(count: number, rowCount: number, defs: Array<[string, string, Float32Array | undefined]>) {
    const attributes = new Map<string, Float32Array>();
    const attributeInfo = new Map<string, { label: string, min: number, max: number }>();
    for (const [key, label, buf] of defs) {
        if (!buf) continue;
        const arr = count === rowCount ? buf : buf.slice(0, count);
        let min = Infinity, max = -Infinity;
        for (let i = 0; i < count; i++) {
            const v = arr[i];
            if (isFinite(v)) {
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        if (!isFinite(min)) continue;
        attributes.set(key, arr);
        attributeInfo.set(key, { label, min, max });
    }
    return attributes.size > 0
        ? { attributes, attributeInfo }
        : { attributes: undefined, attributeInfo: undefined };
}

export function getDynamoTblTomogramIds(data: DynamoTblFile) {
    const tomograms = new Set<number>();
    const tomo = data.fields.tomo;
    for (let i = 0, il = tomo.rowCount; i < il; ++i) {
        if (tomo.valueKind(i) !== Column.ValueKinds.Present) continue;
        const v = tomo.value(i);
        if (Number.isFinite(v)) tomograms.add(v);
    }
    return Array.from(tomograms).sort((a, b) => a - b);
}

export function createParticleListFromDynamoTbl(data: DynamoTblFile, options: DynamoParticleListOptions = {}): ParticleList {
    const { x, y, z, dx, dy, dz, tdrot, tilt, narot, tomo, reg, apix } = data.fields;
    const classCol = data.fields.class;
    const annotationCol = data.fields.annotation;
    const rowCount = data.rowCount;

    const tomoFilter = options.tomos?.length
        ? new Set<number>(options.tomos)
        : void 0;

    const _keys = new Int32Array(rowCount);
    const _coordinates = new Float32Array(rowCount * 3);
    const _rotations = new Float32Array(rowCount * 4);

    const ccCol = data.fields.cc;
    const cc2Col = data.fields.cc2;
    const _ccAttr = ccCol.isDefined ? new Float32Array(rowCount) : undefined;
    const _cc2Attr = cc2Col.isDefined ? new Float32Array(rowCount) : undefined;
    const _classAttr = classCol.isDefined ? new Float32Array(rowCount) : undefined;

    const rotation = Mat4();
    const quaternion = Quat();

    let count = 0;
    for (let row = 0; row < rowCount; ++row) {
        if (tomoFilter !== void 0 && !tomoFilter.has(tomo.value(row))) continue;

        let ps = 1;
        if (options.pixelSize) {
            ps = options.pixelSize;
        } else if (apix.valueKind(row) === Column.ValueKinds.Present) {
            ps = apix.value(row);
        }

        const cOffset = count * 3;
        _coordinates[cOffset + 0] = (x.value(row) + dx.value(row)) * ps;
        _coordinates[cOffset + 1] = (y.value(row) + dy.value(row)) * ps;
        _coordinates[cOffset + 2] = (z.value(row) + dz.value(row)) * ps;

        dynamoEulerToRotation(rotation, tdrot.value(row), tilt.value(row), narot.value(row));
        Quat.normalize(quaternion, Quat.fromMat4(quaternion, rotation));
        const qOffset = count * 4;
        _rotations[qOffset + 0] = quaternion[0];
        _rotations[qOffset + 1] = quaternion[1];
        _rotations[qOffset + 2] = quaternion[2];
        _rotations[qOffset + 3] = quaternion[3];

        if (_ccAttr && ccCol.valueKind(row) === Column.ValueKinds.Present) _ccAttr[count] = ccCol.value(row);
        if (_cc2Attr && cc2Col.valueKind(row) === Column.ValueKinds.Present) _cc2Attr[count] = cc2Col.value(row);
        if (_classAttr && classCol.valueKind(row) === Column.ValueKinds.Present) _classAttr[count] = classCol.value(row);

        _keys[count] = row;
        ++count;
    }

    if (count === 0) {
        throw new Error(tomoFilter !== void 0
            ? `No Dynamo particle rows matched tomos '${options.tomos!.join(', ')}'.`
            : 'No readable Dynamo table rows were found.');
    }

    const keys = count === rowCount ? _keys : _keys.slice(0, count);
    const coordinates = count === rowCount ? _coordinates : _coordinates.slice(0, count * 3);
    const rotations = count === rowCount ? _rotations : _rotations.slice(0, count * 4);

    const { attributes, attributeInfo } = buildTblAttrMaps(count, rowCount, [
        ['cc', 'CC', _ccAttr],
        ['cc2', 'CC2', _cc2Attr],
        ['class', 'Class', _classAttr],
    ]);

    const radii = options.particleRadius && options.particleRadius > 0
        ? new Float32Array(count).fill(options.particleRadius)
        : undefined;

    return {
        label: options.label ?? buildDynamoLabel(options.tomos),
        count,
        keys,
        targets: new Int32Array(count),
        coordinates,
        rotations,
        radii,
        attributes,
        attributeInfo,
        getParticleLabel: (index: number) => {
            const row = keys[index];
            const parts: string[] = [`#${row + 1}`];
            if (tomo.valueKind(row) === Column.ValueKinds.Present) parts.push(`tomo ${tomo.value(row)}`);
            if (reg.valueKind(row) === Column.ValueKinds.Present) parts.push(`reg ${reg.value(row)}`);
            if (classCol.valueKind(row) === Column.ValueKinds.Present) parts.push(`class ${classCol.value(row)}`);
            if (annotationCol.valueKind(row) === Column.ValueKinds.Present) parts.push(`ann ${annotationCol.value(row)}`);
            return parts.join(' | ');
        },
        sourceData: DynamoTblFormat.create(data),
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

//

export { DynamoTblFormat };

type DynamoTblFormat = ModelFormat<DynamoTblFile>

namespace DynamoTblFormat {
    export function is(x?: ModelFormat): x is DynamoTblFormat {
        return x?.kind === 'dynamo-tbl';
    }

    export function create(dynamoTbl: DynamoTblFile): DynamoTblFormat {
        return { kind: 'dynamo-tbl', name: 'dynamo-tbl', data: dynamoTbl };
    }
}
