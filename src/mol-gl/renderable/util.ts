/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from 'mol-math/geometry'
import { Mat4, Vec3 } from 'mol-math/linear-algebra'
import { ValueCell } from 'mol-util';

export function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n / width) : 0
    return { width, height, length: width * height * itemSize }
}

export interface TextureImage {
    array: Uint8Array
    width: number
    height: number
}

export function createTextureImage (n: number, itemSize: number): TextureImage {
    const { length, width, height } = calculateTextureInfo(n, itemSize)
    return { array: new Uint8Array(length), width, height }
}

export interface PositionValues {
    aPosition: ValueCell<Float32Array>
    drawCount: ValueCell<number>,
    aTransform: ValueCell<Float32Array>,
    instanceCount: ValueCell<number>,
}

function getPositionDataFromValues(values: PositionValues) {
    return {
        position: values.aPosition.ref.value,
        positionCount: values.drawCount.ref.value / 3 / 3,
        transform: values.aTransform.ref.value,
        transformCount: values.instanceCount.ref.value
    }
}

export function calculateBoundingSphereFromValues(values: PositionValues) {
    const { position, positionCount, transform, transformCount } = getPositionDataFromValues(values)
    return calculateBoundingSphere(position, positionCount, transform, transformCount)
}

export function calculateBoundingSphere(position: Float32Array, positionCount: number, transform: Float32Array, transformCount: number): Sphere3D {

    const m = Mat4.zero()

    let cx = 0, cy = 0, cz = 0;
    let radiusSq = 0;

    for (let i = 0, _i = positionCount * 3; i < _i; i += 3) {
        cx += position[i];
        cy += position[i + 1];
        cz += position[i + 2];
    }

    if (positionCount > 0) {
        cx /= positionCount;
        cy /= positionCount;
        cz /= positionCount;
    }

    for (let i = 0, _i = positionCount * 3; i < _i; i += 3) {
        const dx = position[i] - cx
        const dy = position[i + 1] - cy
        const dz = position[i + 2] - cz;
        const d = dx * dx + dy * dy + dz * dz;
        if (d > radiusSq) radiusSq = d;
    }

    const c = Vec3.create(cx, cy, cz)
    const ct = Vec3.zero()

    const center = Vec3.zero()
    const centers = new Float32Array(3 * transformCount)

    for (let i = 0, _i = transformCount; i < _i; ++i) {
        Mat4.fromArray(m, transform, i * 16)
        Vec3.transformMat4(ct, c, m)
        Vec3.add(center, center, ct)
        Vec3.toArray(ct, centers, i * 3)
    }

    Vec3.scale(center, center, 1 / transformCount)

    let r = Math.sqrt(radiusSq)
    let radius = r

    for (let i = 0, _i = transformCount; i < _i; ++i) {
        Vec3.fromArray(ct, centers, i * 3)
        radius = Math.max(radius, Vec3.distance(center, ct) + r)
    }

    return { center, radius };
}