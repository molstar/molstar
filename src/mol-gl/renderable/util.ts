/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from 'mol-math/geometry'
import { Mat4, Vec3 } from 'mol-math/linear-algebra'

export function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n * itemSize)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n * itemSize / width) : 0
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

export function fillSerial<T extends Helpers.NumberArray> (array: T) {
    const n = array.length
    for (let i = 0; i < n; ++i) array[ i ] = i
    return array
}

export interface PositionData {
    position: Float32Array
    positionCount: number,
    transform: Float32Array,
    transformCount: number,
}

export function calculateBoundingSphere(data: PositionData): Sphere3D {
    const { position, positionCount, transform, transformCount } = data

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