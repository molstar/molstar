/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from 'mol-math/geometry'
import { Mat4, Vec3 } from 'mol-math/linear-algebra'
import { BoundaryHelper } from 'mol-math/geometry/boundary-helper';

export function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n / width) : 0
    return { width, height, length: width * height * itemSize }
}

export interface TextureImage<T extends Uint8Array | Float32Array> {
    readonly array: T
    readonly width: number
    readonly height: number
}

export interface TextureVolume<T extends Uint8Array | Float32Array> {
    readonly array: T
    readonly width: number
    readonly height: number
    readonly depth: number
}

export function createTextureImage(n: number, itemSize: number): TextureImage<Uint8Array> {
    const { length, width, height } = calculateTextureInfo(n, itemSize)
    return { array: new Uint8Array(length), width, height }
}

//

const m = Mat4.zero()
const v = Vec3.zero()
const boundaryHelper = new BoundaryHelper()

export function calculateInvariantBoundingSphere(position: Float32Array, positionCount: number): Sphere3D {
    boundaryHelper.reset(0)
    for (let i = 0, _i = positionCount * 3; i < _i; i += 3) {
        Vec3.fromArray(v, position, i)
        boundaryHelper.boundaryStep(v, 0)
    }
    boundaryHelper.finishBoundaryStep()
    for (let i = 0, _i = positionCount * 3; i < _i; i += 3) {
        Vec3.fromArray(v, position, i)
        boundaryHelper.extendStep(v, 0)
    }
    return boundaryHelper.getSphere()
}

export function calculateTransformBoundingSphere(invariantBoundingSphere: Sphere3D, transform: Float32Array, transformCount: number): Sphere3D {
    const { center, radius } = invariantBoundingSphere
    boundaryHelper.reset(0)
    for (let i = 0, _i = transformCount; i < _i; ++i) {
        Vec3.transformMat4(v, center, Mat4.fromArray(m, transform, i * 16))
        boundaryHelper.boundaryStep(v, radius)
    }
    boundaryHelper.finishBoundaryStep()
    for (let i = 0, _i = transformCount; i < _i; ++i) {
        Vec3.transformMat4(v, center, Mat4.fromArray(m, transform, i * 16))
        boundaryHelper.extendStep(v, radius)
    }
    return boundaryHelper.getSphere()
}

export function calculateBoundingSphere(position: Float32Array, positionCount: number, transform: Float32Array, transformCount: number): { boundingSphere: Sphere3D, invariantBoundingSphere: Sphere3D } {
    const invariantBoundingSphere = calculateInvariantBoundingSphere(position, positionCount)
    const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform, transformCount)
    return { boundingSphere, invariantBoundingSphere }
}