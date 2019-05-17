/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from 'mol-math/geometry'
import { Vec3 } from 'mol-math/linear-algebra'
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

export function createTextureImage<T extends Uint8Array | Float32Array>(n: number, itemSize: number, arrayCtor: new (length: number) => T, array?: T): TextureImage<T> {
    const { length, width, height } = calculateTextureInfo(n, itemSize)
    array = array && array.length >= length ? array : new arrayCtor(length)
    return { array, width, height }
}

export function printTextureImage(textureImage: TextureImage<any>, scale = 1) {
    const { array, width, height } = textureImage
    const itemSize = array.length / (width * height)
    const data = new Uint8ClampedArray(width * height * 4)
    if (itemSize === 1) {
        for (let y = 0; y < height; ++y) {
            for (let x = 0; x < width; ++x) {
                data[(y * width + x) * 4 + 3] = array[y * width + x]
            }
        }
    } else if (itemSize === 4) {
        data.set(array)
    } else {
        console.warn(`itemSize '${itemSize}' not supported`)
    }
    return printImageData(new ImageData(data, width, height), scale)
}

export function printImageData(imageData: ImageData, scale = 1, pixelated = false) {
    const canvas = document.createElement('canvas')
    canvas.width = imageData.width
    canvas.height = imageData.height
    const ctx = canvas.getContext('2d')
    if (!ctx) throw new Error('Could not create canvas 2d context')
    ctx.putImageData(imageData, 0, 0)
    canvas.toBlob(imgBlob => {
        const objectURL = window.URL.createObjectURL(imgBlob)
        const img = document.createElement('img')
        img.src = objectURL
        img.style.width = imageData.width * scale + 'px'
        img.style.height = imageData.height * scale + 'px';
        if (pixelated) (img.style as any).imageRendering = 'pixelated' // supported only in Chrome
        img.style.position = 'absolute'
        img.style.top = '0px'
        img.style.left = '0px'
        img.style.border = 'solid grey'
        document.body.appendChild(img)
    }, 'image/png')
}

//

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
        Vec3.transformMat4Offset(v, center, transform, 0, 0, i * 16)
        boundaryHelper.boundaryStep(v, radius)
    }
    boundaryHelper.finishBoundaryStep()
    for (let i = 0, _i = transformCount; i < _i; ++i) {
        Vec3.transformMat4Offset(v, center, transform, 0, 0, i * 16)
        boundaryHelper.extendStep(v, radius)
    }
    return boundaryHelper.getSphere()
}

export function calculateBoundingSphere(position: Float32Array, positionCount: number, transform: Float32Array, transformCount: number, padding = 0): { boundingSphere: Sphere3D, invariantBoundingSphere: Sphere3D } {
    const invariantBoundingSphere = calculateInvariantBoundingSphere(position, positionCount)
    const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform, transformCount)
    Sphere3D.expand(boundingSphere, boundingSphere, padding)
    Sphere3D.expand(invariantBoundingSphere, invariantBoundingSphere, padding)
    return { boundingSphere, invariantBoundingSphere }
}