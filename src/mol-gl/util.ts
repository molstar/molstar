/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n * itemSize)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n * itemSize / width) : 0
    return { width, height, length: width * height * itemSize }
}

export interface ColorTexture extends Uint8Array {
    width: number,
    height: number
}

export function createColorTexture (n: number): ColorTexture {
    const colorTexInfo = calculateTextureInfo(n, 3)
    const colorTexture = new Uint8Array(colorTexInfo.length)
    return Object.assign(colorTexture, {
        width: colorTexInfo.width,
        height: colorTexInfo.height
    })
}