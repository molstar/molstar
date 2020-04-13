/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export { PixelData };

interface PixelData {
    readonly array: Uint8Array
    readonly width: number
    readonly height: number
}

namespace PixelData {
    export function create(array: Uint8Array, width: number, height: number): PixelData {
        return { array, width, height };
    }

    /** horizontally flips the pixel data in-place */
    export function flipY(pixelData: PixelData): PixelData {
        const { array, width, height } = pixelData;
        const width4 = width * 4;
        for (let i = 0, maxI = height / 2; i < maxI; ++i) {
            for (let j = 0, maxJ = width4; j < maxJ; ++j) {
                const index1 = i * width4 + j;
                const index2 = (height - i - 1) * width4 + j;
                const tmp = array[index1];
                array[index1] = array[index2];
                array[index2] = tmp;
            }
        }
        return pixelData;
    }
}