/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export { PixelData };

interface PixelData {
    readonly array: Uint8Array | Float32Array
    readonly width: number
    readonly height: number
}

namespace PixelData {
    export function create(array: Uint8Array | Float32Array, width: number, height: number): PixelData {
        return { array, width, height };
    }

    /** horizontally flips the pixel data in-place */
    export function flipY(pixelData: PixelData): PixelData {
        const { array, width, height } = pixelData;
        const itemSize = array.length / (width * height);
        const widthIS = width * itemSize;
        for (let i = 0, maxI = height / 2; i < maxI; ++i) {
            for (let j = 0, maxJ = widthIS; j < maxJ; ++j) {
                const index1 = i * widthIS + j;
                const index2 = (height - i - 1) * widthIS + j;
                const tmp = array[index1];
                array[index1] = array[index2];
                array[index2] = tmp;
            }
        }
        return pixelData;
    }

    /** to undo pre-multiplied alpha */
    export function divideByAlpha(pixelData: PixelData): PixelData {
        const { array } = pixelData;
        const factor = (array instanceof Uint8Array) ? 255 : 1;
        for (let i = 0, il = array.length; i < il; i += 4) {
            const a = array[i + 3] / factor;
            array[i] /= a;
            array[i + 1] /= a;
            array[i + 2] /= a;
        }
        return pixelData;
    }
}