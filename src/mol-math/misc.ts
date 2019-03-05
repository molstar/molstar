import { normalize } from "./interpolate";

/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function degToRad (deg: number) {
    return deg * 0.01745  // deg * Math.PI / 180
}

export function radToDeg (rad: number) {
    return rad * 57.29578  // rad * 180 / Math.PI
}

export interface Histogram {
    bins: ArrayLike<number>
    counts: ArrayLike<number>
}

export function calcHistogram(data: ArrayLike<number>, min: number, max: number, length = 100): Histogram {
    const bins = new Float32Array(length + 1)
    const counts = new Uint32Array(100)
    const widthOfBin = (Math.abs(max) + Math.abs(min))/100;
    let start = min;

    // Putting the values into the bins they belong to
    for(let i = 0; i < data.length; i++) {
        const num = Math.floor(normalize(data[i], min, max)*100);
        counts[num] = counts[num] + 1;
    }
    // Setting the bins
    for(let i = 0; i < bins.length; i++) {
        bins[i] = start;
        start += widthOfBin;
    }

    return { bins, counts }
}