/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { arrayMinMax } from '../mol-util/array';

export interface Histogram {
    min: number,
    max: number,
    binWidth: number,
    counts: Int32Array
}

export function calculateHistogram(data: ArrayLike<number>, binCount: number, options?: { min: number, max: number, }): Histogram {
    if (!options) {
        const [min, max] = arrayMinMax(data);
        return _calcHistogram(data, binCount, min, max);
    } else {
        return _calcHistogram(data, binCount, options.min, options.max);
    }
}

function _calcHistogram(data: ArrayLike<number>, binCount: number, min: number, max: number): Histogram {
    let binWidth = (max - min) / binCount;
    if (binWidth === 0) binWidth = 1;

    const counts = new Int32Array(binCount);

    for (let i = 0, _i = data.length; i < _i; i++) {
        let bin = Math.floor((data[i] - min) / binWidth);
        if (bin >= binCount) bin = binCount - 1;
        else if (bin < 0) bin = 0;
        counts[bin]++;
    }

    return { min, max, binWidth, counts };
}