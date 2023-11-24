/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { sortIfNeeded } from '../mol-util/array';


function randomFloats(n: number) {
    const SCALE = 1000;
    const data = new Array(n);
    for (let i = 0; i < n; i++) {
        data[i] = SCALE * Math.random();
    }
    return data;
}

type SortFunction = (data: number[], compareFn: (a: number, b: number) => number) => any

function benchmarkSortFunction(sortFunction: SortFunction, arrayLength: number, alreadySorted: boolean) {
    const data = randomFloats(arrayLength);
    if (alreadySorted) data.sort((a, b) => a - b);
    const label = `${sortFunction.name} (${arrayLength}, ${alreadySorted ? 'pre-sorted' : 'not sorted'})`;
    console.time(label);
    sortFunction(data, (a, b) => a - b);
    console.timeEnd(label);
}

export function runBenchmarks(arrayLength: number, nRepeats: number) {
    const sort: SortFunction = (data, cmp) => data.sort(cmp);
    for (let i = 0; i < nRepeats; i++) benchmarkSortFunction(sort, arrayLength, true);
    for (let i = 0; i < nRepeats; i++) benchmarkSortFunction(sortIfNeeded, arrayLength, true);
    for (let i = 0; i < nRepeats; i++) benchmarkSortFunction(sort, arrayLength, false);
    for (let i = 0; i < nRepeats; i++) benchmarkSortFunction(sortIfNeeded, arrayLength, false);
    console.log('`sortIfNeeded` should be faster than `sort` on pre-sorted data, same speed on non-sorted data');
}

runBenchmarks(10 ** 6, 3);
