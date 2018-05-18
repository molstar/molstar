/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type Comparer<T = any> = (data: T, i: number, j: number) => number
export type Swapper<T = any> = (data: T, i: number, j: number) => void

type Ctx = { cmp: Comparer, swap: Swapper, parts: number[], data: any }

export function arrayLess(arr: ArrayLike<number>, i: number, j: number) {
    return arr[i] - arr[j];
}

export function arraySwap(arr: ArrayLike<any>, i: number, j: number) {
    const temp = arr[i];
    (arr as any[])[i] = arr[j];
    (arr as any[])[j] = temp;
}

function medianPivotIndex(data: any, cmp: Comparer, l: number, r: number) {
    const m = (l + r) >> 1;
    if (cmp(data, l, r) > 0) return cmp(data, l, m) > 0 ? cmp(data, m, r) > 0 ? m : r : l;
    else return cmp(data, r, m) > 0 ? cmp(data, m, l) > 0 ? m : l : r;
}

function partition(ctx: Ctx, l: number, r: number) {
    const { cmp, swap, data, parts } = ctx;
    let equals = l + 1, tail = r;

    // move the median to the 1st spot
    swap(data, l, medianPivotIndex(data, cmp, l, r));

    while (cmp(data, tail, l) > 0) { --tail; }
    for (let i = l + 1; i <= tail; i++) {
        const c = cmp(data, i, l);
        if (c > 0) {
            swap(data, i, tail);
            --tail;
            while (cmp(data, tail, l) > 0) { --tail; }
            i--;
        } else if (c === 0) {
            swap(data, i, equals);
            equals++;
        }
    }

    // move the medians to the correct spots
    for (let i = l; i < equals; i++) { swap(data, i, l + tail - i); }
    parts[0] = tail - equals + l + 1;
    parts[1] = tail;
}

function insertionSort({ data, cmp, swap }: Ctx, start: number, end: number) {
    for (let i = start + 1; i <= end; i++) {
        let j = i - 1;
        while (j >= start && cmp(data, j, j + 1) > 0) {
            swap(data, j, j + 1);
            j = j - 1;
        }
    }
}

function quickSort(ctx: Ctx, low: number, high: number) {
    const { parts } = ctx;
    while (low < high) {
        if (high - low < 16) {
            insertionSort(ctx, low, high);
            return;
        }

        partition(ctx, low, high);
        const li = parts[0], ri = parts[1];

        if (li - low < high - ri) {
            quickSort(ctx, low, li - 1);
            low = ri + 1;
        } else {
            quickSort(ctx, ri + 1, high);
            high = li - 1;
        }
    }
}

function partitionArrayAsc(data: number[], parts: number[], l: number, r: number) {
    let equals = l + 1, tail = r;

    // move the median to the 1st spot
    arraySwap(data, l, medianPivotIndex(data, arrayLess, l, r));
    const pivot = data[l];

    while (data[tail] > pivot) { --tail; }
    for (let i = l + 1; i <= tail; i++) {
        const v = data[i];
        if (v > pivot) {
            arraySwap(data, i, tail);
            --tail;
            while (data[tail] > pivot) { --tail; }
            i--;
        } else if (v === pivot) {
            arraySwap(data, i, equals);
            ++equals;
        }
    }

    // move all medians to the correct spots
    for (let i = l; i < equals; i++) { arraySwap(data, i, l + tail - i); }
    parts[0] = tail - equals + l + 1;
    parts[1] = tail;
}

function insertionSortArrayAsc(data: number[], start: number, end: number) {
    for (let i = start + 1; i <= end; i++) {
        const key = data[i];
        let j = i - 1;
        while (j >= start && data[j] > key) {
            data[j + 1] = data[j];
            j = j - 1;
        }
        data[j + 1] = key;
    }
}

function quickSortArrayAsc(data: number[], parts: number[], low: number, high: number) {
    while (low < high) {
        if (high - low < 16) {
            insertionSortArrayAsc(data, low, high);
            return;
        }

        partitionArrayAsc(data, parts, low, high);
        const li = parts[0], ri = parts[1];

        if (li - low < high - ri) {
            quickSortArrayAsc(data, parts, low, li - 1);
            low = ri + 1;
        } else {
            quickSortArrayAsc(data, parts, ri + 1, high);
            high = li - 1;
        }
    }
}

export function sortArray(data: ArrayLike<number>, cmp: Comparer<ArrayLike<number>> = arrayLess): ArrayLike<number> {
    return sortArrayRange(data, 0, data.length, cmp);
}

export function sortArrayRange(data: ArrayLike<number>, start: number, end: number, cmp: Comparer<ArrayLike<number>> = arrayLess): ArrayLike<number> {
    if (cmp === arrayLess) quickSortArrayAsc(data as any, [0, 0], start, end - 1);
    else quickSort({ data, cmp, swap: arraySwap, parts: [0, 0] }, start, end - 1);
    return data;
}

export function sort<T>(data: T, start: number, end: number, cmp: Comparer<T>, swap: Swapper<T>): T {
    const ctx: Ctx = { data, cmp, swap, parts: [0, 0] };
    quickSort(ctx, start, end - 1);
    return data;
}