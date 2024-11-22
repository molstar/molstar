import { findPredecessorIndex } from '../mol-data/int/impl/sorted-array';


function createData(n: number) {
    const data = []; // new Int32Array(n);
    let last = (15 * Math.random()) | 0;
    for (let i = 0; i < n; i++) {
        data[i] = last;
        last += (15 * Math.random()) | 0;
    }
    return data;
}

function binarySearchHelper(list: ArrayLike<number>, t: number) {
    let min = 0, max = list.length - 1;
    while (min <= max) {
        if (min + 11 > max) {
            for (let i = min; i <= max; i++) {
                if (t === list[i]) return i;
            }
            return -1;
        }

        const mid = (min + max) >> 1;
        const v = list[mid];
        if (t < v) max = mid - 1;
        else if (t > v) min = mid + 1;
        else return mid;
    }
    return -1;
}

function findPredecessorIndexHelper_old(list: ArrayLike<number>, t: number) {
    const index = SortedArray_old.findPredecessorIndex(list, t);
    if (list[index] === t) return index;
    return -1;
}

function findPredecessorIndexHelper_new(list: ArrayLike<number>, t: number) {
    const index = findPredecessorIndex(list, t);
    if (list[index] === t) return index;
    return -1;
}

/** Old implementation of `findPredecessorIndex` from src/mol-data/int/impl/sorted-array.ts */
namespace SortedArray_old {
    export function findPredecessorIndex(xs: ArrayLike<number>, v: number) {
        const len = xs.length;
        if (v <= xs[0]) return 0;
        if (v > xs[len - 1]) return len;
        return binarySearchPredIndexRange(xs, v, 0, len);
    }

    function binarySearchPredIndexRange(xs: ArrayLike<number>, value: number, start: number, end: number) {
        let min = start, max = end - 1;
        while (min < max) {
            // do a linear search if there are only 10 or less items remaining
            if (min + 11 > max) {
                for (let i = min; i <= max; i++) {
                    if (value <= xs[i]) {
                        return i === min ? normalizePred(xs, start, end, i) : i;
                    }
                }
                return max + 1;
            }
            const mid = (min + max) >> 1;
            const v = xs[mid];
            if (value < v) max = mid - 1;
            else if (value > v) min = mid + 1;
            else return normalizePred(xs, start, end, mid);
        }
        if (min > max) return max + 1;
        return xs[min] >= value ? min : min + 1;
    }

    function normalizePred(xs: ArrayLike<number>, start: number, end: number, idx: number) {
        if (idx <= start) return idx;

        const value = xs[idx];
        let i: number;
        for (i = idx - 1; i >= start; i--) {
            if (xs[i] < value) break;
        }
        return i + 1;
    }
}

function objSearch(obj: any, val: number) {
    return typeof obj[val] !== 'undefined';
}

function setSearch(set: Set<number>, val: number) {
    return set.has(val);
}

type Mask = { min: number, max: number, mask: ArrayLike<number> }
function maskSearch({ min, max, mask }: Mask, val: number) {
    return val >= min && val <= max && !!mask[val - min];
}

function prepareSet(list: ArrayLike<number>) {
    const set = new Set<number>();
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        set.add(v);
    }
    return set;
}

function prepareObj(list: ArrayLike<number>) {
    const obj = Object.create(null);
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        obj[v] = i;
    }

    return obj;
}

function prepareMask(list: ArrayLike<number>): Mask {
    let max = Number.NEGATIVE_INFINITY, min = Number.POSITIVE_INFINITY;
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        if (max < v) max = v;
        if (min > v) min = v;
    }
    const mask = new Uint8Array(max - min + 1);
    for (let i = 0; i < list.length; i++) {
        const v = list[i];
        mask[v - min] = 1;
    }

    return { min, max, mask };
}

function testBinary(list: ArrayLike<number>, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (binarySearchHelper(list, points[i]) >= 0) r += points[i];
    }
    return r;
}

function testObj(obj: any, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (objSearch(obj, points[i])) r += points[i];
    }
    return r;
}

function testSet(set: Set<number>, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (setSearch(set, points[i])) r += points[i];
    }
    return r;
}

function testMask(mask: Mask, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (maskSearch(mask, points[i])) r += points[i];
    }
    return r;
}

function testFindPredecessorIndex_old(list: ArrayLike<number>, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (findPredecessorIndexHelper_old(list, points[i]) >= 0) r += points[i];
    }
    return r;
}

function testFindPredecessorIndex_new(list: ArrayLike<number>, points: ArrayLike<number>) {
    let r = 0;
    for (let i = 0, _i = points.length; i < _i; i++) {
        if (findPredecessorIndexHelper_new(list, points[i]) >= 0) r += points[i];
    }
    return r;
}

function run(f: () => number, n: number) {
    for (let i = 0; i < n; i++) f();
}

(function () {
    const size = 10000;
    const list = createData(size);
    const queryPoints = createData(size);

    const obj = prepareObj(list);
    const set = prepareSet(list);
    const mask = prepareMask(list);

    console.log('list', testBinary(list, queryPoints));
    console.log('obj', testObj(obj, queryPoints));
    console.log('set', testSet(set, queryPoints));
    console.log('mask', testMask(mask, queryPoints));
    console.log('list-findPredecessorIndex-old', testFindPredecessorIndex_old(list, queryPoints));
    console.log('list-findPredecessorIndex-new', testFindPredecessorIndex_new(list, queryPoints));

    console.log('----------------------');

    console.time('obj');
    run(() => testObj(obj, queryPoints), 100);
    console.timeEnd('obj');

    console.time('set');
    run(() => testSet(set, queryPoints), 100);
    console.timeEnd('set');

    console.time('bin-search');
    run(() => testBinary(list, queryPoints), 100);
    console.timeEnd('bin-search');

    console.time('mask-search');
    run(() => testMask(mask, queryPoints), 100);
    console.timeEnd('mask-search');

    console.time('findPredecessorIndex-old');
    run(() => testFindPredecessorIndex_old(list, queryPoints), 100);
    console.timeEnd('findPredecessorIndex-old');

    console.time('findPredecessorIndex-new');
    run(() => testFindPredecessorIndex_new(list, queryPoints), 100);
    console.timeEnd('findPredecessorIndex-new');

    console.log('----------------------');

    console.time('prepare-obj');
    run(() => prepareObj(list), 1);
    console.timeEnd('prepare-obj');

    console.time('prepare-set');
    run(() => prepareSet(list).size, 1);
    console.timeEnd('prepare-set');

    console.time('prepare-mask');
    run(() => prepareMask(list).min, 1);
    console.timeEnd('prepare-mask');
}());
