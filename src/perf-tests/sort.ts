import * as B from 'benchmark'
import * as Sort from '../structure/collections/sort'

function createTestData(n: number) {
    const data = new Int32Array(n); //new Array(n);
    for (let i = 0; i < n; i++) {
        data[i] = (n * Math.random()) | 0;
    }
    return data;
}

export function copyArray(xs: any) {
    const ret = new xs.constructor(xs.length);
    for (let i = 0, _i = xs.length; i < _i; i++) ret[i] = xs[i];
    return ret;
}

export function checkSorted(arr: ArrayLike<number>) {
    for (let i = 0; i < arr.length - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            return false;
        }
    }
    return true;
}

const _data = createTestData(10000);
const data = () => copyArray(_data);

const suite = new B.Suite();

function le(x: number, y: number) { return x - y; }

suite
    .add('native', () => Array.prototype.sort.call(data(), le))
    .add('qsort (array asc)', () => Sort.sortArray(data()))
    .add('qsort (generic)', () => Sort.sort(data(), _data.length, Sort.arrayLess, Sort.arraySwap))
    .add('native sorted', () => Array.prototype.sort.call(_data, le))
    .add('qsort sorted (array asc)', () => Sort.sortArray(_data))
    .add('qsort sorted (generic)', () => Sort.sort(_data, _data.length, Sort.arrayLess, Sort.arraySwap))
    .on('cycle', (e: any) => {
        console.log(String(e.target));
    })
    .run();