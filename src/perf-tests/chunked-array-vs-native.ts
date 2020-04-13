import * as B from 'benchmark';
import { ChunkedArray } from '../mol-data/util';

function testNative(size: number) {
    const xs = new Array(size);
    for (let i = 0; i < size; i++) xs[i] = i * i;
    return xs;
}

function testChunkedTyped(size: number, chunk: number) {
    const xs = ChunkedArray.create(Int32Array, 1, chunk);
    for (let i = 0; i < size; i++) ChunkedArray.add(xs, i * i);
    return ChunkedArray.compact(xs);
}

function testChunkedNative(size: number, chunk: number) {
    const xs = ChunkedArray.create(Array, 1, chunk);
    for (let i = 0; i < size; i++) ChunkedArray.add(xs, i * i);
    return ChunkedArray.compact(xs);
}

const suite = new B.Suite();

const N = 70000;

suite
    .add('native', () => testNative(N))
    // .add('chunkedT 0.1k', () => testChunkedTyped(N, 100, false))
    // .add('chunkedT 4k', () => testChunkedTyped(N, 4096, false))
    .add('chunkedT 4k lin', () => testChunkedTyped(N, 4096))
    // .add('chunkedT N / 2', () => testChunkedTyped(N, N / 2, false))
    // .add('chunkedT N', () => testChunkedTyped(N, N, false))
    // .add('chunkedT 2 * N', () => testChunkedTyped(N, 2 * N, false))

    .add('chunkedN N', () => testChunkedNative(N, N))
    .add('chunkedN 0.1k', () => testChunkedNative(N, 100))
    .add('chunkedN N / 2', () => testChunkedNative(N, N / 2))
    .add('chunkedN 2 * N', () => testChunkedNative(N, 2 * N))
    .on('cycle', (e: any) => {
        console.log(String(e.target));
    })
    .run();

// console.log(testChunkedTyped(10, 16));
