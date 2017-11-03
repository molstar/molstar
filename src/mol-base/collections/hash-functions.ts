/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// from http://burtleburtle.net/bob/hash/integer.html
export function hash1(i: number) {
    let a = i ^ (i >> 4);
    a = (a ^ 0xdeadbeef) + (a << 5);
    a = a ^ (a >> 11);
    return a;
}

export function hash2(i: number, j: number) {
    let a = 23;
    a = (31 * a + i) | 0;
    a = (31 * a + j) | 0;
    a = a ^ (a >> 4)
    a = (a ^ 0xdeadbeef) + (a << 5);
    a = a ^ (a >> 11);
    return a;
}

export function hash3(i: number, j: number, k: number) {
    let a = 23;
    a = (31 * a + i) | 0;
    a = (31 * a + j) | 0;
    a = (31 * a + k) | 0;
    a = a ^ (a >> 4)
    a = (a ^ 0xdeadbeef) + (a << 5);
    a = a ^ (a >> 11);
    return a;
}

export function hash4(i: number, j: number, k: number, l: number) {
    let a = 23;
    a = (31 * a + i) | 0;
    a = (31 * a + j) | 0;
    a = (31 * a + k) | 0;
    a = (31 * a + l) | 0;
    a = a ^ (a >> 4)
    a = (a ^ 0xdeadbeef) + (a << 5);
    a = a ^ (a >> 11);
    return a;
}

/**
 * A unique number for each pair of integers
 * Biggest representable pair is (67108863, 67108863) (limit imposed by Number.MAX_SAFE_INTEGER)
 */
export function cantorPairing(a: number, b: number) {
    return (a + b) * (a + b + 1) / 2 + b;
}