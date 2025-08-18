/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
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
    a = a ^ (a >> 4);
    a = (a ^ 0xdeadbeef) + (a << 5);
    a = a ^ (a >> 11);
    return a;
}

export function hash3(i: number, j: number, k: number) {
    let a = 23;
    a = (31 * a + i) | 0;
    a = (31 * a + j) | 0;
    a = (31 * a + k) | 0;
    a = a ^ (a >> 4);
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
    a = a ^ (a >> 4);
    a = (a ^ 0xdeadbeef) + (a << 5);
    a = a ^ (a >> 11);
    return a;
}

export function hashString(s: string) {
    let h = 0;
    for (let i = 0, l = s.length; i < l; i++) {
        h = (h << 5) - h + s.charCodeAt(i) | 0;
    }
    return h;
}

/**
 * A unique number for each pair of integers
 * Biggest representable pair is (67108863, 67108863) (limit imposed by Number.MAX_SAFE_INTEGER)
 */
export function cantorPairing(a: number, b: number) {
    return (a + b) * (a + b + 1) / 2 + b;
}

/**
 * A unique number for each sorted pair of integers
 * Biggest representable pair is (67108863, 67108863) (limit imposed by Number.MAX_SAFE_INTEGER)
 */
export function sortedCantorPairing(a: number, b: number) {
    return a < b ? cantorPairing(a, b) : cantorPairing(b, a);
}

export function invertCantorPairing(out: [number, number], z: number) {
    const w = Math.floor((Math.sqrt(8 * z + 1) - 1) / 2);
    const t = (w * w + w) / 2;
    const y = z - t;
    out[0] = w - y;
    out[1] = y;
    return out;
}

/**
 * 32 bit FNV-1a hash, see http://isthe.com/chongo/tech/comp/fnv/
 */
export function hashFnv32a(array: ArrayLike<number>) {
    let hval = 0x811c9dc5;
    for (let i = 0, il = array.length; i < il; ++i) {
        hval ^= array[i];
        hval += (hval << 1) + (hval << 4) + (hval << 7) + (hval << 8) + (hval << 24);
    }
    return hval >>> 0;
}

/**
 * 256 bit FNV-1a hash, returns 8 32-bit words
 * Based on the FNV-1a algorithm extended to 256 bits
 */
export function hashFnv256a(array: ArrayLike<number>, out: Uint32Array) {
    out.set(Fnv256Base);

    for (let i = 0, il = array.length; i < il; ++i) {
        // XOR with input byte
        out[0] ^= array[i] & 0xff;

        // Multiply by FNV prime (256-bit multiplication)
        multiplyBy256BitPrime(out);
    }

    return out;
}

/**
 * 256-bit object hash function using FNV-1a
 */
export function hashFnv256o(obj: any): string {
    return _Hasher256.hash(obj);
}

class ObjectHasher256 {
    private hashTarget: Uint32Array = new Uint32Array(8);
    private numberBytes = new Uint8Array(8);
    private numberView = new DataView(this.numberBytes.buffer);

    hash(obj: any): string {
        this.hashTarget.set(Fnv256Base);
        this.hashValue(obj, 0);
        return hashFnv256aToHex(this.hashTarget);
    }

    private hashValue(value: any, depth: number): void {
        if (depth > 50) return;

        const type = typeof value;
        this.addByte(type.charCodeAt(0));

        switch (type) {
            case 'string':
                this.addString(value);
                break;
            case 'number':
                this.addNumber(value);
                break;
            case 'boolean':
                this.addByte(value ? 1 : 0);
                break;
            case 'object':
                if (value === null) {
                    this.addByte(0);
                } else if (Array.isArray(value)) {
                    this.addArray(value, depth);
                } else {
                    this.addObject(value, depth);
                }
                break;
            case 'undefined':
                this.addByte(255);
                break;
        }
    }

    private addByte(byte: number): void {
        // XOR with input byte
        this.hashTarget[0] ^= byte & 0xff;
        // Multiply by FNV prime (256-bit multiplication)
        multiplyBy256BitPrime(this.hashTarget);
    }

    private addString(str: string): void {
        for (let i = 0; i < str.length; i++) {
            const code = str.charCodeAt(i);
            if (code < 128) {
                this.addByte(code);
            } else if (code < 2048) {
                this.addByte(0xc0 | (code >> 6));
                this.addByte(0x80 | (code & 0x3f));
            } else {
                this.addByte(0xe0 | (code >> 12));
                this.addByte(0x80 | ((code >> 6) & 0x3f));
                this.addByte(0x80 | (code & 0x3f));
            }
        }
    }

    private addNumber(num: number): void {
        if (Number.isNaN(num)) {
            this.addByte(0x7f); this.addByte(0xc0); this.addByte(0x00); this.addByte(0x00);
            this.addByte(0x00); this.addByte(0x00); this.addByte(0x00); this.addByte(0x00);
        } else if (!Number.isFinite(num)) {
            if (num > 0) {
                this.addByte(0x7f); this.addByte(0x80); this.addByte(0x00); this.addByte(0x00);
            } else {
                this.addByte(0xff); this.addByte(0x80); this.addByte(0x00); this.addByte(0x00);
            }
            this.addByte(0x00); this.addByte(0x00); this.addByte(0x00); this.addByte(0x00);
        } else {
            this.numberView.setFloat64(0, num, false);
            for (let i = 0; i < 8; i++) {
                this.addByte(this.numberBytes[i]);
            }
        }
    }

    private addArray(arr: any[], depth: number): void {
        this.addNumber(arr.length);
        for (let i = 0; i < arr.length; i++) {
            this.addNumber(i);
            this.hashValue(arr[i], depth + 1);
        }
    }

    private addObject(obj: any, depth: number): void {
        const keys = Object.keys(obj).sort();
        this.addNumber(keys.length);

        for (const key of keys) {
            this.addString(key);
            this.hashValue(obj[key], depth + 1);
        }
    }
}

const _Hasher256 = new ObjectHasher256();

const Fnv256Base = new Uint32Array([
    0x6c62272e, 0x07bb0142, 0x62b82175, 0x6295c58d,
    0x16d67530, 0xdd7121e3, 0xb3174000, 0x00000100
]);

const MultTmp1 = new Uint32Array(8);
const MultTmp2 = new Uint32Array(8);

/**
 * Helper function to multiply 256-bit number by FNV prime
 */
function multiplyBy256BitPrime(hash: Uint32Array): void {
    // Since FNV 256-bit prime is 2^88 + 2^8 + 0x3b, we can optimize:
    // hash * prime = hash * (2^88 + 2^8 + 0x3b) = (hash << 88) + (hash << 8) + hash * 0x3b

    // hash << 88 (shift left by 88 bits = 2 full 32-bit words + 24 bits)
    MultTmp1[0] = 0;
    MultTmp1[1] = 0;
    MultTmp1[2] = hash[0] << 24;
    MultTmp1[3] = (hash[0] >>> 8) | (hash[1] << 24);
    MultTmp1[4] = (hash[1] >>> 8) | (hash[2] << 24);
    MultTmp1[5] = (hash[2] >>> 8) | (hash[3] << 24);
    MultTmp1[6] = (hash[3] >>> 8) | (hash[4] << 24);
    MultTmp1[7] = (hash[4] >>> 8) | (hash[5] << 24);

    // hash << 8
    MultTmp2[0] = hash[0] << 8;
    MultTmp2[1] = (hash[0] >>> 24) | (hash[1] << 8);
    MultTmp2[2] = (hash[1] >>> 24) | (hash[2] << 8);
    MultTmp2[3] = (hash[2] >>> 24) | (hash[3] << 8);
    MultTmp2[4] = (hash[3] >>> 24) | (hash[4] << 8);
    MultTmp2[5] = (hash[4] >>> 24) | (hash[5] << 8);
    MultTmp2[6] = (hash[5] >>> 24) | (hash[6] << 8);
    MultTmp2[7] = (hash[6] >>> 24) | (hash[7] << 8);

    // hash * 0x3b (simple multiplication by small constant)
    let carry = 0;
    for (let i = 0; i < 8; i++) {
        const product = hash[i] * 0x3b + carry;
        hash[i] = product >>> 0;
        carry = Math.floor(product / 0x100000000);
    }

    // Add all three components: (hash << 88) + (hash << 8) + hash * 0x3b
    carry = 0;
    for (let i = 0; i < 8; i++) {
        const sum = hash[i] + MultTmp1[i] + MultTmp2[i] + carry;
        hash[i] = sum >>> 0;
        carry = sum >= 0x100000000 ? 1 : 0;
    }
}

const _8digit_padding = [
    '00000000',
    '0000000',
    '000000',
    '00000',
    '0000',
    '000',
    '00',
    '0'
];


function padHexNumber(num: number): string {
    const base = num.toString(16);
    if (base.length >= 8) return base; // No padding needed
    return _8digit_padding[base.length] + base;
}

/**
 * Convert 256-bit hash to hex string
 */
function hashFnv256aToHex(hash: Uint32Array): string {
    let result = '';
    for (let i = 7; i >= 0; i--) {
        result += padHexNumber(hash[i]);
    }
    return result;
}

/**
 * 32-bit Murmur hash
 */
export function hashMurmur32o(obj: any, seed: number = 42): number {
    const jsonString = JSON.stringify(obj);
    return murmurHash3_32(jsonString, seed);
}

/**
 * 128-bit Murmur hash
 */
export function hashMurmur128o(obj: any, seed: number = 42): string {
    const jsonString = JSON.stringify(obj);
    return murmurHash3_128(jsonString, seed);
}

/**
 * MurmurHash3 32-bit implementation
 * @param key - The input string to hash
 * @param seed - The seed value (default: 0)
 * @returns The 32-bit hash as a number
 */
export function murmurHash3_32(key: string, seed: number): number {
    let h = seed >>> 0;
    const remainder = key.length % 4;
    const bytes = key.length - remainder;

    for (let i = 0; i < bytes; i += 4) {
        let k = (key.charCodeAt(i) & 0xff) |
            ((key.charCodeAt(i + 1) & 0xff) << 8) |
            ((key.charCodeAt(i + 2) & 0xff) << 16) |
            ((key.charCodeAt(i + 3) & 0xff) << 24);

        k = Math.imul(k, 0xcc9e2d51);
        k = (k << 15) | (k >>> 17);
        k = Math.imul(k, 0x1b873593);

        h ^= k;
        h = (h << 13) | (h >>> 19);
        h = Math.imul(h, 5) + 0xe6546b64;
    }

    let k = 0;
    switch (remainder) {
        case 3: k ^= (key.charCodeAt(bytes + 2) & 0xff) << 16;
        case 2: k ^= (key.charCodeAt(bytes + 1) & 0xff) << 8;
        case 1: k ^= (key.charCodeAt(bytes) & 0xff);
            k = Math.imul(k, 0xcc9e2d51);
            k = (k << 15) | (k >>> 17);
            k = Math.imul(k, 0x1b873593);
            h ^= k;
    }

    h ^= key.length;
    h ^= h >>> 16;
    h = Math.imul(h, 0x85ebca6b);
    h ^= h >>> 13;
    h = Math.imul(h, 0xc2b2ae35);
    h ^= h >>> 16;

    return h >>> 0;
}

/**
 * MurmurHash3 128-bit implementation
 * @param key - The input data to hash
 * @param seed - The seed value (default: 0)
 * @returns The 128-bit hash as a hexadecimal string
 */
export function murmurHash3_128_fromBytes(key: Uint8Array, seed: number): string {
    // This fakeString approach is much faster than `new TextDecoder('ascii').decode(key)`
    const fakeString = {
        length: key.length,
        charCodeAt(i: number) { return key[i]; },
    };
    return murmurHash3_128(fakeString as string, seed);
}

/**
 * MurmurHash3 128-bit implementation
 * @param key - The input string to hash
 * @param seed - The seed value (default: 0)
 * @returns The 128-bit hash as a hexadecimal string
 */
export function murmurHash3_128(key: string, seed: number): string {
    let h1 = seed >>> 0;
    let h2 = seed >>> 0;
    let h3 = seed >>> 0;
    let h4 = seed >>> 0;

    const remainder = key.length % 16;
    const bytes = key.length - remainder;

    for (let i = 0; i < bytes; i += 16) {
        let k1 = (key.charCodeAt(i) & 0xff) |
            ((key.charCodeAt(i + 1) & 0xff) << 8) |
            ((key.charCodeAt(i + 2) & 0xff) << 16) |
            ((key.charCodeAt(i + 3) & 0xff) << 24);

        let k2 = (key.charCodeAt(i + 4) & 0xff) |
            ((key.charCodeAt(i + 5) & 0xff) << 8) |
            ((key.charCodeAt(i + 6) & 0xff) << 16) |
            ((key.charCodeAt(i + 7) & 0xff) << 24);

        let k3 = (key.charCodeAt(i + 8) & 0xff) |
            ((key.charCodeAt(i + 9) & 0xff) << 8) |
            ((key.charCodeAt(i + 10) & 0xff) << 16) |
            ((key.charCodeAt(i + 11) & 0xff) << 24);

        let k4 = (key.charCodeAt(i + 12) & 0xff) |
            ((key.charCodeAt(i + 13) & 0xff) << 8) |
            ((key.charCodeAt(i + 14) & 0xff) << 16) |
            ((key.charCodeAt(i + 15) & 0xff) << 24);

        k1 = Math.imul(k1, 0x239b961b);
        k1 = (k1 << 15) | (k1 >>> 17);
        k1 = Math.imul(k1, 0xab0e9789);
        h1 ^= k1;

        h1 = (h1 << 19) | (h1 >>> 13);
        h1 += h2;
        h1 = Math.imul(h1, 5) + 0x561ccd1b;

        k2 = Math.imul(k2, 0xab0e9789);
        k2 = (k2 << 16) | (k2 >>> 16);
        k2 = Math.imul(k2, 0x38b34ae5);
        h2 ^= k2;

        h2 = (h2 << 17) | (h2 >>> 15);
        h2 += h3;
        h2 = Math.imul(h2, 5) + 0x0bcaa747;

        k3 = Math.imul(k3, 0x38b34ae5);
        k3 = (k3 << 17) | (k3 >>> 15);
        k3 = Math.imul(k3, 0xa1e38b93);
        h3 ^= k3;

        h3 = (h3 << 15) | (h3 >>> 17);
        h3 += h4;
        h3 = Math.imul(h3, 5) + 0x96cd1c35;

        k4 = Math.imul(k4, 0xa1e38b93);
        k4 = (k4 << 13) | (k4 >>> 19);
        k4 = Math.imul(k4, 0x239b961b);
        h4 ^= k4;

        h4 = (h4 << 13) | (h4 >>> 19);
        h4 += h1;
        h4 = Math.imul(h4, 5) + 0x32ac3b17;
    }

    let k1 = 0, k2 = 0, k3 = 0, k4 = 0;

    switch (remainder) {
        case 15: k4 ^= key.charCodeAt(bytes + 14) << 16;
        case 14: k4 ^= key.charCodeAt(bytes + 13) << 8;
        case 13: k4 ^= key.charCodeAt(bytes + 12);
            k4 = Math.imul(k4, 0xa1e38b93);
            k4 = (k4 << 13) | (k4 >>> 19);
            k4 = Math.imul(k4, 0x239b961b);
            h4 ^= k4;

        case 12: k3 ^= key.charCodeAt(bytes + 11) << 24;
        case 11: k3 ^= key.charCodeAt(bytes + 10) << 16;
        case 10: k3 ^= key.charCodeAt(bytes + 9) << 8;
        case 9: k3 ^= key.charCodeAt(bytes + 8);
            k3 = Math.imul(k3, 0x38b34ae5);
            k3 = (k3 << 17) | (k3 >>> 15);
            k3 = Math.imul(k3, 0xa1e38b93);
            h3 ^= k3;

        case 8: k2 ^= key.charCodeAt(bytes + 7) << 24;
        case 7: k2 ^= key.charCodeAt(bytes + 6) << 16;
        case 6: k2 ^= key.charCodeAt(bytes + 5) << 8;
        case 5: k2 ^= key.charCodeAt(bytes + 4);
            k2 = Math.imul(k2, 0xab0e9789);
            k2 = (k2 << 16) | (k2 >>> 16);
            k2 = Math.imul(k2, 0x38b34ae5);
            h2 ^= k2;

        case 4: k1 ^= key.charCodeAt(bytes + 3) << 24;
        case 3: k1 ^= key.charCodeAt(bytes + 2) << 16;
        case 2: k1 ^= key.charCodeAt(bytes + 1) << 8;
        case 1: k1 ^= key.charCodeAt(bytes);
            k1 = Math.imul(k1, 0x239b961b);
            k1 = (k1 << 15) | (k1 >>> 17);
            k1 = Math.imul(k1, 0xab0e9789);
            h1 ^= k1;
    }

    h1 ^= key.length;
    h2 ^= key.length;
    h3 ^= key.length;
    h4 ^= key.length;

    h1 += h2;
    h1 += h3;
    h1 += h4;
    h2 += h1;
    h3 += h1;
    h4 += h1;

    h1 ^= h1 >>> 16;
    h1 = Math.imul(h1, 0x85ebca6b);
    h1 ^= h1 >>> 13;
    h1 = Math.imul(h1, 0xc2b2ae35);
    h1 ^= h1 >>> 16;

    h2 ^= h2 >>> 16;
    h2 = Math.imul(h2, 0x85ebca6b);
    h2 ^= h2 >>> 13;
    h2 = Math.imul(h2, 0xc2b2ae35);
    h2 ^= h2 >>> 16;

    h3 ^= h3 >>> 16;
    h3 = Math.imul(h3, 0x85ebca6b);
    h3 ^= h3 >>> 13;
    h3 = Math.imul(h3, 0xc2b2ae35);
    h3 ^= h3 >>> 16;

    h4 ^= h4 >>> 16;
    h4 = Math.imul(h4, 0x85ebca6b);
    h4 ^= h4 >>> 13;
    h4 = Math.imul(h4, 0xc2b2ae35);
    h4 ^= h4 >>> 16;

    h1 += h2;
    h1 += h3;
    h1 += h4;
    h2 += h1;
    h3 += h1;
    h4 += h1;

    return (
        (h1 >>> 0).toString(16).padStart(8, '0') +
        (h2 >>> 0).toString(16).padStart(8, '0') +
        (h3 >>> 0).toString(16).padStart(8, '0') +
        (h4 >>> 0).toString(16).padStart(8, '0')
    );
}