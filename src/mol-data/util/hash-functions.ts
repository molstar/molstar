/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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