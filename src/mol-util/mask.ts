/**
 * Copyright (c) 2017 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO check if the removal of FastSet and the removal of the context object for forEach
// have any performance implications

function _ascSort(a: number, b: number) {
    return a - b;
}

export function sortAsc<T extends ArrayLike<number>>(array: T): T {
    Array.prototype.sort.call(array, _ascSort);
    return array;
}

interface Mask {
    '@type': 'mask'
    size: number;
    has(i: number): boolean;
    /** in-order iteration of all "masked elements". */
    forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx?: Ctx): Ctx | undefined;
}

namespace Mask {
    class EmptyMask implements Mask {
        '@type': 'mask'
        size = 0;
        has(i: number) { return false; }
        forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx: Ctx) { return ctx; }
        constructor() { }
    }

    class SingletonMask implements Mask {
        '@type': 'mask'
        size = 1;
        has(i: number) { return i === this.idx; }
        forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx?: Ctx) { f(this.idx, ctx); return ctx; }
        constructor(private idx: number) { }
    }

    class BitMask implements Mask {
        '@type': 'mask'
        private length: number;
        has(i: number) { return i < this.length && !!this.mask[i] as any; }

        private _forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx: Ctx | undefined) {
            for (let i = 0; i < this.length; i++) {
                if (this.mask[i]) f(i, ctx);
            }
        }
        forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx?: Ctx) {
            this._forEach(f, ctx);
            return ctx;
        }
        constructor(private mask: boolean[], public size: number) { this.length = mask.length;  }
    }

    class AllMask implements Mask {
        '@type': 'mask'
        has(i: number) { return true; }
        private _forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx: Ctx | undefined) {
            for (let i = 0; i < this.size; i++) {
                f(i, ctx);
            }
        }
        forEach<Ctx>(f: (i: number, ctx?: Ctx) => void, ctx?: Ctx) {
            this._forEach(f, ctx);
            return ctx;
        }
        constructor(public size: number) { }
    }

    class SetMask implements Mask {
        '@type': 'mask'
        private _flat: number[] | undefined = void 0;
        size: number;
        has(i: number) { return this.set.has(i); }

        private _forEach<Ctx>(f: (i: number, ctx: Ctx) => void, ctx: Ctx) {
            for (const idx of this.flatten()) {
                f(idx, ctx);
            }
        }
        private flatten() {
            if (this._flat) return this._flat;
            const indices = new Int32Array(this.size);
            let offset = 0;
            this.set.forEach(i => indices[offset++] = i);
            sortAsc(indices);
            this._flat = indices as any as number[];
            return this._flat;
        }
        forEach<Ctx>(f: (i: number, ctx: Ctx) => void, ctx: Ctx) {
            this._forEach(f, ctx);
            return ctx;
        }
        constructor(private set: Set<number>) {
            this.size = set.size;
        }
    }

    export function always(size: number) { return new AllMask(size); }
    export const never = new EmptyMask();

    export function ofSet(set: Set<number>): Mask {
        return new SetMask(set);
    }

    export function singleton(i: number) {
        return new SingletonMask(i);
    }

    export function ofUniqueIndices(indices: ArrayLike<number>): Mask {
        const len = indices.length;
        if (len === 0) return new EmptyMask();
        if (len === 1) return new SingletonMask(indices[0]);

        let max = 0;
        for (const i of (indices as number[])) {
            if (i > max) max = i;
        }
        if (len === max) return new AllMask(len);

        const f = len / max;
        if (f < 1 / 12) {
            const set = new Set<number>();
            for (const i of (indices as number[])) set.add(i);
            return new SetMask(set);
        }

        const mask = new Int8Array(max + 1);
        for (const i of (indices as number[])) {
            mask[i] = 1;
        }
        return new BitMask(mask as any as boolean[], indices.length);
    }

    export function ofMask(mask: boolean[], size: number): Mask {
        return new BitMask(mask, size);
    }

    export function hasAny(mask: Mask, xs: number[]) {
        for (const x of xs) {
            if (mask.has(x)) return true;
        }
        return false;
    }

    export function complement(mask: Mask, against: Mask) {
        let count = 0;
        let max = 0;
        against.forEach(i => {
            if (!mask.has(i)) {
                count++;
                if (i > max) max = i;
            }
        });

        if (count / max < 1 / 12) {
            // set based
            const set = new Set<number>();
            against.forEach(i => {
                if (!mask.has(i)) {
                    set.add(i);
                }
            });
            return ofSet(set);
        } else {
            // mask based
            const target = new Uint8Array(max + 1);
            against.forEach(i => {
                if (!mask.has(i)) {
                    target[i] = 1;
                }
            });
            return ofMask(target as any as boolean[], count);
        }
    }
}

export default Mask;