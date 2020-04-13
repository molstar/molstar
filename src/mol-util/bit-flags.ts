/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

type BitFlags<Flags> = number & Flags

namespace BitFlags {
    export function create<F>(flags: F): BitFlags<F> { return flags as any; }

    export function has<F>(flags: BitFlags<F>, flag: F) { return ((flags as any) & (flag as any)) !== 0; }
    /** toCheck must be non-zero */
    export function hasAll<F>(flags: BitFlags<F>, toCheck: BitFlags<F>) { return !!toCheck && ((flags as any) & (toCheck as any)) === (toCheck as any); }
}

export default BitFlags;