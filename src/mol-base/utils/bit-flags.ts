/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface BitFlags<Flags> { '@type': Flags }

namespace BitFlags {
    export function has<F>(flags: BitFlags<F>, flag: F) { return ((flags as any) & (flag as any)) !== 0; }
}

export default BitFlags