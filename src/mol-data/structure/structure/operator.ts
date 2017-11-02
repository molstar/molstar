/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from 'mol-base/math/linear-algebra'

interface Operator extends Readonly<{
    name: string,
    hkl: number[], // defaults to [0, 0, 0] for non symmetry entries
    transform: Mat4,
    // cache the inverse of the transform
    inverse: Mat4,
    // optimize the identity case
    isIdentity: boolean
}> { }

namespace Operator {
    export const Identity: Operator = { name: '1_555', hkl: [0, 0, 0], transform: Mat4.identity(), inverse: Mat4.identity(), isIdentity: true };
}

export default Operator