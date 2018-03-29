/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import BitFlags from './bit-flags'
import StringBuilder from './string-builder'
import UUID from './uuid'
import Mask from './mask'

export { BitFlags, StringBuilder, UUID, Mask }

export function arrayEqual<T>(arr1: T[], arr2: T[]) {
    const length = arr1.length
    if (length !== arr2.length) return false
    for (let i = 0; i < length; i++) {
        if (arr1[i] !== arr2[i]) {
            return false
        }
    }
    return true
}