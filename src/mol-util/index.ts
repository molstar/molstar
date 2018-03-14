/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import BitFlags from './bit-flags'
import Computation from './computation'
import Scheduler from './scheduler'
import StringBuilder from './string-builder'
import Time from './time'
import UUID from './uuid'

export { BitFlags, Computation, Scheduler, StringBuilder, Time, UUID }

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