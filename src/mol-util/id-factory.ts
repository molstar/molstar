/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function idFactory(firstId = 0) {
    let _nextId = firstId
    return () => _nextId++
}