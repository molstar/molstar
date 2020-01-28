/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { deflate, inflate, parse, encode } from '../zip/zip'

describe('zip', () => {
    it('roundtrip deflate/inflate', () => {
        const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7])
        const deflated = deflate(data)
        console.log(deflated)
        const inflated = inflate(deflated)
        console.log(inflated)
    })

    it('roundtrip zip', () => {
        const zipped = encode({
            'test.foo': new Uint8Array([1, 2, 3, 4, 5, 6, 7])
        })
        console.log(zipped)
        const unzipped = parse(zipped)
        console.log(unzipped)
    })
})