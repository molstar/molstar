/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');

interface Attribute<T extends Helpers.TypedArray> {
    readonly buffer: REGL.Buffer
    readonly length: number
    readonly data: T
    set(value: number, index: number): void
    growIfNeeded(size: number): void
    update(mutator: Attribute.Mutator<T>): void
    reload(): void
}

namespace Attribute {
    export type Mutator<T extends Helpers.TypedArray> = (data: T) => (boolean | void)

    export function create<T extends Helpers.TypedArray>(regl: REGL.Regl, data: T): Attribute<T> {
        let _data = data
        let _length = _data.length
        const buffer = regl.buffer(_data)
        const growIfNeeded = function(length: number) {
            if (length > _data.length) {
                _data = new (_data as any).constructor(_data)
                _length = _data.length
                buffer(_data)
            }
        }
        return {
            buffer,
            get length() { return _length },
            get data() { return _data },
            set: (value: number, index: number) => {
                growIfNeeded(index)
                _data[index] = value
                buffer.subdata([value], index * data.BYTES_PER_ELEMENT)
            },
            growIfNeeded(size: number) {
                growIfNeeded(size)
            },
            update: (mutator: Mutator<T>) => {
                mutator(_data)
                buffer(_data)
            },
            reload: () => {
                buffer(_data)
            }
        }
    }
}

export default Attribute