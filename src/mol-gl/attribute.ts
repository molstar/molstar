/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'

// export type AttributeGroupMutator<T extends AttributesData> = (data: T) => (boolean | void)
// export type AttributeGroupData = { [k: string]: Helpers.TypedArray }
// export type AttributesBuffers<T extends AttributesData> = { [K in keyof T]: REGL.Buffer }

// interface AttributeGroup<T extends AttributeGroup.Data> {
//     readonly buffer: REGL.Buffer
//     readonly length: number
//     readonly data: T
//     setCount(size: number): void
//     update(mutator: AttributeGroup.Mutator<T>): void
// }

// namespace AttributeGroup {
//     export type Data = { [k: string]: Helpers.TypedArray }
//     export type Mutator<T extends Data> = (data: T) => (UpdateInfo<T> | void)
//     export type UpdateInfo<T extends Data> = boolean | { [k in keyof T]: Attribute.UpdateInfo }
//     export type Attributes<T extends Data> = { [K in keyof T]: Attribute<T[K]> }

//     export function create<T extends Data>(regl: REGL.Regl, data: T): AttributeGroup<T> {
//         const attributes: Attributes<any> = {}
//         for (const k of Object.keys(data)) {
//             attributes[k] = Attribute.create(regl, data[k])
//         }
//         return {
//             update: (mutator: Mutator<T>) => {

//             }
//         }
//     }
// }


interface Attribute<T extends Helpers.TypedArray> {
    readonly buffer: REGL.AttributeConfig
    getCount(): number
    setCount(count: number): void
    getArray(): T
    set(index: number, ...values: number[]): void
    update(mutator: Attribute.Mutator<T>): void
    reload(): void
    destroy(): void
}

interface AttributeProps {
    size: 1 | 2 | 3 | 4,
    divisor?: number,
    offset?: number,
    stride?: number
}

namespace Attribute {
    export type Mutator<T extends Helpers.TypedArray> = (data: T) => (UpdateInfo | void)
    export type UpdateInfo = boolean | { offset: number, count: number }
    export type ArrayCell<T> = { array: ReferenceCell<T> }
    export type ReferenceCell<T> = { readonly version: number, readonly value: T }

    export function create<T extends Float32Array>(regl: REGL.Regl, array: ValueCell<T>, count: number, props: AttributeProps): Attribute<T> {
        const itemSize = props.size
        let _array = array.ref.value
        let _count = count
        // if (props.stride) _count = _array.length / (props.stride / _array.BYTES_PER_ELEMENT)
        // console.log(_array.length, props.stride)
        // console.log('buffer', {
        //     data: _array,
        //     length: count * itemSize,
        //     usage: 'dynamic',
        //     type: 'float32'
        // })
        const buffer = regl.buffer({
            data: _array,
            // length: count * itemSize * _array.BYTES_PER_ELEMENT,
            usage: 'dynamic',
            type: 'float32',
            dimension: itemSize
        } as any)
        // console.log(buffer)
        const attribute: REGL.AttributeConfig = { ...props, buffer }
        const growIfNeeded = function(count: number) {
            if (count * itemSize > _array.length) {
                const newArray: T = new (_array as any).constructor(count * itemSize)
                newArray.set(_array)
                _array = newArray
                buffer(_array)
            }
            _count = count
        }
        return {
            buffer: attribute,
            getCount: () =>  _count,
            setCount: (count: number) => growIfNeeded(count),
            getArray: () => _array,
            set: (index: number, ...values: number[]) => {
                if (values.length !== itemSize) throw new Error('wrong number of values given')
                growIfNeeded(index)
                for (let i = 0; i < itemSize; ++i) {
                    _array[index * itemSize + i] = values[i]
                }
                buffer.subdata(values, index * itemSize * _array.BYTES_PER_ELEMENT)
            },
            update: (mutator: Mutator<T>, offset?: number, count?: number) => {
                if (offset && count) growIfNeeded(offset + count)
                mutator(_array)
                buffer(_array)
            },
            reload: () => buffer(_array),
            destroy: () => buffer.destroy()
        }
    }
}

export default Attribute