/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';

export type UniformSize = { type: 'uniform', value: number }
export type AttributeSize = { type: 'attribute', value: ValueCell<Float32Array> }
export type SizeData = UniformSize | AttributeSize

export interface OffsetData {
    primitiveCount: number,
    offsetCount: number,
    offsets: Uint32Array
}

export interface UniformSizeProps {
    value: number
}

/** Creates size uniform */
export function createUniformSize(props: UniformSizeProps): UniformSize {
    return { type: 'uniform', value: props.value }
}

export interface AttributeSizeProps {
    sizeFn: (elementIdx: number) => number
    offsetData: OffsetData
}

/** Creates size attribute with size for each element (i.e. shared across indtances/units) */
export function createAttributeSize(props: AttributeSizeProps): AttributeSize {
    const { sizeFn, offsetData } = props
    const { primitiveCount, offsetCount, offsets } = offsetData
    const sizes = new Float32Array(primitiveCount);
    const _offsets = offsets // .ref.value
    for (let i = 0, il = offsetCount - 1; i < il; ++i) {
        const start = _offsets[i]
        const end = _offsets[i + 1]
        const size = sizeFn(i)
        for (let i = start, il = end; i < il; ++i) {
            sizes[i] = size
        }
    }
    return { type: 'attribute', value: ValueCell.create(sizes) }
}