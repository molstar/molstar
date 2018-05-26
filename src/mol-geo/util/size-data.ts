/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import VertexMap from '../shape/vertex-map';

export type SizeData = {
    uSize: ValueCell<number>,
    aSize: ValueCell<Float32Array>,
    dSizeType: ValueCell<string>,
}

export interface UniformSizeProps {
    value: number
}

/** Creates size uniform */
export function createUniformSize(props: UniformSizeProps): SizeData {
    return {
        uSize: ValueCell.create(props.value),
        aSize: ValueCell.create(new Float32Array(0)),
        dSizeType: ValueCell.create('uniform'),
    }
}

export interface AttributeSizeProps {
    sizeFn: (elementIdx: number) => number
    vertexMap: VertexMap
}

/** Creates size attribute with size for each element (i.e. shared across indtances/units) */
export function createAttributeSize(props: AttributeSizeProps): SizeData {
    const { sizeFn, vertexMap } = props
    const { idCount, offsetCount, offsets } = vertexMap
    const sizes = new Float32Array(idCount);
    for (let i = 0, il = offsetCount - 1; i < il; ++i) {
        const start = offsets[i]
        const end = offsets[i + 1]
        const size = sizeFn(i)
        for (let i = start, il = end; i < il; ++i) {
            sizes[i] = size
        }
    }
    return {
        uSize: ValueCell.create(0),
        aSize: ValueCell.create(sizes),
        dSizeType: ValueCell.create('attribute'),
    }
}