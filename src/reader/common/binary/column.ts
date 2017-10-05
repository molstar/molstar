/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getArrayBounds, ToArrayParams } from '../column'

export function isTypedArray(data: any) {
    return data.buffer && typeof data.byteLength === 'number' && data.BYTES_PER_ELEMENT;
}

export function typedArrayWindow(data: any, params?: ToArrayParams): ArrayLike<number> {
    const { constructor, buffer, length, byteOffset, BYTES_PER_ELEMENT } = data;
    const { start, end } = getArrayBounds(length, params);
    if (start === 0 && end === length) return data;
    return new constructor(buffer, byteOffset + BYTES_PER_ELEMENT * start, Math.min(length, end - start));
}