/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from './context'
import { ValueCell } from 'mol-util';
import { RenderableSchema } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { ValueOf } from 'mol-util/type-helpers';

const getNextBufferId = idFactory()

export type UsageHint = 'static' | 'dynamic' | 'stream'
export type DataType = 'uint8' | 'int8' | 'uint16' | 'int16' | 'uint32' | 'int32' | 'float32'
export type BufferType = 'attribute' | 'elements'

export type DataTypeArrayType = {
    'uint8': Uint8Array
    'int8': Int8Array
    'uint16': Uint16Array
    'int16': Int16Array
    'uint32': Uint32Array
    'int32': Int32Array
    'float32': Float32Array
}
export type ArrayType = ValueOf<DataTypeArrayType>
export type ArrayKind = keyof DataTypeArrayType

export type BufferItemSize = 1 | 2 | 3 | 4 | 16

export function getUsageHint(ctx: WebGLContext, usageHint: UsageHint) {
    const { gl } = ctx
    switch (usageHint) {
        case 'static': return gl.STATIC_DRAW
        case 'dynamic': return gl.DYNAMIC_DRAW
        case 'stream': return gl.STREAM_DRAW
    }
}

export function getDataType(ctx: WebGLContext, dataType: DataType) {
    const { gl } = ctx
    switch (dataType) {
        case 'uint8': return gl.UNSIGNED_BYTE
        case 'int8': return gl.BYTE
        case 'uint16': return gl.UNSIGNED_SHORT
        case 'int16': return gl.SHORT
        case 'uint32': return gl.UNSIGNED_INT
        case 'int32': return gl.INT
        case 'float32': return gl.FLOAT
    }
}

function dataTypeFromArray(ctx: WebGLContext, array: ArrayType) {
    const { gl } = ctx
    if (array instanceof Uint8Array) {
        return gl.UNSIGNED_BYTE
    } else if (array instanceof Int8Array) {
        return gl.BYTE
    } else if (array instanceof Uint16Array) {
        return gl.UNSIGNED_SHORT
    } else if (array instanceof Int16Array) {
        return gl.SHORT
    } else if (array instanceof Uint32Array) {
        return gl.UNSIGNED_INT
    } else if (array instanceof Int32Array) {
        return gl.INT
    } else if (array instanceof Float32Array) {
        return gl.FLOAT
    } else {
        throw new Error('Should nevver happen')
    }
}

export function getBufferType(ctx: WebGLContext, bufferType: BufferType) {
    const { gl } = ctx
    switch (bufferType) {
        case 'attribute': return gl.ARRAY_BUFFER
        case 'elements': return gl.ELEMENT_ARRAY_BUFFER
    }
}

export interface Buffer {
    readonly id: number

    readonly _buffer: WebGLBuffer
    readonly _usageHint: number
    readonly _bufferType: number
    readonly _dataType: number
    readonly _bpe: number

    readonly itemSize: number
    readonly itemCount: number
    readonly length: number

    updateData: (array: ArrayType) => void
    updateSubData: (array: ArrayType, offset: number, count: number) => void
    destroy: () => void
}

export function createBuffer(ctx: WebGLContext, array: ArrayType, itemSize: BufferItemSize, usageHint: UsageHint, bufferType: BufferType): Buffer {
    const { gl } = ctx
    const _buffer = gl.createBuffer()
    if (_buffer === null) {
        throw new Error('Could not create WebGL buffer')
    }

    const _usageHint = getUsageHint(ctx, usageHint)
    const _bufferType = getBufferType(ctx, bufferType)
    const _dataType = dataTypeFromArray(ctx, array)
    const _bpe = array.BYTES_PER_ELEMENT
    const _length = array.length
    const _itemCount = Math.floor(_length / itemSize)

    function updateData(array: ArrayType) {
        gl.bindBuffer(_bufferType, _buffer);
        (gl as WebGLRenderingContext).bufferData(_bufferType, array, _usageHint) // TODO remove cast when webgl2 types are fixed
    }
    updateData(array)

    let destroyed = false
    ctx.bufferCount += 1

    return {
        id: getNextBufferId(),

        _buffer,
        _usageHint,
        _bufferType,
        _dataType,
        _bpe,

        get itemSize () { return itemSize },
        get itemCount () { return _itemCount },
        get length () { return _length },

        updateData,
        updateSubData: (array: ArrayType, offset: number, count: number) => {
            gl.bindBuffer(_bufferType, _buffer);
            (gl as WebGLRenderingContext).bufferSubData(_bufferType, offset * _bpe, array.subarray(offset, offset + count)) // TODO remove cast when webgl2 types are fixed
        },

        destroy: () => {
            if (destroyed) return
            gl.deleteBuffer(_buffer)
            destroyed = true
            ctx.bufferCount -= 1
        }
    }
}

export type AttributeDefs = {
    [k: string]: { kind: ArrayKind, itemSize: BufferItemSize, divisor: number }
}
export type AttributeValues = { [k: string]: ValueCell<ArrayType> }
export type AttributeBuffers = { [k: string]: AttributeBuffer }

export interface AttributeBuffer extends Buffer {
    bind: (location: number) => void
}

export function createAttributeBuffer<T extends ArrayType, S extends BufferItemSize>(ctx: WebGLContext, array: ArrayType, itemSize: S, divisor: number, usageHint: UsageHint = 'dynamic'): AttributeBuffer {
    const { gl } = ctx
    const { instancedArrays } = ctx.extensions

    const buffer = createBuffer(ctx, array, itemSize, usageHint, 'attribute')
    const { _buffer, _bufferType, _dataType, _bpe } = buffer

    return {
        ...buffer,
        bind: (location: number) => {
            gl.bindBuffer(_bufferType, _buffer)
            if (itemSize === 16) {
                for (let i = 0; i < 4; ++i) {
                    gl.enableVertexAttribArray(location + i)
                    gl.vertexAttribPointer(location + i, 4, _dataType, false, 4 * 4 * _bpe, i * 4 * _bpe)
                    instancedArrays.vertexAttribDivisor(location + i, divisor)
                }
            } else {
                gl.enableVertexAttribArray(location)
                gl.vertexAttribPointer(location, itemSize, _dataType, false, 0, 0)
                instancedArrays.vertexAttribDivisor(location, divisor)
            }
        }
    }
}

export function createAttributeBuffers(ctx: WebGLContext, schema: RenderableSchema, values: AttributeValues) {
    const buffers: AttributeBuffers = {}
    Object.keys(schema).forEach(k => {
        const spec = schema[k]
        if (spec.type === 'attribute') {
            buffers[k] = createAttributeBuffer(ctx, values[k].ref.value, spec.itemSize, spec.divisor)
        }
    })
    return buffers as AttributeBuffers
}

export type ElementsType = Uint16Array | Uint32Array
export type ElementsKind = 'uint16' | 'uint32'

export interface ElementsBuffer extends Buffer {
    bind: () => void
}

export function createElementsBuffer(ctx: WebGLContext, array: ElementsType, usageHint: UsageHint = 'static'): ElementsBuffer {
    const { gl } = ctx
    const buffer = createBuffer(ctx, array, 1, usageHint, 'elements')
    const { _buffer } = buffer

    return {
        ...buffer,
        bind: () => {
            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, _buffer);
        }
    }
}