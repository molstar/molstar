/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'

export type UsageHint = 'static' | 'dynamic' | 'stream'
export type DataType = 'uint8' | 'int8' | 'uint16' | 'int16' | 'uint32' | 'int32' | 'float32'
export type BufferType = 'attribute' | 'element'

export type DataTypeArrayType = {
    'uint8': Uint8Array
    'int8': Int8Array
    'uint16': Uint16Array
    'int16': Int16Array
    'uint32': Uint32Array
    'int32': Int32Array
    'float32': Float32Array
}
export type ArrayType = Helpers.ValueOf<DataTypeArrayType>
export type ArrayKind = keyof DataTypeArrayType

export type BufferItemSize = 1 | 2 | 3 | 4

export function getUsageHint(gl: WebGLRenderingContext, usageHint: UsageHint) {
    switch (usageHint) {
        case 'static': return gl.STATIC_DRAW
        case 'dynamic': return gl.DYNAMIC_DRAW
        case 'stream': return gl.STREAM_DRAW
    }
}

export function getDataType(gl: WebGLRenderingContext, dataType: DataType) {
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

function dataTypeFromArray(gl: WebGLRenderingContext, array: ArrayType) {
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

export function getBufferType(gl: WebGLRenderingContext, bufferType: BufferType) {
    switch (bufferType) {
        case 'attribute': return gl.ARRAY_BUFFER
        case 'element': return gl.ELEMENT_ARRAY_BUFFER
    }
}

export interface Buffer<T extends ArrayType, S extends BufferItemSize, B extends BufferType> {
    updateData: (array: T) => void
    updateSubData: (array: T, offset: number, count: number) => void
    bind: (location: number, stride: number, offset: number) => void
    destroy: () => void
}

export function createBuffer<T extends ArrayType, S extends BufferItemSize, B extends BufferType>(ctx: Context, array: T, itemSize: S, usageHint: UsageHint, bufferType: B): Buffer<T, S, B> {
    const { gl } = ctx
    const buffer = gl.createBuffer()
    if (buffer === null) {
        throw new Error('Could not create WebGL buffer')
    }

    const _usageHint = getUsageHint(gl, usageHint)
    const _bufferType = getBufferType(gl, bufferType)
    const _dataType = dataTypeFromArray(gl, array)

    function updateData(array: T) {
        gl.bindBuffer(_bufferType, buffer)
        gl.bufferData(_bufferType, array, _usageHint)
    }
    updateData(array)

    return {
        updateData,
        updateSubData: (array: T, offset: number, count: number) => {
            gl.bindBuffer(_bufferType, buffer)
            gl.bufferSubData(_bufferType, offset * array.BYTES_PER_ELEMENT, array.subarray(offset, offset + count))
        },
        bind: (location: number, stride: number, offset: number) => {
            gl.bindBuffer(_bufferType, buffer);
            gl.enableVertexAttribArray(location);
            gl.vertexAttribPointer(location, itemSize, _dataType, false, stride, offset);
        },
        destroy: () => {
            gl.deleteBuffer(buffer)
        }
    }
}

export type AttributeDefs = { [k: string]: { kind: ArrayKind, itemSize: BufferItemSize, divisor: number } }
export type AttributeValues<T extends AttributeDefs> = { [K in keyof T]: ArrayType }
export type AttributeBuffers<T extends AttributeDefs> = {
    [K in keyof T]: AttributeBuffer<DataTypeArrayType[T[K]['kind']], T[K]['itemSize']>
}

export interface AttributeBuffer<T extends ArrayType, S extends BufferItemSize> extends Buffer<T, S, 'attribute'> {}

export function createAttributeBuffer<T extends ArrayType, S extends BufferItemSize>(ctx: Context, array: T, itemSize: S, divisor: number, usageHint: UsageHint = 'dynamic'): AttributeBuffer<T, S> {
    const buffer = createBuffer(ctx, array, itemSize, usageHint, 'attribute')
    const { angleInstancedArrays } = ctx.extensions

    return {
        ...buffer,
        bind: (location: number, stride: number, offset: number) => {
            buffer.bind(location, stride, offset)
            angleInstancedArrays.vertexAttribDivisorANGLE(location, divisor)
        }
    }
}

export function createAttributeBuffers<T extends AttributeDefs>(ctx: Context, props: T, state: AttributeValues<T>) {
    const buffers: Partial<AttributeBuffers<T>> = {}
    Object.keys(props).forEach(k => {
        buffers[k] = createAttributeBuffer(ctx, state[k], props[k].itemSize, props[k].divisor)
    })
    return buffers as AttributeBuffers<T>
}

export type ElementType = Uint16Array | Uint32Array

export interface ElementBuffer<T extends ElementType> extends Buffer<T, 3, 'element'> {}

export function createElementBuffer<T extends ElementType>(ctx: Context, array: T, usageHint: UsageHint = 'static'): ElementBuffer<T> {
    const buffer = createBuffer(ctx, array, 3, usageHint, 'element')

    return {
        ...buffer
    }
}