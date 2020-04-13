/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from './context';
import { ValueCell } from '../../mol-util';
import { RenderableSchema } from '../renderable/schema';
import { idFactory } from '../../mol-util/id-factory';
import { ValueOf } from '../../mol-util/type-helpers';
import { GLRenderingContext } from './compat';
import { WebGLExtensions } from './extensions';

const getNextBufferId = idFactory();

export type UsageHint = 'static' | 'dynamic' | 'stream'
export type DataType = 'uint8' | 'int8' | 'uint16' | 'int16' | 'uint32' | 'int32' | 'float32'
export type BufferType = 'attribute' | 'elements' | 'uniform'

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

export function getUsageHint(gl: GLRenderingContext, usageHint: UsageHint) {
    switch (usageHint) {
        case 'static': return gl.STATIC_DRAW;
        case 'dynamic': return gl.DYNAMIC_DRAW;
        case 'stream': return gl.STREAM_DRAW;
    }
}

export function getDataType(gl: GLRenderingContext, dataType: DataType) {
    switch (dataType) {
        case 'uint8': return gl.UNSIGNED_BYTE;
        case 'int8': return gl.BYTE;
        case 'uint16': return gl.UNSIGNED_SHORT;
        case 'int16': return gl.SHORT;
        case 'uint32': return gl.UNSIGNED_INT;
        case 'int32': return gl.INT;
        case 'float32': return gl.FLOAT;
    }
}

function dataTypeFromArray(gl: GLRenderingContext, array: ArrayType) {
    if (array instanceof Uint8Array) {
        return gl.UNSIGNED_BYTE;
    } else if (array instanceof Int8Array) {
        return gl.BYTE;
    } else if (array instanceof Uint16Array) {
        return gl.UNSIGNED_SHORT;
    } else if (array instanceof Int16Array) {
        return gl.SHORT;
    } else if (array instanceof Uint32Array) {
        return gl.UNSIGNED_INT;
    } else if (array instanceof Int32Array) {
        return gl.INT;
    } else if (array instanceof Float32Array) {
        return gl.FLOAT;
    } else {
        throw new Error('Should nevver happen');
    }
}

export function getBufferType(gl: GLRenderingContext, bufferType: BufferType) {
    switch (bufferType) {
        case 'attribute': return gl.ARRAY_BUFFER;
        case 'elements': return gl.ELEMENT_ARRAY_BUFFER;
        case 'uniform': return (gl as WebGL2RenderingContext).UNIFORM_BUFFER;
    }
}

export interface Buffer {
    readonly id: number

    readonly _usageHint: number
    readonly _bufferType: number
    readonly _dataType: number
    readonly _bpe: number

    readonly length: number

    getBuffer: () => WebGLBuffer
    updateData: (array: ArrayType) => void
    updateSubData: (array: ArrayType, offset: number, count: number) => void

    reset: () => void
    destroy: () => void
}

function getBuffer(gl: GLRenderingContext) {
    const buffer = gl.createBuffer();
    if (buffer === null) {
        throw new Error('Could not create WebGL buffer');
    }
    return buffer;
}

function createBuffer(gl: GLRenderingContext, array: ArrayType, usageHint: UsageHint, bufferType: BufferType): Buffer {
    let _buffer = getBuffer(gl);

    const _usageHint = getUsageHint(gl, usageHint);
    const _bufferType = getBufferType(gl, bufferType);
    const _dataType = dataTypeFromArray(gl, array);
    const _bpe = array.BYTES_PER_ELEMENT;
    const _length = array.length;

    function updateData(array: ArrayType) {
        gl.bindBuffer(_bufferType, _buffer);
        gl.bufferData(_bufferType, array, _usageHint);
    }
    updateData(array);

    let destroyed = false;

    return {
        id: getNextBufferId(),

        _usageHint,
        _bufferType,
        _dataType,
        _bpe,

        length: _length,
        getBuffer: () => _buffer,

        updateData,
        updateSubData: (array: ArrayType, offset: number, count: number) => {
            gl.bindBuffer(_bufferType, _buffer);
            gl.bufferSubData(_bufferType, offset * _bpe, array.subarray(offset, offset + count));
        },

        reset: () => {
            _buffer = getBuffer(gl);
            updateData(array);
        },
        destroy: () => {
            if (destroyed) return;
            gl.deleteBuffer(_buffer);
            destroyed = true;
        }
    };
}

//

export type AttributeItemSize = 1 | 2 | 3 | 4 | 16
export type AttributeKind = 'float32' | 'int32'

export function getAttribType(gl: GLRenderingContext, kind: AttributeKind, itemSize: AttributeItemSize) {
    switch (kind) {
        case 'int32':
            switch (itemSize) {
                case 1: return gl.INT;
                case 2: return gl.INT_VEC2;
                case 3: return gl.INT_VEC3;
                case 4: return gl.INT_VEC4;
            }
            break;
        case 'float32':
            switch (itemSize) {
                case 1: return gl.FLOAT;
                case 2: return gl.FLOAT_VEC2;
                case 3: return gl.FLOAT_VEC3;
                case 4: return gl.FLOAT_VEC4;
                case 16: return gl.FLOAT_MAT4;
            }
            break;
    }
    throw new Error(`unknown attribute type for kind '${kind}' and itemSize '${itemSize}'`);
}

export type AttributeDefs = {
    [k: string]: { kind: AttributeKind, itemSize: AttributeItemSize, divisor: number }
}
export type AttributeValues = { [k: string]: ValueCell<ArrayType> }
export type AttributeBuffers = [string, AttributeBuffer][]

export interface AttributeBuffer extends Buffer {
    bind: (location: number) => void
}

export function createAttributeBuffer<T extends ArrayType, S extends AttributeItemSize>(gl: GLRenderingContext, extensions: WebGLExtensions, array: T, itemSize: S, divisor: number, usageHint: UsageHint = 'dynamic'): AttributeBuffer {
    const { instancedArrays } = extensions;

    const buffer = createBuffer(gl, array, usageHint, 'attribute');
    const { _bufferType, _dataType, _bpe } = buffer;

    return {
        ...buffer,
        bind: (location: number) => {
            gl.bindBuffer(_bufferType, buffer.getBuffer());
            if (itemSize === 16) {
                for (let i = 0; i < 4; ++i) {
                    gl.enableVertexAttribArray(location + i);
                    gl.vertexAttribPointer(location + i, 4, _dataType, false, 4 * 4 * _bpe, i * 4 * _bpe);
                    instancedArrays.vertexAttribDivisor(location + i, divisor);
                }
            } else {
                gl.enableVertexAttribArray(location);
                gl.vertexAttribPointer(location, itemSize, _dataType, false, 0, 0);
                instancedArrays.vertexAttribDivisor(location, divisor);
            }
        }
    };
}

export function createAttributeBuffers(ctx: WebGLContext, schema: RenderableSchema, values: AttributeValues) {
    const buffers: AttributeBuffers = [];
    Object.keys(schema).forEach(k => {
        const spec = schema[k];
        if (spec.type === 'attribute') {
            buffers[buffers.length] = [k, ctx.resources.attribute(values[k].ref.value, spec.itemSize, spec.divisor)];
        }
    });
    return buffers;
}

//

export type ElementsType = Uint16Array | Uint32Array
export type ElementsKind = 'uint16' | 'uint32'

export interface ElementsBuffer extends Buffer {
    bind: () => void
}

export function createElementsBuffer(gl: GLRenderingContext, array: ElementsType, usageHint: UsageHint = 'static'): ElementsBuffer {
    const buffer = createBuffer(gl, array, usageHint, 'elements');

    return {
        ...buffer,
        bind: () => {
            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buffer.getBuffer());
        }
    };
}