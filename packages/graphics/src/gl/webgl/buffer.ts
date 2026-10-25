/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { WebGLContext } from './context.js';
import { type DefineLookup, isUsed, type RenderableSchema } from '../renderable/schema.js';
import { idFactory } from '@molstar/core/util/id-factory';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { WebGLExtensions } from './extensions.js';
import type { WebGLState } from './state.js';
import { getFormat, getType } from './texture.js';
import { getBytesPerElement, type TextureFormat, type TextureType } from '../texture.js';
import type {
  ArrayType,
  AttributeItemSize,
  AttributeKind,
  AttributeValues,
  DataType,
  ElementsType,
  UsageHint,
} from '../buffer.js';
import type { AttributeLayout } from '../attribute-layout.js';
import { SegmentedArray } from '../segmented-array.js';
import { isDebugMode } from '@molstar/core/util/debug';

const getNextBufferId = idFactory();

export type BufferType = 'attribute' | 'elements' | 'uniform' | 'pixel-pack';

export function getUsageHint(gl: GLRenderingContext, usageHint: UsageHint) {
  switch (usageHint) {
    case 'static':
      return gl.STATIC_DRAW;
    case 'dynamic':
      return gl.DYNAMIC_DRAW;
    case 'stream':
      return gl.STREAM_DRAW;
    default:
      assertUnreachable(usageHint);
  }
}

export function getDataType(gl: GLRenderingContext, dataType: DataType) {
  switch (dataType) {
    case 'uint8':
      return gl.UNSIGNED_BYTE;
    case 'int8':
      return gl.BYTE;
    case 'uint16':
      return gl.UNSIGNED_SHORT;
    case 'int16':
      return gl.SHORT;
    case 'uint32':
      return gl.UNSIGNED_INT;
    case 'int32':
      return gl.INT;
    case 'float32':
      return gl.FLOAT;
    default:
      assertUnreachable(dataType);
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
  }
  assertUnreachable(array);
}

export function getBufferType(gl: GLRenderingContext, bufferType: BufferType) {
  switch (bufferType) {
    case 'attribute':
      return gl.ARRAY_BUFFER;
    case 'elements':
      return gl.ELEMENT_ARRAY_BUFFER;
    case 'uniform':
      if (isWebGL2(gl)) {
        return gl.UNIFORM_BUFFER;
      } else {
        throw new Error('WebGL2 is required for uniform buffers');
      }
    case 'pixel-pack':
      if (isWebGL2(gl)) {
        return gl.PIXEL_PACK_BUFFER;
      } else {
        throw new Error('WebGL2 is required for pixel-pack buffers');
      }
  }
}

export interface Buffer {
  readonly id: number;

  readonly _usageHint: number;
  readonly _bufferType: number;
  readonly _dataType: number;
  readonly _bpe: number;

  readonly length: number;

  getByteCount: () => number;

  getBuffer: () => WebGLBuffer;
  updateData: (array: ArrayType | SegmentedArray) => void;
  /** of a segmented array, the parts changed since it was written (`offset` and `count` unused) */
  updateSubData: (array: ArrayType | SegmentedArray, offset: number, count: number) => void;

  destroy: () => void;
}

export function getBuffer(gl: GLRenderingContext) {
  const buffer = gl.createBuffer();
  if (buffer === null) {
    throw new Error('Could not create WebGL buffer');
  }
  return buffer;
}

function createBuffer(
  gl: GLRenderingContext,
  array: ArrayType | SegmentedArray,
  usageHint: UsageHint,
  bufferType: BufferType,
): Buffer {
  const _buffer = getBuffer(gl);

  const _usageHint = getUsageHint(gl, usageHint);
  const _bufferType = getBufferType(gl, bufferType);
  // segmented arrays: of merged render objects, float attributes and uint32 elements
  const _dataType = SegmentedArray.is(array)
    ? bufferType === 'elements'
      ? gl.UNSIGNED_INT
      : gl.FLOAT
    : dataTypeFromArray(gl, array);
  const _bpe = array.BYTES_PER_ELEMENT;
  const _length = array.length;
  /** the parts of a segmented array written */
  let segmented: SegmentedArray.Uploaded | undefined;

  const writePart = (data: ArrayType, offset: number) => gl.bufferSubData(_bufferType, offset * _bpe, data);

  function updateData(array: ArrayType | SegmentedArray) {
    gl.bindBuffer(_bufferType, _buffer);
    if (SegmentedArray.is(array)) {
      gl.bufferData(_bufferType, array.byteLength, _usageHint);
      segmented = SegmentedArray.writeChanged(array, undefined, writePart);
    } else {
      segmented = undefined;
      gl.bufferData(_bufferType, array, _usageHint);
    }
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

    getByteCount: () => {
      return _bpe * _length;
    },

    getBuffer: () => _buffer,

    updateData,
    updateSubData: (array: ArrayType | SegmentedArray, offset: number, count: number) => {
      gl.bindBuffer(_bufferType, _buffer);
      if (SegmentedArray.is(array)) {
        segmented = SegmentedArray.writeChanged(array, segmented, writePart);
      } else if (count - offset === array.length) {
        gl.bufferSubData(_bufferType, 0, array);
      } else {
        gl.bufferSubData(_bufferType, offset * _bpe, array.subarray(offset, offset + count));
      }
    },

    destroy: () => {
      if (destroyed) return;
      gl.deleteBuffer(_buffer);
      destroyed = true;
    },
  };
}

//

export function getAttribType(gl: GLRenderingContext, kind: AttributeKind, itemSize: AttributeItemSize) {
  switch (kind) {
    case 'float32':
      switch (itemSize) {
        case 1:
          return gl.FLOAT;
        case 2:
          return gl.FLOAT_VEC2;
        case 3:
          return gl.FLOAT_VEC3;
        case 4:
          return gl.FLOAT_VEC4;
        case 16:
          return gl.FLOAT_MAT4;
      }
    default:
      assertUnreachable(kind);
  }
}

export type AttributeBuffers = [string, AttributeBuffer][];

export interface AttributeBuffer extends Buffer {
  readonly itemSize: AttributeItemSize;
  readonly divisor: number;
  bind: (location: number) => void;
  changeOffset: (location: number, offset: number) => void;
}

export function createAttributeBuffer<T extends ArrayType | SegmentedArray, S extends AttributeItemSize>(
  gl: GLRenderingContext,
  state: WebGLState,
  extensions: WebGLExtensions,
  array: T,
  itemSize: S,
  divisor: number,
  usageHint: UsageHint = 'static',
): AttributeBuffer {
  const { instancedArrays } = extensions;

  const buffer = createBuffer(gl, array, usageHint, 'attribute');
  const { _bufferType, _dataType, _bpe } = buffer;

  return {
    ...buffer,
    itemSize,
    divisor,
    bind: (location: number) => {
      gl.bindBuffer(_bufferType, buffer.getBuffer());
      if (itemSize === 16) {
        for (let i = 0; i < 4; ++i) {
          state.enableVertexAttrib(location + i);
          gl.vertexAttribPointer(location + i, 4, _dataType, false, 4 * 4 * _bpe, i * 4 * _bpe);
          instancedArrays.vertexAttribDivisor(location + i, divisor);
        }
      } else {
        state.enableVertexAttrib(location);
        gl.vertexAttribPointer(location, itemSize, _dataType, false, 0, 0);
        instancedArrays.vertexAttribDivisor(location, divisor);
      }
    },
    changeOffset: (location: number, offset: number) => {
      const o = offset * _bpe * itemSize;
      gl.bindBuffer(_bufferType, buffer.getBuffer());
      if (itemSize === 16) {
        for (let i = 0; i < 4; ++i) {
          gl.vertexAttribPointer(location + i, 4, _dataType, false, 4 * 4 * _bpe, i * 4 * _bpe + o);
        }
      } else {
        gl.vertexAttribPointer(location, itemSize, _dataType, false, 0, o);
      }
    },
  };
}

/**
 * Binds each buffer at its fixed location, whether a program uses it or not,
 * so the state (or VAO) is the same for all programs of a render item.
 */
export function bindAttributeBuffers(state: WebGLState, layout: AttributeLayout, attributeBuffers: AttributeBuffers) {
  state.clearVertexAttribsState();
  for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
    const [k, buffer] = attributeBuffers[i];
    const entry = layout.byName[k];
    if (isDebugMode && (buffer.itemSize !== entry?.itemSize || buffer.divisor !== entry.divisor)) {
      throw new Error(
        `attribute '${k}' buffer (itemSize ${buffer.itemSize}, divisor ${buffer.divisor}) does not match the layout (itemSize ${entry?.itemSize}, divisor ${entry?.divisor})`,
      );
    }
    buffer.bind(entry.location);
  }
  state.disableUnusedVertexAttribs();
}

export function createAttributeBuffers(
  ctx: WebGLContext,
  schema: RenderableSchema,
  values: AttributeValues,
  defines: DefineLookup,
) {
  const buffers: AttributeBuffers = [];
  Object.keys(schema).forEach((k) => {
    const spec = schema[k];
    if (spec.type === 'attribute' && isUsed(spec, defines)) {
      buffers[buffers.length] = [k, ctx.resources.attribute(values[k].ref.value, spec.itemSize, spec.divisor)];
    }
  });
  return buffers;
}

//

export interface ElementsBuffer extends Buffer {
  bind: () => void;
}

export function createElementsBuffer(
  gl: GLRenderingContext,
  array: ElementsType | SegmentedArray,
  usageHint: UsageHint = 'static',
): ElementsBuffer {
  const buffer = createBuffer(gl, array, usageHint, 'elements');

  return {
    ...buffer,
    bind: () => {
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buffer.getBuffer());
    },
  };
}

//

export interface PixelPackBuffer {
  readonly id: number;

  readonly _type: number;
  readonly _format: number;
  readonly _bpe: number;

  getByteCount: () => number;

  read: (x: number, y: number, width: number, height: number) => void;
  getSubData: (array: ArrayType) => void;

  destroy: () => void;
}

export function createPixelPackBuffer(
  gl: WebGL2RenderingContext,
  extensions: WebGLExtensions,
  format: TextureFormat,
  type: TextureType,
): PixelPackBuffer {
  const _buffer = getBuffer(gl);

  const _type = getType(gl, extensions, type);
  const _format = getFormat(gl, format, type);
  const _bpe = getBytesPerElement(format, type);

  let _width = 0;
  let _height = 0;

  function read(x: number, y: number, width: number, height: number) {
    _width = width;
    _height = height;
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, _buffer);
    gl.bufferData(gl.PIXEL_PACK_BUFFER, width * height * _bpe, gl.STREAM_READ);
    gl.readPixels(x, y, width, height, _format, _type, 0);
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
  }

  function getSubData(array: ArrayType) {
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, _buffer);
    gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, array);
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
  }

  let destroyed = false;

  return {
    id: getNextBufferId(),

    _type,
    _format,
    _bpe,

    getByteCount: () => {
      return _bpe * _width * _height;
    },

    read,
    getSubData,

    destroy: () => {
      if (destroyed) return;
      gl.deleteBuffer(_buffer);
      destroyed = true;
    },
  };
}
