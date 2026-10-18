/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat3, Mat4, Quat, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { getWindowRectCell } from '../window-rect.js';
import { ValueCell } from '@molstar/core/util';
import { GlobalTextureSchema, GlobalUniformSchema, TextureSpec, UniformSpec } from '../renderable/schema.js';
import { stripUniformDeclarations } from '../webgl/shader-code.js';
import {
  createUniformBlockData,
  createUniformBlockWriter,
  getBindingLayout,
  getUniformBlockLayout,
} from '../uniform-layout.js';
import type { WebGLState } from '../webgl/state.js';
import type { UniformsList } from '../uniform.js';
import { createUniformArena, createUniformBuffer, createUniformRingBuffer } from '../webgl/uniform-buffer.js';
import { getShaderUniformDeclarations, getUniformBlockDeclaration } from '../webgl/shader-declarations.js';

const TestSchema = {
  a: UniformSpec('f', undefined, 'object'),
  b: UniformSpec('v3', undefined, 'object'),
  c: UniformSpec('f', undefined, 'object'),
  d: UniformSpec('m3', undefined, 'object'),
  e: UniformSpec('v2', undefined, 'object'),
  g: UniformSpec('iv3', undefined, 'object'),
  h: UniformSpec('b', undefined, 'object'),
  i: UniformSpec('m4', undefined, 'object'),
  j: UniformSpec('q', undefined, 'object'),
  classic: UniformSpec('f', undefined, 'draw'),
  material: UniformSpec('f', undefined, 'material'),
};

const ArraySchema = {
  s: UniformSpec('f'),
  arrF: UniformSpec('f[]', undefined, undefined, 3),
  arrV3: UniformSpec('v3[]', undefined, undefined, 2),
  arrM3: UniformSpec('m3[]', undefined, undefined, 2),
  after: UniformSpec('f'),
  arrB: UniformSpec('b[]', undefined, undefined, 2),
};

/** simulates GPU buffer contents and uniform buffer range bindings */
function createMockGl() {
  const buffers = new Map<object, Uint8Array>();
  let current: object | undefined;
  const counts = { uploads: 0 };
  const gl = {
    UNIFORM_BUFFER: 0x8a11,
    DYNAMIC_DRAW: 0x88e8,
    UNIFORM_BUFFER_OFFSET_ALIGNMENT: 0x8a34,
    getParameter: () => 256,
    createBuffer: () => ({}),
    bindBuffer: (_target: number, buffer: object) => {
      current = buffer;
    },
    bufferData: (_target: number, size: number) => {
      buffers.set(current!, new Uint8Array(size));
    },
    bufferSubData: (
      _target: number,
      dstOffset: number,
      src: Uint8Array,
      srcOffset = 0,
      length = src.byteLength - srcOffset,
    ) => {
      buffers.get(current!)!.set(new Uint8Array(src.buffer, src.byteOffset + srcOffset, length), dstOffset);
      counts.uploads += 1;
    },
    deleteBuffer: (buffer: object) => {
      buffers.delete(buffer);
    },
  } as unknown as WebGL2RenderingContext;

  const bound: { buffer: object; offset: number; size: number }[] = [];
  const state = {
    bindUniformBufferRange: (index: number, _id: number, buffer: object, offset: number, size: number) => {
      bound[index] = { buffer, offset, size };
    },
  } as unknown as WebGLState;

  /** GPU data bound at `index` */
  function read(index: number) {
    const { buffer, offset, size } = bound[index];
    return new Float32Array(buffers.get(buffer)!.buffer, offset, size / 4);
  }

  return { gl, state, arena: createUniformArena(gl, state), counts, bound, read };
}

const GlobalIndex = 0;
const ObjectIndex = 2;

describe('uniform-buffer', () => {
  it('std140 layout', () => {
    const layout = getUniformBlockLayout('object', TestSchema)!;
    expect(layout.members.map((m) => [m.name, m.offset])).toEqual([
      ['a', 0],
      ['b', 16], // vec3 aligned to 16
      ['c', 28], // scalar packed after vec3
      ['d', 32],
      ['e', 80], // mat3 is 3 x vec4
      ['g', 96],
      ['h', 108],
      ['i', 112],
      ['j', 176],
    ]);
    expect(layout.size).toBe(192);
  });

  it('std140 array layout', () => {
    const layout = getUniformBlockLayout('object', ArraySchema)!;
    expect(layout.members.map((m) => [m.name, m.offset, m.stride])).toEqual([
      ['s', 0, 0],
      ['arrF', 16, 16], // elements padded to vec4
      ['arrV3', 64, 16],
      ['arrM3', 96, 48],
      ['after', 192, 0],
      ['arrB', 208, 16],
    ]);
    expect(layout.size).toBe(240);
  });

  it('layout per block', () => {
    expect(getUniformBlockLayout('material', TestSchema)!.members.map((m) => m.name)).toEqual(['material']);
    expect(getUniformBlockLayout('global', TestSchema)).toBeUndefined();
    expect(getUniformBlockLayout('object', { classic: UniformSpec('f', undefined, 'draw') })).toBeUndefined();
    expect(getShaderUniformDeclarations({ arr: UniformSpec('f[]') }, true)).toBeUndefined();
    expect(
      getUniformBlockLayout('object', { ...ArraySchema, arrF: UniformSpec('f[]', undefined, undefined, 4) })!.hash,
    ).not.toBe(getUniformBlockLayout('object', ArraySchema)!.hash);
  });

  it('global layout', () => {
    const layout = getUniformBlockLayout('global', GlobalUniformSchema)!;
    // any change here changes the 'Global' block of all programs
    expect(layout.members.map((m) => [m.name, m.offset])).toEqual([
      ['uModel', 0],
      ['uView', 64],
      ['uInvView', 128],
      ['uModelView', 192],
      ['uInvModelView', 256],
      ['uProjection', 320],
      ['uInvProjection', 384],
      ['uModelViewProjection', 448],
      ['uInvModelViewProjection', 512],
      ['uHasHeadRotation', 576],
      ['uInvHeadRotation', 592],
      ['uIsAsymmetricProjection', 656],
      ['uHasEyeCamera', 660],
      ['uModelViewEye', 672],
      ['uInvModelViewEye', 736],
      ['uIsOrtho', 800],
      ['uPixelRatio', 804],
      ['uViewport', 816],
      ['uViewOffset', 832],
      ['uModelScale', 840],
      ['uDrawingBufferSize', 848],
      ['uCameraPosition', 864],
      ['uCameraDir', 880],
      ['uCameraPlane', 896],
      ['uNear', 912],
      ['uFar', 916],
      ['uFog', 920],
      ['uFogNear', 924],
      ['uFogFar', 928],
      ['uFogColor', 944],
      ['uTransparentBackground', 956],
      ['uLightDirection', 960],
      ['uLightColor', 1216],
      ['uAmbientColor', 1472],
      ['uPickingAlphaThreshold', 1484],
      ['uHighlightColor', 1488],
      ['uSelectColor', 1504],
      ['uDimColor', 1520],
      ['uHighlightStrength', 1532],
      ['uSelectStrength', 1536],
      ['uDimStrength', 1540],
      ['uMarkerPriority', 1544],
      ['uMarkerAverage', 1548],
      ['uXrayEdgeFalloff', 1552],
      ['uCelSteps', 1556],
      ['uExposure', 1560],
      ['uRenderMask', 1564],
      ['uMarkingDepthTest', 1568],
      ['uDepthBack', 1572],
      ['uMarkingType', 1576],
      ['uPickType', 1580],
      ['uTime', 1584],
      ['uEnableAnimation', 1588],
    ]);
    expect(layout.size).toBe(1600);

    // any schema that includes global uniforms gets the full global layout
    const schema = { ...GlobalUniformSchema, x: UniformSpec('f', undefined, 'object') };
    expect(getUniformBlockLayout('global', schema)).toEqual(layout);
    expect(getUniformBlockLayout('global', { uTime: GlobalUniformSchema.uTime })).toEqual(layout);
  });

  it('declaration', () => {
    const layout = getUniformBlockLayout('material', {
      uAlpha: UniformSpec('f', undefined, 'material'),
      uColor: UniformSpec('v3', undefined, 'material'),
    })!;
    expect(getUniformBlockDeclaration(layout)).toBe(
      ['layout(std140) uniform Material {', '    float uAlpha;', '    vec3 uColor;', '};'].join('\n'),
    );
    expect(getUniformBlockDeclaration(getUniformBlockLayout('object', ArraySchema)!)).toBe(
      [
        'layout(std140) uniform Object {',
        '    float s;',
        '    float arrF[3];',
        '    vec3 arrV3[2];',
        '    mat3 arrM3[2];',
        '    float after;',
        '    bool arrB[2];',
        '};',
      ].join('\n'),
    );
  });

  it('spec validation', () => {
    expect(UniformSpec('f', undefined, 'material').variant).toBe('material');
    expect(UniformSpec('f').block).toBe('object');
    expect(UniformSpec('f', 'material').block).toBe('material');
    expect(UniformSpec('f', 'buffered').block).toBe('object');
    expect(UniformSpec('f[]').block).toBeUndefined();
    expect(UniformSpec('t').block).toBeUndefined();
    expect(() => UniformSpec('v3[]', undefined, 'object')).toThrow();
    expect(() => UniformSpec('f', undefined, undefined, 2)).toThrow();
    expect(() => UniformSpec('t', undefined, 'object')).toThrow();
    expect(() => UniformSpec('f', 'buffered', 'material')).toThrow();
    expect(() => UniformSpec('f', 'material', 'object')).toThrow();
    expect(() => UniformSpec('f', 'material', 'global')).toThrow();
  });

  it('strip declarations', () => {
    const source = [
      'uniform float uA;',
      'uniform mat4 uB, uC, uD;',
      '  uniform highp vec3 uE[2];',
      'uniform sampler2D uA;',
      'uniform vec2 uKeep;',
      '// uniform float uA;',
    ].join('\n');
    const members = new Map([
      ['uA', 'float'],
      ['uC', 'mat4'],
      ['uE', 'vec3'],
    ]);
    expect(stripUniformDeclarations(source, members)).toBe(
      ['', 'uniform mat4 uB, uD;', '', 'uniform sampler2D uA;', 'uniform vec2 uKeep;', '// uniform float uA;'].join(
        '\n',
      ),
    );
    expect(() => stripUniformDeclarations('uniform int uA;', members)).toThrow();
  });

  it('generated declarations', () => {
    const schema = {
      uA: UniformSpec('f'),
      uM: UniformSpec('v3', 'material'),
      uD: UniformSpec('i', undefined, 'draw'),
      uArr: UniformSpec('v3[]', undefined, undefined, 2),
      uDef: UniformSpec('f[]', undefined, undefined, 4, 'dN'),
      uFree: UniformSpec('f[]'),
    };
    const withBlocks = getShaderUniformDeclarations(schema, true)!;
    expect(withBlocks.code).toContain('layout(std140) uniform Material {\n    vec3 uM;\n};');
    expect(withBlocks.code).toContain(
      'layout(std140) uniform Object {\n    float uA;\n    vec3 uArr[2];\n    float uDef[4];\n};',
    );
    expect(withBlocks.code).toContain('uniform int uD;');
    expect(withBlocks.code).not.toContain('uFree');
    expect([...withBlocks.types.keys()].sort()).toEqual(['uA', 'uArr', 'uD', 'uDef', 'uM']);

    const classic = getShaderUniformDeclarations(schema, false)!;
    expect(classic.code).not.toContain('layout(std140)');
    expect(classic.code).toContain(
      'uniform float uA;\nuniform vec3 uM;\nuniform int uD;\nuniform vec3 uArr[2];\n#if dN > 0\nuniform float uDef[dN];\n#endif',
    );
    // arrays without `length` or `lengthDefine` stay in the shader code
    expect(classic.code).not.toContain('uFree');
    expect([...classic.types.keys()]).toEqual(['uA', 'uM', 'uD', 'uArr', 'uDef']);
    expect(getShaderUniformDeclarations(schema, false)).toBe(classic);
    expect(() => UniformSpec('f', undefined, undefined, undefined, 'dN')).toThrow();
  });

  it('packing', () => {
    const { arena, counts, read } = createMockGl();
    const layout = getUniformBlockLayout('object', TestSchema)!;
    const buffer = createUniformBuffer(arena, layout);

    const values = {
      a: ValueCell.create(1.5),
      b: ValueCell.create(Vec3.create(1, 2, 3)),
      c: ValueCell.create(4),
      d: ValueCell.create(Mat3.fromArray(Mat3(), [1, 2, 3, 4, 5, 6, 7, 8, 9], 0)),
      e: ValueCell.create(Vec2.create(5, 6)),
      g: ValueCell.create(Vec3.create(-1, 2, 3)),
      h: ValueCell.create(true),
      i: ValueCell.create(Mat4.identity()),
      j: ValueCell.create(Quat.create(1, 2, 3, 4)),
      classic: ValueCell.create(99),
    };
    const write = buffer.createWriter(Object.entries(values) as UniformsList);
    write();
    buffer.sync();
    buffer.bind();
    expect(counts.uploads).toBe(1);

    const f32 = read(ObjectIndex);
    const i32 = new Int32Array(f32.buffer, f32.byteOffset, f32.length);
    expect(f32.length).toBeGreaterThanOrEqual(layout.size / 4);
    expect(f32[0]).toBe(1.5);
    expect(Array.from(f32.subarray(4, 8))).toEqual([1, 2, 3, 4]);
    expect(Array.from(f32.subarray(8, 20))).toEqual([1, 2, 3, 0, 4, 5, 6, 0, 7, 8, 9, 0]);
    expect(Array.from(f32.subarray(20, 22))).toEqual([5, 6]);
    expect(Array.from(i32.subarray(24, 28))).toEqual([-1, 2, 3, 1]);
    expect(Array.from(f32.subarray(28, 44))).toEqual(Array.from(Mat4.identity()));
    expect(Array.from(f32.subarray(44, 48))).toEqual([1, 2, 3, 4]);
    // not a block member, ignored
    expect(Array.from(f32)).not.toContain(99);

    // no upload without changes
    write();
    buffer.sync();
    expect(counts.uploads).toBe(1);

    ValueCell.update(values.h, false);
    write();
    buffer.sync();
    expect(counts.uploads).toBe(2);
    expect(i32[27]).toBe(0);
  });

  it('array packing', () => {
    const { arena, read } = createMockGl();
    const layout = getUniformBlockLayout('object', ArraySchema)!;
    const buffer = createUniformBuffer(arena, layout);

    const values = {
      arrF: ValueCell.create([1, 2, 3, 4]), // longer than length
      arrV3: ValueCell.create([1, 2, 3, 4, 5, 6]),
      arrM3: ValueCell.create([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19]),
      arrB: ValueCell.create([true]), // shorter than length
    };
    buffer.createWriter(Object.entries(values) as UniformsList)();
    buffer.sync();
    buffer.bind();

    const f32 = read(ObjectIndex);
    const i32 = new Int32Array(f32.buffer, f32.byteOffset, f32.length);
    expect(Array.from(f32.subarray(4, 16))).toEqual([1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0]);
    expect(Array.from(f32.subarray(16, 24))).toEqual([1, 2, 3, 0, 4, 5, 6, 0]);
    expect(Array.from(f32.subarray(24, 48))).toEqual([
      1, 2, 3, 0, 4, 5, 6, 0, 7, 8, 9, 0, 11, 12, 13, 0, 14, 15, 16, 0, 17, 18, 19, 0,
    ]);
    expect(Array.from(i32.subarray(52, 60))).toEqual([1, 0, 0, 0, 0, 0, 0, 0]);
  });

  it('arena sub-allocation', () => {
    const { arena, counts, bound, read } = createMockGl();
    const layout = getUniformBlockLayout('object', TestSchema)!;
    const values = { a: ValueCell.create(1) };

    const a = createUniformBuffer(arena, layout);
    const b = createUniformBuffer(arena, layout);
    a.bind();
    const aBinding = bound[ObjectIndex];
    b.bind();
    const bBinding = bound[ObjectIndex];
    expect(aBinding.buffer).toBe(bBinding.buffer);
    expect(aBinding.offset % 256).toBe(0);
    expect(bBinding.offset % 256).toBe(0);
    expect(aBinding.offset).not.toBe(bBinding.offset);
    expect(aBinding.size).toBeGreaterThanOrEqual(layout.size);

    // pending writes of both are coalesced into one upload
    const uploads = counts.uploads;
    a.createWriter(Object.entries(values) as UniformsList)();
    b.createWriter(Object.entries(values) as UniformsList)();
    a.sync();
    b.sync();
    expect(counts.uploads).toBe(uploads + 1);

    // freed slices are reused and zeroed
    a.destroy();
    const c = createUniformBuffer(arena, layout);
    c.sync();
    c.bind();
    expect(bound[ObjectIndex].offset).toBe(aBinding.offset);
    expect(read(ObjectIndex)[0]).toBe(0);

    expect(arena.getByteCount()).toBe(64 * 1024);
  });

  it('ring buffer', () => {
    const { arena, bound, read } = createMockGl();
    const layout = getUniformBlockLayout('global', GlobalUniformSchema)!;
    const buffer = createUniformRingBuffer(arena, layout, 2);
    const uTime = ValueCell.create(1);
    const write = buffer.createWriter([['uTime', uTime]]);
    const timeIndex = layout.members.find((m) => m.name === 'uTime')!.offset / 4;

    write();
    buffer.sync();
    buffer.bind();
    const first = bound[GlobalIndex];
    expect(read(GlobalIndex)[timeIndex]).toBe(1);

    // unchanged, same slot
    write();
    buffer.sync();
    buffer.bind();
    expect(bound[GlobalIndex]).toEqual(first);

    // changed, next slot, previous snapshot is kept
    ValueCell.update(uTime, 2);
    write();
    buffer.sync();
    buffer.bind();
    const second = bound[GlobalIndex];
    expect(second.offset).toBe(first.size);
    expect(read(GlobalIndex)[timeIndex]).toBe(2);
    bound[GlobalIndex] = first;
    expect(read(GlobalIndex)[timeIndex]).toBe(1);

    // wraps around
    ValueCell.update(uTime, 3);
    write();
    buffer.sync();
    buffer.bind();
    expect(bound[GlobalIndex].offset).toBe(0);
  });

  it('binding layout', () => {
    const schema = {
      ...GlobalUniformSchema,
      tDepth: GlobalTextureSchema.tDepth,
      tMaterial: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest', 'material'),
      tColor: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest'),
      uAlpha: UniformSpec('f', 'material'),
      uObjectId: UniformSpec('i'),
    };
    const layout = getBindingLayout(schema);
    const names = (tier: keyof typeof layout) => layout[tier].textures.map((t) => t.name);
    expect(names('global')).toEqual(['tDepth']);
    // global by spec, not by name: a compute schema's own `tDepth` is per object
    const compute = getBindingLayout({ tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'linear') });
    expect(compute.global.textures).toEqual([]);
    expect(compute.object.textures.map((t) => t.name)).toEqual(['tDepth']);
    expect(layout.global.uniforms).toEqual(getUniformBlockLayout('global', GlobalUniformSchema));
    expect(names('material')).toEqual(['tMaterial']);
    expect(layout.material.uniforms!.members.map((m) => m.name)).toEqual(['uAlpha']);
    expect(names('object')).toEqual(['tColor']);
    expect(layout.object.uniforms!.members.map((m) => m.name)).toEqual(['uObjectId']);
  });
});

describe('uniform block writer', () => {
  it('checks the cells only after one changed, until disposed', () => {
    const layout = getUniformBlockLayout('global', GlobalUniformSchema)!;
    const data = createUniformBlockData(new ArrayBuffer(layout.size), 0, layout.size);
    const timeIndex = layout.members.find((m) => m.name === 'uTime')!.offset / 4;
    const uTime = ValueCell.create(1);
    let writes = 0;
    const write = createUniformBlockWriter(layout, data, [['uTime', uTime]], () => writes++);

    write();
    expect(writes).toBe(1);
    expect(data.f32[timeIndex]).toBe(1);

    // nothing changed: no write
    write();
    expect(writes).toBe(1);

    ValueCell.update(uTime, 2);
    write();
    expect(writes).toBe(2);
    expect(data.f32[timeIndex]).toBe(2);

    // the same value again: the cell changed, the data does not need writing but is checked
    ValueCell.update(uTime, 2);
    write();
    expect(data.f32[timeIndex]).toBe(2);

    write.dispose();
    expect(uTime.observers?.size).toBe(0);
  });

  it('writes flipped window rectangles when their source or the target height changes', () => {
    const schema = { uViewport: { ...UniformSpec('v4'), windowRect: 'pixels' as const } };
    const layout = getUniformBlockLayout('object', schema)!;
    const data = createUniformBlockData(new ArrayBuffer(layout.size), 0, layout.size);
    const index = layout.members.find((m) => m.name === 'uViewport')!.offset / 4;
    const viewport = () => Array.from(data.f32.subarray(index, index + 4));
    const source = ValueCell.create(Vec4.create(0, 0, 10, 20));
    let height = 100;
    const flipped = getWindowRectCell(schema, 'uViewport', source, () => height);
    const write = createUniformBlockWriter(layout, data, [['uViewport', flipped]], () => {});

    write();
    expect(viewport()).toEqual([0, 80, 10, 20]);

    // the source is updated, the writer observes the flipped cell
    ValueCell.update(source, Vec4.create(5, 0, 10, 20));
    write();
    expect(viewport()).toEqual([5, 80, 10, 20]);

    // another target, without a notification
    height = 50;
    write();
    expect(viewport()).toEqual([5, 30, 10, 20]);

    write.dispose();
    expect(source.observers?.size).toBe(0);
  });
});
