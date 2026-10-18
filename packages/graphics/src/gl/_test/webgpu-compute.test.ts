/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { KernelCode, type KernelSchema, StorageBufferSpec, StorageTextureSpec } from '../compute.js';
import {
  getDispatchSize,
  ScanAddKernel,
  ScanAddSchema,
  ScanBlockKernel,
  ScanBlockSchema,
} from '../compute/prefix-sum.js';
import {
  ColorAccumulateKernel,
  ColorAccumulateSchema,
  ColorNormalizeKernel,
  ColorNormalizeSchema,
} from '../compute/color-smoothing.js';
import {
  GdAccumulateKernel,
  GdAccumulateSchema,
  GdGroupKernel,
  GdGroupSchema,
  GdWriteKernel,
  GdWriteSchema,
} from '../compute/gaussian-density.js';
import {
  McBufferInput,
  McChunkAllocateKernel,
  McChunkAllocateSchema,
  McChunkCountKernel,
  McChunkCountOutput,
  McChunkEmitKernel,
  McChunkEmitOutput,
  McChunkFreeKernel,
  McChunkFreeSchema,
  McChunkPushKernel,
  McChunkPushSchema,
  McTextureInput,
} from '../compute/marching-cubes/compute.js';
import { DefineSpec, TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { DefineValues } from '../shader-code.js';
import { preprocessShader } from '../shader-preprocessor.js';
import { getWgslKernelDeclarations } from '../webgpu/compute.js';
import { getWgslDefines } from '../webgpu/program.js';
import { getWgslKernel } from '../webgpu/shader/kernels.js';

/** WebGPU core limits, the minimum spec */
const Limits = {
  maxComputeInvocationsPerWorkgroup: 256,
  maxComputeWorkgroupSizeX: 256,
  maxComputeWorkgroupSizeY: 256,
  maxComputeWorkgroupSizeZ: 64,
  maxStorageBuffersPerShaderStage: 8,
  maxStorageTexturesPerShaderStage: 4,
};

/** kernels with the define values of their variants */
const BuiltInKernels: [KernelCode, KernelSchema, { [k: string]: string | boolean }][] = [
  [ScanBlockKernel, ScanBlockSchema, {}],
  [ScanAddKernel, ScanAddSchema, {}],
  [McChunkCountKernel, { ...McBufferInput, ...McChunkCountOutput }, { dInput: 'buffer' }],
  [McChunkCountKernel, { ...McTextureInput, ...McChunkCountOutput }, { dInput: 'texture', dChannel: 'alpha' }],
  [McChunkAllocateKernel, McChunkAllocateSchema, {}],
  [McChunkEmitKernel, { ...McBufferInput, ...McChunkEmitOutput }, { dInput: 'buffer' }],
  [
    McChunkEmitKernel,
    { ...McTextureInput, ...McChunkEmitOutput },
    { dInput: 'texture', dChannel: 'alpha', dPackedGroup: true },
  ],
  [McChunkFreeKernel, McChunkFreeSchema, {}],
  [McChunkPushKernel, McChunkPushSchema, {}],
  [ColorAccumulateKernel, ColorAccumulateSchema, { dColorType: 'group' }],
  [ColorAccumulateKernel, ColorAccumulateSchema, { dColorType: 'groupInstance' }],
  [ColorNormalizeKernel, ColorNormalizeSchema, {}],
  [GdAccumulateKernel, GdAccumulateSchema, {}],
  [GdGroupKernel, GdGroupSchema, {}],
  [GdWriteKernel, GdWriteSchema, {}],
];

describe('webgpu compute', () => {
  it('declares the bindings of a kernel schema in group 0', () => {
    const code = KernelCode('test', [8, 4]);
    const schema = {
      uCount: UniformSpec('i'),
      uScale: UniformSpec('v3'),
      dMode: DefineSpec('string', ['a', 'b']),
      input: StorageBufferSpec('vec4<f32>', 'read'),
      counts: StorageBufferSpec('atomic<u32>', 'read-write'),
      tGrid: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, { dimension: '3d' }),
      outTexture: StorageTextureSpec('rgba', 'float', '2d', 'write'),
    };
    const defines = { dMode: ValueCell.create('a') } as DefineValues;
    const code_ = getWgslKernelDeclarations(code, schema, defines);
    expect(code_).toContain('const WorkgroupSizeX = 8u;');
    expect(code_).toContain('const WorkgroupSizeY = 4u;');
    expect(code_).toContain('const WorkgroupSizeZ = 1u;');
    expect(code_).toMatch(/struct KernelParams \{[^}]*uCount: i32,[^}]*uScale: vec3<f32>,[^}]*\};/);
    expect(code_).toContain('@group(0) @binding(0) var<uniform> uParams: KernelParams;');
    expect(code_).toContain('@group(0) @binding(1) var<storage, read> input: array<vec4<f32>>;');
    expect(code_).toContain('@group(0) @binding(2) var<storage, read_write> counts: array<atomic<u32>>;');
    expect(code_).toContain('@group(0) @binding(3) var tGrid: texture_3d<f32>;');
    expect(code_).toContain('@group(0) @binding(4) var outTexture: texture_storage_2d<rgba32float, write>;');
  });

  it('declares no uniform block without uniforms', () => {
    const code_ = getWgslKernelDeclarations(KernelCode('test', [64]), { data: StorageBufferSpec('u32') }, {});
    expect(code_).not.toContain('KernelParams');
    expect(code_).toContain('@group(0) @binding(0) var<storage, read> data: array<u32>;');
  });

  for (const [code, schema, defineValues] of BuiltInKernels) {
    it(`kernel '${code.name}' ${JSON.stringify(defineValues)} preprocesses and is within the core limits`, () => {
      const source = getWgslKernel(code.name);
      expect(source).toBeDefined();
      const defines = Object.fromEntries(
        Object.entries(defineValues).map(([k, v]) => [k, ValueCell.create(v)]),
      ) as DefineValues;
      const body = preprocessShader(source!, getWgslDefines(defines), { expandMacros: true });
      expect(() => getWgslKernelDeclarations(code, schema, defines)).not.toThrow();
      expect(body).toContain('fn main(');
      expect(body).not.toMatch(/^\s*#/m);

      const [x, y, z] = code.workgroupSize;
      expect(x * y * z).toBeLessThanOrEqual(Limits.maxComputeInvocationsPerWorkgroup);
      expect(x).toBeLessThanOrEqual(Limits.maxComputeWorkgroupSizeX);
      expect(y).toBeLessThanOrEqual(Limits.maxComputeWorkgroupSizeY);
      expect(z).toBeLessThanOrEqual(Limits.maxComputeWorkgroupSizeZ);
      const specs = Object.values(schema);
      const buffers = specs.filter((s) => s.type === 'storage-buffer').length;
      const textures = specs.filter((s) => s.type === 'storage-texture').length;
      expect(buffers).toBeLessThanOrEqual(Limits.maxStorageBuffersPerShaderStage);
      expect(textures).toBeLessThanOrEqual(Limits.maxStorageTexturesPerShaderStage);
    });
  }

  it('splits dispatches beyond the per-dimension limit into 2d', () => {
    expect(getDispatchSize(10, 65535)).toEqual({ x: 10, y: 1 });
    expect(getDispatchSize(65535, 65535)).toEqual({ x: 65535, y: 1 });
    expect(getDispatchSize(65536, 65535)).toEqual({ x: 65535, y: 2 });
    expect(getDispatchSize(9, 4)).toEqual({ x: 4, y: 3 });
  });
});
