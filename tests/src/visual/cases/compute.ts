/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { type ComputeContext, KernelCode, StorageBufferSpec, StorageTextureSpec } from '@molstar/graphics/gl/compute';
import { createPrefixSum } from '@molstar/graphics/gl/compute/prefix-sum';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { UniformSpec } from '@molstar/graphics/gl/renderable/schema';
import { registerWgslKernel } from '@molstar/graphics/gl/webgpu/shader/kernels';
import { type VisualTestCase } from '../framework/test-case.js';

/** Kernels of the conformance checks, registered with the WebGPU backend */
export const FillKernel = KernelCode('test-fill', [64]);
registerWgslKernel(
  'test-fill',
  /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if (id.x < u32(uParams.uCount)) {
        data[u32(uParams.uOffset) + id.x] = u32(uParams.uValue) + id.x;
    }
}
`,
);
export const FillSchema = {
  uCount: UniformSpec('i'),
  uOffset: UniformSpec('i'),
  uValue: UniformSpec('i'),
  data: StorageBufferSpec('u32', 'read-write'),
};

const AtomicsKernel = KernelCode('test-atomics', [64]);
registerWgslKernel(
  'test-atomics',
  /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if (id.x >= arrayLength(&values)) {
        return;
    }
    let v = values[id.x];
    if (v > u32(uParams.uThreshold)) {
        atomicAdd(&result[0], 1u);
    }
    atomicMin(&result[1], v);
    atomicMax(&result[2], v);
}
`,
);
const AtomicsSchema = {
  uThreshold: UniformSpec('i'),
  values: StorageBufferSpec('u32', 'read'),
  result: StorageBufferSpec('atomic<u32>', 'read-write'),
};

const TextureKernel = KernelCode('test-storage-texture', [8, 8]);
registerWgslKernel(
  'test-storage-texture',
  /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let size = textureDimensions(outTexture);
    if (id.x >= size.x || id.y >= size.y) {
        return;
    }
    textureStore(outTexture, id.xy, vec4<f32>(f32(id.x) / 255.0, f32(id.y) / 255.0, 128.0 / 255.0, 1.0));
}
`,
);
const TextureSchema = {
  outTexture: StorageTextureSpec('rgba', 'ubyte', '2d', 'write'),
};

function u32(buffer: ArrayBuffer) {
  return new Uint32Array(buffer);
}

function getCompute(graphics: GraphicsContext): ComputeContext {
  if (!graphics.compute) throw new Error('no compute context');
  return graphics.compute;
}

/** Up to 5 mismatches of `actual` against `expected` */
function compareArrays(label: string, actual: ArrayLike<number>, expected: ArrayLike<number>) {
  const errors: string[] = [];
  if (actual.length < expected.length) return [`${label}: ${actual.length} values, expected ${expected.length}`];
  for (let i = 0; i < expected.length && errors.length < 5; ++i) {
    if (actual[i] !== expected[i]) errors.push(`${label}[${i}] = ${actual[i]}, expected ${expected[i]}`);
  }
  return errors;
}

/** Deterministic pseudo-random values in [0, max) */
function randomValues(count: number, max: number, seed = 1) {
  const out = new Uint32Array(count);
  let s = seed;
  for (let i = 0; i < count; ++i) {
    s = (s * 1664525 + 1013904223) >>> 0;
    out[i] = s % max;
  }
  return out;
}

async function checkStorageBuffer(graphics: GraphicsContext) {
  const compute = getCompute(graphics);
  const a = compute.createBuffer(64);
  const b = compute.createBuffer(64);
  a.write(new Uint32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]));
  a.copyTo(b);
  b.clear(8, 16);
  const errors = [
    ...compareArrays('a', u32(await a.readAsync()), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]),
    ...compareArrays('b', u32(await b.readAsync()), [1, 2, 0, 0, 0, 0, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]),
    ...compareArrays('b[4..8]', u32(await b.readAsync(16, 16)), [0, 0, 7, 8]),
  ];
  a.destroy();
  b.destroy();
  return errors;
}

/** Two dispatches of one kernel in one pass with different uniforms, and an indirect dispatch */
async function checkDispatches(graphics: GraphicsContext) {
  const compute = getCompute(graphics);
  const data = compute.createBuffer(4 * 256);
  const args = compute.createBuffer(12);
  data.clear();
  const values = {
    uCount: ValueCell.create(100),
    uOffset: ValueCell.create(0),
    uValue: ValueCell.create(1000),
    data: ValueCell.create(data),
  };
  const kernel = compute.createKernel(FillKernel, FillSchema, values);
  compute.beginComputePass({ label: 'test dispatches' });
  kernel.dispatch(Math.ceil(100 / 64));
  ValueCell.update(values.uCount, 50);
  ValueCell.update(values.uOffset, 100);
  ValueCell.update(values.uValue, 2000);
  kernel.dispatch(1);
  // indirect: 2 workgroups of 64, the count limits them to 100
  args.write(new Uint32Array([2, 1, 1]));
  ValueCell.update(values.uCount, 100);
  ValueCell.update(values.uOffset, 150);
  ValueCell.update(values.uValue, 3000);
  kernel.dispatchIndirect(args);
  compute.endComputePass();

  const expected = new Uint32Array(256);
  for (let i = 0; i < 100; ++i) expected[i] = 1000 + i;
  for (let i = 0; i < 50; ++i) expected[100 + i] = 2000 + i;
  for (let i = 0; i < 100 && 150 + i < 256; ++i) expected[150 + i] = 3000 + i;
  const errors = compareArrays('data', u32(await data.readAsync()), expected);
  kernel.destroy();
  data.destroy();
  args.destroy();
  return errors;
}

async function checkAtomics(graphics: GraphicsContext) {
  const compute = getCompute(graphics);
  const input = randomValues(10000, 1000000, 7);
  const values = compute.createBuffer(input.byteLength);
  values.write(input);
  const result = compute.createBuffer(12);
  result.write(new Uint32Array([0, 0xffffffff, 0]));
  const kernel = compute.createKernel(AtomicsKernel, AtomicsSchema, {
    uThreshold: ValueCell.create(500000),
    values: ValueCell.create(values),
    result: ValueCell.create(result),
  });
  compute.beginComputePass({ label: 'test atomics' });
  kernel.dispatch(Math.ceil(input.length / 64));
  compute.endComputePass();

  let count = 0;
  let min = 0xffffffff;
  let max = 0;
  for (const v of input) {
    if (v > 500000) count += 1;
    min = Math.min(min, v);
    max = Math.max(max, v);
  }
  const errors = compareArrays('[count, min, max]', u32(await result.readAsync()), [count, min, max]);
  kernel.destroy();
  values.destroy();
  result.destroy();
  return errors;
}

async function checkPrefixSum(graphics: GraphicsContext) {
  const compute = getCompute(graphics);
  const errors: string[] = [];
  // the default dispatch limit, and a limit of 4 workgroups per dimension to exercise 2d dispatches
  for (const maxWorkgroups of [65535, 4]) {
    const prefixSum = createPrefixSum(compute, maxWorkgroups);
    for (const count of [1, 7, 256, 257, 1000, 70000, 300000]) {
      const input = randomValues(count, 100, count);
      const inputBuffer = compute.createBuffer(input.byteLength);
      const outputBuffer = compute.createBuffer(input.byteLength);
      inputBuffer.write(input);
      compute.beginComputePass({ label: 'test prefix sum' });
      const total = prefixSum.run(inputBuffer, outputBuffer, count);
      compute.endComputePass();

      const expected = new Uint32Array(count);
      let sum = 0;
      for (let i = 0; i < count; ++i) {
        expected[i] = sum;
        sum += input[i];
      }
      const label = `scan of ${count} (max ${maxWorkgroups} workgroups)`;
      errors.push(...compareArrays(label, u32(await outputBuffer.readAsync()), expected));
      errors.push(...compareArrays(`${label} total`, u32(await total.readAsync(0, 4)), [sum]));
      inputBuffer.destroy();
      outputBuffer.destroy();
    }
    prefixSum.destroy();
  }
  return errors;
}

async function checkStorageTexture(graphics: GraphicsContext) {
  const compute = getCompute(graphics);
  const width = 37;
  const height = 21;
  const texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
  texture.define(width, height);
  const kernel = compute.createKernel(TextureKernel, TextureSchema, { outTexture: ValueCell.create(texture) });
  compute.beginComputePass({ label: 'test storage texture' });
  kernel.dispatch(Math.ceil(width / 8), Math.ceil(height / 8));
  compute.endComputePass();

  // rows top to bottom: row `y` of the read is the row the kernel stored at `y`
  const data = await graphics.readTextureAsync(texture, 0, 0, width, height);
  const expected = new Uint8Array(width * height * 4);
  for (let y = 0; y < height; ++y) {
    for (let x = 0; x < width; ++x) {
      expected.set([x, y, 128, 255], (y * width + x) * 4);
    }
  }
  const errors = compareArrays('texels', data, expected);
  kernel.destroy();
  texture.destroy();
  return errors;
}

function computeCase(id: string, description: string, check: VisualTestCase['check']): VisualTestCase {
  return { id: `compute/${id}`, description, requires: ['compute'], scene: async () => [], check };
}

export const ComputeCases: VisualTestCase[] = [
  computeCase('storage-buffer', 'Storage buffer write, copy, clear and read back.', checkStorageBuffer),
  computeCase(
    'dispatches',
    'Dispatches of one kernel with different uniforms in one pass, and an indirect dispatch.',
    checkDispatches,
  ),
  computeCase('atomics', 'Atomic add, min and max on storage buffers.', checkAtomics),
  computeCase('prefix-sum', 'Exclusive prefix sum of u32 against the CPU, 1d and 2d dispatches.', checkPrefixSum),
  computeCase('storage-texture', 'Kernel writes into a storage texture, read back.', checkStorageTexture),
];
