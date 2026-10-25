/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isTimingMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import { idFactory } from '@molstar/core/util/id-factory';
import type {
  ComputeContext,
  ComputeKernel,
  KernelCode,
  KernelSchema,
  KernelValues,
  StorageBuffer,
  StorageTextureAccess,
  KernelRequest,
} from '../compute.js';
import { ValueCell } from '@molstar/core/util';
import type { PendingWork } from '../pending-work.js';
import {
  type DefineLookup,
  getTextureDimension,
  type RenderableSchema,
  type TextureSpec,
} from '../renderable/schema.js';
import type { DefineType, DefineValues } from '../shader-code.js';
import { preprocessShader } from '../shader-preprocessor.js';
import type { GraphicsStats } from '../stats.js';
import type { TextureKind } from '../texture.js';
import type { UniformsList } from '../uniform.js';
import { createUniformBlockData, createUniformBlockWriter, getUniformBlockLayout } from '../uniform-layout.js';
import { BufferUsage, MapMode, ShaderStage } from './constants.js';
import { getOverrideConstants, getWgslDefineLookup, getWgslDefines } from './program.js';
import { getWgslKernel } from './shader/kernels.js';
import { asWebGPUTexture, getGPUTextureFormat, type Release, type WebGPUByteCounts } from './texture.js';
import { getWgslType } from './uniform-declarations.js';
import { createUniformSnapshot, type UniformRing } from './uniform-ring.js';
import { getWgslTextureType } from './wgsl-declarations.js';

/** What the compute context needs of the WebGPU context */
export interface WebGPUComputeEnv {
  readonly device: GPUDevice;
  readonly ring: UniformRing;
  readonly stats: GraphicsStats;
  readonly pending: PendingWork;
  release: Release;
  /** Records commands outside of passes (ends the open pass encoder) */
  encode: (commands: (encoder: GPUCommandEncoder) => void) => void;
  /** Submits the recorded commands */
  flush: () => void;
  /** The encoder of the open compute pass, begun on first use */
  getComputePassEncoder: () => GPUComputePassEncoder;
  beginComputePass: (options: { label: string }) => void;
  endComputePass: () => void;
  onError: (message: string) => void;
  /** bytes of the storage buffers created and not destroyed yet */
  storageBufferBytes: number;
  /** bytes of the other live resources, here the staging buffers of reads (`pixelPack`) */
  readonly byteCounts: WebGPUByteCounts;
}

const getNextStorageBufferId = idFactory();
const getNextKernelPipelineId = idFactory();

/**
 * Storage buffers can also be copied, cleared, read back, used for indirect dispatches and draws, and bound as vertex
 * buffers (instanced attributes written by kernels, GPU culling)
 */
const StorageBufferUsage =
  BufferUsage.STORAGE | BufferUsage.COPY_SRC | BufferUsage.COPY_DST | BufferUsage.INDIRECT | BufferUsage.VERTEX;

export interface WebGPUStorageBuffer extends StorageBuffer {
  getGPUBuffer: () => GPUBuffer;
  /** bumped when the GPU buffer is recreated (device loss) */
  readonly version: number;
}

const StorageBuffers = new WeakSet<StorageBuffer>();

export function isWebGPUStorageBuffer(x: unknown): x is WebGPUStorageBuffer {
  return StorageBuffers.has(x as StorageBuffer);
}

export function asWebGPUStorageBuffer(buffer: StorageBuffer): WebGPUStorageBuffer {
  if (!StorageBuffers.has(buffer)) throw new Error('expected a WebGPU storage buffer');
  return buffer as WebGPUStorageBuffer;
}

function createStorageBuffer(env: WebGPUComputeEnv, byteLength: number, label?: string): WebGPUStorageBuffer {
  if (byteLength <= 0 || byteLength % 4 !== 0) {
    throw new Error(`storage buffer size must be a positive multiple of 4, got ${byteLength}`);
  }
  const create = () => env.device.createBuffer({ size: byteLength, usage: StorageBufferUsage, label });
  const gpuBuffer = create();
  const version = 0;
  let destroyed = false;
  env.stats.resourceCounts.storageBuffer += 1;
  env.storageBufferBytes += byteLength;

  const range = (byteOffset = 0, length = byteLength - byteOffset) => {
    if (byteOffset < 0 || length < 0 || byteOffset + length > byteLength || byteOffset % 4 || length % 4) {
      throw new Error(`invalid storage buffer range ${byteOffset}, ${length} of ${byteLength} bytes`);
    }
    return length;
  };

  const buffer: WebGPUStorageBuffer = {
    id: getNextStorageBufferId(),
    byteLength,
    get version() {
      return version;
    },
    getGPUBuffer: () => gpuBuffer,
    write: (data, byteOffset = 0) => {
      range(byteOffset, data.byteLength);
      env.device.queue.writeBuffer(gpuBuffer, byteOffset, data.buffer, data.byteOffset, data.byteLength);
    },
    clear: (byteOffset = 0, length?: number) => {
      const size = range(byteOffset, length);
      env.encode((encoder) => encoder.clearBuffer(gpuBuffer, byteOffset, size));
    },
    copyTo: (target, srcByteOffset = 0, dstByteOffset = 0, length?: number) => {
      const size = range(
        srcByteOffset,
        length ?? Math.min(byteLength - srcByteOffset, target.byteLength - dstByteOffset),
      );
      const dst = asWebGPUStorageBuffer(target).getGPUBuffer();
      env.encode((encoder) => encoder.copyBufferToBuffer(gpuBuffer, srcByteOffset, dst, dstByteOffset, size));
    },
    readAsync: (byteOffset = 0, length?: number) => {
      const size = range(byteOffset, length);
      const { device } = env;
      const staging = device.createBuffer({ size, usage: BufferUsage.MAP_READ | BufferUsage.COPY_DST });
      env.byteCounts.pixelPack += size;
      env.stats.resourceCounts.pixelPack += 1;
      env.encode((encoder) => encoder.copyBufferToBuffer(gpuBuffer, byteOffset, staging, 0, size));
      // a buffer can only be mapped once the commands using it are submitted
      env.flush();
      const read = staging
        .mapAsync(MapMode.READ)
        .then(() => {
          const out = staging.getMappedRange().slice(0);
          staging.unmap();
          return out;
        })
        .finally(() => {
          staging.destroy();
          env.byteCounts.pixelPack -= size;
          env.stats.resourceCounts.pixelPack -= 1;
        });
      return env.pending.track(read, 'storage buffer read');
    },
    getByteCount: () => byteLength,
    destroy: () => {
      if (destroyed) return;
      destroyed = true;
      env.release(gpuBuffer);
      env.stats.resourceCounts.storageBuffer -= 1;
      env.storageBufferBytes -= byteLength;
    },
  };
  StorageBuffers.add(buffer);
  return buffer;
}

//

type KernelBinding =
  | { kind: 'uniforms'; binding: number }
  | { kind: 'storage-buffer'; name: string; binding: number; readOnly: boolean }
  | { kind: 'texture'; name: string; binding: number; spec: TextureSpec<TextureKind> }
  | {
      kind: 'storage-texture';
      name: string;
      binding: number;
      format: GPUTextureFormat;
      access: GPUStorageTextureAccess;
      dimension: '2d' | '3d';
    };

/** Bindings of a kernel schema in group 0: the uniform block (if any) first, then the schema order */
function getKernelBindings(schema: KernelSchema, defines: DefineValues) {
  const lookup = getWgslDefineLookup(getWgslDefines(defines));
  const uniforms = getUniformBlockLayout('object', schema as unknown as RenderableSchema);
  const bindings: KernelBinding[] = [];
  if (uniforms) bindings.push({ kind: 'uniforms', binding: 0 });
  for (const name of Object.keys(schema)) {
    const spec = schema[name];
    const binding = bindings.length;
    if (spec.type === 'storage-buffer') {
      bindings.push({ kind: 'storage-buffer', name, binding, readOnly: spec.access === 'read' });
    } else if (spec.type === 'texture') {
      bindings.push({ kind: 'texture', name, binding, spec });
    } else if (spec.type === 'storage-texture') {
      const format = getGPUTextureFormat(spec.format, spec.dataType);
      const access = GPUStorageAccess[spec.access];
      bindings.push({ kind: 'storage-texture', name, binding, format, access, dimension: spec.dimension });
    }
  }
  return { uniforms, bindings, lookup };
}

const GPUStorageAccess: { [k in StorageTextureAccess]: GPUStorageTextureAccess } = {
  write: 'write-only',
  read: 'read-only',
  'read-write': 'read-write',
};

const WgslStorageAccess: { [k in GPUStorageTextureAccess]: string } = {
  'write-only': 'write',
  'read-only': 'read',
  'read-write': 'read_write',
};

/** WGSL declarations of the bindings of `schema` and the `WorkgroupSize*` constants of `code` */
/** The WGSL of a kernel: its own source, or the registered one */
function getKernelSource(code: KernelCode) {
  return code.sources?.wgsl ?? getWgslKernel(code.name);
}

export function getWgslKernelDeclarations(code: KernelCode, schema: KernelSchema, defines: DefineValues) {
  const { uniforms, bindings, lookup } = getKernelBindings(schema, defines);
  const lines: string[] = [];
  const [x, y, z] = code.workgroupSize;
  lines.push(`const WorkgroupSizeX = ${x}u;`, `const WorkgroupSizeY = ${y}u;`, `const WorkgroupSizeZ = ${z}u;`);
  if (uniforms) {
    lines.push('struct KernelParams {');
    for (const m of uniforms.members) lines.push(`    ${m.name}: ${getWgslType(m)},`);
    lines.push('};');
  }
  for (const b of bindings) {
    const at = `@group(0) @binding(${b.binding})`;
    switch (b.kind) {
      case 'uniforms':
        lines.push(`${at} var<uniform> uParams: KernelParams;`);
        break;
      case 'storage-buffer': {
        const spec = schema[b.name] as { element: string };
        lines.push(`${at} var<storage, ${b.readOnly ? 'read' : 'read_write'}> ${b.name}: array<${spec.element}>;`);
        break;
      }
      case 'texture': {
        const type = getWgslTextureType({
          dimension: getTextureDimension(b.spec, lookup),
          sampleType: b.spec.sampleType,
        });
        lines.push(`${at} var ${b.name}: ${type};`);
        break;
      }
      case 'storage-texture':
        lines.push(`${at} var ${b.name}: texture_storage_${b.dimension}<${b.format}, ${WgslStorageAccess[b.access]}>;`);
        break;
    }
  }
  return lines.join('\n');
}

function getBindGroupLayoutEntry(b: KernelBinding, lookup: DefineLookup): GPUBindGroupLayoutEntry {
  const visibility = ShaderStage.COMPUTE;
  switch (b.kind) {
    case 'uniforms':
      return { binding: b.binding, visibility, buffer: { type: 'uniform', hasDynamicOffset: true } };
    case 'storage-buffer':
      return { binding: b.binding, visibility, buffer: { type: b.readOnly ? 'read-only-storage' : 'storage' } };
    case 'texture': {
      const dimension = getTextureDimension(b.spec, lookup);
      // kernels read textures with `textureLoad`
      const sampleType = b.spec.sampleType === 'float' ? 'unfilterable-float' : b.spec.sampleType;
      return {
        binding: b.binding,
        visibility,
        texture: { sampleType, viewDimension: dimension === 'cube' ? 'cube' : dimension },
      };
    }
    case 'storage-texture':
      return {
        binding: b.binding,
        visibility,
        storageTexture: { access: b.access, format: b.format, viewDimension: b.dimension },
      };
  }
}

type KernelPipeline = {
  readonly id: number;
  readonly pipeline: GPUComputePipeline;
  readonly layout: GPUBindGroupLayout;
  readonly bindings: KernelBinding[];
  readonly uniforms: ReturnType<typeof getUniformBlockLayout>;
};

function getDefineKey(defines: DefineValues) {
  return [...getWgslDefines(defines)]
    .map(([k, v]) => (v ? `${k}=${v}` : k))
    .sort()
    .join(',');
}

function getKernelKey(code: KernelCode, schema: KernelSchema, defines: DefineValues) {
  return `${code.name}|${getDefineKey(defines)}|${Object.keys(schema).join(',')}`;
}

/** The pipeline descriptor of a kernel, and what its kernels bind */
function getKernelPipelineDescriptor(
  env: WebGPUComputeEnv,
  code: KernelCode,
  schema: KernelSchema,
  defines: DefineValues,
) {
  const source = getKernelSource(code);
  if (!source) throw new Error(`no WGSL source for kernel '${code.name}'`);
  const wgslDefines = getWgslDefines(defines);
  const body = preprocessShader(source, wgslDefines, { expandMacros: true });
  const wgsl = `${getWgslKernelDeclarations(code, schema, defines)}\n\n${body}`;
  const { device } = env;
  const module = device.createShaderModule({ code: wgsl, label: code.name });
  env.stats.resourceCounts.shader += 1;
  module.getCompilationInfo().then((info) => {
    const lines = wgsl.split('\n');
    for (const m of info.messages) {
      if (m.type !== 'error') continue;
      env.onError(
        `WGSL error in kernel '${code.name}' line ${m.lineNum}: ${m.message}\n  ${lines[m.lineNum - 1]?.trim()}`,
      );
    }
  });
  const { uniforms, bindings, lookup } = getKernelBindings(schema, defines);
  const layout = device.createBindGroupLayout({
    entries: bindings.map((b) => getBindGroupLayoutEntry(b, lookup)),
    label: code.name,
  });
  const descriptor: GPUComputePipelineDescriptor = {
    layout: device.createPipelineLayout({ bindGroupLayouts: [layout], label: code.name }),
    compute: { module, entryPoint: 'main', constants: getOverrideConstants(body, defines, 'main') },
    label: code.name,
  };
  return { descriptor, layout, bindings, uniforms };
}

/** Synchronously, for kernels not prepared: their dispatches must not be skipped, the queue waits for the compile */
function createKernelPipeline(
  env: WebGPUComputeEnv,
  code: KernelCode,
  schema: KernelSchema,
  defines: DefineValues,
): KernelPipeline {
  const { descriptor, layout, bindings, uniforms } = getKernelPipelineDescriptor(env, code, schema, defines);
  const start = isTimingMode ? now() : 0;
  const pipeline = env.device.createComputePipeline(descriptor);
  if (isTimingMode)
    console.log(`Created kernel pipeline ${code.name} synchronously in ${(now() - start).toFixed(2)} ms`);
  env.stats.resourceCounts.kernel += 1;
  return { id: getNextKernelPipelineId(), pipeline, layout, bindings, uniforms };
}

function getKernelRequestKey({ code, schema, defines: values }: KernelRequest) {
  const defines: DefineValues = {};
  for (const k of Object.keys(values ?? {})) defines[k] = ValueCell.create(values![k] as DefineType);
  return { key: getKernelKey(code, schema, defines), defines };
}

/** In the background, see `ComputeContext.prepareKernels` */
function prepareKernelPipelines(
  env: WebGPUComputeEnv,
  pipelines: Map<string, KernelPipeline>,
  compiling: Map<string, Promise<void>>,
  requests: Iterable<KernelRequest>,
): Promise<void> {
  const waits: Promise<void>[] = [];
  for (const request of requests) {
    const { code, schema } = request;
    const { key, defines } = getKernelRequestKey(request);
    if (pipelines.has(key)) continue;
    let compiled = compiling.get(key);
    if (!compiled) {
      const { descriptor, layout, bindings, uniforms } = getKernelPipelineDescriptor(env, code, schema, defines);
      const start = now();
      compiled = env.device.createComputePipelineAsync(descriptor).then(
        (pipeline) => {
          compiling.delete(key);
          if (isTimingMode) console.log(`Compiled kernel pipeline ${code.name} in ${(now() - start).toFixed(2)} ms`);
          // a kernel created meanwhile compiled its own synchronously
          if (pipelines.has(key)) return;
          env.stats.resourceCounts.kernel += 1;
          pipelines.set(key, { id: getNextKernelPipelineId(), pipeline, layout, bindings, uniforms });
        },
        // reported by the synchronous creation of the first kernel using it
        () => {
          compiling.delete(key);
        },
      );
      compiling.set(key, compiled);
      env.pending.track(compiled, `kernel pipeline ${code.name}`);
    }
    waits.push(compiled);
  }
  return waits.length ? Promise.all(waits).then(() => {}) : Promise.resolve();
}

function createKernel<S extends KernelSchema>(
  env: WebGPUComputeEnv,
  pipelines: Map<string, KernelPipeline>,
  code: KernelCode,
  schema: S,
  values: KernelValues<S>,
): ComputeKernel {
  const defineValues = () => {
    const d: { [k: string]: DefineValues[string] } = {};
    for (const k of Object.keys(schema)) if (schema[k].type === 'define') d[k] = values[k] as DefineValues[string];
    return d as DefineValues;
  };
  let current: KernelPipeline;
  function getPipeline() {
    const defines = defineValues();
    const key = getKernelKey(code, schema, defines);
    let p = pipelines.get(key);
    if (!p) pipelines.set(key, (p = createKernelPipeline(env, code, schema, defines)));
    return p;
  }
  current = getPipeline();

  const uniformEntries: UniformsList = Object.keys(schema)
    .filter((k) => schema[k].type === 'uniform')
    .map((k) => [k, values[k] as UniformsList[number][1]]);
  const snapshot = current.uniforms ? createUniformSnapshot(env.ring, current.uniforms.size) : undefined;
  const writeUniforms =
    current.uniforms && snapshot
      ? createUniformBlockWriter(
          current.uniforms,
          createUniformBlockData(snapshot.bytes.buffer as ArrayBuffer, 0, current.uniforms.size),
          uniformEntries,
          snapshot.markChanged,
        )
      : undefined;

  let bindGroup: GPUBindGroup | undefined;
  let bindKey = '';
  function getBindGroup(uniformBuffer: GPUBuffer | undefined) {
    const entries: GPUBindGroupEntry[] = [];
    const keys: (string | number)[] = [current.id];
    for (const b of current.bindings) {
      if (b.kind === 'uniforms') {
        entries.push({ binding: b.binding, resource: { buffer: uniformBuffer!, size: current.uniforms!.size } });
        keys.push(uniformBuffers.get(uniformBuffer!)!);
      } else if (b.kind === 'storage-buffer') {
        const buffer = asWebGPUStorageBuffer(values[b.name].ref.value as StorageBuffer);
        entries.push({ binding: b.binding, resource: { buffer: buffer.getGPUBuffer() } });
        keys.push(buffer.id, buffer.version);
      } else {
        const texture = asWebGPUTexture(values[b.name].ref.value as never);
        // read before written: allocated (zeros), as a kernel input of undefined contents
        if (b.kind === 'texture') texture.getGPUTexture();
        const view = b.kind === 'texture' ? texture.getView() : texture.getAttachmentView();
        if (!view) throw new Error(`texture '${b.name}' of kernel '${code.name}' has no storage`);
        entries.push({ binding: b.binding, resource: view });
        keys.push(texture.id, texture.version);
      }
    }
    const key = keys.join('|');
    if (!bindGroup || key !== bindKey) {
      bindGroup = env.device.createBindGroup({ layout: current.layout, entries, label: code.name });
      bindKey = key;
    }
    return bindGroup;
  }

  function prepare(pass: GPUComputePassEncoder) {
    writeUniforms?.();
    const slot = snapshot?.getSlot();
    if (slot && !uniformBuffers.has(slot.buffer)) uniformBuffers.set(slot.buffer, nextUniformBufferId++);
    pass.setPipeline(current.pipeline);
    pass.setBindGroup(0, getBindGroup(slot?.buffer), slot ? [slot.offset] : []);
  }

  let destroyed = false;
  return {
    update: () => {
      current = getPipeline();
    },
    dispatch: (x: number, y = 1, z = 1) => {
      const max = env.device.limits.maxComputeWorkgroupsPerDimension;
      if (x > max || y > max || z > max) {
        throw new Error(`kernel '${code.name}': ${x} × ${y} × ${z} workgroups exceed ${max} per dimension`);
      }
      if (x === 0 || y === 0 || z === 0) return;
      const pass = env.getComputePassEncoder();
      prepare(pass);
      pass.dispatchWorkgroups(x, y, z);
    },
    dispatchIndirect: (buffer: StorageBuffer, byteOffset = 0) => {
      const pass = env.getComputePassEncoder();
      prepare(pass);
      pass.dispatchWorkgroupsIndirect(asWebGPUStorageBuffer(buffer).getGPUBuffer(), byteOffset);
    },
    getByteCount: () => 0,
    destroy: () => {
      if (destroyed) return;
      destroyed = true;
      writeUniforms?.dispose();
      bindGroup = undefined;
    },
  };
}

/** ids of the uniform ring's buffers, for bind group keys */
const uniformBuffers = new WeakMap<GPUBuffer, number>();
let nextUniformBufferId = 0;

/** The storage buffers a kernel of `schema` binds (one under a condition counted: it may bind it) */
function getStorageBufferCount(schema: KernelSchema) {
  let count = 0;
  for (const k in schema) if (schema[k].type === 'storage-buffer') count += 1;
  return count;
}

export function createWebGPUComputeContext(env: WebGPUComputeEnv): ComputeContext {
  /** pipelines per kernel name, defines and schema */
  const pipelines = new Map<string, KernelPipeline>();
  /** pipelines compiling in the background, by the same key */
  const compiling = new Map<string, Promise<void>>();
  const { limits } = env.device;
  const ownBuffers = new WeakSet<StorageBuffer>();
  return {
    maxStorageBufferSize: Math.min(limits.maxStorageBufferBindingSize, limits.maxBufferSize),
    getStorageBufferByteCount: () => env.storageBufferBytes,
    isKernelSupported: (code, schema) =>
      getKernelSource(code) !== undefined &&
      (!schema || getStorageBufferCount(schema) <= env.device.limits.maxStorageBuffersPerShaderStage),
    createKernel: (code, schema, values) => createKernel(env, pipelines, code, schema, values),
    prepareKernels: (requests) => prepareKernelPipelines(env, pipelines, compiling, requests),
    hasKernels: (requests) => {
      for (const r of requests) if (!pipelines.has(getKernelRequestKey(r).key)) return false;
      return true;
    },
    createBuffer: (byteLength, label) => {
      const buffer = createStorageBuffer(env, byteLength, label);
      ownBuffers.add(buffer);
      return buffer;
    },
    tryCreateBuffer: async (byteLength, label) => {
      env.device.pushErrorScope('out-of-memory');
      const buffer = createStorageBuffer(env, byteLength, label);
      const error = await env.device.popErrorScope();
      if (error) {
        buffer.destroy();
        return undefined;
      }
      ownBuffers.add(buffer);
      return buffer;
    },
    ownsBuffer: (buffer) => ownBuffers.has(buffer),
    beginComputePass: (options) => env.beginComputePass(options),
    endComputePass: () => env.endComputePass(),
  };
}
