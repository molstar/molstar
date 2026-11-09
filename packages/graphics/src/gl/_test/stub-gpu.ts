/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// A stub WebGPU device, for what the backend's JavaScript does (the allocation guards): every call accepted, nothing
// drawn or validated. Objects a call makes (`create*`, `begin*`, `finish`, views) are new, as a device's are; other calls
// return one object per method, made once. The values the backend reads are real: the spec's default limits, no
// features, buffers and textures with their descriptors' sizes, and promises that resolve (a device never lost).

/** WebGPU's default limits (what an adapter guarantees) */
const Limits = {
  maxTextureDimension1D: 8192,
  maxTextureDimension2D: 8192,
  maxTextureDimension3D: 2048,
  maxTextureArrayLayers: 256,
  maxBindGroups: 4,
  maxBindGroupsPlusVertexBuffers: 24,
  maxBindingsPerBindGroup: 1000,
  maxDynamicUniformBuffersPerPipelineLayout: 8,
  maxDynamicStorageBuffersPerPipelineLayout: 4,
  maxSampledTexturesPerShaderStage: 16,
  maxSamplersPerShaderStage: 16,
  maxStorageBuffersPerShaderStage: 8,
  maxStorageTexturesPerShaderStage: 4,
  maxUniformBuffersPerShaderStage: 12,
  maxUniformBufferBindingSize: 65536,
  maxStorageBufferBindingSize: 134217728,
  minUniformBufferOffsetAlignment: 256,
  minStorageBufferOffsetAlignment: 256,
  maxVertexBuffers: 8,
  maxBufferSize: 268435456,
  maxVertexAttributes: 16,
  maxVertexBufferArrayStride: 2048,
  maxInterStageShaderVariables: 16,
  maxColorAttachments: 8,
  maxColorAttachmentBytesPerSample: 32,
  maxComputeWorkgroupStorageSize: 16384,
  maxComputeInvocationsPerWorkgroup: 256,
  maxComputeWorkgroupSizeX: 256,
  maxComputeWorkgroupSizeY: 256,
  maxComputeWorkgroupSizeZ: 64,
  maxComputeWorkgroupsPerDimension: 65535,
};

const Makes = /^(create|begin)|^(finish|getCurrentTexture)$/;

/** The calls of the stubs' own methods by name (`draw`, `drawIndexed`, `setBindGroup`, ...), for a test to read */
export const stubCalls = new Map<string, number>();

/**
 * The stubs' methods by name, shared by all of them (a method's result doesn't depend on its object): a new object a
 * frame (a pass encoder) costs its own allocation only
 */
const methods = new Map<PropertyKey, unknown>();
const handler: ProxyHandler<object> = {
  get(target, key) {
    if (key in target) return (target as any)[key];
    if (key === 'then' || typeof key === 'symbol') return undefined;
    let m = methods.get(key);
    if (m === undefined) {
      const name = String(key);
      stubCalls.set(name, 0);
      if (Makes.test(name)) {
        m = () => {
          stubCalls.set(name, stubCalls.get(name)! + 1);
          return stub();
        };
      } else {
        const result = stub();
        m = () => {
          stubCalls.set(name, stubCalls.get(name)! + 1);
          return result;
        };
      }
      methods.set(key, m);
    }
    return m;
  },
};

/** An object whose members are stubs: properties of `own` first, then methods (`Makes`: a new object a call) */
function stub(own: object = {}): any {
  return new Proxy(own, handler);
}

function stubBuffer(desc: GPUBufferDescriptor) {
  return stub({
    size: desc.size,
    usage: desc.usage,
    mapState: 'unmapped',
    label: desc.label ?? '',
    mapAsync: () => Promise.resolve(),
    getMappedRange: (offset = 0, size?: number) => new ArrayBuffer(size ?? desc.size - offset),
    unmap: () => {},
    destroy: () => {},
  });
}

function stubTexture(desc: GPUTextureDescriptor) {
  const size = desc.size as number[] | GPUExtent3DDict;
  const [width, height = 1, depthOrArrayLayers = 1] = Array.isArray(size)
    ? size
    : [size.width, size.height, size.depthOrArrayLayers];
  return stub({
    width,
    height,
    depthOrArrayLayers,
    format: desc.format,
    usage: desc.usage,
    dimension: desc.dimension ?? '2d',
    mipLevelCount: desc.mipLevelCount ?? 1,
    sampleCount: desc.sampleCount ?? 1,
    label: desc.label ?? '',
    createView: () => stub(),
    destroy: () => {},
  });
}

/** A device that accepts every call (`Limits`, no features) */
export function createStubDevice(): GPUDevice {
  const queue = stub({
    submit: () => {},
    writeBuffer: () => {},
    writeTexture: () => {},
    onSubmittedWorkDone: () => Promise.resolve(),
  });
  return stub({
    label: 'stub',
    limits: Limits,
    features: new Set<string>(),
    adapterInfo: { vendor: 'stub', architecture: '', device: '', description: 'stub' },
    queue,
    // never lost
    lost: new Promise(() => {}),
    createBuffer: stubBuffer,
    createTexture: stubTexture,
    createShaderModule: () => stub({ getCompilationInfo: () => Promise.resolve({ messages: [] }) }),
    createRenderPipelineAsync: () => Promise.resolve(stub()),
    createComputePipelineAsync: () => Promise.resolve(stub()),
    pushErrorScope: () => {},
    popErrorScope: () => Promise.resolve(null),
    addEventListener: () => {},
    removeEventListener: () => {},
    destroy: () => {},
  }) as GPUDevice;
}
