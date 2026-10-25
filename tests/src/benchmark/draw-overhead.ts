/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { BufferUsage, ShaderStage, TextureUsage } from '@molstar/graphics/gl/webgpu/constants';
import { isWebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import type { BenchmarkCase, PathRun } from './framework.js';

/** Draws per frame: the objects of `draw-table/boxes-10k` */
const DrawCount = 10_000;
/** Frames per run, issued back to back */
const Frames = 30;
const Size = 512;
/** Byte stride of per-object uniforms with dynamic offsets (the minimum uniform offset alignment) */
const Stride = 256;
const VertexAttributes = 5;

/**
 * How a frame's draws bind their objects:
 * - `per-object`: as Mol* draws boxes today (see `draw-table/boxes-10k`): per draw an object bind group, a draw block
 *   bind group with two dynamic offsets, five vertex buffers and an index buffer of its own, then the draw
 * - `per-object-groups`, `per-object-buffers`: the same commands, only the bind groups or only the buffers of their
 *   own (the others shared): which of them costs
 * - `buffer-offsets`: the same commands, shared bind groups, the buffers shared with an offset per draw (geometry
 *   sub-allocated in shared buffers, set per draw)
 * - `shared-offsets`: one shared bind group set per draw with the object's dynamic offset, shared buffers
 * - `shared-offsets-indirect`: the same, each draw an indexed indirect draw of its arguments in one buffer (what GPU
 *   culling writes); `-culled`: half of them with an instance count of 0 (culled)
 * - `shared-first-instance`: shared bind groups and buffers set once, per draw only the draw, the object selected by
 *   `firstInstance` from a storage array (layer B of `webgpu-draw-submission.md`); `-indirect`: each draw an indexed
 *   indirect draw with its first instance in the arguments (`indirect-first-instance`; what GPU culling of merged
 *   render objects draws: a fixed list, one per segment and level), `-culled` half of them, `-empty` all with an
 *   instance count of 0
 * - `instanced`: one draw of all objects (batched)
 * - `merged`: the objects' geometry in one buffer, indices rebased to it, an object id per vertex selecting the object
 *   from a storage array: one indexed draw of all (batched without instancing, objects of different geometry);
 *   `-culled`: half the objects culled at random, the visible ones drawn as runs of consecutive objects
 */
type Variant =
  | 'per-object'
  | 'per-object-groups'
  | 'per-object-buffers'
  | 'buffer-offsets'
  | 'shared-offsets'
  | 'shared-offsets-indirect'
  | 'shared-offsets-indirect-culled'
  | 'shared-first-instance'
  | 'shared-first-instance-indirect'
  | 'shared-first-instance-indirect-culled'
  | 'shared-first-instance-indirect-empty'
  | 'instanced'
  | 'merged'
  | 'merged-culled';

const CubeIndices = new Uint32Array([
  0, 1, 2, 0, 2, 3, 4, 6, 5, 4, 7, 6, 0, 4, 5, 0, 5, 1, 1, 5, 6, 1, 6, 2, 2, 6, 7, 2, 7, 3, 3, 7, 4, 3, 4, 0,
]);
const CubeCorners = [
  [-1, -1, -1],
  [1, -1, -1],
  [1, 1, -1],
  [-1, 1, -1],
  [-1, -1, 1],
  [1, -1, 1],
  [1, 1, 1],
  [-1, 1, 1],
];

/** the per-object vertex data: the cube's corners, the other attributes zero */
function cubeVertices(attribute: number) {
  const data = new Float32Array(8 * 4);
  if (attribute === 0) {
    for (let i = 0; i < 8; ++i) data.set([...CubeCorners[i], 1], i * 4);
  }
  return data;
}

/** where object `i` sits: a grid filling clip space, small cubes */
function objectOffset(i: number) {
  const n = Math.ceil(Math.sqrt(DrawCount));
  const x = ((i % n) + 0.5) / n;
  const y = (Math.floor(i / n) + 0.5) / n;
  return [x * 2 - 1, y * 2 - 1, 0.5, 0.3 / n];
}

const PerObjectShader = /* wgsl */ `
@group(0) @binding(0) var<uniform> globals: vec4f;
@group(1) @binding(0) var<uniform> object: vec4f;
@group(2) @binding(0) var<uniform> drawA: vec4f;
@group(2) @binding(1) var<uniform> drawB: vec4f;
struct In { @location(0) p: vec4f, @location(1) a: vec4f, @location(2) b: vec4f, @location(3) c: vec4f, @location(4) d: vec4f }
@vertex fn vs(v: In) -> @builtin(position) vec4f {
  let extra = (v.a + v.b + v.c + v.d + drawA + drawB) * 0.0;
  return vec4f(object.xyz + v.p.xyz * object.w + extra.xyz, 1.0) * globals.x;
}
@fragment fn fs() -> @location(0) vec4f { return vec4f(1.0, 0.5, 0.2, 1.0); }
`;

const SharedOffsetShader = /* wgsl */ `
@group(0) @binding(0) var<uniform> globals: vec4f;
@group(1) @binding(0) var<uniform> object: vec4f;
@group(1) @binding(1) var<storage, read> corners: array<vec4f>;
@vertex fn vs(@builtin(vertex_index) v: u32) -> @builtin(position) vec4f {
  return vec4f(object.xyz + corners[v].xyz * object.w, 1.0) * globals.x;
}
@fragment fn fs() -> @location(0) vec4f { return vec4f(1.0, 0.5, 0.2, 1.0); }
`;

const MergedShader = /* wgsl */ `
@group(0) @binding(0) var<uniform> globals: vec4f;
@group(1) @binding(0) var<storage, read> objects: array<vec4f>;
struct In { @location(0) p: vec4f, @location(1) object: u32 }
@vertex fn vs(v: In) -> @builtin(position) vec4f {
  let object = objects[v.object];
  return vec4f(object.xyz + v.p.xyz * object.w, 1.0) * globals.x;
}
@fragment fn fs() -> @location(0) vec4f { return vec4f(1.0, 0.5, 0.2, 1.0); }
`;

const FirstInstanceShader = /* wgsl */ `
@group(0) @binding(0) var<uniform> globals: vec4f;
@group(1) @binding(0) var<storage, read> objects: array<vec4f>;
@group(1) @binding(1) var<storage, read> corners: array<vec4f>;
@vertex fn vs(@builtin(vertex_index) v: u32, @builtin(instance_index) i: u32) -> @builtin(position) vec4f {
  let object = objects[i];
  return vec4f(object.xyz + corners[v].xyz * object.w, 1.0) * globals.x;
}
@fragment fn fs() -> @location(0) vec4f { return vec4f(1.0, 0.5, 0.2, 1.0); }
`;

type Scene = {
  readonly device: GPUDevice;
  readonly target: GPUTexture;
  readonly pipeline: GPURenderPipeline;
  /** draws per frame where not `DrawCount` */
  readonly runs?: number;
  /** encodes the frame's draws into `encoder` */
  readonly encode: (encoder: GPURenderPassEncoder | GPURenderBundleEncoder) => void;
};

function createBuffer(device: GPUDevice, data: ArrayBufferView, usage: number) {
  const buffer = device.createBuffer({ size: Math.max(16, data.byteLength), usage: usage | BufferUsage.COPY_DST });
  device.queue.writeBuffer(buffer, 0, data.buffer, data.byteOffset, data.byteLength);
  return buffer;
}

function createScene(device: GPUDevice, variant: Variant): Scene {
  const target = device.createTexture({
    size: { width: Size, height: Size },
    format: 'rgba8unorm',
    usage: TextureUsage.RENDER_ATTACHMENT,
  });
  const U = BufferUsage.UNIFORM;
  const S = BufferUsage.STORAGE;
  const globalsBuffer = createBuffer(device, new Float32Array([1, 0, 0, 0]), U);
  const indexBuffer = createBuffer(device, CubeIndices, BufferUsage.INDEX);
  const corners = new Float32Array(8 * 4);
  for (let i = 0; i < 8; ++i) corners.set([...CubeCorners[i], 1], i * 4);
  const cornerBuffer = createBuffer(device, corners, S);

  const perObject =
    variant === 'per-object' ||
    variant === 'per-object-groups' ||
    variant === 'per-object-buffers' ||
    variant === 'buffer-offsets';
  const sharedOffsets = variant.startsWith('shared-offsets');
  const merged = variant.startsWith('merged');
  const shader = perObject
    ? PerObjectShader
    : sharedOffsets
      ? SharedOffsetShader
      : merged
        ? MergedShader
        : FirstInstanceShader;
  const module = device.createShaderModule({ code: shader });
  const vertexBuffers: GPUVertexBufferLayout[] = perObject
    ? Array.from({ length: VertexAttributes }, (_, k) => ({
        arrayStride: 16,
        attributes: [{ shaderLocation: k, offset: 0, format: 'float32x4' as const }],
      }))
    : merged
      ? [
          { arrayStride: 16, attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x4' as const }] },
          { arrayStride: 4, attributes: [{ shaderLocation: 1, offset: 0, format: 'uint32' as const }] },
        ]
      : [];
  const createPipeline = (layout: GPUPipelineLayout | 'auto') =>
    device.createRenderPipeline({
      layout,
      vertex: { module, entryPoint: 'vs', buffers: vertexBuffers },
      fragment: { module, entryPoint: 'fs', targets: [{ format: 'rgba8unorm' }] },
      primitive: { topology: 'triangle-list', cullMode: 'back' },
    });
  const V = ShaderStage.VERTEX;
  const uniform = (binding: number, hasDynamicOffset = false): GPUBindGroupLayoutEntry => ({
    binding,
    visibility: V,
    buffer: { type: 'uniform', hasDynamicOffset },
  });
  const storage = (binding: number): GPUBindGroupLayoutEntry => ({
    binding,
    visibility: V,
    buffer: { type: 'read-only-storage' },
  });
  const globalsLayout = device.createBindGroupLayout({ entries: [uniform(0)] });
  const globals = device.createBindGroup({
    layout: globalsLayout,
    entries: [{ binding: 0, resource: { buffer: globalsBuffer } }],
  });

  if (perObject) {
    const distinctGroups = variant === 'per-object' || variant === 'per-object-groups';
    const distinctBuffers = variant === 'per-object' || variant === 'per-object-buffers';
    // `buffer-offsets`: each attribute and the indices in one shared buffer, the objects' ranges one after another
    const sharedVertices = Array.from({ length: VertexAttributes }, (_, k) => {
      const one = cubeVertices(k);
      const all = new Float32Array(one.length * DrawCount);
      for (let i = 0; i < DrawCount; ++i) all.set(one, i * one.length);
      return createBuffer(device, all, BufferUsage.VERTEX);
    });
    const allIndices = new Uint32Array(CubeIndices.length * DrawCount);
    for (let i = 0; i < DrawCount; ++i) allIndices.set(CubeIndices, i * CubeIndices.length);
    const sharedIndices = createBuffer(device, allIndices, BufferUsage.INDEX);
    const objectLayout = device.createBindGroupLayout({ entries: [uniform(0)] });
    const drawLayout = device.createBindGroupLayout({ entries: [uniform(0, true), uniform(1, true)] });
    const pipeline = createPipeline(
      device.createPipelineLayout({ bindGroupLayouts: [globalsLayout, objectLayout, drawLayout] }),
    );
    const createObject = (i: number) => {
      const objectBuffer = createBuffer(device, new Float32Array(objectOffset(i)), U);
      const drawBuffer = createBuffer(device, new Float32Array(Stride / 2), U);
      return {
        object: device.createBindGroup({
          layout: objectLayout,
          entries: [{ binding: 0, resource: { buffer: objectBuffer } }],
        }),
        draw: device.createBindGroup({
          layout: drawLayout,
          entries: [
            { binding: 0, resource: { buffer: drawBuffer, size: 16 } },
            { binding: 1, resource: { buffer: drawBuffer, size: 16 } },
          ],
        }),
        vertices: Array.from({ length: VertexAttributes }, (_, k) =>
          createBuffer(device, cubeVertices(k), BufferUsage.VERTEX),
        ),
        indices: createBuffer(device, CubeIndices, BufferUsage.INDEX),
      };
    };
    const first = createObject(0);
    const objects = Array.from({ length: DrawCount }, (_, i) => {
      if (i === 0) return first;
      const o = createObject(i);
      return {
        object: distinctGroups ? o.object : first.object,
        draw: distinctGroups ? o.draw : first.draw,
        vertices: distinctBuffers ? o.vertices : first.vertices,
        indices: distinctBuffers ? o.indices : first.indices,
      };
    });
    const offsets = [0, 0];
    return {
      device,
      target,
      pipeline,
      encode: (e) => {
        e.setPipeline(pipeline);
        e.setBindGroup(0, globals);
        if (variant === 'buffer-offsets') {
          for (let i = 0; i < DrawCount; ++i) {
            e.setBindGroup(1, first.object);
            for (let k = 0; k < VertexAttributes; ++k) e.setVertexBuffer(k, sharedVertices[k], i * 8 * 16, 8 * 16);
            e.setIndexBuffer(sharedIndices, 'uint32', i * CubeIndices.byteLength, CubeIndices.byteLength);
            e.setBindGroup(2, first.draw, offsets);
            e.drawIndexed(36, 1);
          }
          return;
        }
        for (const o of objects) {
          e.setBindGroup(1, o.object);
          for (let k = 0; k < VertexAttributes; ++k) e.setVertexBuffer(k, o.vertices[k]);
          e.setIndexBuffer(o.indices, 'uint32');
          e.setBindGroup(2, o.draw, offsets);
          e.drawIndexed(36, 1);
        }
      },
    };
  }

  if (sharedOffsets) {
    // indexed indirect arguments per draw: index count, instance count, first index, base vertex, first instance
    const args = new Uint32Array(DrawCount * 5);
    const culled = variant === 'shared-offsets-indirect-culled';
    for (let i = 0; i < DrawCount; ++i) args.set([36, culled && i % 2 ? 0 : 1, 0, 0, 0], i * 5);
    const argsBuffer = createBuffer(device, args, BufferUsage.INDIRECT);
    const indirect = variant !== 'shared-offsets';
    const objectData = new Float32Array((DrawCount * Stride) / 4);
    for (let i = 0; i < DrawCount; ++i) objectData.set(objectOffset(i), (i * Stride) / 4);
    const objectBuffer = createBuffer(device, objectData, U);
    const objectLayout = device.createBindGroupLayout({ entries: [uniform(0, true), storage(1)] });
    const pipeline = createPipeline(device.createPipelineLayout({ bindGroupLayouts: [globalsLayout, objectLayout] }));
    const objects = device.createBindGroup({
      layout: objectLayout,
      entries: [
        { binding: 0, resource: { buffer: objectBuffer, size: 16 } },
        { binding: 1, resource: { buffer: cornerBuffer } },
      ],
    });
    const offset = [0];
    return {
      device,
      target,
      pipeline,
      encode: (e) => {
        e.setPipeline(pipeline);
        e.setBindGroup(0, globals);
        e.setIndexBuffer(indexBuffer, 'uint32');
        for (let i = 0; i < DrawCount; ++i) {
          offset[0] = i * Stride;
          e.setBindGroup(1, objects, offset);
          if (indirect) e.drawIndexedIndirect(argsBuffer, i * 20);
          else e.drawIndexed(36, 1);
        }
      },
    };
  }

  const objectData = new Float32Array(DrawCount * 4);
  for (let i = 0; i < DrawCount; ++i) objectData.set(objectOffset(i), i * 4);
  const objectBuffer = createBuffer(device, objectData, S);

  if (merged) {
    const positions = new Float32Array(DrawCount * 8 * 4);
    const ids = new Uint32Array(DrawCount * 8);
    const indices = new Uint32Array(DrawCount * CubeIndices.length);
    for (let i = 0; i < DrawCount; ++i) {
      positions.set(corners, i * 8 * 4);
      ids.fill(i, i * 8, (i + 1) * 8);
      for (let k = 0; k < CubeIndices.length; ++k) indices[i * CubeIndices.length + k] = CubeIndices[k] + i * 8;
    }
    const positionBuffer = createBuffer(device, positions, BufferUsage.VERTEX);
    const idBuffer = createBuffer(device, ids, BufferUsage.VERTEX);
    const mergedIndices = createBuffer(device, indices, BufferUsage.INDEX);
    const mergedLayout = device.createBindGroupLayout({ entries: [storage(0)] });
    const pipeline = createPipeline(device.createPipelineLayout({ bindGroupLayouts: [globalsLayout, mergedLayout] }));
    const objects = device.createBindGroup({
      layout: mergedLayout,
      entries: [{ binding: 0, resource: { buffer: objectBuffer } }],
    });
    // runs of consecutive visible objects: [first object, count]
    const runs: number[] = [];
    let seed = 1;
    const random = () => {
      seed = (seed * 16807) % 2147483647;
      return seed / 2147483647;
    };
    for (let i = 0; i < DrawCount; ++i) {
      if (variant === 'merged-culled' && random() < 0.5) continue;
      const last = runs.length - 2;
      if (last >= 0 && runs[last] + runs[last + 1] === i) runs[last + 1] += 1;
      else runs.push(i, 1);
    }
    const n = CubeIndices.length;
    return {
      device,
      target,
      pipeline,
      runs: runs.length / 2,
      encode: (e) => {
        e.setPipeline(pipeline);
        e.setBindGroup(0, globals);
        e.setBindGroup(1, objects);
        e.setVertexBuffer(0, positionBuffer);
        e.setVertexBuffer(1, idBuffer);
        e.setIndexBuffer(mergedIndices, 'uint32');
        for (let r = 0; r < runs.length; r += 2) e.drawIndexed(runs[r + 1] * n, 1, runs[r] * n);
      },
    };
  }
  // indexed indirect arguments per draw: index count, instance count, first index, base vertex, first instance
  let firstInstanceArgs: GPUBuffer | undefined;
  if (variant.startsWith('shared-first-instance-indirect')) {
    if (!device.features.has('indirect-first-instance')) throw new Error("needs 'indirect-first-instance'");
    const args = new Uint32Array(DrawCount * 5);
    for (let i = 0; i < DrawCount; ++i) {
      const empty = variant.endsWith('-empty') || (variant.endsWith('-culled') && i % 2 === 1);
      args.set([36, empty ? 0 : 1, 0, 0, i], i * 5);
    }
    firstInstanceArgs = createBuffer(device, args, BufferUsage.INDIRECT);
  }
  const objectLayout = device.createBindGroupLayout({ entries: [storage(0), storage(1)] });
  const pipeline = createPipeline(device.createPipelineLayout({ bindGroupLayouts: [globalsLayout, objectLayout] }));
  const objects = device.createBindGroup({
    layout: objectLayout,
    entries: [
      { binding: 0, resource: { buffer: objectBuffer } },
      { binding: 1, resource: { buffer: cornerBuffer } },
    ],
  });
  return {
    device,
    target,
    pipeline,
    encode: (e) => {
      e.setPipeline(pipeline);
      e.setBindGroup(0, globals);
      e.setBindGroup(1, objects);
      e.setIndexBuffer(indexBuffer, 'uint32');
      if (variant === 'instanced') {
        e.drawIndexed(36, DrawCount);
      } else if (firstInstanceArgs) {
        for (let i = 0; i < DrawCount; ++i) e.drawIndexedIndirect(firstInstanceArgs, i * 20);
      } else {
        for (let i = 0; i < DrawCount; ++i) e.drawIndexed(36, 1, 0, 0, i);
      }
    },
  };
}

/** `Frames` frames replaying the scene's draws (recorded once into a render bundle, or encoded per frame) */
async function run(graphics: GraphicsContext, variant: Variant, bundled: boolean): Promise<PathRun | string> {
  if (!isWebGPUContext(graphics)) return 'WebGPU only';
  const { device } = graphics;
  const scene = createScene(device, variant);
  const view = scene.target.createView();
  let bundle: GPURenderBundle | undefined;
  let recordMs = 0;
  if (bundled) {
    const t0 = performance.now();
    const encoder = device.createRenderBundleEncoder({ colorFormats: ['rgba8unorm'] });
    scene.encode(encoder);
    bundle = encoder.finish();
    recordMs = performance.now() - t0;
  }
  const frame = () => {
    const encoder = device.createCommandEncoder();
    const pass = encoder.beginRenderPass({
      colorAttachments: [{ view, loadOp: 'clear', storeOp: 'store', clearValue: [0, 0, 0, 1] }],
    });
    if (bundle) pass.executeBundles([bundle]);
    else scene.encode(pass);
    pass.end();
    device.queue.submit([encoder.finish()]);
  };
  // warm-up: pipelines and the first replays
  for (let i = 0; i < 3; ++i) frame();
  await device.queue.onSubmittedWorkDone();

  let submit = 0;
  const start = performance.now();
  for (let i = 0; i < Frames; ++i) {
    const t0 = performance.now();
    frame();
    submit += performance.now() - t0;
  }
  await device.queue.onSubmittedWorkDone();
  const total = performance.now() - start;
  scene.target.destroy();
  return {
    stages: [
      { label: 'submit', ms: submit / Frames },
      { label: 'wait', ms: (total - submit) / Frames },
    ],
    triangles: 0,
    drawCalls: variant === 'instanced' ? 1 : (scene.runs ?? DrawCount),
    notes: bundled ? [`bundle recorded in ${recordMs.toFixed(1)} ms`] : undefined,
  };
}

function overheadCase(variant: Variant, bundled: boolean, description: string): BenchmarkCase<undefined> {
  const id = `draw-overhead/10k-${variant}${bundled ? '-bundle' : '-direct'}`;
  return {
    id,
    description,
    setup: () => undefined,
    paths: { webgpu: (graphics) => run(graphics, variant, bundled) },
  };
}

/**
 * What replaying 10k draws costs on WebGPU by how they bind their objects (`draw-overhead/*`, not in the default
 * runs): the upper bound of shared object records (layer B of `webgpu-draw-submission.md`) and of batching, against
 * draws binding as Mol* does today; tiny cubes, so that the times are the submission's, not the GPU's drawing
 */
export const DrawOverheadCases: BenchmarkCase<any>[] = [
  overheadCase('per-object', false, '10k draws, per-object bind groups and buffers (as Mol*), encoded per frame'),
  overheadCase('per-object', true, '10k draws, per-object bind groups and buffers (as Mol*), one render bundle'),
  overheadCase('per-object-groups', true, '10k draws, per-object bind groups, shared buffers, one render bundle'),
  overheadCase('per-object-buffers', true, '10k draws, shared bind groups, per-object buffers, one render bundle'),
  overheadCase(
    'buffer-offsets',
    true,
    '10k draws, shared bind groups, shared buffers set per draw at an offset, bundle',
  ),
  overheadCase('shared-offsets', true, '10k draws, a shared bind group with a dynamic offset per draw, one bundle'),
  overheadCase('shared-offsets-indirect', true, '10k indexed indirect draws, shared bind group with dynamic offsets'),
  overheadCase('shared-offsets-indirect-culled', true, 'as -indirect, half the draws with instance count 0'),
  overheadCase('shared-first-instance', true, '10k draws, shared bindings set once, objects by firstInstance, bundle'),
  overheadCase('shared-first-instance-indirect', true, 'as shared-first-instance, indexed indirect draws'),
  overheadCase('shared-first-instance-indirect-culled', true, 'as -indirect, half the draws with instance count 0'),
  overheadCase('shared-first-instance-indirect-empty', true, 'as -indirect, all draws with instance count 0'),
  overheadCase('instanced', true, 'one instanced draw of 10k objects, one bundle'),
  overheadCase('merged', true, 'one indexed draw of 10k objects merged (rebased indices, object id per vertex)'),
  overheadCase('merged-culled', true, 'as merged, half the objects culled at random, drawn as runs'),
];
