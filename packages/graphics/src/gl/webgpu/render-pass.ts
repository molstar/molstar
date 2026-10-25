/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isTimingMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import type { TargetDescriptor } from '../context.js';
import type { ClearValues } from '../render-state.js';
import { BufferUsage, ShaderStage } from './constants.js';
import { notSupported } from './errors.js';
import { createUniformSnapshot, type UniformRing } from './uniform-ring.js';

/** The attachment formats of a render pass target, what its pipelines depend on */
export interface WebGPUTargetDescriptor extends TargetDescriptor {
  readonly colorFormats: readonly GPUTextureFormat[];
  readonly depthFormat?: GPUTextureFormat;
  readonly hasStencil: boolean;
}

/** The pipeline cache key of the attachment formats */
export function getTargetKey(colorFormats: readonly GPUTextureFormat[], depthFormat?: GPUTextureFormat) {
  return `${colorFormats.join(',')}|${depthFormat ?? ''}`;
}

/**
 * `target` with its depth (and stencil) attachment read-only: its own pipelines (key), whose depth and stencil writes
 * are masked; a target without depth attachment stays as it is
 */
export function withDepthReadOnly<T extends WebGPUTargetDescriptor>(target: T): T {
  if (!target.depthFormat || target.depthReadOnly) return target;
  return { ...target, key: `${target.key}|depth-read-only`, depthReadOnly: true };
}

/** Attachments of a render pass */
export interface PassTarget extends WebGPUTargetDescriptor {
  readonly colorViews: readonly GPUTextureView[];
  /** the slices rendered into, for attachments of 3d textures */
  readonly colorDepthSlices?: readonly (number | undefined)[];
  readonly depthView?: GPUTextureView;
  readonly width: number;
  readonly height: number;
}

function getClearModule(colorCount: number) {
  const outputs = Array.from({ length: colorCount }, (_, i) => `    @location(${i}) c${i}: vec4<f32>,`);
  const assign = Array.from({ length: colorCount }, (_, i) => `    out.c${i} = uClear.color;`);
  return /* wgsl */ `
struct Clear {
    color: vec4<f32>,
    depth: f32,
};
@group(0) @binding(0) var<uniform> uClear: Clear;

@vertex
fn vertexMain(@builtin(vertex_index) i: u32) -> @builtin(position) vec4<f32> {
    // one triangle covering the viewport
    let p = vec2<f32>(f32((i << 1u) & 2u), f32(i & 2u)) * 2.0 - 1.0;
    return vec4<f32>(p, uClear.depth, 1.0);
}
${
  colorCount > 0
    ? `
struct FragmentOutput {
${outputs.join('\n')}
};

@fragment
fn fragmentMain() -> FragmentOutput {
    var out: FragmentOutput;
${assign.join('\n')}
    return out;
}`
    : ''
}
`;
}

const UniformSize = 32;

/**
 * Clears within the scissor rectangle by drawing: WebGPU only clears whole attachments, as load operations
 * at the start of a pass. Like GL clears, this ignores the write masks of the current render state.
 */
/**
 * Clears inside a pass as a draw: `draw` with the target's viewport and the stencil reference (pass state, for a pass);
 * `drawInBundle` without them, for a render bundle: under the pass's viewport (and scissor), a stencil clear to 0 with
 * the `zero` operation, the values in a constant buffer per values (the ring's slots are rewritten every frame, a
 * replayed bundle reads them later); false for other stencil values (the reference is the pass's)
 */
export interface ClearDraw {
  draw: (pass: GPURenderPassEncoder, target: PassTarget, values: ClearValues) => void;
  drawInBundle: (bundle: GPURenderBundleEncoder, target: PassTarget, values: ClearValues) => boolean;
}

export function createClearDraw(device: GPUDevice, ring: UniformRing): ClearDraw {
  const bindGroupLayout = device.createBindGroupLayout({
    entries: [
      {
        binding: 0,
        visibility: ShaderStage.VERTEX | ShaderStage.FRAGMENT,
        buffer: { type: 'uniform', hasDynamicOffset: true, minBindingSize: UniformSize },
      },
    ],
  });
  const pipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [bindGroupLayout] });
  const modules = new Map<number, GPUShaderModule>();
  const pipelines = new Map<string, GPURenderPipeline>();
  const bindGroups = new Map<GPUBuffer, GPUBindGroup>();
  const snapshot = createUniformSnapshot(ring, UniformSize);
  const data = new Float32Array(snapshot.bytes.buffer as ArrayBuffer);

  function getModule(colorCount: number) {
    let module = modules.get(colorCount);
    if (!module) {
      module = device.createShaderModule({ code: getClearModule(colorCount), label: 'clear' });
      modules.set(colorCount, module);
    }
    return module;
  }

  function getPipeline(target: PassTarget, values: ClearValues, zeroStencil = false) {
    const hasColor = values.color !== undefined;
    const hasDepth = values.depth !== undefined && !!target.depthFormat;
    const hasStencil = values.stencil !== undefined && target.hasStencil;
    const key = `${target.key}|${hasColor}|${hasDepth}|${hasStencil}|${zeroStencil}`;
    let pipeline = pipelines.get(key);
    if (!pipeline) {
      const module = getModule(target.colorFormats.length);
      const stencilFace: GPUStencilFaceState = {
        compare: 'always',
        failOp: 'keep',
        depthFailOp: 'keep',
        passOp: hasStencil ? (zeroStencil ? 'zero' : 'replace') : 'keep',
      };
      const start = isTimingMode ? now() : 0;
      pipeline = device.createRenderPipeline({
        label: `clear ${key}`,
        layout: pipelineLayout,
        vertex: { module, entryPoint: 'vertexMain' },
        fragment:
          target.colorFormats.length > 0
            ? {
                module,
                entryPoint: 'fragmentMain',
                targets: target.colorFormats.map((format) => {
                  if (format.endsWith('sint') || format.endsWith('uint'))
                    notSupported('partial clears of integer targets');
                  return { format, writeMask: hasColor ? 0xf : 0 };
                }),
              }
            : undefined,
        primitive: { topology: 'triangle-list' },
        depthStencil: target.depthFormat
          ? {
              format: target.depthFormat,
              depthWriteEnabled: hasDepth,
              depthCompare: 'always',
              stencilFront: stencilFace,
              stencilBack: stencilFace,
              stencilReadMask: 0xff,
              stencilWriteMask: hasStencil ? 0xff : 0,
            }
          : undefined,
      });
      if (isTimingMode) console.log(`Created pipeline clear ${key} synchronously in ${(now() - start).toFixed(2)} ms`);
      pipelines.set(key, pipeline);
    }
    return pipeline;
  }

  function setData(values: ClearValues) {
    const { color, depth = 0 } = values;
    data[0] = color?.[0] ?? 0;
    data[1] = color?.[1] ?? 0;
    data[2] = color?.[2] ?? 0;
    data[3] = color?.[3] ?? 0;
    data[4] = depth;
  }

  /** bind groups of constant buffers, per values (a few: the clears of solid interior caps) */
  const constantGroups = new Map<string, GPUBindGroup>();
  function getConstantGroup(values: ClearValues) {
    setData(values);
    const key = Array.from(data.subarray(0, 5)).join(',');
    let group = constantGroups.get(key);
    if (!group) {
      const buffer = device.createBuffer({
        size: UniformSize,
        usage: BufferUsage.UNIFORM | BufferUsage.COPY_DST,
        label: 'clear values',
      });
      device.queue.writeBuffer(buffer, 0, data.buffer, data.byteOffset, UniformSize);
      group = device.createBindGroup({
        layout: bindGroupLayout,
        entries: [{ binding: 0, resource: { buffer, offset: 0, size: UniformSize } }],
      });
      constantGroups.set(key, group);
    }
    return group;
  }

  const drawInBundle = (bundle: GPURenderBundleEncoder, target: PassTarget, values: ClearValues) => {
    if (values.stencil !== undefined && values.stencil !== 0 && target.hasStencil) return false;
    bundle.setPipeline(getPipeline(target, values, true));
    bundle.setBindGroup(0, getConstantGroup(values), [0]);
    bundle.draw(3);
    return true;
  };

  const draw = (pass: GPURenderPassEncoder, target: PassTarget, values: ClearValues) => {
    const { stencil = 0 } = values;
    setData(values);
    snapshot.markChanged();
    const slot = snapshot.getSlot();
    let bindGroup = bindGroups.get(slot.buffer);
    if (!bindGroup) {
      bindGroup = device.createBindGroup({
        layout: bindGroupLayout,
        entries: [{ binding: 0, resource: { buffer: slot.buffer, offset: 0, size: UniformSize } }],
      });
      bindGroups.set(slot.buffer, bindGroup);
    }
    pass.setPipeline(getPipeline(target, values));
    pass.setBindGroup(0, bindGroup, [slot.offset]);
    pass.setViewport(0, 0, target.width, target.height, 0, 1);
    pass.setStencilReference(stencil);
    pass.draw(3);
  };

  return { draw, drawInBundle };
}
