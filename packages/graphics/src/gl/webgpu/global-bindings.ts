/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import type { GlobalBindings } from '../render-item.js';
import type { RenderableSchema, TextureDimension, TextureSampleType, UniformSpec } from '../renderable/schema.js';
import type { Textures } from '../texture.js';
import type { UniformKind, UniformsList } from '../uniform.js';
import {
  createUniformBlockData,
  createUniformBlockWriter,
  getUniformBlockLayout,
  type TextureBinding,
} from '../uniform-layout.js';
import { TextureUsage } from './constants.js';
import type { WebGPUProgram } from './program.js';
import { asWebGPUTexture, type Release } from './texture.js';
import { createUniformSnapshot, type UniformRing, type UniformSnapshot } from './uniform-ring.js';
import { WgslSamplers, type WgslSampler } from './wgsl-declarations.js';

/** clip z to `[0, 1]`: `z' = (z + w) / 2` */
const DepthRange = Mat4.fromArray(Mat4(), [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 1], 0);
const InvDepthRange = Mat4.invert(Mat4(), DepthRange);

/** Matrices that include the projection, for WebGPU mapped to the `[0, 1]` depth range */
const ProjectionConversions: { [k: string]: (out: Mat4, m: Mat4) => Mat4 } = {
  uProjection: (out, m) => Mat4.mul(out, DepthRange, m),
  uModelViewProjection: (out, m) => Mat4.mul(out, DepthRange, m),
  uInvProjection: (out, m) => Mat4.mul(out, m, InvDepthRange),
  uInvModelViewProjection: (out, m) => Mat4.mul(out, m, InvDepthRange),
};

/** A cell with the converted value of `cell`, recomputed when its version changes */
function convertedCell(cell: ValueCell<Mat4>, convert: (out: Mat4, m: Mat4) => Mat4): ValueCell<Mat4> {
  const ref = { version: -1, value: Mat4() };
  return {
    get ref() {
      if (cell.ref.version !== ref.version) {
        convert(ref.value, cell.ref.value);
        ref.version = cell.ref.version;
      }
      return ref;
    },
  } as unknown as ValueCell<Mat4>;
}

export interface WebGPUBindingContext {
  readonly device: GPUDevice;
  readonly release: Release;
  readonly ring: UniformRing;
  readonly samplers: { readonly [k in WgslSampler]: GPUSampler };
}

export interface WebGPUGlobalBindings extends GlobalBindings {
  /** 'draw' uniforms of the global schema, written per draw by the render items */
  readonly drawUniforms: UniformsList;
  /** bind group 0 for `program` and its dynamic offset */
  getBindGroup: (program: WebGPUProgram) => { group: GPUBindGroup; offset: number };
}

const GlobalBindingsSet = new WeakSet<GlobalBindings>();

export function asWebGPUGlobalBindings(globals: GlobalBindings): WebGPUGlobalBindings {
  if (!GlobalBindingsSet.has(globals)) throw new Error('expected WebGPU global bindings');
  return globals as WebGPUGlobalBindings;
}

const PlaceholderFormats: { [k in TextureSampleType]: GPUTextureFormat } = {
  float: 'rgba8unorm',
  'unfilterable-float': 'rgba8unorm',
  sint: 'rgba8sint',
  uint: 'rgba8uint',
  depth: 'depth32float',
};

const placeholders = new WeakMap<GPUDevice, Map<string, GPUTextureView>>();

/**
 * A 1x1 texture of zeros for bindings whose texture has no storage yet (e.g. images still loading), as sampling an
 * undefined texture reads zero on WebGL
 */
function getPlaceholderView(device: GPUDevice, sampleType: TextureSampleType, dimension: TextureDimension) {
  let views = placeholders.get(device);
  if (!views) placeholders.set(device, (views = new Map()));
  const key = `${sampleType}:${dimension}`;
  let view = views.get(key);
  if (!view) {
    const texture = device.createTexture({
      size: { width: 1, height: 1, depthOrArrayLayers: dimension === 'cube' ? 6 : 1 },
      dimension: dimension === '3d' ? '3d' : '2d',
      format: PlaceholderFormats[sampleType],
      usage: TextureUsage.TEXTURE_BINDING,
      label: `placeholder ${key}`,
    });
    view = texture.createView({ dimension });
    views.set(key, view);
  }
  return view;
}

/**
 * Entries for the textures of `bindings` (from the layout); textures that are not set or have no storage get a
 * placeholder
 */
export function getTextureEntries(
  device: GPUDevice,
  textures: Textures,
  bindings: readonly TextureBinding[],
  referenced: ReadonlySet<string>,
): { entries: GPUBindGroupEntry[]; key: string } {
  const entries: GPUBindGroupEntry[] = [];
  let key = '';
  for (const { name, binding, sampleType, dimension } of bindings) {
    if (!referenced.has(name)) continue;
    const t = textures.find((e) => e[0] === name)?.[1];
    if (!t) {
      // not bound yet (e.g. the dpoit textures in the opaque pass), reads zero as an unbound sampler on WebGL
      entries.push({ binding, resource: getPlaceholderView(device, sampleType, dimension) });
      key += `${name}:unset,`;
      continue;
    }
    const texture = asWebGPUTexture(t);
    const view = texture.getView();
    if (view) {
      entries.push({ binding, resource: view });
      key += `${texture.id}:${texture.version},`;
    } else {
      entries.push({ binding, resource: getPlaceholderView(device, sampleType, dimension) });
      key += `${texture.id}:none,`;
    }
  }
  return { entries, key };
}

export function createWebGPUGlobalBindings(
  ctx: WebGPUBindingContext,
  schema: RenderableSchema,
  uniforms: UniformsList,
  textures: Textures,
): WebGPUGlobalBindings {
  const layout = getUniformBlockLayout('global', schema);
  if (!layout) throw new Error('global bindings need uniforms of the global block');

  const blockUniforms: UniformsList = [];
  const drawUniforms: UniformsList = [];
  for (const [k, v] of uniforms) {
    const spec = schema[k] as UniformSpec<UniformKind> | undefined;
    if (spec?.block === 'draw') drawUniforms.push([k, v]);
    else if (spec?.block === 'global') {
      const convert = ProjectionConversions[k];
      blockUniforms.push([k, convert ? convertedCell(v as ValueCell<Mat4>, convert) : v]);
    }
  }

  const snapshot: UniformSnapshot = createUniformSnapshot(ctx.ring, layout.size);
  const data = createUniformBlockData(snapshot.bytes.buffer as ArrayBuffer, 0, layout.size);
  const write = createUniformBlockWriter(layout, data, blockUniforms, snapshot.markChanged);
  write();

  type Cached = { group: GPUBindGroup; buffer: GPUBuffer; key: string };
  const groups = new Map<GPUBindGroupLayout, Cached>();

  const bindings: WebGPUGlobalBindings = {
    textures,
    drawUniforms,
    invalidate: () => {
      write();
    },
    getBindGroup: (program: WebGPUProgram) => {
      write();
      const slot = snapshot.getSlot();
      const groupLayout = program.bindGroupLayouts[0];
      const tex = getTextureEntries(ctx.device, textures, program.bindingLayout.global.textures, program.referenced);
      let cached = groups.get(groupLayout);
      if (!cached || cached.buffer !== slot.buffer || cached.key !== tex.key) {
        const entries: GPUBindGroupEntry[] = [
          { binding: 0, resource: { buffer: slot.buffer, offset: 0, size: layout.size } },
          ...tex.entries,
        ];
        const samplerBinding = program.bindingLayout.global.textures.length + 1;
        for (let i = 0; i < WgslSamplers.length; ++i) {
          const name = WgslSamplers[i];
          if (program.referenced.has(name)) {
            entries.push({ binding: samplerBinding + i, resource: ctx.samplers[name] });
          }
        }
        cached = {
          group: ctx.device.createBindGroup({ layout: groupLayout, entries }),
          buffer: slot.buffer,
          key: tex.key,
        };
        groups.set(groupLayout, cached);
      }
      return { group: cached.group, offset: slot.offset };
    },
    destroy: () => {
      groups.clear();
    },
  };
  GlobalBindingsSet.add(bindings);
  return bindings;
}
