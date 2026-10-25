/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GlobalBindings } from '../render-item.js';
import type { RenderableSchema, TextureDimension, TextureSampleType, UniformSpec } from '../renderable/schema.js';
import type { Texture, Textures } from '../texture.js';
import type { UniformKind, UniformsList } from '../uniform.js';
import {
  createUniformBlockData,
  createUniformBlockWriter,
  getUniformBlockLayout,
  type TextureBinding,
} from '../uniform-layout.js';
import { TextureUsage } from './constants.js';
import { getDepthRangeCell, ReversedZeroToOne } from '../depth-range.js';
import type { WebGPUProgram } from './program.js';
import { asWebGPUTexture, type Release, type WebGPUTexture } from './texture.js';
import { createUniformSnapshot, type UniformRing, type UniformSnapshot } from './uniform-ring.js';
import { WgslSamplers, type WgslSampler } from './wgsl-declarations.js';

export interface WebGPUBindingContext {
  readonly device: GPUDevice;
  readonly release: Release;
  readonly ring: UniformRing;
  /**
   * The ring of the global blocks only: their snapshots of a submission take offsets in the order the frame's passes
   * change the globals, so a frame of the same passes binds the same offsets (what a recorded pass relies on), whatever
   * other blocks and compute kernels write in between
   */
  readonly globalRing: UniformRing;
  readonly samplers: { readonly [k in WgslSampler]: GPUSampler };
}

export interface WebGPUGlobalBindings extends GlobalBindings {
  /** 'draw' uniforms of the global schema, written per draw by the render items */
  readonly drawUniforms: UniformsList;
  /** the cap the following draws are of (`setSolidInteriorCap`), their draw blocks' values in place of the globals' */
  readonly solidInteriorCap: { readonly pass: number; readonly clip: number };
  /** bind group 0 for `program` and its dynamic offset */
  getBindGroup: (program: WebGPUProgram) => { group: GPUBindGroup; offset: number };
  /**
   * What the draws of the current pass bind, the block's slot taken as by a draw (so that the slots of later passes
   * stay where a frame of the same passes has them): compared by a draw list before replaying
   */
  getState: () => { readonly buffer: GPUBuffer; readonly offset: number; readonly textures: string };
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

/** global bind groups kept per layout, for passes alternating between texture sets (see `createWebGPUGlobalBindings`) */
const MaxGlobalGroupsPerLayout = 8;

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
      // the projections mapped to reversed `[0, 1]` depth
      blockUniforms.push([k, getDepthRangeCell(schema, k, v, ReversedZeroToOne)]);
    }
  }

  const snapshot: UniformSnapshot = createUniformSnapshot(ctx.globalRing, layout.size);
  const data = createUniformBlockData(snapshot.bytes.buffer as ArrayBuffer, 0, layout.size);
  const write = createUniformBlockWriter(layout, data, blockUniforms, snapshot.markChanged);
  write();

  type Cached = { group: GPUBindGroup; buffer: GPUBuffer };
  /**
   * per layout, the groups by their textures (`getTextureEntries`): passes alternating between texture sets (the
   * marking passes, depth and mask) keep binding the same groups, which recordings compare by identity
   * (`GraphicsContext.captureDraws`); the most recent `MaxGlobalGroupsPerLayout`
   */
  const groups = new Map<GPUBindGroupLayout, Map<string, Cached>>();

  /**
   * per layout the group last returned, with what it was resolved from: a draw compares those (the uniform buffer, the
   * textures list's generation, the textures' views) instead of building the textures' entries and key
   */
  type Last = {
    group: GPUBindGroup;
    buffer: GPUBuffer;
    generation: number;
    textures: (WebGPUTexture | undefined)[];
    views: (GPUTextureView | undefined)[];
  };
  const lastGroups = new Map<GPUBindGroupLayout, Last>();
  /** the textures list's entries as last seen; `generation` moves when one changes (another texture, name or count) */
  const seenNames: string[] = [];
  const seenTextures: Texture[] = [];
  let generation = 0;
  const syncTextures = () => {
    let changed = seenTextures.length !== textures.length;
    for (let i = 0, il = textures.length; i < il; ++i) {
      const [name, texture] = textures[i];
      if (seenNames[i] !== name || seenTextures[i] !== texture) {
        seenNames[i] = name;
        seenTextures[i] = texture;
        changed = true;
      }
    }
    if (changed) {
      seenNames.length = seenTextures.length = textures.length;
      generation += 1;
    }
  };
  const isLastCurrent = (last: Last, buffer: GPUBuffer) => {
    if (last.buffer !== buffer || last.generation !== generation) return false;
    const { textures, views } = last;
    for (let i = 0, il = textures.length; i < il; ++i) {
      if (textures[i]?.getView() !== views[i]) return false;
    }
    return true;
  };
  /** the result of `getBindGroup`, one object (read before the next call) */
  const result = { group: undefined as unknown as GPUBindGroup, offset: 0 };

  const solidInteriorCap = { pass: 0, clip: -1 };

  const bindings: WebGPUGlobalBindings = {
    textures,
    drawUniforms,
    solidInteriorCap,
    invalidate: () => {
      write();
    },
    setSolidInteriorCap: (pass: number, clip: number) => {
      solidInteriorCap.pass = pass;
      solidInteriorCap.clip = clip;
    },
    getState: () => {
      write();
      const slot = snapshot.getSlot();
      let key = '';
      for (const [name, t] of textures) {
        const texture = asWebGPUTexture(t);
        key += `${name}:${texture.id}:${texture.version},`;
      }
      return { buffer: slot.buffer, offset: slot.offset, textures: key };
    },
    getBindGroup: (program: WebGPUProgram) => {
      write();
      const slot = snapshot.getSlot();
      const groupLayout = program.bindGroupLayouts[0];
      syncTextures();
      const last = lastGroups.get(groupLayout);
      if (last && isLastCurrent(last, slot.buffer)) {
        result.group = last.group;
        result.offset = slot.offset;
        return result;
      }
      const tex = getTextureEntries(ctx.device, textures, program.bindingLayout.global.textures, program.referenced);
      let byTextures = groups.get(groupLayout);
      if (!byTextures) groups.set(groupLayout, (byTextures = new Map()));
      let cached = byTextures.get(tex.key);
      if (!cached || cached.buffer !== slot.buffer) {
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
        cached = { group: ctx.device.createBindGroup({ layout: groupLayout, entries }), buffer: slot.buffer };
        byTextures.delete(tex.key);
        if (byTextures.size >= MaxGlobalGroupsPerLayout) byTextures.delete(byTextures.keys().next().value!);
        byTextures.set(tex.key, cached);
      }
      const resolved: (WebGPUTexture | undefined)[] = [];
      for (const { name } of program.bindingLayout.global.textures) {
        if (!program.referenced.has(name)) continue;
        const t = textures.find((e) => e[0] === name)?.[1];
        resolved.push(t ? asWebGPUTexture(t) : undefined);
      }
      lastGroups.set(groupLayout, {
        group: cached.group,
        buffer: slot.buffer,
        generation,
        textures: resolved,
        views: resolved.map((t) => t?.getView()),
      });
      result.group = cached.group;
      result.offset = slot.offset;
      return result;
    },
    destroy: () => {
      write.dispose();
      groups.clear();
      lastGroups.clear();
    },
  };
  GlobalBindingsSet.add(bindings);
  return bindings;
}
