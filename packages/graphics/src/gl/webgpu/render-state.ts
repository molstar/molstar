/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { DrawMode } from '../render-item.js';
import type { ColorWriteMask, RenderState, StencilFaceState } from '../render-state.js';

/** ids per state object (fast path) and per state content */
const stateIds = new WeakMap<RenderState, number>();
const contentIds = new Map<string, number>();

/** The fields of `value` in sorted key order, without the stencil reference (dynamic state of the pass) */
function getStateContent(value: unknown, key = ''): string {
  if (value === null || typeof value !== 'object') return JSON.stringify(value) ?? '';
  if (Array.isArray(value)) return `[${value.map((v) => getStateContent(v)).join(',')}]`;
  const keys = Object.keys(value)
    .filter((k) => !(key === 'stencil' && k === 'reference'))
    .sort();
  return `{${keys.map((k) => `${k}:${getStateContent((value as Record<string, unknown>)[k], k)}`).join(',')}}`;
}

/**
 * Id of a render state for pipeline keys, equal for states with equal content (the renderer also builds states per
 * draw, e.g. for the stencil passes of solid interiors)
 */
export function getRenderStateId(state: RenderState) {
  let id = stateIds.get(state);
  if (id === undefined) {
    const content = getStateContent(state);
    id = contentIds.get(content);
    if (id === undefined) {
      id = contentIds.size;
      contentIds.set(content, id);
    }
    stateIds.set(state, id);
  }
  return id;
}

export function getPrimitiveTopology(drawMode: DrawMode): GPUPrimitiveTopology {
  switch (drawMode) {
    case 'points':
      return 'point-list';
    case 'lines':
      return 'line-list';
    case 'line-strip':
      return 'line-strip';
    case 'triangles':
      return 'triangle-list';
    case 'triangle-strip':
      return 'triangle-strip';
    default:
      throw new Error(`draw mode '${drawMode}' is not supported by WebGPU`);
  }
}

function getWriteMask(mask: ColorWriteMask | undefined): GPUColorWriteFlags {
  if (!mask) return 0xf;
  return (mask[0] ? 0x1 : 0) | (mask[1] ? 0x2 : 0) | (mask[2] ? 0x4 : 0) | (mask[3] ? 0x8 : 0);
}

export function getColorTargetState(state: RenderState, format: GPUTextureFormat): GPUColorTargetState {
  const { blend } = state;
  return {
    format,
    blend: blend ? { color: { ...blend.color }, alpha: { ...blend.alpha } } : undefined,
    writeMask: getWriteMask(state.colorWrite),
  };
}

function getStencilFaceState(face: StencilFaceState): GPUStencilFaceState {
  return { compare: face.compare, failOp: face.failOp, depthFailOp: face.depthFailOp, passOp: face.passOp };
}

/** A pass with a depth attachment needs a depth state; an unset depth state means no test and no write */
export function getDepthStencilState(
  state: RenderState,
  format: GPUTextureFormat,
  hasStencil: boolean,
): GPUDepthStencilState {
  const { depth, stencil } = state;
  const result: GPUDepthStencilState = {
    format,
    depthWriteEnabled: depth?.write ?? false,
    depthCompare: depth?.compare ?? 'always',
  };
  if (stencil && hasStencil) {
    result.stencilFront = getStencilFaceState(stencil.front);
    result.stencilBack = getStencilFaceState(stencil.back);
    result.stencilReadMask = stencil.readMask;
    result.stencilWriteMask = stencil.writeMask;
  }
  return result;
}

export function getPrimitiveState(state: RenderState, drawMode: DrawMode): GPUPrimitiveState {
  return {
    topology: getPrimitiveTopology(drawMode),
    cullMode: state.cull ?? 'none',
    frontFace: state.frontFace ?? 'ccw',
  };
}
