/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Textures } from './texture.js';
import type { UniformsList } from './uniform.js';
import type { StorageBuffer } from './compute.js';

export type DrawMode =
  | 'points'
  | 'lines'
  | 'line-strip'
  | 'line-loop'
  | 'triangles'
  | 'triangle-strip'
  | 'triangle-fan';

export type Transparency = 'blended' | 'wboit' | 'dpoit';

const GraphicsRenderVariant = { color: '', pick: '', depth: '', marking: '', emissive: '', tracing: '' };
export type GraphicsRenderVariant = keyof typeof GraphicsRenderVariant;
export const GraphicsRenderVariants = Object.keys(GraphicsRenderVariant) as GraphicsRenderVariant[];

const ComputeRenderVariant = { compute: '' };
export type ComputeRenderVariant = keyof typeof ComputeRenderVariant;
export const ComputeRenderVariants = Object.keys(ComputeRenderVariant) as ComputeRenderVariant[];

export type ProgramVariant = GraphicsRenderVariant | ComputeRenderVariant;

/** Program (pipeline) of a render item variant, compiled and linked asynchronously */
export interface GraphicsProgram {
  readonly id: number;
  readonly variant: ProgramVariant;

  isReady(): boolean;
  link(): void;
  finalize(force?: boolean): boolean;
}

export type MultiDrawBaseData = {
  /** Only used for `multiDrawArraysInstancedBaseInstance` */
  firsts: Int32Array;
  counts: Int32Array;
  /** Only used for `multiDrawElementsInstancedBaseVertexBaseInstance` */
  offsets: Int32Array;
  instanceCounts: Int32Array;
  /** Only used for `multiDrawElementsInstancedBaseVertexBaseInstance` */
  baseVertices: Int32Array;
  baseInstances: Uint32Array;
  count: number;
  uniforms: UniformsList;
  /**
   * The draws' arguments in a buffer, written on the GPU (GPU culling), in place of the arrays: per draw (`count`) a
   * record of `IndirectRecordBytes`; first instances other than 0 need `indirectInstanceRanges` (WebGPU only)
   */
  indirect?: MultiDrawIndirect;
};

/**
 * Indirect arguments of a list's draws: u32 index count, instance count, first index, base vertex, first instance
 * (indexed draws), or vertex count, instance count, first vertex, first instance
 */
export type MultiDrawIndirect = {
  readonly buffer: StorageBuffer;
  readonly byteOffset: number;
  /** whether the records are indexed (`drawIndexedIndirect`); undefined: where the item has elements */
  readonly indexed?: boolean;
};

/** Bytes of an indirect argument record, indexed (`drawIndexedIndirect`) or not (`drawIndirect`) */
export const IndirectRecordBytes = { indexed: 20, vertex: 16 } as const;

/**
 * Uniforms and textures shared by the render items of a pass (the renderer globals), applied by the
 * render items before drawing; a uniform buffer plus classic uniforms on WebGL, a bind group on WebGPU
 */
export interface GlobalBindings {
  /** Bound before the textures of the render items, may be changed followed by `invalidate` */
  readonly textures: Textures;
  /** Marks uniform values or textures as changed, they are applied again before the next draw */
  invalidate: () => void;
  /**
   * The solid interior cap the following draws are of: `pass` 2 marks, 1 fills, 0 none (the object itself); `clip` the
   * clip object cut, -1 the near plane. A draw state, not a change of the globals: WebGPU binds a draw block per cap
   * (constant, so that recorded draws replay), WebGL sets the uniforms
   */
  setSolidInteriorCap: (pass: number, clip: number) => void;
  destroy: () => void;
}

/**
 * The draws of a pass recorded for replay in later frames (WebGPU render bundles, `GraphicsContext.createDrawList`):
 * whether the pass would draw the recorded draws the caller decides, the list whether they bind the same globals
 */
export interface DrawList {
  /** Replays the recording into the open pass; false without one or when `globals` bind other than at recording */
  replay: (globals: GlobalBindings) => boolean;
  /**
   * Draws `draw` into the open pass, recorded for `replay` where it can be (else drawn directly, nothing recorded); a
   * list may keep a recording per target layout of the pass
   */
  record: (globals: GlobalBindings, draw: () => void) => void;
  /** Forgets the recordings (the draws changed) */
  reset: () => void;
}

export interface RenderItem<T extends string> {
  readonly id: number;
  readonly materialId: number;
  /** GPU bytes of the item, a shared texture by the item's share (`.v6/plans/gl/shared-textures.md` §5) */
  getByteCount(): number;
  getProgram: (variant: T) => GraphicsProgram;
  /**
   * Names of the bound textures with linear filtering (WebGPU programs read them filtered, `dLinear_*`), for the
   * signature of the item's object (`ObjectSignature`)
   */
  getLinearTextures?: () => readonly string[];
  setTransparency: (transparency: Transparency) => void;

  render: (variant: T, globals?: GlobalBindings, mdbDataList?: MultiDrawBaseData[]) => void;
  update: () => void;
  destroy: () => void;
}

export type GraphicsRenderItem = RenderItem<GraphicsRenderVariant>;

/**
 * The programs of a render object without its data (e.g. from an `ObjectSignature`), for preparing pipelines: requests
 * (`PipelineRequest`) take it in place of a render item; its programs come from the same caches as render items'
 * (`GraphicsContext.createProgramItem`)
 */
export interface ProgramItem {
  readonly id: number;
  getProgram: (variant: GraphicsRenderVariant) => GraphicsProgram;
  /** releases its programs (render items using them keep them) */
  destroy: () => void;
}
export type ComputeRenderItem = RenderItem<ComputeRenderVariant>;
