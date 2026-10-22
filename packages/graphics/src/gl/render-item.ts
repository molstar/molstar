/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Textures } from './texture.js';
import type { UniformsList } from './uniform.js';

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
};

/**
 * Uniforms and textures shared by the render items of a pass (the renderer globals), applied by the
 * render items before drawing; a uniform buffer plus classic uniforms on WebGL, a bind group on WebGPU
 */
export interface GlobalBindings {
  /** Bound before the textures of the render items, may be changed followed by `invalidate` */
  readonly textures: Textures;
  /** Marks uniform values or textures as changed, they are applied again before the next draw */
  invalidate: () => void;
  destroy: () => void;
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
