/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { TextureReadArray, TextureReader, TextureReaderOptions, TextureReadFormat } from './texture-reader.js';
import { isDebugMode } from '@molstar/core/util/debug';
import type { ComputeContext } from './compute.js';
import type { Subject } from 'rxjs';
import type { Mat4 } from '@molstar/core/math/linear-algebra';
import type { Viewport } from '@molstar/core/util/viewport';
import type { ValueCell } from '@molstar/core/util';
import type { PendingWork } from './pending-work.js';
import type { now } from '@molstar/core/util/now';
import type { GraphicsCapabilities } from './capabilities.js';
import type { Framebuffer } from './framebuffer.js';
import type { RenderTarget } from './render-target.js';
import type {
  ComputeRenderItem,
  DrawMode,
  DrawList,
  GlobalBindings,
  GraphicsRenderItem,
  ProgramItem,
  ProgramVariant,
  Transparency,
} from './render-item.js';
import type { ComputeRenderable } from './renderable.js';
import type { RenderableSchema, RenderableValues } from './renderable/schema.js';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from './renderbuffer.js';
import type { ClearValues, RenderState } from './render-state.js';
import type { DefineValues, ShaderCode } from './shader-code.js';
import type { GraphicsStats, ResourceByteCounts } from './stats.js';
import type {
  CubeFaces,
  Texture,
  TextureFilter,
  TextureFormat,
  TextureKind,
  Textures,
  TextureType,
} from './texture.js';
import type { GraphicsTimer } from './timer.js';
import type { UniformBuffer, UniformsList } from './uniform.js';
import type { UniformBlockLayout } from './uniform-layout.js';

export interface GraphicsLimits {
  readonly maxTextureSize: number;
  /** 0 without `texture3d` */
  readonly max3dTextureSize: number;
  readonly maxRenderbufferSize: number;
  /** 0 without `multipleRenderTargets` */
  readonly maxDrawBuffers: number;
  readonly maxTextureImageUnits: number;
  /** The most quad items a draw of the shared quad indices reaches (`quad-indices.ts`): past it quads draw unindexed */
  readonly maxQuadItemCount: number;
}

export interface GraphicsResources {
  readonly limits: GraphicsLimits;

  texture: (kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter) => Texture;
  /** Whether `texture` was created by `texture` of these resources (not by another context, e.g. before a backend switch) */
  ownsTexture: (texture: Texture) => boolean;
  /** Loads the faces, tracked in `pending`; `onload` reports when they are loaded or failed to */
  cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: (errored?: boolean) => void) => Texture;
  framebuffer: () => Framebuffer;
  renderbuffer: (
    format: RenderbufferFormat,
    attachment: RenderbufferAttachment,
    width: number,
    height: number,
  ) => Renderbuffer;
  uniformBuffer: (layout: UniformBlockLayout) => UniformBuffer;
  /** writes a new snapshot on every change instead of updating in place */
  uniformRingBuffer: (layout: UniformBlockLayout) => UniformBuffer;
  /** cached per `materialId`, shared by render items with the same material */
  materialUniformBuffer: (materialId: number, layout: UniformBlockLayout) => UniformBuffer;

  getByteCounts: () => ResourceByteCounts;

  linkPrograms: (variants?: ProgramVariant[]) => void;
  finalizePrograms: (variants?: ProgramVariant[], isSynchronous?: boolean) => boolean;
  /**
   * Makes the pipelines of `requests` ready without drawing (WebGPU: compiles them in the background; WebGL: links the
   * programs, state and target don't select programs there) and resolves when they and other pending work are ready.
   * The program of a request is the one its item selects at the call (its current defines and textures).
   */
  preparePipelines: (requests: Iterable<PipelineRequest>) => Promise<void>;

  /** Asynchronous work frames may depend on; `ImagePass` waits for it, `Canvas3D` redraws when it settles */
  readonly pending: PendingWork;
  /**
   * Work that did not get the GPU memory it needed, with what it left out (e.g. surfaces drawn incomplete, WebGPU
   * out-of-memory errors), for the user; `Canvas3D.outOfMemory` forwards it
   */
  readonly outOfMemory: Subject<string>;
  /**
   * Loads the image at `source` (a URL or data URL) into the rgba `ubyte` `texture` with the raw pixel values (no
   * color space conversion), then updates `cell`; tracked in `pending`
   */
  loadImage: (texture: Texture, source: string, cell?: ValueCell<Texture>, label?: string) => Promise<void>;
}

export type RenderTargetDepthStencil = 'none' | 'depth' | 'depth-stencil';
export type RenderTargetType = 'uint8' | 'float32' | 'fp16';
export type RenderTargetFormat = 'rgba' | 'r';

/** Target of a render pass, `drawing-buffer` is the canvas or XR layer */
export type RenderPassTarget = RenderTarget | Framebuffer | 'drawing-buffer';

/**
 * What pipelines depend on of a render pass target: the formats of its attachments (`GraphicsContext.describeTarget`);
 * targets with equal keys draw with the same pipelines
 */
export interface TargetDescriptor {
  readonly key: string;
  /** whether the target has a stencil attachment (solid interior caps draw only into those) */
  readonly hasStencil: boolean;
  /** depth and stencil only read in the pass (`RenderPassOptions.depthReadOnly`) */
  readonly depthReadOnly?: boolean;
}

/**
 * A pipeline a draw needs: the program of `item` for `variant`, drawn with `state` into a target like `target`; render
 * items of full-screen passes are compute items drawn with the `compute` variant
 */
export interface PipelineRequest {
  readonly item: GraphicsRenderItem | ComputeRenderItem | ProgramItem;
  readonly variant: ProgramVariant;
  readonly state: RenderState;
  readonly target: TargetDescriptor;
}

export interface RenderPassOptions {
  /** Clears the attachments at the start of the pass, like `clear` right after `beginRenderPass` (current scissor) */
  readonly clear?: ClearValues;
  /** Name of the pass in timer results, derived from the target if not given */
  readonly label?: string;
  /**
   * The depth (and stencil) attachment is only read in the pass: draws test against it but don't write it (their
   * render states' writes are masked) and it is not cleared; on WebGPU its texture can also be bound for reading in the
   * pass (a writable attachment can't), on WebGL it must not be (a feedback loop)
   */
  readonly depthReadOnly?: boolean;
}

/** How a view of an XR viewer pose is rendered, see `GraphicsXR.getView` */
export interface GraphicsXRView {
  /** the rectangle of the view in the drawing buffer, top-left origin */
  readonly viewport: Viewport;
  /** the projection of the view with GL's clip depth (`[-1, 1]`), as the camera builds its projections */
  readonly projection: Mat4;
}

/**
 * XR sessions rendered by the context into a layer of the session (WebGL: `XRWebGLLayer`; WebGPU: a projection layer
 * of `XRGPUBinding`), whose views share the drawing buffer
 */
export interface GraphicsXR {
  readonly session?: XRSession;
  readonly changed: Subject<void>;
  set: (session: XRSession | undefined, options?: { resolutionScale?: number }) => Promise<void>;
  end: () => Promise<void>;
  /** Whether the session has the layer the context renders into */
  isPresenting: () => boolean;
  /** The viewport and projection of `view` (of a pose of the session) into `out`, false without the layer */
  getView: (view: XRView, out: GraphicsXRView) => boolean;
}

/** Backend-neutral graphics context, implemented by `WebGLContext` */
/** Texture contents read back, four values per texel, see `GraphicsContext.readTextureData` */
export interface TextureData {
  readonly array: Uint8Array | Float32Array;
  readonly width: number;
  readonly height: number;
}

export interface GraphicsContext {
  readonly capabilities: GraphicsCapabilities;
  readonly resources: GraphicsResources;
  readonly stats: GraphicsStats;
  readonly timer: GraphicsTimer;
  /** Compute shaders, where available (`capabilities.compute`) */
  readonly compute?: ComputeContext;

  readonly pixelRatio: number;
  setPixelScale: (value: number) => void;

  readonly isContextLost: boolean;
  /** Emits when the context is lost (`setContextLost`), e.g. on a WebGL context loss event or a lost WebGPU device */
  readonly contextLost: Subject<now.Timestamp>;
  setContextLost: () => void;
  /**
   * After a loss: a new context of the same backend and options on the same canvas, once one can be created (WebGL when
   * the platform restored the context of the canvas, WebGPU with a new device), to replace this one (then destroyed with
   * `lost`); undefined if none can be created or this one was destroyed meanwhile
   */
  recreate: () => Promise<GraphicsContext | undefined>;
  /**
   * Loses the context as the platform would, for debugging; replaced with `recreate` (on WebGL the context of the canvas
   * is restored after a second). Returns false where it can't (WebGL without `WEBGL_lose_context`, a lost context).
   */
  simulateContextLoss: () => boolean;
  /**
   * Loses the context on purpose, which frees its GPU resources (e.g. a canvas out of view), reported as a loss
   * (`contextLost`); `recreate` brings it back, not the platform (WebGL `WEBGL_lose_context`, WebGPU the device
   * destroyed). Returns false where it can't (WebGL without `WEBGL_lose_context`, a lost context).
   */
  release: () => boolean;

  readonly xr: GraphicsXR;

  /** Cache for frambuffers, managed by consumers */
  readonly namedFramebuffers: { [name: string]: Framebuffer };
  /** Cache for textures, managed by consumers */
  readonly namedTextures: { [name: string]: Texture };
  /** Cache for compute renderables, managed by consumers */
  readonly namedComputeRenderables: { [name: string]: ComputeRenderable<any> };

  /** Whether the backend implements `shaderCode`; passes whose programs it doesn't are unsupported and fall back */
  isProgramSupported: (shaderCode: ShaderCode) => boolean;
  createGraphicsRenderItem: (
    drawMode: DrawMode,
    shaderCode: ShaderCode,
    schema: RenderableSchema,
    values: RenderableValues,
    materialId: number,
    transparency: Transparency,
  ) => GraphicsRenderItem;
  createComputeRenderItem: (
    drawMode: DrawMode,
    shaderCode: ShaderCode,
    schema: RenderableSchema,
    values: RenderableValues,
    materialId?: number,
  ) => ComputeRenderItem;
  /**
   * The programs of a render item of `shaderCode` and `schema` with `defineValues` (the item's, without the render
   * variant) and `linearTextures` (`RenderItem.getLinearTextures`), drawn with `transparency`, without its data
   */
  createProgramItem: (
    drawMode: DrawMode,
    shaderCode: ShaderCode,
    schema: RenderableSchema,
    defineValues: DefineValues,
    linearTextures: readonly string[],
    transparency: Transparency,
  ) => ProgramItem;
  /** Uniforms (`schema` and `uniforms`) and textures shared by the render items of a pass */
  createGlobalBindings: (schema: RenderableSchema, uniforms: UniformsList, textures: Textures) => GlobalBindings;
  /** Draw lists for replaying a pass's draws (WebGPU render bundles); not where draws can't be recorded (WebGL) */
  createDrawList?: () => DrawList;
  /**
   * Advances with each submission of the queued commands, where the per-submission uniform rings restart: the n-th
   * global block written after a submission gets the same slot each time, what a draw list's recording binds (with
   * draw lists, `createDrawList`)
   */
  readonly submitCount?: number;
  /**
   * Advances when a resource that recorded draws bind is replaced (the shared quad indices grown, a shared texture's
   * storage replaced by a grown image): recordings made before are not current (with draw lists, `createDrawList`)
   */
  readonly drawResourcesVersion?: number;
  /**
   * Runs `draw` in the open pass without encoding it, its data written as by drawing (uniform slots, vertex data), and
   * returns what its commands bind (pipelines, bind groups with offsets, buffers, draw arguments): equal for the
   * draws of a recording means the recording draws the current data (with draw lists, `createDrawList`)
   */
  captureDraws?: (draw: () => void) => string;

  createRenderTarget: (
    width: number,
    height: number,
    depthStencil?: RenderTargetDepthStencil,
    type?: RenderTargetType,
    filter?: TextureFilter,
    format?: RenderTargetFormat,
  ) => RenderTarget;
  /** Render target of the drawing buffer (canvas or XR layer) */
  createDrawTarget: () => RenderTarget;
  getDrawingBufferSize: () => { width: number; height: number };

  /**
   * Starts a render pass on `target`, ending the open one. Draws and clears need an open pass.
   * Copies, readback and attachment changes end the open pass (WebGPU runs them outside of passes).
   */
  beginRenderPass: (target: RenderPassTarget, options?: RenderPassOptions) => void;
  /** The attachment formats of `target` as currently attached, for `resources.preparePipelines` */
  describeTarget: (target: RenderPassTarget, options?: Pick<RenderPassOptions, 'depthReadOnly'>) => TargetDescriptor;
  endRenderPass: () => void;
  /** Whether the target of the open render pass has a stencil buffer */
  hasRenderPassStencil: () => boolean;

  /** Fixed-function state for the following draws, unset fields are reset to their defaults */
  setRenderState: (state: RenderState) => void;
  /** Used by the `constant` blend factors */
  setBlendConstant: (red: number, green: number, blue: number, alpha: number) => void;
  /**
   * Rectangles of the context (viewport, scissor, readback) are in pixels of the target with a top-left origin (y down),
   * as in the browser; a viewport partly outside the target is cropped
   */
  setViewport: (x: number, y: number, width: number, height: number) => void;
  /** Limits draws and clears to the rectangle (top-left origin), until `disableScissor` */
  setScissor: (x: number, y: number, width: number, height: number) => void;
  /** Draws and clears cover the whole target */
  disableScissor: () => void;
  /**
   * Clears the attachments of the target of the open render pass, within the scissor rectangle and ignoring the
   * write masks. On WebGPU, clears of the whole target before the first draw of a pass are load operations, all
   * others are draws.
   */
  clear: (values: ClearValues) => void;
  /** Ends the open render pass and submits pending commands */
  flush: () => void;
  /** Resolves when the GPU finished the submitted commands */
  waitForGpuCommandsComplete: () => Promise<void>;
  /**
   * Reads the color attachment of `target` as rgba8 in a rectangle with a top-left origin, rows top to bottom.
   * Ends the open render pass.
   */
  readRenderTarget: (target: RenderTarget, x: number, y: number, width: number, height: number) => Promise<Uint8Array>;
  /**
   * Reads an rgba8 2d texture in a rectangle with a top-left origin, rows top to bottom; the read is ordered after
   * the commands recorded so far. Ends the open render pass.
   */
  readTextureAsync: (texture: Texture, x: number, y: number, width: number, height: number) => Promise<Uint8Array>;
  /**
   * Reads a whole 2d color texture in data order (texel row 0 first, as `load` uploads), four values per texel
   * (single-channel formats in `r`): bytes for `ubyte` textures, floats for `float` and `fp16` textures. The read is
   * ordered after the commands recorded so far. Ends the open render pass. For data written by `load` or compute
   * kernels; rendered content is stored with the backend's framebuffer origin (rows flipped between WebGL and
   * WebGPU), read it with `readTextureAsync`.
   */
  readTextureData: (texture: Texture) => Promise<TextureData>;
  /** A reader for repeated reads of the same rectangle of textures (e.g. picking, Hi-Z), see `TextureReader` */
  createTextureReader: <F extends TextureReadFormat = 'rgba8'>(
    options?: TextureReaderOptions<F>,
  ) => TextureReader<TextureReadArray<F>>;

  destroy: (options?: GraphicsDestroyOptions) => void;
}

export interface GraphicsDestroyOptions {
  /** WebGL: keeps the context of the canvas usable instead of forcing its loss (which helps garbage collection) */
  readonly doNotForceWebGLContextLoss?: boolean;
  /**
   * The context was lost and is replaced by a new one on the same canvas: released without GPU calls (its objects
   * belong to the lost context) and the canvas kept usable
   */
  readonly lost?: boolean;
}

/**
 * In debug mode, logs what a new graphics context supports: the neutral part (backend, capabilities, limits) and the
 * backend's own (e.g. WebGL parameters and extensions, WebGPU adapter info and device limits)
 */
export function logGraphicsContextInfo(ctx: GraphicsContext, backendInfo: Record<string, unknown>) {
  if (!isDebugMode) return;
  const { capabilities } = ctx;
  console.log({
    graphics: { backend: capabilities.backend, capabilities, limits: ctx.resources.limits },
    [capabilities.backend]: backendInfo,
  });
}
