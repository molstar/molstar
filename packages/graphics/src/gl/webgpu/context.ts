/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createWebGPUComputeContext } from './compute.js';
import { testBundleFragDepth } from './bundle-frag-depth.js';
import { testViewportOutsideTarget } from './viewport-outside-target.js';
import { Subject } from 'rxjs';
import { isDebugMode } from '@molstar/core/util/debug';
import { fromHalfFloat } from '@molstar/core/util/number-conversion';
import type { TextureReadFormat } from '../texture-reader.js';
import { createWebGPUTextureReader } from './texture-reader.js';
import { now } from '@molstar/core/util/now';
import { isTimingMode } from '@molstar/core/util/debug';
import type { GraphicsCapabilities } from '../capabilities.js';
import {
  type GraphicsContext,
  type GraphicsDestroyOptions,
  type GraphicsResources,
  type GraphicsXR,
  logGraphicsContextInfo,
  type PipelineRequest,
  type RenderPassOptions,
  type RenderPassTarget,
} from '../context.js';
import {
  ComputeRenderVariants,
  type ComputeRenderVariant,
  type DrawList,
  type DrawMode,
  type GlobalBindings,
  type ComputeRenderItem,
  type GraphicsRenderItem,
  GraphicsRenderVariants,
  type GraphicsRenderVariant,
  type ProgramItem,
  type ProgramVariant,
} from '../render-item.js';
import { type ClearValues, intersectTargetRect, type RenderState, type TargetRect } from '../render-state.js';
import { createGraphicsStats, type GraphicsStats } from '../stats.js';
import { createGraphicsTimer, type TimerSource } from '../timer.js';
import { createWebGPUTimerSource } from './timer.js';
import { BufferUsage, MapMode, TextureUsage } from './constants.js';
import { notSupported } from './errors.js';
import { createWebGPUFramebuffer, isWebGPUFramebuffer, type WebGPUFramebuffer } from './framebuffer.js';
import { asWebGPUGlobalBindings, createWebGPUGlobalBindings, type WebGPUBindingContext } from './global-bindings.js';
import { createWebGPUProgramCache, type WebGPUProgram, type WebGPUProgramCache } from './program.js';
import { createWebGPUProgramItem, createWebGPURenderItem } from './render-item.js';
import {
  createClearDraw,
  getTargetKey,
  type PassTarget,
  type WebGPUTargetDescriptor,
  withDepthReadOnly,
} from './render-pass.js';
import { getColorTargetState, getDepthStencilState, getPrimitiveState, getRenderStateId } from './render-state.js';
import { createMipmapGenerator } from './mipmap.js';
import {
  asWebGPUTexture,
  CubeLayerSides,
  countResource,
  createWebGPUByteCounts,
  createWebGPUCubeTexture,
  createWebGPUTexture,
  type Release,
  type WebGPUResourceEnv,
} from './texture.js';
import { createWebGPURenderbuffer } from './renderbuffer.js';
import { createRenderTargetFromResources, type RenderTarget } from '../render-target.js';
import { type CubeFaces, createNullTexture, decodeImage, type Texture } from '../texture.js';
import { createSharedTextures, type SharedTextures } from '../shared-texture.js';
import {
  BoxItems,
  createQuadIndices,
  LineItems,
  getMaxQuadItemCount,
  type ItemIndexPattern,
  type QuadIndices,
  QuadItems,
} from '../quad-indices.js';
import { createPendingWork } from '../pending-work.js';
import { ValueCell } from '@molstar/core/util';
import { createUniformRing } from './uniform-ring.js';
import { type BufferPages, createBufferPages } from './buffer-pages.js';
import { createBundlePassState, createCaptureEncoder, createPassState, type PassState } from './pass-state.js';
import { getWgslProgram } from './shader/programs.js';
import type { WgslSampler } from './wgsl-declarations.js';

const AllProgramVariants: ProgramVariant[] = [...GraphicsRenderVariants, ...ComputeRenderVariants];
const DepthStencilFormat: GPUTextureFormat = 'depth24plus-stencil8';

/**
 * The color format of a canvas: the browser's preferred one (`bgra8unorm` on Windows and macOS, `rgba8unorm` on
 * Android), which the compositor takes without a copy per frame
 */
function getCanvasFormat(): GPUTextureFormat {
  const preferred = typeof navigator !== 'undefined' ? navigator.gpu?.getPreferredCanvasFormat?.() : undefined;
  return preferred === 'bgra8unorm' ? 'bgra8unorm' : 'rgba8unorm';
}

/** Swaps red and blue of rgba8 texels read from a `bgra8unorm` texture */
function swapRedBlue(data: Uint8Array) {
  for (let i = 0, il = data.length; i < il; i += 4) {
    const r = data[i];
    data[i] = data[i + 2];
    data[i + 2] = r;
  }
}

/** Internals the WebGPU render items draw with */
export interface WebGPUDrawContext extends WebGPUBindingContext, WebGPUResourceEnv {
  readonly stats: GraphicsStats;
  /** The state of the open render pass (its encoder begun on first use), with the dynamic state applied */
  getPassState: () => PassState;
  /** Pipeline of `program` for the current render state and pass target, `undefined` while it compiles */
  getPipeline: (program: WebGPUProgram, drawMode: DrawMode) => GPURenderPipeline | undefined;
  /** Starts compiling the pipeline of `program` for `state` and `target` unless it exists (pending work) */
  preparePipeline: (
    program: WebGPUProgram,
    drawMode: DrawMode,
    state: RenderState,
    target: WebGPUTargetDescriptor,
  ) => void;
  /** How render items request their pipelines for `resources.preparePipelines` */
  readonly itemPreparers: WeakMap<GraphicsRenderItem | ComputeRenderItem | ProgramItem, WebGPUItemPreparer>;
  /** Programs shared by the render items of the context */
  readonly programs: WebGPUProgramCache;
  /** Current render state, for dynamic state (stencil reference) */
  readonly renderState: RenderState;
  onError: (message: string) => void;
  /** Whether render items of programs that are not ported throw instead of drawing nothing */
  readonly strict: boolean;
  /** Warns once per message about a render item that draws nothing */
  warnUnsupported: (message: string) => void;
  /** Render items whose programs `linkPrograms` creates ahead of their first draw */
  readonly linkables: WebGPULinkables;
  /**
   * Shared buffers the vertex and index data of render items is sub-allocated in: draws binding ranges of a few buffers
   * replay far faster than draws binding buffers of their own (`shared-draw-resources.md`)
   */
  readonly geometryPages: BufferPages;
  /** Shared storage buffers the theme layers of render items are sub-allocated in (`data-layers.ts`) */
  readonly dataPages: BufferPages;
  /** The bytes a data page binding may span: an item's layers beyond it go into more pages (`DataPageNames`) */
  readonly dataBindingSize: number;
  /**
   * Shared uniform buffers the object and draw blocks of render items are sub-allocated in: items binding the same
   * share their object bind group (`render-item.ts`)
   */
  readonly uniformPages: BufferPages;
  /** Textures shared by render items (`.v6/plans/gl/shared-textures.md` §2), gone with the device */
  readonly sharedTextures: SharedTextures;
  /** The index buffer of quad items (spheres, points), shared by all (`quad-indices.ts`), a buffer of its own */
  readonly quadIndices: QuadIndices<GPUBuffer>;
  /** The index buffer of box items (cylinders), shared by all */
  readonly boxIndices: QuadIndices<GPUBuffer>;
  /** The index buffer of line items, shared by all */
  readonly lineIndices: QuadIndices<GPUBuffer>;
}

/**
 * Render items with programs to link per variant: added when created or when their programs changed, `linkPrograms`
 * links the pending ones of the requested variants
 */
export interface WebGPULinkables {
  /** marks `linkable` pending for all variants */
  add: (linkable: WebGPULinkable) => void;
  delete: (linkable: WebGPULinkable) => void;
}

/** Requests the pipeline of a render item's program for `variant`, `state` and `target` */
export type WebGPUItemPreparer = (
  variant: ProgramVariant,
  state: RenderState,
  target: WebGPUTargetDescriptor,
  define?: string,
) => void;

/** Creates the programs of a render item for a render variant (shader module, pipelines of known keys) */
export interface WebGPULinkable {
  link: (variant: ProgramVariant) => void;
}

export interface WebGPUContext extends GraphicsContext {
  readonly device: GPUDevice;
  /**
   * Whether a viewport partly outside the target is set as given (cropped); false where the device rejects one, which
   * is then clipped to the target; undefined until the probe answers (`testViewportOutsideTarget`)
   */
  readonly viewportOutsideTarget: boolean | undefined;
  /** The index buffer of quad items, shared by all (`quad-indices.ts`) */
  readonly quadIndices: QuadIndices<GPUBuffer>;
  /** The index buffer of box items (cylinders), shared by all */
  readonly boxIndices: QuadIndices<GPUBuffer>;
  /** The index buffer of line items, shared by all */
  readonly lineIndices: QuadIndices<GPUBuffer>;
  /** Textures shared by render items (`.v6/plans/gl/shared-textures.md`), as WebGL's `resources.sharedTextures` */
  readonly sharedTextures: SharedTextures;
  /**
   * Records the draws of `draw` into a render bundle for the open pass's target, instead of the pass; undefined when
   * the draws can't be replayed from one (a clear between them, a pipeline still compiling, the viewport, scissor or
   * stencil reference changing between them): the caller draws them directly then
   */
  recordRenderBundle: (draw: () => void) => GPURenderBundle | undefined;
  /** Replays `bundles` in the open pass, with the pass's viewport, scissor and stencil reference */
  executeRenderBundles: (bundles: readonly GPURenderBundle[]) => void;
  /**
   * Reads the drawing buffer, rows top to bottom; ends the open pass and submits pending commands. With a canvas, call it
   * in the task that rendered the frame (the canvas texture changes per frame).
   */
  readDrawingBuffer: () => Promise<ImageData>;
  /** Submits the recorded commands now (otherwise at the end of the task); the browser presents the canvas texture */
  present: () => void;
  /** Returns and clears device and shader compilation errors */
  takeErrors: () => string[];
  setDrawingBufferSize: (width: number, height: number) => void;
}

/**
 * Bytes of a page of `WebGPUDrawContext.geometryPages`: many objects' geometry in one buffer, few pages for typical
 * scenes; ranges beyond a quarter get a buffer of their own
 */
const GeometryPageSize = 32 * 1024 * 1024;

/** Bytes of a page of `WebGPUDrawContext.dataPages`, bound whole as storage (at most the binding size limit) */
const DataPageSize = 32 * 1024 * 1024;

/** Bytes of a page of `WebGPUDrawContext.uniformPages`: the blocks of thousands of objects */
const UniformPageSize = 8 * 1024 * 1024;

/** Features used where available, see `capabilities` */
const OptionalFeatures: GPUFeatureName[] = [
  'float32-filterable',
  'float32-blendable',
  'depth-clip-control',
  'depth32float-stencil8',
  // GPU timer
  'timestamp-query',
  // indirect draws of instance ranges (GPU culling of merged render objects)
  'indirect-first-instance',
  // 64-bit atomic min and max (`atomicStoreMax` on `atomic<vec2u>`), not in the WebGPU types yet
  'atomic-vec2u-min-max' as GPUFeatureName,
  // draws of indirect records by one call, on render passes (`--enable-unsafe-webgpu`), not in the WebGPU types
  'chromium-experimental-multi-draw-indirect' as GPUFeatureName,
];

/** Image bitmap of `source` with the raw pixel values, as WebGL uploads with `UNPACK_COLORSPACE_CONVERSION_WEBGL` `NONE` */
async function decodeImageBitmap(source: string | File | Promise<Blob>): Promise<ImageBitmap> {
  const options: ImageBitmapOptions = { colorSpaceConversion: 'none', premultiplyAlpha: 'none' };
  if (typeof source === 'string') return createImageBitmap(await decodeImage(source), options);
  return createImageBitmap(await source, options);
}

export function isWebGPUContext(ctx: GraphicsContext): ctx is WebGPUContext {
  return ctx.capabilities.backend === 'webgpu';
}

/**
 * Limits of a device: the maxima of its adapter (`adapter`), or the core defaults of the spec (`core`, the minimum
 * spec)
 */
export type WebGPULimits = 'adapter' | 'core';

/** Optional features of a device: those of `OptionalFeatures` its adapter has (`optional`), or none (the minimum spec) */
export type WebGPUFeatures = 'optional' | 'none';

/**
 * Every `max*` limit of `adapter` at its value; alignment limits (`min*`) keep their defaults (smaller is "better"
 * there, and the uniform ring assumes 256 bytes)
 */
export function getAdapterMaxLimits(adapter: Pick<GPUAdapter, 'limits'>): Record<string, number> {
  const limits: Record<string, number> = {};
  // limits are getters of the prototype in browsers
  for (const key in adapter.limits) {
    const value = (adapter.limits as unknown as Record<string, unknown>)[key];
    if (key.startsWith('max') && typeof value === 'number') limits[key] = value;
  }
  return limits;
}

/**
 * The device of an adapter of `gpu` (by default the browser's `navigator.gpu`, e.g. Dawn's in Node.js); `undefined`
 * where the adapter lacks one of `requiredFeatures`, or is a fallback (software) adapter and `rejectFallbackAdapter` is
 * set (callers fall back, e.g. to WebGL)
 */
export async function requestWebGPUDevice(
  powerPreference: 'default' | 'high-performance' | 'low-power' = 'high-performance',
  gpu: GPU | undefined = typeof navigator === 'undefined' ? undefined : navigator.gpu,
  limits: WebGPULimits = 'adapter',
  features: WebGPUFeatures = 'optional',
  requiredFeatures: readonly GPUFeatureName[] = [],
  rejectFallbackAdapter = false,
): Promise<GPUDevice | undefined> {
  if (!gpu) return;
  const adapter = await gpu.requestAdapter(powerPreference === 'default' ? undefined : { powerPreference });
  if (!adapter) return;
  // `isFallbackAdapter` on the adapter itself in older browsers
  if (
    rejectFallbackAdapter &&
    (adapter.info?.isFallbackAdapter || (adapter as { isFallbackAdapter?: boolean }).isFallbackAdapter)
  ) {
    if (isDebugMode) console.log('WebGPU adapter is a fallback adapter');
    return;
  }
  const lacking = requiredFeatures.filter((f) => !adapter.features.has(f));
  if (lacking.length) {
    if (isDebugMode) console.log(`WebGPU adapter without required features: ${lacking.join(', ')}`);
    return;
  }
  const optional = features === 'optional' ? OptionalFeatures.filter((f) => adapter.features.has(f)) : [];
  const requiredLimits = limits === 'adapter' ? getAdapterMaxLimits(adapter) : undefined;
  try {
    return await adapter.requestDevice({
      requiredFeatures: [...new Set([...requiredFeatures, ...optional])],
      requiredLimits,
    });
  } catch (e) {
    // e.g. the D3D12 shader compiler failing to load; callers fall back (other backends)
    console.warn('WebGPU device request failed:', e);
    return undefined;
  }
}

/** The properties of a WebGPU info object (getters of its prototype) as a plain object, for logging */
function toPlainObject(o: object | undefined) {
  const plain: Record<string, unknown> = {};
  for (const key in o) plain[key] = (o as Record<string, unknown>)[key];
  return plain;
}

type Rect = { x: number; y: number; width: number; height: number };

/**
 * WebGPU graphics context. The drawing buffer is the canvas texture of the current frame (the preferred canvas format,
 * see `getCanvasFormat`, with a `depth24plus-stencil8` attachment), or an offscreen `rgba8unorm` texture without canvas. Recorded commands are submitted at the
 * end of the task that recorded them, as WebGL flushes implicitly.
 */
export interface WebGPUContextOptions {
  width: number;
  height: number;
  canvas?: HTMLCanvasElement;
  pixelScale?: number;
  /**
   * Unported programs are reported as supported and using them throws a `WebGPUUnsupportedError`, instead of passes
   * falling back and render items drawing nothing (with a warning)
   */
  strict?: boolean;
  /** Requests a device for a new context replacing a lost one (`recreate`); without it a lost context is not replaced */
  requestDevice?: () => Promise<GPUDevice | undefined>;
  /**
   * The context owns its device (requested for it): `destroy` destroys it, and contexts replacing it (`recreate`) own
   * theirs; not for a device shared with others (e.g. tests creating contexts on one device)
   */
  ownsDevice?: boolean;
  /** The bytes a data page binding may span, at most the device's `maxStorageBufferBindingSize` (lower in tests) */
  maxDataBindingSize?: number;
  /**
   * Whether draws writing `frag_depth` may be recorded into render bundles, instead of probing the device
   * (`testBundleFragDepth`); for tests
   */
  fragDepthInBundles?: boolean;
  /**
   * Whether viewports partly outside the target may be set as given, instead of probing the device
   * (`testViewportOutsideTarget`); for tests
   */
  viewportOutsideTarget?: boolean;
}

/** The least time between reports of out-of-memory errors (`GraphicsResources.outOfMemory`) */
const OutOfMemoryReportMs = 5000;

export function createWebGPUContext(device: GPUDevice, options: WebGPUContextOptions): WebGPUContext {
  /** the options of the context, where `destroy`'s parameter shadows `options` */
  const contextOptions = options;
  let destroyed = false;
  /** pending errors until `takeErrors`, the oldest are dropped when nobody takes them */
  const errors: string[] = [];
  const onError = (message: string) => {
    if (isDebugMode) console.error(message);
    if (errors.length >= 100) errors.shift();
    errors.push(message);
  };
  const contextLost = new Subject<now.Timestamp>();
  const outOfMemory = new Subject<string>();
  /** when out-of-memory errors were last reported, see `OutOfMemoryReportMs` */
  let outOfMemoryReported = -Infinity;
  let isContextLost = false;
  /** the device destroyed on purpose (`release`): its loss is expected */
  let released = false;
  const setContextLost = () => {
    if (isContextLost) return;
    isContextLost = true;
    contextLost.next(now());
  };

  /**
   * Removes the context's error listener from its device: the device can outlive the context (Dawn's Node.js bindings
   * keep it by a handle of their own, also destroyed), and the listener's scope is the whole context, its resources
   * and their memory (in RAM on a software adapter: a test run's contexts held to the end)
   */
  let detachDevice: (() => void) | undefined;
  function attachDevice(d: GPUDevice) {
    const onUncapturedError = (e: Event) => {
      const { error } = e as GPUUncapturedErrorEvent;
      onError(`WebGPU error: ${error.message}`);
      // one allocation failing is usually followed by many: reported once in a while
      if (typeof GPUOutOfMemoryError !== 'undefined' && error instanceof GPUOutOfMemoryError) {
        const time = now();
        if (time - outOfMemoryReported < OutOfMemoryReportMs) return;
        outOfMemoryReported = time;
        outOfMemory.next(`Not enough GPU memory: ${error.message}`);
      }
    };
    d.addEventListener('uncapturederror', onUncapturedError);
    detachDevice = () => d.removeEventListener('uncapturederror', onUncapturedError);
    d.lost.then((info) => {
      // losses of replaced devices and of destroyed contexts are expected
      if (d !== device || destroyed) return;
      if (!released) onError(`WebGPU device lost (${info.reason}): ${info.message}`);
      setContextLost();
    });
  }
  attachDevice(device);

  const canvasContext = options.canvas?.getContext('webgpu') as GPUCanvasContext | null | undefined;
  const drawingBufferFormat: GPUTextureFormat = canvasContext ? getCanvasFormat() : 'rgba8unorm';
  /** The drawing buffer's formats, the canvas' (and the zero-size buffer's) or the offscreen buffer's */
  const drawingBufferTarget: WebGPUTargetDescriptor = {
    key: getTargetKey([drawingBufferFormat], DepthStencilFormat),
    colorFormats: [drawingBufferFormat],
    depthFormat: DepthStencilFormat,
    hasStencil: true,
  };
  let canvasConfigured = false;
  /**
   * Configures the canvas context for the current device once the canvas has a size: Firefox validates the size in
   * `configure` (zero-sized canvases fail, e.g. before the page is laid out), Chrome when the frame's texture is taken
   */
  function configureCanvas() {
    canvasConfigured = false;
    if (!canvasContext || !options.canvas || options.canvas.width === 0 || options.canvas.height === 0) return;
    canvasContext.configure({
      device,
      format: drawingBufferFormat,
      usage: TextureUsage.RENDER_ATTACHMENT | TextureUsage.COPY_SRC | TextureUsage.COPY_DST,
      alphaMode: 'premultiplied',
    });
    canvasConfigured = true;
  }
  configureCanvas();
  /** drawn into while the canvas has no size (nothing to show) */
  let zeroSizeBuffer: { color: GPUTexture; depth: GPUTexture; target: PassTarget } | undefined;

  const stats = createGraphicsStats();
  const pending = createPendingWork();

  /**
   * Whether draws writing `frag_depth` may be recorded into render bundles: not where they lose their depth there
   * (`testBundleFragDepth`, Safari); unknown, and so not recorded, until the probe answers
   */
  let fragDepthInBundles = options.fragDepthInBundles;
  if (fragDepthInBundles === undefined) {
    pending.track(
      testBundleFragDepth(device).then((kept) => {
        fragDepthInBundles = kept;
        if (!kept && isDebugMode) console.log('WebGPU: draws writing frag_depth are not recorded into render bundles');
      }),
      'render bundle frag_depth probe',
    );
  }
  /**
   * Whether a viewport partly outside the target is set as given (cropped, as in WebGL): not where the device rejects
   * it (`testViewportOutsideTarget`, Safari), there it is clipped to the target (squashed into it, but its pass drawn);
   * clipped until the probe answers, and beyond the spec's bounds (`toTargetViewport`). Canvas3D crops its frame itself
   * (`Camera.Crop`), so this is a guard
   */
  let viewportOutsideTarget = options.viewportOutsideTarget;
  if (viewportOutsideTarget === undefined) {
    pending.track(
      testViewportOutsideTarget(device).then((taken) => {
        viewportOutsideTarget = taken;
        if (!taken && isDebugMode) console.log('WebGPU: viewports outside the target are clipped to it');
      }),
      'viewport outside the target probe',
    );
  }
  let loggedClippedViewport = false;
  /** the spec's bounds of `setViewport` (2 × `maxTextureDimension2D` either way, its size at most one) */
  const maxViewportSize = device.limits.maxTextureDimension2D;
  const maxViewportRange = 2 * maxViewportSize;
  /**
   * the viewport as set: as given where the device takes one outside the target and the spec's bounds hold (an
   * oversized one, e.g. a frame in stereo, would invalidate the pass everywhere), else clipped to the target
   */
  function toTargetViewport(rect: Rect, target: PassTarget): Rect {
    if (
      viewportOutsideTarget === true &&
      rect.x >= -maxViewportRange &&
      rect.y >= -maxViewportRange &&
      rect.width <= maxViewportSize &&
      rect.height <= maxViewportSize &&
      rect.x + rect.width <= maxViewportRange - 1 &&
      rect.y + rect.height <= maxViewportRange - 1
    ) {
      return rect;
    }
    const x0 = Math.min(Math.max(rect.x, 0), target.width);
    const y0 = Math.min(Math.max(rect.y, 0), target.height);
    const x1 = Math.min(Math.max(rect.x + rect.width, 0), target.width);
    const y1 = Math.min(Math.max(rect.y + rect.height, 0), target.height);
    if (x0 === rect.x && y0 === rect.y && x1 === rect.x + rect.width && y1 === rect.y + rect.height) return rect;
    if (isDebugMode && !loggedClippedViewport) {
      loggedClippedViewport = true;
      console.log(
        `WebGPU: viewport ${rect.x},${rect.y},${rect.width},${rect.height} outside the ${target.width}x${target.height} target clipped to it`,
      );
    }
    return { x: x0, y: y0, width: x1 - x0, height: y1 - y0 };
  }

  /** per program, whether its fragment stage writes `frag_depth` */
  const writesFragDepth = new WeakMap<WebGPUProgram, boolean>();
  function isWritingFragDepth(program: WebGPUProgram) {
    let writes = writesFragDepth.get(program);
    if (writes === undefined) writesFragDepth.set(program, (writes = program.code.includes('@builtin(frag_depth)')));
    return writes;
  }
  const timerSource = createWebGPUTimerSource(device);
  /** the timer keeps its source, which changes with the device */
  const timerSourceRef: TimerSource | undefined = timerSource && {
    begin: () => timerSource?.begin() ?? -1,
    end: () => timerSource?.end(),
    poll: () => timerSource?.poll() ?? [],
    clear: () => timerSource?.clear(),
  };
  const timer = createGraphicsTimer(timerSourceRef, stats);
  /** bytes of the live resources, see `resources.getByteCounts` */
  const byteCounts = createWebGPUByteCounts();
  const ring = createUniformRing(device, byteCounts, stats);
  const globalRing = createUniformRing(device, byteCounts, stats);
  const samplers = {} as { [k in WgslSampler]: GPUSampler };
  function createSamplers() {
    const sampler = (filter: GPUFilterMode, mipmapFilter: GPUMipmapFilterMode = 'nearest') =>
      device.createSampler({ magFilter: filter, minFilter: filter, mipmapFilter });
    samplers.sNearest = sampler('nearest');
    samplers.sLinear = sampler('linear');
    samplers.sLinearMipmap = sampler('linear', 'linear');
  }
  createSamplers();

  // drawing buffer

  let drawingBuffer: { color: GPUTexture; depth: GPUTexture; target: PassTarget } | undefined;
  /** `setPreserveDrawingBuffer`: the canvas's image, copied into each frame's new texture */
  let preserveDrawingBuffer = false;
  let preservedColor: GPUTexture | undefined;
  /** the canvas texture to copy from at a submit: valid in the task that got it, presented at its end */
  let colorToPreserve: GPUTexture | undefined;
  let drawingBufferWidth = options.width;
  let drawingBufferHeight = options.height;
  let pixelScale = options.pixelScale ?? 1;

  function createDrawingBufferTarget(color: GPUTexture, depth: GPUTexture): PassTarget {
    return {
      key: getTargetKey([color.format], DepthStencilFormat),
      colorViews: [color.createView()],
      colorFormats: [color.format],
      depthView: depth.createView(),
      depthFormat: DepthStencilFormat,
      hasStencil: true,
      width: color.width,
      height: color.height,
    };
  }

  function createDepthTexture(width: number, height: number) {
    return device.createTexture({
      size: { width, height },
      format: DepthStencilFormat,
      usage: TextureUsage.RENDER_ATTACHMENT,
      label: 'drawing buffer depth',
    });
  }

  /** A 1×1 target for frames while the canvas has no size */
  function getZeroSizeBuffer() {
    if (!zeroSizeBuffer) {
      const size = { width: 1, height: 1 };
      const color = device.createTexture({
        size,
        format: drawingBufferFormat,
        usage: TextureUsage.RENDER_ATTACHMENT | TextureUsage.COPY_SRC | TextureUsage.TEXTURE_BINDING,
        label: 'zero size drawing buffer',
      });
      const depth = createDepthTexture(1, 1);
      zeroSizeBuffer = { color, depth, target: createDrawingBufferTarget(color, depth) };
    }
    return zeroSizeBuffer;
  }

  /** Into a new canvas texture, in a submit of its own: the current encoder may hold an open pass */
  function restorePreservedColor(color: GPUTexture) {
    if (!preservedColor || preservedColor.width !== color.width || preservedColor.height !== color.height) return;
    const encoder = device.createCommandEncoder();
    encoder.copyTextureToTexture({ texture: preservedColor }, { texture: color }, [color.width, color.height]);
    device.queue.submit([encoder.finish()]);
  }

  function savePreservedColor(encoder: GPUCommandEncoder) {
    const color = colorToPreserve;
    if (!color) return;
    if (preservedColor && (preservedColor.width !== color.width || preservedColor.height !== color.height)) {
      preservedColor.destroy();
      preservedColor = undefined;
    }
    preservedColor ??= device.createTexture({
      size: [color.width, color.height],
      format: drawingBufferFormat,
      usage: TextureUsage.COPY_SRC | TextureUsage.COPY_DST,
      label: 'preserved drawing buffer',
    });
    encoder.copyTextureToTexture({ texture: color }, { texture: preservedColor }, [color.width, color.height]);
  }

  /** The canvas texture of the current frame, or an offscreen texture without canvas */
  function getDrawingBuffer() {
    if (canvasContext) {
      // a canvas resized to zero after it was configured has no texture either
      if (options.canvas && (options.canvas.width === 0 || options.canvas.height === 0)) return getZeroSizeBuffer();
      if (!canvasConfigured) configureCanvas();
      if (!canvasConfigured) return getZeroSizeBuffer();
      // a new texture per frame, presented when the task that rendered into it ends
      const color = canvasContext.getCurrentTexture();
      if (preserveDrawingBuffer && colorToPreserve !== color) {
        colorToPreserve = color;
        setTimeout(() => {
          if (colorToPreserve === color) colorToPreserve = undefined;
        });
        restorePreservedColor(color);
      }
      if (drawingBuffer?.color !== color) {
        let depth = drawingBuffer?.depth;
        if (!depth || depth.width !== color.width || depth.height !== color.height) {
          if (depth) release(depth);
          depth = createDepthTexture(color.width, color.height);
        }
        drawingBuffer = { color, depth, target: createDrawingBufferTarget(color, depth) };
      }
      return drawingBuffer;
    }
    if (!drawingBuffer) {
      const size = { width: Math.max(1, drawingBufferWidth), height: Math.max(1, drawingBufferHeight) };
      const color = device.createTexture({
        size,
        format: drawingBufferFormat,
        usage: TextureUsage.RENDER_ATTACHMENT | TextureUsage.COPY_SRC | TextureUsage.TEXTURE_BINDING,
        label: 'drawing buffer',
      });
      const depth = device.createTexture({
        size,
        format: DepthStencilFormat,
        usage: TextureUsage.RENDER_ATTACHMENT,
        label: 'drawing buffer depth',
      });
      drawingBuffer = {
        color,
        depth,
        target: {
          key: drawingBufferTarget.key,
          colorViews: [color.createView()],
          colorFormats: [drawingBufferFormat],
          depthView: depth.createView(),
          depthFormat: DepthStencilFormat,
          hasStencil: true,
          ...size,
        },
      };
    }
    return drawingBuffer;
  }

  function describeFramebuffer(framebuffer: WebGPUFramebuffer): WebGPUTargetDescriptor {
    const colorFormats = framebuffer.getColorAttachments().map((c) => {
      if (!c) throw new Error('color attachment of the target is missing');
      return c.format;
    });
    const depth = framebuffer.getDepthAttachment();
    return {
      key: getTargetKey(colorFormats, depth?.format),
      colorFormats,
      depthFormat: depth?.format,
      hasStencil: !!depth?.hasStencil,
    };
  }

  /** Depth and stencil of a read-only depth attachment are not cleared, see `RenderPassOptions.depthReadOnly` */
  function checkReadOnlyClear(values: ClearValues) {
    if (passTarget?.depthReadOnly && (values.depth !== undefined || values.stencil !== undefined)) {
      throw new Error('clearing the depth or stencil of a render pass with a read-only depth attachment');
    }
  }

  function describeTarget(
    target: RenderPassTarget,
    options?: Pick<RenderPassOptions, 'depthReadOnly'>,
  ): WebGPUTargetDescriptor {
    const descriptor = describeWritableTarget(target);
    return options?.depthReadOnly ? withDepthReadOnly(descriptor) : descriptor;
  }

  function describeWritableTarget(target: RenderPassTarget): WebGPUTargetDescriptor {
    if (target === 'drawing-buffer') return drawingBufferTarget;
    const framebuffer = 'framebuffer' in target ? target.framebuffer : target;
    if ('framebuffer' in target && drawTargets.has(target)) return drawingBufferTarget;
    if (isWebGPUFramebuffer(framebuffer)) return describeFramebuffer(framebuffer);
    throw new Error('expected a WebGPU framebuffer or render target');
  }

  function getFramebufferTarget(framebuffer: WebGPUFramebuffer): PassTarget {
    const colors = framebuffer.getColorAttachments();
    const depth = framebuffer.getDepthAttachment();
    const first = colors[0] ?? depth;
    if (!first) throw new Error('render pass on a framebuffer without attachments');
    const colorViews: GPUTextureView[] = [];
    for (let i = 0; i < colors.length; ++i) {
      const view = colors[i]?.getAttachmentView();
      if (!view) throw new Error(`color attachment ${i} of the render pass is missing or undefined`);
      colorViews.push(view);
    }
    const colorFormats = colors.map((c) => c!.format);
    return {
      key: getTargetKey(colorFormats, depth?.format),
      colorViews,
      colorDepthSlices: colors.map((c) => c!.depthSlice),
      colorFormats,
      depthView: depth?.getAttachmentView(),
      depthFormat: depth?.format,
      hasStencil: !!depth?.hasStencil,
      width: first.getWidth(),
      height: first.getHeight(),
    };
  }

  /** render targets of the drawing buffer, see `createDrawTarget` */
  const drawTargets = new WeakSet<RenderTarget>();

  function getDrawingBufferTarget() {
    return getDrawingBuffer().target;
  }

  function getPassTarget(target: RenderPassTarget): PassTarget {
    if (target === 'drawing-buffer') return getDrawingBufferTarget();
    const framebuffer = 'framebuffer' in target ? target.framebuffer : target;
    if ('framebuffer' in target && drawTargets.has(target)) return getDrawingBufferTarget();
    if (isWebGPUFramebuffer(framebuffer)) return getFramebufferTarget(framebuffer);
    throw new Error('expected a WebGPU framebuffer or render target');
  }

  /** Rows of a 2d texture (or of the rectangle with a top-left origin), tightly packed, top to bottom */
  async function readGPUTexture(
    texture: GPUTexture,
    bytesPerTexel: number,
    rect = { x: 0, y: 0, width: texture.width, height: texture.height },
  ) {
    submit();
    const { width, height } = rect;
    const rowSize = width * bytesPerTexel;
    const bytesPerRow = Math.ceil(rowSize / 256) * 256;
    const size = bytesPerRow * height;
    const buffer = device.createBuffer({ size, usage: BufferUsage.COPY_DST | BufferUsage.MAP_READ });
    // the staging buffer of the read, as WebGL's pixel pack buffers
    byteCounts.pixelPack += size;
    stats.resourceCounts.pixelPack += 1;
    try {
      const encoder = device.createCommandEncoder();
      encoder.copyTextureToBuffer(
        { texture, origin: { x: rect.x, y: rect.y } },
        { buffer, bytesPerRow },
        { width, height },
      );
      device.queue.submit([encoder.finish()]);
      await buffer.mapAsync(MapMode.READ);
      const src = new Uint8Array(buffer.getMappedRange());
      const out = new Uint8Array(rowSize * height);
      for (let y = 0; y < height; ++y) out.set(src.subarray(y * bytesPerRow, y * bytesPerRow + rowSize), y * rowSize);
      buffer.unmap();
      return out;
    } finally {
      buffer.destroy();
      byteCounts.pixelPack -= size;
      stats.resourceCounts.pixelPack -= 1;
    }
  }

  /** texture formats of the `TextureReader` formats, four bytes per texel */
  const ReaderFormats: { [F in TextureReadFormat]: GPUTextureFormat } = { rgba8: 'rgba8unorm', r32f: 'r32float' };

  /** bytes per texel and its channels, of formats `readTextureData` reads */
  const ReadableFormats: { [k: string]: { bytes: number; channels: number; type: 'ubyte' | 'float' | 'half' } } = {
    r8unorm: { bytes: 1, channels: 1, type: 'ubyte' },
    rg8unorm: { bytes: 2, channels: 2, type: 'ubyte' },
    rgba8unorm: { bytes: 4, channels: 4, type: 'ubyte' },
    r32float: { bytes: 4, channels: 1, type: 'float' },
    rg32float: { bytes: 8, channels: 2, type: 'float' },
    rgba32float: { bytes: 16, channels: 4, type: 'float' },
    r16float: { bytes: 2, channels: 1, type: 'half' },
    rgba16float: { bytes: 8, channels: 4, type: 'half' },
  };

  // passes

  let commandEncoder: GPUCommandEncoder | undefined;
  const encode = (commands: (encoder: GPUCommandEncoder) => void) => {
    // copies run outside of passes: they end the open render pass, and the encoder of the open compute pass (the
    // next dispatch continues the compute pass in a new encoder)
    endPassEncoder();
    endComputePassEncoder();
    if (passTarget) timer.passEnd();
    passTarget = undefined;
    commands(getCommandEncoder());
  };
  const pendingReleases: { destroy(): void }[] = [];
  const release: Release = (resource) => {
    if (commandEncoder) {
      pendingReleases.push(resource);
    } else {
      resource.destroy();
      scheduleReclaim();
    }
  };
  /**
   * Dawn frees the memory of destroyed resources when work is submitted after it: without frames (an idle canvas, a
   * released pass target) a destroy keeps the memory until the next one, so an empty submission follows (once per task)
   */
  let reclaimScheduled = false;
  function scheduleReclaim() {
    if (reclaimScheduled) return;
    reclaimScheduled = true;
    setTimeout(() => {
      reclaimScheduled = false;
      if (!destroyed && !isContextLost) device.queue.submit([]);
    }, 0);
  }
  let computePassEncoder: GPUComputePassEncoder | undefined;
  /** label of the open compute pass, the encoder is begun on the first dispatch */
  let computePassLabel: string | undefined;

  function endComputePassEncoder() {
    computePassEncoder?.end();
    computePassEncoder = undefined;
  }

  function getComputePassEncoder() {
    if (computePassLabel === undefined) throw new Error('dispatch outside of a compute pass');
    if (!computePassEncoder) {
      const descriptor: GPUComputePassDescriptor = { label: computePassLabel };
      const timestampWrites = timerSource?.getTimestampWrites();
      if (timestampWrites) descriptor.timestampWrites = timestampWrites;
      computePassEncoder = getCommandEncoder().beginComputePass(descriptor);
    }
    return computePassEncoder;
  }

  function endComputePass() {
    endComputePassEncoder();
    if (computePassLabel !== undefined) timer.passEnd();
    computePassLabel = undefined;
  }

  function beginComputePass(options: { label: string }) {
    endPassEncoder();
    if (passTarget) timer.passEnd();
    passTarget = undefined;
    endComputePass();
    timer.passBegin('compute', options.label);
    computePassLabel = options.label;
  }

  let passEncoder: GPURenderPassEncoder | undefined;
  let passState: PassState | undefined;
  /** the passes' state, rebound to each pass (`beginPassEncoder`) */
  let directPassState: PassState | undefined;
  /** target of the open pass, the encoder is begun on first draw */
  let passTarget: PassTarget | undefined;
  type PassClear = { readonly color?: ClearValues['color']; readonly depth?: number; readonly stencil?: number };
  /** a pass clearing nothing: shared, not made a pass */
  const NoClear: PassClear = {};
  let passClear: PassClear = NoClear;
  /**
   * The descriptor of each target's passes, its attachments' views, load operations and clear values set anew each
   * pass: a stable target (the drawing buffer's) begins its passes without allocating
   */
  const passDescriptors = new WeakMap<
    PassTarget,
    { descriptor: GPURenderPassDescriptor; colors: GPURenderPassColorAttachment[]; clears: GPUColorDict[] }
  >();
  function getPassDescriptor(target: PassTarget) {
    let d = passDescriptors.get(target);
    if (!d) {
      const colors: GPURenderPassColorAttachment[] = [];
      const clears: GPUColorDict[] = [];
      for (let i = 0, il = target.colorViews.length; i < il; ++i) {
        colors.push({ view: target.colorViews[i], loadOp: 'load', storeOp: 'store' });
        clears.push({ r: 0, g: 0, b: 0, a: 0 });
      }
      const descriptor: GPURenderPassDescriptor = { colorAttachments: colors };
      if (target.depthView && target.depthReadOnly) {
        // read in the pass (and bindable as a texture), so neither loaded nor stored
        descriptor.depthStencilAttachment = {
          view: target.depthView,
          depthReadOnly: true,
          stencilReadOnly: target.hasStencil ? true : undefined,
        };
      } else if (target.depthView) {
        descriptor.depthStencilAttachment = {
          view: target.depthView,
          depthLoadOp: 'load',
          depthStoreOp: 'store',
          depthClearValue: 0,
          stencilLoadOp: target.hasStencil ? 'load' : undefined,
          stencilStoreOp: target.hasStencil ? 'store' : undefined,
          stencilClearValue: 0,
        };
      }
      d = { descriptor, colors, clears };
      passDescriptors.set(target, d);
    }
    return d;
  }

  let renderState: RenderState = {};
  let blendConstant: GPUColorDict = { r: 0, g: 0, b: 0, a: 0 };
  let viewport: Rect | undefined;
  let scissor: Rect | undefined;

  function getCommandEncoder() {
    if (!commandEncoder) {
      commandEncoder = device.createCommandEncoder();
      // WebGL flushes implicitly: submit what the current task recorded (e.g. a frame), before the canvas is presented
      queueMicrotask(submit);
    }
    return commandEncoder;
  }

  function beginPassEncoder(target: PassTarget) {
    const { color, depth, stencil } = passClear;
    const { descriptor, colors, clears } = getPassDescriptor(target);
    for (let i = 0, il = colors.length; i < il; ++i) {
      const a = colors[i];
      a.view = target.colorViews[i];
      a.depthSlice = target.colorDepthSlices?.[i];
      a.loadOp = color ? 'clear' : 'load';
      if (color) {
        const c = clears[i];
        c.r = color[0];
        c.g = color[1];
        c.b = color[2];
        c.a = color[3];
        a.clearValue = c;
      } else {
        a.clearValue = undefined;
      }
    }
    const ds = descriptor.depthStencilAttachment;
    if (ds && target.depthView) {
      ds.view = target.depthView;
      if (!target.depthReadOnly) {
        ds.depthLoadOp = depth !== undefined ? 'clear' : 'load';
        ds.depthClearValue = depth ?? 0;
        if (target.hasStencil) ds.stencilLoadOp = stencil !== undefined ? 'clear' : 'load';
        ds.stencilClearValue = stencil ?? 0;
      }
    }
    passClear = NoClear;
    descriptor.timestampWrites = timerSource?.getTimestampWrites();
    const pass = getCommandEncoder().beginRenderPass(descriptor);
    pass.setBlendConstant(blendConstant);
    if (directPassState) directPassState.rebind(pass);
    else directPassState = createPassState(pass);
    passState = directPassState;
    return pass;
  }

  function endPassEncoder() {
    if (passEncoder) flushBundles();
    if (passTarget && !passEncoder && hasClear(passClear)) {
      // a pass with clears only
      passEncoder = beginPassEncoder(passTarget);
    }
    passEncoder?.end();
    passEncoder = undefined;
    passState = undefined;
  }

  function hasClear(c: typeof passClear) {
    return c.color !== undefined || c.depth !== undefined || c.stencil !== undefined;
  }

  function coversTarget(rect: Rect | undefined, target: PassTarget) {
    return (
      !rect ||
      (rect.x <= 0 && rect.y <= 0 && rect.x + rect.width >= target.width && rect.y + rect.height >= target.height)
    );
  }

  const scissorRect: TargetRect = { x: 0, y: 0, width: 0, height: 0 };

  /** Whether the scissor (the target without one) lies inside the viewport (the target without one) */
  function isScissorInViewport(target: PassTarget) {
    const s = scissor ? toTargetScissor(scissor, target) : { x: 0, y: 0, width: target.width, height: target.height };
    const v = viewport
      ? toTargetViewport(viewport, target)
      : { x: 0, y: 0, width: target.width, height: target.height };
    return s.x >= v.x && s.y >= v.y && s.x + s.width <= v.x + v.width && s.y + s.height <= v.y + v.height;
  }

  /** the part of the scissor inside the target (`intersectTargetRect`) */
  function toTargetScissor(rect: Rect, target: PassTarget): TargetRect {
    return intersectTargetRect(scissorRect, rect.x, rect.y, rect.width, rect.height, target.width, target.height);
  }

  /** the render bundle being recorded (`recordRenderBundle`), the draws' state meanwhile */
  let recording: { state: PassState; failed?: string } | undefined;
  /** the draws being captured (`captureDraws`): their state, encoding nothing */
  let capturing: PassState | undefined;

  /**
   * Bundles replayed one after another, executed together before the next command of the pass (`getPassState`, the
   * pass's end, a blend constant): one `executeBundles` for a pass replayed in chunks; `pendingDynamic` the dynamic
   * state they run with (`readDynamicState`). Both reused: a replay allocates nothing
   */
  const pendingBundles: GPURenderBundle[] = [];
  const pendingDynamic = new Float64Array(9);
  const currentDynamic = new Float64Array(9);
  function flushBundles() {
    if (pendingBundles.length === 0) return;
    // the encoder copies the list
    passEncoder!.executeBundles(pendingBundles);
    pendingBundles.length = 0;
    // a bundle leaves the pipeline, bind groups and buffers unset
    passState!.invalidate();
  }
  /** The viewport, scissor (`-Infinity` where unset) and stencil reference into `out` */
  function readDynamicState(out: Float64Array) {
    const v = viewport;
    const s = scissor;
    out[0] = v ? v.x : -Infinity;
    out[1] = v ? v.y : -Infinity;
    out[2] = v ? v.width : -Infinity;
    out[3] = v ? v.height : -Infinity;
    out[4] = s ? s.x : -Infinity;
    out[5] = s ? s.y : -Infinity;
    out[6] = s ? s.width : -Infinity;
    out[7] = s ? s.height : -Infinity;
    out[8] = renderState.stencil?.reference ?? 0;
  }
  function isPendingDynamicState() {
    for (let i = 0; i < 9; ++i) if (currentDynamic[i] !== pendingDynamic[i]) return false;
    return true;
  }

  function getPassState() {
    if (!passTarget) throw new Error('draw outside of a render pass');
    if (!passEncoder) passEncoder = beginPassEncoder(passTarget);
    flushBundles();
    const state = capturing ?? recording?.state ?? passState!;
    const target = passTarget;
    if (viewport) {
      // as given where the device takes it: a viewport partly outside the target is cropped as in WebGL
      const { x, y, width, height } = toTargetViewport(viewport, target);
      state.setViewport(x, y, width, height);
    }
    if (scissor) {
      const r = toTargetScissor(scissor, target);
      state.setScissorRect(r.x, r.y, r.width, r.height);
    } else {
      state.setScissorRect(0, 0, target.width, target.height);
    }
    state.setStencilReference(renderState.stencil?.reference ?? 0);
    return state;
  }

  /** pipeline keys compiling in the background, per program */
  const compilingPipelines = new WeakMap<WebGPUProgram, Set<string>>();
  /** pipeline keys whose background compilation failed: created synchronously, which reports the error once */
  const failedPipelines = new WeakMap<WebGPUProgram, Set<string>>();

  /** The program independent part of a pipeline: render state, pass target formats and draw mode */
  type PipelineConfig = Pick<GPURenderPipelineDescriptor, 'primitive' | 'depthStencil'> & {
    readonly targets: GPUColorTargetState[];
  };
  /** Pipeline configurations by key (render state, target formats, draw mode), shared by the programs */
  const pipelineConfigs = new Map<string, PipelineConfig>();

  function onProgramCreated(_program: WebGPUProgram) {
    stats.resourceCounts.program += 1;
    // its shader module, holding both stages (programs and modules live as long as the context)
    stats.resourceCounts.shader += 1;
  }

  /**
   * The pipeline of `program` for the current render state and pass target, `undefined` while it compiles. Missing
   * pipelines compile in the background (`createRenderPipelineAsync`, off the GPU process' main thread) and
   * their draws are skipped until ready, so first uses (e.g. the first hover: picking and marking) don't stall frames;
   * the compilation is pending work, `Canvas3D` redraws when it settles and `ImagePass` renders again
   */
  function getPipeline(program: WebGPUProgram, drawMode: DrawMode): GPURenderPipeline | undefined {
    const target = passTarget;
    if (!target) throw new Error('draw outside of a render pass');
    // drawn directly where a render bundle loses the depth it writes
    if (recording && fragDepthInBundles !== true && isWritingFragDepth(program)) {
      recording.failed ??= 'a draw writing frag_depth';
    }
    return getPipelineFor(program, drawMode, renderState, target);
  }

  function preparePipeline(
    program: WebGPUProgram,
    drawMode: DrawMode,
    state: RenderState,
    target: WebGPUTargetDescriptor,
  ) {
    getPipelineFor(program, drawMode, state, target);
  }

  const itemPreparers = new WeakMap<GraphicsRenderItem | ComputeRenderItem | ProgramItem, WebGPUItemPreparer>();

  function getPipelineFor(
    program: WebGPUProgram,
    drawMode: DrawMode,
    state: RenderState,
    target: WebGPUTargetDescriptor,
  ): GPURenderPipeline | undefined {
    const last = program.lastPipeline;
    if (last && last.state === state && last.target === target && last.drawMode === drawMode) return last.pipeline;
    const key = `${getRenderStateId(state)}|${target.key}|${drawMode}`;
    const existing = program.pipelines.get(key);
    if (existing) {
      program.lastPipeline = { state, target, drawMode, pipeline: existing };
      return existing;
    }
    let config = pipelineConfigs.get(key);
    if (!config) {
      config = {
        targets: target.colorFormats.map((f) => getColorTargetState(state, f)),
        primitive: getPrimitiveState(state, drawMode),
        depthStencil: target.depthFormat
          ? getDepthStencilState(state, target.depthFormat, target.hasStencil, target.depthReadOnly)
          : undefined,
      };
      pipelineConfigs.set(key, config);
    }
    return requestPipeline(program, key, config);
  }

  /**
   * Pipelines compiling in the background together (timing mode): from the first request while none compiled to the
   * last of them ready, wall time (Dawn compiles on worker threads, in parallel)
   */
  let compileBatch: { start: number; count: number; longest: number; longestLabel: string } | undefined;
  let compilingCount = 0;

  function onPipelineCompiled(label: string, start: number) {
    compilingCount -= 1;
    if (!compileBatch) return;
    const time = now() - start;
    if (time > compileBatch.longest) {
      compileBatch.longest = time;
      compileBatch.longestLabel = label;
    }
    if (compilingCount > 0) return;
    const { count, longest, longestLabel } = compileBatch;
    console.log(
      `Compiled ${count} pipelines in ${(now() - compileBatch.start).toFixed(2)} ms (longest ${longest.toFixed(2)} ms: ${longestLabel})`,
    );
    compileBatch = undefined;
  }

  /** The pipeline of `program` for `config`, `undefined` while it compiles in the background */
  function requestPipeline(program: WebGPUProgram, key: string, config: PipelineConfig) {
    const existing = program.pipelines.get(key);
    if (existing) return existing;
    const descriptor: GPURenderPipelineDescriptor = {
      label: `${program.name} ${program.variant} ${key}`,
      layout: program.pipelineLayout,
      vertex: {
        module: program.module,
        entryPoint: 'vertexMain',
        buffers: [...program.vertexBuffers],
        constants: program.vertexConstants,
      },
      fragment: {
        module: program.module,
        entryPoint: 'fragmentMain',
        targets: config.targets,
        constants: program.fragmentConstants,
      },
      primitive: config.primitive,
      depthStencil: config.depthStencil,
    };
    if (failedPipelines.get(program)?.has(key)) {
      const start = isTimingMode ? now() : 0;
      const pipeline = device.createRenderPipeline(descriptor);
      if (isTimingMode)
        console.log(`Created pipeline ${descriptor.label} synchronously in ${(now() - start).toFixed(2)} ms`);
      program.pipelines.set(key, pipeline);
      return pipeline;
    }
    let compiling = compilingPipelines.get(program);
    if (!compiling) compilingPipelines.set(program, (compiling = new Set()));
    if (!compiling.has(key)) {
      compiling.add(key);
      const pipelineDevice = device;
      const start = now();
      compilingCount += 1;
      if (isTimingMode) {
        if (!compileBatch) compileBatch = { start, count: 0, longest: 0, longestLabel: '' };
        compileBatch.count += 1;
      }
      const label = descriptor.label!;
      const compiled = pipelineDevice.createRenderPipelineAsync(descriptor).then(
        (p) => {
          // pipelines of a replaced device are gone with it
          if (pipelineDevice === device) program.pipelines.set(key, p);
        },
        (e) => {
          let failed = failedPipelines.get(program);
          if (!failed) failedPipelines.set(program, (failed = new Set()));
          failed.add(key);
          throw e;
        },
      );
      const keys = compiling;
      // settling requests a redraw, which draws with the pipeline
      pending.track(
        compiled.finally(() => {
          keys.delete(key);
          onPipelineCompiled(label, start);
        }),
        `pipeline ${program.name} ${program.variant} ${key}`,
      );
    }
    // the draw is skipped: a recording without it would replay without it
    if (recording) recording.failed ??= 'a pipeline compiling';
    return undefined;
  }

  /** Advances when one of the shared index buffers is replaced (grown): their ids only increase */
  function getSharedIndicesVersion() {
    return drawContext.quadIndices.id + drawContext.boxIndices.id + drawContext.lineIndices.id;
  }

  function recordRenderBundle(draw: () => void): GPURenderBundle | undefined {
    const target = passTarget;
    if (!target) throw new Error('recording a render bundle outside of a render pass');
    if (recording) throw new Error('recording a render bundle while recording one');
    const pass = getPassState();
    const encoder = device.createRenderBundleEncoder({
      colorFormats: [...target.colorFormats],
      depthStencilFormat: target.depthFormat,
      depthReadOnly: !!target.depthReadOnly,
      stencilReadOnly: !!target.depthReadOnly && target.hasStencil,
    });
    const current: { state: PassState; failed?: string } = { state: pass };
    current.state = createBundlePassState(encoder, pass, (what) => {
      current.failed ??= `the ${what} changing`;
    });
    const arenaReplacements = ring.arenaReplacements;
    const sharedIndices = getSharedIndicesVersion();
    recording = current;
    try {
      draw();
    } finally {
      recording = undefined;
    }
    // draws recorded before bind the released buffer
    if (ring.arenaReplacements !== arenaReplacements) current.failed ??= 'a uniform arena growing';
    // likewise the shared quad or box indices grown by a later draw (`quad-indices.ts`): drawn directly, recorded again
    // next frame with the grown buffer
    if (getSharedIndicesVersion() !== sharedIndices) {
      current.failed ??= 'the shared indices growing';
    }
    const bundle = encoder.finish();
    return current.failed ? undefined : bundle;
  }

  const captureIds = new WeakMap<object, number>();
  let nextCaptureId = 0;
  const getCaptureId = (o: object) => {
    let id = captureIds.get(o);
    if (id === undefined) captureIds.set(o, (id = nextCaptureId++));
    return id;
  };

  function captureDraws(draw: () => void): string {
    if (!passTarget) throw new Error('capturing draws outside of a render pass');
    if (recording || capturing) throw new Error('capturing draws while recording or capturing');
    const tokens: string[] = [];
    // a state of its own: every draw sets what it binds, whatever the draws before it bound
    capturing = createPassState(createCaptureEncoder(tokens, getCaptureId));
    try {
      draw();
    } finally {
      capturing = undefined;
    }
    return tokens.join(' ');
  }

  /** Recordings a draw list keeps, one per target layout its pass draws into (e.g. the canvas and image passes) */
  const MaxDrawListTargets = 4;

  /**
   * Render bundles of a pass's draws, one per target layout (`PassTarget.key`: formats, read-only depth), each replayed
   * while the pass binds the globals it was recorded with
   */
  function createDrawList(): DrawList {
    type Recording = {
      bundle: GPURenderBundle;
      state: ReturnType<ReturnType<typeof asWebGPUGlobalBindings>['getState']>;
      /** `UniformRing.arenaReplacements` when recorded: a replaced arena buffer may be bound by the bundle */
      arenaReplacements: number;
    };
    const recordings = new Map<string, Recording>();
    return {
      replay: (globals) => {
        if (!passTarget) return false;
        const recording = recordings.get(passTarget.key);
        if (!recording) return false;
        const state = asWebGPUGlobalBindings(globals).getState();
        const recorded = recording.state;
        if (
          state.buffer !== recorded.buffer ||
          state.offset !== recorded.offset ||
          state.textures !== recorded.textures ||
          recording.arenaReplacements !== ring.arenaReplacements
        ) {
          recordings.delete(passTarget.key);
          return false;
        }
        executeRenderBundle(recording.bundle);
        return true;
      },
      record: (globals, draw) => {
        const target = passTarget;
        if (!target) throw new Error('recording a draw list outside of a render pass');
        // copied: the bindings update their state in place
        const { buffer, offset, textures } = asWebGPUGlobalBindings(globals).getState();
        const state = { buffer, offset, textures };
        const bundle = recordRenderBundle(draw);
        recordings.delete(target.key);
        if (bundle) {
          if (recordings.size >= MaxDrawListTargets) recordings.delete(recordings.keys().next().value!);
          recordings.set(target.key, { bundle, state, arenaReplacements: ring.arenaReplacements });
          executeRenderBundle(bundle);
        } else {
          draw();
        }
      },
      reset: () => {
        recordings.clear();
      },
    };
  }

  function executeRenderBundle(bundle: GPURenderBundle) {
    if (recording) throw new Error('replaying render bundles while recording one');
    readDynamicState(currentDynamic);
    // the dynamic state applied (and earlier bundles executed) unless it is the one of the bundles pending
    if (pendingBundles.length === 0 || !isPendingDynamicState()) {
      getPassState();
      pendingDynamic.set(currentDynamic);
    }
    pendingBundles.push(bundle);
  }
  function executeRenderBundles(bundles: readonly GPURenderBundle[]) {
    for (let i = 0, il = bundles.length; i < il; ++i) executeRenderBundle(bundles[i]);
  }

  let submitCount = 0;
  function submit() {
    endPassEncoder();
    endComputePassEncoder();
    if (!commandEncoder) return;
    savePreservedColor(commandEncoder);
    ring.upload();
    globalRing.upload();
    const afterSubmit = timerSource?.beforeSubmit(commandEncoder);
    device.queue.submit([commandEncoder.finish()]);
    afterSubmit?.();
    commandEncoder = undefined;
    ring.reset();
    globalRing.reset();
    submitCount += 1;
    reclaimPages(false);
    if (pendingReleases.length) {
      for (const r of pendingReleases) r.destroy();
      pendingReleases.length = 0;
      scheduleReclaim();
    }
  }

  const clearDraw = createClearDraw(device, ring);

  const programs = createWebGPUProgramCache(device, onError, onProgramCreated);
  /** per variant, the render items not linked for it */
  const pendingLinks = new Map<ProgramVariant, Set<WebGPULinkable>>(AllProgramVariants.map((v) => [v, new Set()]));
  const linkables: WebGPULinkables = {
    add: (linkable) => {
      for (const pending of pendingLinks.values()) pending.add(linkable);
    },
    delete: (linkable) => {
      for (const pending of pendingLinks.values()) pending.delete(linkable);
    },
  };

  const generateMipmaps = createMipmapGenerator(device);
  const ownTextures = new WeakSet<Texture>();
  const env: WebGPUResourceEnv = {
    device,
    stats,
    byteCounts,
    release,
    encode,
    flush: submit,
    get generateMipmaps() {
      return generateMipmaps;
    },
  };

  const compute = createWebGPUComputeContext({
    device,
    ring,
    stats,
    pending,
    release,
    encode,
    flush: submit,
    getComputePassEncoder,
    beginComputePass,
    endComputePass,
    onError,
    storageBufferBytes: 0,
    byteCounts,
  });

  const strict = options.strict ?? false;
  const warned = new Set<string>();

  const geometryPages = createBufferPages(device, {
    usage: BufferUsage.VERTEX | BufferUsage.INDEX,
    pageSize: Math.min(GeometryPageSize, device.limits.maxBufferSize),
    label: 'geometry page',
    getSubmitCount: () => submitCount,
    onFree: () => schedulePageReclaim(),
  });
  const dataBindingSize = Math.min(
    device.limits.maxStorageBufferBindingSize,
    contextOptions.maxDataBindingSize ?? Number.POSITIVE_INFINITY,
  );
  const dataPages = createBufferPages(device, {
    usage: BufferUsage.STORAGE,
    pageSize: Math.min(DataPageSize, dataBindingSize),
    label: 'data page',
    getSubmitCount: () => submitCount,
    onFree: () => schedulePageReclaim(),
  });
  const uniformPages = createBufferPages(device, {
    usage: BufferUsage.UNIFORM,
    pageSize: Math.min(UniformPageSize, device.limits.maxBufferSize),
    alignment: device.limits.minUniformBufferOffsetAlignment,
    label: 'uniform page',
    getSubmitCount: () => submitCount,
    onFree: () => schedulePageReclaim(),
  });
  /**
   * The pages' freed ranges reclaimed (and empty pages dropped) after each submission, and once the canvas is idle
   * after a free (no frame follows a removal): else only by a later allocation, its pages kept meanwhile
   */
  const reclaimPages = (idle: boolean) => {
    geometryPages.reclaim(idle);
    dataPages.reclaim(idle);
    uniformPages.reclaim(idle);
  };
  let pageReclaimScheduled = false;
  function schedulePageReclaim() {
    if (pageReclaimScheduled) return;
    pageReclaimScheduled = true;
    setTimeout(() => {
      pageReclaimScheduled = false;
      if (destroyed || isContextLost) return;
      // commands recorded but not submitted may still read a freed range: those wait for their submission
      reclaimPages(!commandEncoder);
      // the memory of the pages dropped is freed on a submission
      scheduleReclaim();
    }, 0);
  }

  /** The shared index buffer of the items of `items` (`quad-indices.ts`), a buffer of its own */
  function createItemIndices(items: ItemIndexPattern) {
    return createQuadIndices(
      (indices) => {
        const buffer = device.createBuffer({
          size: indices.byteLength,
          usage: BufferUsage.INDEX | BufferUsage.COPY_DST,
          label: `${items.name} indices`,
        });
        device.queue.writeBuffer(buffer, 0, indices);
        byteCounts.elements += indices.byteLength;
        stats.resourceCounts.elements += 1;
        return buffer;
      },
      (buffer) => {
        byteCounts.elements -= buffer.size;
        stats.resourceCounts.elements -= 1;
        // recorded work may still bind it
        release(buffer);
      },
      device.limits.maxBufferSize,
      items,
    );
  }

  const drawContext: WebGPUDrawContext = {
    byteCounts,
    get generateMipmaps() {
      return generateMipmaps;
    },
    encode,
    flush: submit,
    device,
    programs,
    linkables,
    release,
    ring,
    globalRing,
    geometryPages,
    dataPages,
    dataBindingSize,
    uniformPages,
    samplers,
    stats,
    getPassState,
    getPipeline,
    preparePipeline,
    itemPreparers,
    get renderState() {
      return renderState;
    },
    onError,
    strict,
    warnUnsupported: (message: string) => {
      if (warned.has(message)) return;
      warned.add(message);
      console.warn(`${message}, not drawn`);
    },
    sharedTextures: createSharedTextures(
      (spec) => createWebGPUTexture(drawContext, spec.kind, spec.format, spec.dataType, spec.filter),
      (texture) => asWebGPUTexture(texture).version,
    ),
    quadIndices: createItemIndices(QuadItems),
    boxIndices: createItemIndices(BoxItems),
    lineIndices: createItemIndices(LineItems),
  };

  // stubs

  // XR with `XRGPUBinding` (sessions with the 'webgpu' feature, an `xrCompatible` adapter): a projection layer, the
  // drawing buffer a view's sub-image (`getViewSubImage`, its viewport with a top-left origin as given), projections
  // of the views with a `[0, 1]` clip depth converted to GL's
  const xr: GraphicsXR = {
    session: undefined,
    changed: new Subject<void>(),
    set: async () => notSupported('XR'),
    end: async () => {},
    isPresenting: () => false,
    getView: () => false,
  };

  const capabilities: GraphicsCapabilities = {
    backend: 'webgpu',
    texture3d: true,
    flatInterpolation: true,
    integerTextures: true,
    mixedAttachmentFormats: true,
    floatTextureLinear: device.features.has('float32-filterable'),
    halfFloatTexture: true,
    halfFloatTextureLinear: true,
    floatRenderTarget: true,
    // blending into float32 targets (wboit, dpoit) needs the feature, as WebGL needs EXT_float_blend
    floatBlend: device.features.has('float32-blendable'),
    halfFloatRenderTarget: true,
    depthTexture: true,
    fragDepth: true,
    multipleRenderTargets: true,
    textureLod: true,
    multiDraw: false,
    clipControl: true,
    depthClamp: device.features.has('depth-clip-control'),
    compute: true,
    indirectInstanceRanges: device.features.has('indirect-first-instance'),
    multiDrawIndirect: device.features.has('chromium-experimental-multi-draw-indirect' as GPUFeatureName),
    atomicMinMax64: device.features.has('atomic-vec2u-min-max' as GPUFeatureName),
    parallelCompile: true,
  };

  const resources: GraphicsResources = {
    limits: {
      maxTextureSize: device.limits.maxTextureDimension2D,
      max3dTextureSize: device.limits.maxTextureDimension3D,
      maxRenderbufferSize: device.limits.maxTextureDimension2D,
      maxDrawBuffers: device.limits.maxColorAttachments,
      maxTextureImageUnits: device.limits.maxSampledTexturesPerShaderStage,
      maxQuadItemCount: getMaxQuadItemCount(device.limits.maxBufferSize),
      maxBoxItemCount: getMaxQuadItemCount(device.limits.maxBufferSize, BoxItems),
      maxLineItemCount: getMaxQuadItemCount(device.limits.maxBufferSize, LineItems),
    },
    texture: (kind, format, type, filter) => {
      const texture = createWebGPUTexture(env, kind, format, type, filter);
      ownTextures.add(texture);
      return texture;
    },
    ownsTexture: (texture) => ownTextures.has(texture),
    cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: (errored?: boolean) => void) => {
      const texture = createWebGPUCubeTexture(env, mipmaps);
      let destroyed = false;
      const load = async () => {
        const images = await Promise.all(CubeLayerSides.map((side) => decodeImageBitmap(faces[side])));
        if (destroyed) return;
        texture.loadFaces(images);
      };
      // `onload` as part of the pending work, so that work waiting for the faces sees what `onload` did
      pending.track(
        load().then(
          () => onload?.(destroyed),
          () => onload?.(true),
        ),
        'cube texture',
      );
      const destroy = texture.destroy;
      texture.destroy = () => {
        destroyed = true;
        destroy();
      };
      return texture;
    },
    framebuffer: () => {
      const framebuffer = createWebGPUFramebuffer();
      const uncount = countResource(stats, 'framebuffer');
      const destroy = framebuffer.destroy;
      framebuffer.destroy = () => {
        destroy();
        uncount();
      };
      return framebuffer;
    },
    renderbuffer: (format, attachment, width, height) =>
      createWebGPURenderbuffer(env, format, attachment, width, height),
    uniformBuffer: () => notSupported('uniformBuffer'),
    uniformRingBuffer: () => notSupported('uniformRingBuffer'),
    materialUniformBuffer: () => notSupported('materialUniformBuffer'),
    getByteCounts: () => ({ ...byteCounts, storageBuffer: compute.getStorageBufferByteCount() }),
    pending,
    outOfMemory,
    loadImage: (texture: Texture, source: string, cell?: ValueCell<Texture>, label = 'image') => {
      const load = async () => {
        const element = await decodeImage(source);
        // the raw pixel values, as WebGL uploads with `UNPACK_COLORSPACE_CONVERSION_WEBGL` `NONE`
        const image =
          typeof createImageBitmap === 'undefined'
            ? element
            : await createImageBitmap(element, { colorSpaceConversion: 'none', premultiplyAlpha: 'none' });
        asWebGPUTexture(texture).load(image);
        if (cell) ValueCell.update(cell, texture);
      };
      return pending.track(load(), label);
    },
    linkPrograms: (variants?: ProgramVariant[]) => {
      for (const v of variants ?? AllProgramVariants) {
        const pending = pendingLinks.get(v);
        if (!pending || pending.size === 0) continue;
        for (const l of pending) l.link(v);
        pending.clear();
      }
    },
    // frames don't wait for programs: draws of pipelines compiling in the background are skipped, see `getPipeline`
    finalizePrograms: () => true,
    preparePipelines: (requests: Iterable<PipelineRequest>) => {
      for (const { item, variant, state, target, define } of requests) {
        itemPreparers.get(item)?.(variant, state, target as WebGPUTargetDescriptor, define);
      }
      return pending.whenReady();
    },
  };

  const ctx: WebGPUContext = {
    device,
    recordRenderBundle,
    executeRenderBundles,
    createDrawList,
    captureDraws,
    get submitCount() {
      return submitCount;
    },
    get viewportOutsideTarget() {
      return viewportOutsideTarget;
    },
    get drawResourcesVersion() {
      // both only advance: their sum advances with either
      return getSharedIndicesVersion() + drawContext.sharedTextures.version;
    },
    quadIndices: drawContext.quadIndices,
    boxIndices: drawContext.boxIndices,
    lineIndices: drawContext.lineIndices,
    sharedTextures: drawContext.sharedTextures,
    capabilities,
    resources,
    stats,
    timer,
    compute,
    get pixelRatio() {
      const dpr = typeof window !== 'undefined' ? window.devicePixelRatio || 1 : 1;
      return dpr * pixelScale;
    },
    setPixelScale: (value: number) => {
      pixelScale = value;
    },

    get isContextLost() {
      return isContextLost;
    },
    contextLost,
    setContextLost,
    simulateContextLoss: () => {
      if (isContextLost || destroyed) return false;
      // reported as a loss by `device.lost`
      device.destroy();
      return true;
    },
    release: () => {
      if (isContextLost || destroyed) return false;
      // the device's resources go with it, reported as a loss by `device.lost`; `recreate` requests a new one
      released = true;
      // reads in flight fail with the device: forgotten, not reported
      pending.clear();
      device.destroy();
      return true;
    },

    xr,

    namedFramebuffers: Object.create(null),
    namedTextures: Object.create(null),
    namedComputeRenderables: Object.create(null),

    isProgramSupported: (shaderCode) => strict || getWgslProgram(shaderCode.name) !== undefined,
    createProgramItem: (drawMode, shaderCode, schema, defineValues, linearTextures, transparency) =>
      createWebGPUProgramItem(drawContext, drawMode, shaderCode, schema, defineValues, linearTextures, transparency),
    createGraphicsRenderItem: (drawMode, shaderCode, schema, values, materialId, transparency) =>
      createWebGPURenderItem<GraphicsRenderVariant>(
        drawContext,
        drawMode,
        shaderCode,
        schema,
        values,
        materialId,
        transparency,
      ),
    createComputeRenderItem: (drawMode, shaderCode, schema, values, materialId = -1) =>
      createWebGPURenderItem<ComputeRenderVariant>(
        drawContext,
        drawMode,
        shaderCode,
        schema,
        values,
        materialId,
        undefined,
      ),
    createGlobalBindings: (schema, uniforms, textures): GlobalBindings =>
      createWebGPUGlobalBindings(drawContext, schema, uniforms, textures),

    createRenderTarget: (width, height, depthStencil, type, filter, format) =>
      createRenderTargetFromResources(resources, width, height, depthStencil, type, filter, format),
    createDrawTarget: () => {
      const target: RenderTarget = {
        id: -1,
        texture: createNullTexture(),
        framebuffer: createWebGPUFramebuffer(),
        depthRenderbuffer: null,
        getByteCount: () => 0,
        getWidth: () => drawingBufferWidth,
        getHeight: () => drawingBufferHeight,
        bind: () => {},
        setSize: () => {},
        destroy: () => {},
      };
      drawTargets.add(target);
      return target;
    },
    getDrawingBufferSize: () =>
      options.canvas
        ? { width: options.canvas.width, height: options.canvas.height }
        : { width: drawingBufferWidth, height: drawingBufferHeight },
    setPreserveDrawingBuffer: (value: boolean) => {
      if (!canvasContext || value === preserveDrawingBuffer) return;
      preserveDrawingBuffer = value;
      if (value) return;
      colorToPreserve = undefined;
      preservedColor?.destroy();
      preservedColor = undefined;
    },
    setDrawingBufferSize: (width: number, height: number) => {
      // with a canvas the drawing buffer follows the canvas size
      if (canvasContext) return;
      if (width === drawingBufferWidth && height === drawingBufferHeight) return;
      if (drawingBuffer) {
        release(drawingBuffer.color);
        release(drawingBuffer.depth);
      }
      drawingBuffer = undefined;
      drawingBufferWidth = width;
      drawingBufferHeight = height;
    },

    describeTarget,
    beginRenderPass: (target: RenderPassTarget, options?: RenderPassOptions) => {
      endComputePass();
      endPassEncoder();
      timer.passBegin(target, options?.label);
      passTarget = getPassTarget(target);
      if (options?.depthReadOnly) passTarget = withDepthReadOnly(passTarget);
      passClear = NoClear;
      if (options?.clear) {
        checkReadOnlyClear(options.clear);
        passClear = { ...options.clear };
      }
    },
    endRenderPass: () => {
      endPassEncoder();
      if (passTarget) timer.passEnd();
      passTarget = undefined;
    },
    hasRenderPassStencil: () => !!passTarget?.hasStencil,

    setRenderState: (state: RenderState) => {
      renderState = state;
    },
    setBlendConstant: (r: number, g: number, b: number, a: number) => {
      blendConstant = { r, g, b, a };
      if (passEncoder) flushBundles();
      passEncoder?.setBlendConstant(blendConstant);
    },
    setViewport: (x: number, y: number, width: number, height: number) => {
      viewport = { x, y, width, height };
    },
    setScissor: (x: number, y: number, width: number, height: number) => {
      scissor = { x, y, width, height };
    },
    disableScissor: () => {
      scissor = undefined;
    },
    clear: (values: ClearValues) => {
      const target = passTarget;
      if (!target) throw new Error('clear outside of a render pass');
      checkReadOnlyClear(values);
      if (capturing) {
        capturing.encoder.draw(0, 0, 0, 0);
        return;
      }
      if (recording) {
        // a clear draw of the bundle covers the pass's viewport: the same as the pass's clear (the scissor) where the
        // scissor lies inside the viewport (as for the clears between solid interior caps)
        const recorded =
          isScissorInViewport(target) &&
          clearDraw.drawInBundle(recording.state.encoder as GPURenderBundleEncoder, target, values);
        if (recorded) recording.state.invalidate();
        else recording.failed ??= 'a clear';
        return;
      }
      if (!passEncoder && coversTarget(scissor, target)) {
        // before the first draw, whole-target clears are load operations of the pass
        passClear = { ...passClear, ...values };
      } else {
        const state = getPassState();
        clearDraw.draw(state.encoder as GPURenderPassEncoder, target, values);
        // the clear draws with pipelines and bind groups of its own
        state.invalidate();
      }
    },
    flush: submit,

    readDrawingBuffer: async () => {
      const { color } = getDrawingBuffer();
      const data = await readGPUTexture(color, 4);
      if (color.format === 'bgra8unorm') swapRedBlue(data);
      return new ImageData(new Uint8ClampedArray(data.buffer), color.width, color.height);
    },
    waitForGpuCommandsComplete: async () => {
      submit();
      await device.queue.onSubmittedWorkDone();
    },
    readRenderTarget: async (target: RenderTarget, x: number, y: number, width: number, height: number) => {
      const texture = drawTargets.has(target)
        ? getDrawingBuffer().color
        : asWebGPUTexture(target.texture).getGPUTexture();
      if (!texture) throw new Error('render target without texture');
      if (texture.format !== 'rgba8unorm' && texture.format !== 'bgra8unorm') {
        notSupported(`reading render targets of format '${texture.format}'`);
      }
      const data = await readGPUTexture(texture, 4, { x, y, width, height });
      if (texture.format === 'bgra8unorm') swapRedBlue(data);
      return data;
    },
    readTextureAsync: async (texture: Texture, x: number, y: number, width: number, height: number) => {
      const t = asWebGPUTexture(texture);
      const gpuTexture = t.getGPUTexture();
      if (!gpuTexture) throw new Error('texture without storage');
      if (t.format !== 'rgba8unorm') notSupported(`reading textures of format '${t.format}'`);
      return readGPUTexture(gpuTexture, 4, { x, y, width, height });
    },
    createTextureReader: (options) =>
      createWebGPUTextureReader((texture, format, x, y, width, height) => {
        const t = asWebGPUTexture(texture);
        const gpuTexture = t.getGPUTexture();
        if (!gpuTexture) throw new Error('texture without storage');
        if (t.format !== ReaderFormats[format]) notSupported(`reading textures of format '${t.format}' as '${format}'`);
        return readGPUTexture(gpuTexture, 4, { x, y, width, height });
      }, options),
    readTextureData: async (texture: Texture) => {
      const t = asWebGPUTexture(texture);
      const gpuTexture = t.getGPUTexture();
      const info = ReadableFormats[t.format];
      if (!gpuTexture || !info) notSupported(`reading textures of format '${t.format}'`);
      const { width, height } = gpuTexture;
      // rows in texture order, which is data order (as uploaded)
      const data = await readGPUTexture(gpuTexture, info.bytes);
      const count = width * height;
      const source =
        info.type === 'ubyte'
          ? data
          : info.type === 'float'
            ? new Float32Array(data.buffer, data.byteOffset, count * info.channels)
            : new Uint16Array(data.buffer, data.byteOffset, count * info.channels);
      const array = info.type === 'ubyte' ? new Uint8Array(count * 4) : new Float32Array(count * 4);
      // four values per texel, single-channel formats in `r`
      for (let i = 0; i < count; ++i) {
        for (let c = 0; c < info.channels; ++c) {
          const v = source[i * info.channels + c];
          array[i * 4 + c] = info.type === 'half' ? fromHalfFloat(v) : v;
        }
      }
      return { array, width, height };
    },
    present: () => {
      // rendering goes into the canvas texture, which the browser presents
      submit();
    },
    takeErrors: () => errors.splice(0, errors.length),

    recreate: async () => {
      if (!isContextLost || destroyed || !options.requestDevice) return undefined;
      // right after a real loss there may be no device yet (the GPU process restarting): asked again a few times, about
      // 3 s in all (as long as `Canvas3D` waits for a restore before falling back)
      let next = await options.requestDevice();
      for (let delay = 100; !next && !destroyed && delay <= 1600; delay *= 2) {
        await new Promise((resolve) => setTimeout(resolve, delay));
        if (!destroyed) next = await options.requestDevice();
      }
      if (!next || destroyed) {
        next?.destroy();
        return undefined;
      }
      // on the same canvas, with the options of this one
      return createWebGPUContext(next, options);
    },

    destroy: (options?: GraphicsDestroyOptions) => {
      preservedColor?.destroy();
      preservedColor = undefined;
      destroyed = true;
      pending.clear();
      contextLost.complete();
      // the resources of a lost device are gone with it; its canvas is the one a restored context configures
      detachDevice?.();
      detachDevice = undefined;
      if (options?.lost) {
        if (contextOptions.ownsDevice) device.destroy();
        return;
      }
      submit();
      timer.destroy();
      timerSource?.destroy();
      ring.destroy();
      globalRing.destroy();
      geometryPages.destroy();
      dataPages.destroy();
      uniformPages.destroy();
      if (!canvasContext) drawingBuffer?.color.destroy();
      zeroSizeBuffer?.color.destroy();
      zeroSizeBuffer?.depth.destroy();
      drawingBuffer?.depth.destroy();
      if (canvasConfigured) canvasContext?.unconfigure();
      // the whole device's memory, which garbage collection would otherwise keep until it runs
      if (contextOptions.ownsDevice) device.destroy();
    },
  };
  logGraphicsContextInfo(ctx, {
    adapterInfo: toPlainObject(device.adapterInfo),
    limits: toPlainObject(device.limits),
    features: Array.from(device.features).sort(),
  });
  return ctx;
}
