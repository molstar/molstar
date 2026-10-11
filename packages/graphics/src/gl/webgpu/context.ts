/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';
import { isDebugMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import type { GraphicsCapabilities } from '../capabilities.js';
import type {
  GraphicsContext,
  GraphicsResources,
  GraphicsXR,
  RenderPassOptions,
  RenderPassTarget,
} from '../context.js';
import type { ComputeRenderVariant, DrawMode, GlobalBindings, GraphicsRenderVariant } from '../render-item.js';
import type { ClearValues, RenderState } from '../render-state.js';
import { createGraphicsStats, type GraphicsStats } from '../stats.js';
import { createGraphicsTimer, type TimerSource } from '../timer.js';
import { createWebGPUTimerSource } from './timer.js';
import { BufferUsage, MapMode, TextureUsage } from './constants.js';
import { notSupported } from './errors.js';
import { createWebGPUFramebuffer, isWebGPUFramebuffer, type WebGPUFramebuffer } from './framebuffer.js';
import { createWebGPUGlobalBindings, type WebGPUBindingContext } from './global-bindings.js';
import { createWebGPUProgramCache, type WebGPUProgram, type WebGPUProgramCache } from './program.js';
import { createWebGPURenderItem } from './render-item.js';
import { createClearDraw, type PassTarget } from './render-pass.js';
import { getColorTargetState, getDepthStencilState, getPrimitiveState, getRenderStateId } from './render-state.js';
import { createMipmapGenerator } from './mipmap.js';
import {
  asWebGPUTexture,
  CubeLayerSides,
  createWebGPUCubeTexture,
  createWebGPUTexture,
  type Release,
  type WebGPUResourceEnv,
  type WebGPURestorable,
} from './texture.js';
import { createWebGPURenderbuffer } from './renderbuffer.js';
import { createRenderTargetFromResources, type RenderTarget } from '../render-target.js';
import { type CubeFaces, createNullTexture, decodeImage, type Texture } from '../texture.js';
import { createPendingWork } from '../pending-work.js';
import { ValueCell } from '@molstar/core/util';
import { createUniformRing } from './uniform-ring.js';
import { getWgslProgram } from './shader/programs.js';
import type { WgslSampler } from './wgsl-declarations.js';

const DrawingBufferFormat: GPUTextureFormat = 'rgba8unorm';
const DepthStencilFormat: GPUTextureFormat = 'depth24plus-stencil8';

/** Internals the WebGPU render items draw with */
export interface WebGPUDrawContext extends WebGPUBindingContext, WebGPUResourceEnv {
  readonly stats: GraphicsStats;
  /** The encoder of the open render pass, begun on first use */
  getPassEncoder: () => GPURenderPassEncoder;
  /** Pipeline of `program` for the current render state and pass target, `undefined` while it compiles */
  getPipeline: (program: WebGPUProgram, drawMode: DrawMode) => GPURenderPipeline | undefined;
  /** Programs shared by the render items of the context */
  readonly programs: WebGPUProgramCache;
  /** Current render state, for dynamic state (stencil reference) */
  readonly renderState: RenderState;
  onError: (message: string) => void;
  /** Whether render items of programs that are not ported throw instead of drawing nothing */
  readonly strict: boolean;
  /** Warns once per message about a render item that draws nothing */
  warnUnsupported: (message: string) => void;
}

export interface WebGPUContext extends GraphicsContext {
  readonly device: GPUDevice;
  /**
   * Reads the drawing buffer, rows top to bottom; ends the open pass and submits pending commands. With a canvas, call it
   * in the task that rendered the frame (the canvas texture changes per frame).
   */
  readDrawingBuffer: () => Promise<ImageData>;
  /**
   * Reads a texture, rows top to bottom (WebGPU texture order); `Uint8Array` for 8-bit formats, `Float32Array`
   * for 32-bit float formats. Ends the open pass and submits pending commands.
   */
  readTexture: (texture: Texture) => Promise<Uint8Array | Float32Array>;
  /** Submits the recorded commands now (otherwise at the end of the task); the browser presents the canvas texture */
  present: () => void;
  /** Returns and clears device and shader compilation errors */
  takeErrors: () => string[];
  setDrawingBufferSize: (width: number, height: number) => void;
}

/** Features used where available, see `capabilities` */
const OptionalFeatures: GPUFeatureName[] = [
  'float32-filterable',
  'float32-blendable',
  'depth-clip-control',
  'depth32float-stencil8',
  // GPU timer
  'timestamp-query',
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

/** The device of an adapter of `gpu` (by default the browser's `navigator.gpu`, e.g. Dawn's in Node.js) */
export async function requestWebGPUDevice(
  powerPreference: 'default' | 'high-performance' | 'low-power' = 'high-performance',
  gpu: GPU | undefined = typeof navigator === 'undefined' ? undefined : navigator.gpu,
): Promise<GPUDevice | undefined> {
  if (!gpu) return;
  const adapter = await gpu.requestAdapter(powerPreference === 'default' ? undefined : { powerPreference });
  if (!adapter) return;
  const requiredFeatures = OptionalFeatures.filter((f) => adapter.features.has(f));
  return adapter.requestDevice({ requiredFeatures });
}

type Rect = { x: number; y: number; width: number; height: number };

/**
 * WebGPU graphics context. The drawing buffer is the canvas texture of the current frame (`rgba8unorm` with a
 * `depth24plus-stencil8` attachment), or an offscreen texture without canvas. Recorded commands are submitted at the
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
  /** Requests a device to replace a lost one, see `restore`; without it a lost device is not replaced */
  requestDevice?: () => Promise<GPUDevice | undefined>;
}

export function createWebGPUContext(initialDevice: GPUDevice, options: WebGPUContextOptions): WebGPUContext {
  /** the current device, replaced by `restore` after a loss */
  let device = initialDevice;
  /** `device` for components that read it on use */
  const deviceRef = {
    get device() {
      return device;
    },
  };
  let destroyed = false;
  /** pending errors until `takeErrors`, the oldest are dropped when nobody takes them */
  const errors: string[] = [];
  const onError = (message: string) => {
    if (isDebugMode) console.error(message);
    if (errors.length >= 100) errors.shift();
    errors.push(message);
  };
  const contextLost = new Subject<now.Timestamp>();
  const contextRestored = new Subject<now.Timestamp>();
  let isContextLost = false;
  const setContextLost = () => {
    if (isContextLost) return;
    isContextLost = true;
    contextLost.next(now());
  };

  function attachDevice(d: GPUDevice) {
    d.addEventListener('uncapturederror', (e) => {
      onError(`WebGPU error: ${(e as GPUUncapturedErrorEvent).error.message}`);
    });
    d.lost.then((info) => {
      // losses of replaced devices and of destroyed contexts are expected
      if (d !== device || destroyed) return;
      onError(`WebGPU device lost (${info.reason}): ${info.message}`);
      setContextLost();
    });
  }
  attachDevice(device);

  const canvasContext = options.canvas?.getContext('webgpu') as GPUCanvasContext | null | undefined;
  function configureCanvas() {
    canvasContext?.configure({
      device,
      format: DrawingBufferFormat,
      usage: TextureUsage.RENDER_ATTACHMENT | TextureUsage.COPY_SRC,
      alphaMode: 'premultiplied',
    });
  }
  configureCanvas();

  const stats = createGraphicsStats();
  const pending = createPendingWork();
  let timerSource = createWebGPUTimerSource(device);
  /** the timer keeps its source, which changes with the device */
  const timerSourceRef: TimerSource | undefined = timerSource && {
    begin: () => timerSource?.begin() ?? -1,
    end: () => timerSource?.end(),
    poll: () => timerSource?.poll() ?? [],
    clear: () => timerSource?.clear(),
  };
  const timer = createGraphicsTimer(timerSourceRef, stats);
  const ring = createUniformRing(deviceRef);
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
  let drawingBufferWidth = options.width;
  let drawingBufferHeight = options.height;
  let pixelScale = options.pixelScale ?? 1;

  function createDrawingBufferTarget(color: GPUTexture, depth: GPUTexture): PassTarget {
    return {
      key: `${color.format}|${DepthStencilFormat}`,
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

  /** The canvas texture of the current frame, or an offscreen texture without canvas */
  function getDrawingBuffer() {
    if (canvasContext) {
      // a new texture per frame, presented when the task that rendered into it ends
      const color = canvasContext.getCurrentTexture();
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
        format: DrawingBufferFormat,
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
          key: `${DrawingBufferFormat}|${DepthStencilFormat}`,
          colorViews: [color.createView()],
          colorFormats: [DrawingBufferFormat],
          depthView: depth.createView(),
          depthFormat: DepthStencilFormat,
          hasStencil: true,
          ...size,
        },
      };
    }
    return drawingBuffer;
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
      key: `${colorFormats.join(',')}|${depth?.format ?? ''}`,
      colorViews,
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

  function getPassTarget(target: RenderPassTarget): PassTarget {
    if (target === 'drawing-buffer') return getDrawingBuffer().target;
    const framebuffer = 'framebuffer' in target ? target.framebuffer : target;
    if ('framebuffer' in target && drawTargets.has(target)) return getDrawingBuffer().target;
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
    const buffer = device.createBuffer({
      size: bytesPerRow * height,
      usage: BufferUsage.COPY_DST | BufferUsage.MAP_READ,
    });
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
    buffer.destroy();
    return out;
  }

  const ReadableFormats: { [k: string]: { bytes: number; float: boolean } } = {
    r8unorm: { bytes: 1, float: false },
    rg8unorm: { bytes: 2, float: false },
    rgba8unorm: { bytes: 4, float: false },
    r32float: { bytes: 4, float: true },
    rg32float: { bytes: 8, float: true },
    rgba32float: { bytes: 16, float: true },
  };

  // passes

  let commandEncoder: GPUCommandEncoder | undefined;
  const encode = (commands: (encoder: GPUCommandEncoder) => void) => {
    // copies run outside of render passes and end the open one
    endPassEncoder();
    if (passTarget) timer.passEnd();
    passTarget = undefined;
    commands(getCommandEncoder());
  };
  const pendingReleases: { destroy(): void }[] = [];
  const release: Release = (resource) => {
    if (commandEncoder) pendingReleases.push(resource);
    else resource.destroy();
  };
  let passEncoder: GPURenderPassEncoder | undefined;
  /** target of the open pass, the encoder is begun on first draw */
  let passTarget: PassTarget | undefined;
  let passClear: { color?: ClearValues['color']; depth?: number; stencil?: number } = {};

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
    const descriptor: GPURenderPassDescriptor = {
      colorAttachments: target.colorViews.map((view) => ({
        view,
        loadOp: color ? 'clear' : 'load',
        storeOp: 'store',
        clearValue: color ? { r: color[0], g: color[1], b: color[2], a: color[3] } : undefined,
      })),
    };
    if (target.depthView) {
      descriptor.depthStencilAttachment = {
        view: target.depthView,
        depthLoadOp: depth !== undefined ? 'clear' : 'load',
        depthStoreOp: 'store',
        depthClearValue: depth ?? 1,
        stencilLoadOp: target.hasStencil ? (stencil !== undefined ? 'clear' : 'load') : undefined,
        stencilStoreOp: target.hasStencil ? 'store' : undefined,
        stencilClearValue: stencil ?? 0,
      };
    }
    passClear = {};
    const timestampWrites = timerSource?.getTimestampWrites();
    if (timestampWrites) descriptor.timestampWrites = timestampWrites;
    const pass = getCommandEncoder().beginRenderPass(descriptor);
    pass.setBlendConstant(blendConstant);
    return pass;
  }

  function endPassEncoder() {
    if (passTarget && !passEncoder && hasClear(passClear)) {
      // a pass with clears only
      passEncoder = beginPassEncoder(passTarget);
    }
    passEncoder?.end();
    passEncoder = undefined;
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

  /** bottom-left origin (GL convention of the neutral API) to top-left, clamped to the target */
  function toTargetRect(rect: Rect, target: PassTarget): Rect {
    const x = Math.max(0, Math.min(rect.x, target.width));
    const y = Math.max(0, Math.min(target.height - rect.y - rect.height, target.height));
    const width = Math.max(0, Math.min(rect.width, target.width - x));
    const height = Math.max(0, Math.min(rect.height, target.height - y));
    return { x, y, width, height };
  }

  function getPassEncoder() {
    if (!passTarget) throw new Error('draw outside of a render pass');
    if (!passEncoder) passEncoder = beginPassEncoder(passTarget);
    const target = passTarget;
    if (viewport) {
      const r = toTargetRect(viewport, target);
      passEncoder.setViewport(r.x, r.y, r.width, r.height, 0, 1);
    }
    if (scissor) {
      const r = toTargetRect(scissor, target);
      passEncoder.setScissorRect(r.x, r.y, r.width, r.height);
    } else {
      passEncoder.setScissorRect(0, 0, target.width, target.height);
    }
    passEncoder.setStencilReference(renderState.stencil?.reference ?? 0);
    return passEncoder;
  }

  /** pipeline keys compiling in the background, per program */
  const compilingPipelines = new WeakMap<WebGPUProgram, Set<string>>();
  /** pipeline keys whose background compilation failed: created synchronously, which reports the error once */
  const failedPipelines = new WeakMap<WebGPUProgram, Set<string>>();

  /**
   * The pipeline of `program` for the current render state and pass target, `undefined` while it compiles. Missing
   * pipelines compile in the background (`createRenderPipelineAsync`, off the GPU process' main thread) and
   * their draws are skipped until ready, so first uses (e.g. the first hover: picking and marking) don't stall frames;
   * the compilation is pending work, `Canvas3D` redraws when it settles and `ImagePass` renders again
   */
  function getPipeline(program: WebGPUProgram, drawMode: DrawMode): GPURenderPipeline | undefined {
    const target = passTarget;
    if (!target) throw new Error('draw outside of a render pass');
    const state = renderState;
    const key = `${getRenderStateId(state)}|${target.key}|${drawMode}`;
    const existing = program.pipelines.get(key);
    if (existing) return existing;
    const descriptor: GPURenderPipelineDescriptor = {
      label: `${program.name} ${key}`,
      layout: program.pipelineLayout,
      vertex: {
        module: program.module,
        entryPoint: 'vertexMain',
        buffers: [...program.vertexBuffers],
        constants: program.constants,
      },
      fragment: {
        module: program.module,
        entryPoint: 'fragmentMain',
        targets: target.colorFormats.map((f) => getColorTargetState(state, f)),
        constants: program.constants,
      },
      primitive: getPrimitiveState(state, drawMode),
      depthStencil: target.depthFormat ? getDepthStencilState(state, target.depthFormat, target.hasStencil) : undefined,
    };
    if (failedPipelines.get(program)?.has(key)) {
      const pipeline = device.createRenderPipeline(descriptor);
      program.pipelines.set(key, pipeline);
      return pipeline;
    }
    let compiling = compilingPipelines.get(program);
    if (!compiling) compilingPipelines.set(program, (compiling = new Set()));
    if (!compiling.has(key)) {
      compiling.add(key);
      const pipelineDevice = device;
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
        compiled.finally(() => keys.delete(key)),
        'pipeline',
      );
    }
    return undefined;
  }

  function submit() {
    endPassEncoder();
    if (!commandEncoder) return;
    ring.upload();
    const afterSubmit = timerSource?.beforeSubmit(commandEncoder);
    device.queue.submit([commandEncoder.finish()]);
    afterSubmit?.();
    commandEncoder = undefined;
    ring.reset();
    for (const r of pendingReleases) r.destroy();
    pendingReleases.length = 0;
  }

  let clearDraw = createClearDraw(device, ring);

  const programs = createWebGPUProgramCache(deviceRef, onError, () => {
    stats.resourceCounts.program += 1;
  });

  let generateMipmaps = createMipmapGenerator(device);
  /** resources recreated by `restore` */
  const restorables = new Set<WebGPURestorable>();
  const env: WebGPUResourceEnv = {
    get device() {
      return device;
    },
    restorables,
    release,
    encode,
    get generateMipmaps() {
      return generateMipmaps;
    },
  };

  const strict = options.strict ?? false;
  const warned = new Set<string>();

  const drawContext: WebGPUDrawContext = {
    get generateMipmaps() {
      return generateMipmaps;
    },
    encode,
    get device() {
      return device;
    },
    restorables,
    programs,
    release,
    ring,
    samplers,
    stats,
    getPassEncoder,
    getPipeline,
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
  };

  // stubs

  const xr: GraphicsXR = {
    session: undefined,
    changed: new Subject<void>(),
    set: async () => notSupported('XR'),
    end: async () => {},
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
    // blending into float32 targets (wboit, dpoit) needs the feature, as WebGL needs EXT_float_blend
    floatRenderTarget: device.features.has('float32-blendable'),
    halfFloatRenderTarget: true,
    depthTexture: true,
    fragDepth: true,
    multipleRenderTargets: true,
    textureLod: true,
    multiDraw: false,
    clipControl: true,
    depthClamp: device.features.has('depth-clip-control'),
  };

  const resources: GraphicsResources = {
    limits: {
      maxTextureSize: device.limits.maxTextureDimension2D,
      max3dTextureSize: device.limits.maxTextureDimension3D,
      maxRenderbufferSize: device.limits.maxTextureDimension2D,
      maxDrawBuffers: device.limits.maxColorAttachments,
      maxTextureImageUnits: device.limits.maxSampledTexturesPerShaderStage,
    },
    texture: (kind, format, type, filter) => {
      stats.resourceCounts.texture += 1;
      return createWebGPUTexture(env, kind, format, type, filter);
    },
    cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: (errored?: boolean) => void) => {
      stats.resourceCounts.cubeTexture += 1;
      const texture = createWebGPUCubeTexture(env, mipmaps);
      let destroyed = false;
      const load = async () => {
        const images = await Promise.all(CubeLayerSides.map((side) => decodeImageBitmap(faces[side])));
        if (destroyed) return;
        texture.loadFaces(images);
      };
      pending.track(load(), 'cube texture').then(
        () => onload?.(destroyed),
        () => onload?.(true),
      );
      const destroy = texture.destroy;
      texture.destroy = () => {
        destroyed = true;
        destroy();
      };
      return texture;
    },
    framebuffer: () => {
      stats.resourceCounts.framebuffer += 1;
      return createWebGPUFramebuffer();
    },
    renderbuffer: (format, attachment, width, height) => {
      stats.resourceCounts.renderbuffer += 1;
      return createWebGPURenderbuffer(env, format, attachment, width, height);
    },
    uniformBuffer: () => notSupported('uniformBuffer'),
    uniformRingBuffer: () => notSupported('uniformRingBuffer'),
    materialUniformBuffer: () => notSupported('materialUniformBuffer'),
    getByteCounts: () => ({
      texture: 0,
      cubeTexture: 0,
      attribute: 0,
      elements: 0,
      pixelPack: 0,
      uniform: 0,
      renderbuffer: 0,
    }),
    pending,
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
    linkPrograms: () => {},
    // frames don't wait for programs: draws of pipelines compiling in the background are skipped, see `getPipeline`
    finalizePrograms: () => true,
  };

  /** recreates the per-device state and the resources on `next`, which replaces the lost device */
  function replaceDevice(next: GPUDevice) {
    device = next;
    attachDevice(next);
    configureCanvas();
    // commands and attachments of the lost device are gone
    commandEncoder = undefined;
    passEncoder = undefined;
    passTarget = undefined;
    passClear = {};
    pendingReleases.length = 0;
    drawingBuffer = undefined;
    timerSource?.destroy();
    timerSource = createWebGPUTimerSource(next);
    timer.clear();
    createSamplers();
    generateMipmaps = createMipmapGenerator(next);
    clearDraw = createClearDraw(next, ring);
    ring.restore();
    programs.clear();
    for (const r of [...restorables]) r.reset();
  }

  return {
    get device() {
      return device;
    },
    capabilities,
    resources,
    stats,
    timer,
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
    contextRestored,
    setContextLost,
    handleContextRestored: (extraResets?: () => void) => {
      extraResets?.();
      isContextLost = false;
      contextRestored.next(now());
    },
    restore: async (extraResets?: () => void) => {
      if (!isContextLost || destroyed || !options.requestDevice) return false;
      const next = await options.requestDevice();
      if (!next || destroyed) return false;
      replaceDevice(next);
      extraResets?.();
      isContextLost = false;
      contextRestored.next(now());
      return true;
    },

    xr,

    namedFramebuffers: Object.create(null),
    namedTextures: Object.create(null),
    namedComputeRenderables: Object.create(null),

    isProgramSupported: (shaderCode) => strict || getWgslProgram(shaderCode.name) !== undefined,
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
        reset: () => {},
        destroy: () => {},
      };
      drawTargets.add(target);
      return target;
    },
    getDrawingBufferSize: () =>
      options.canvas
        ? { width: options.canvas.width, height: options.canvas.height }
        : { width: drawingBufferWidth, height: drawingBufferHeight },
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

    beginRenderPass: (target: RenderPassTarget, options?: RenderPassOptions) => {
      endPassEncoder();
      timer.passBegin(target, options?.label);
      passTarget = getPassTarget(target);
      passClear = {};
      if (options?.clear) passClear = { ...options.clear };
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
      if (!passEncoder && coversTarget(scissor, target)) {
        // before the first draw, whole-target clears are load operations of the pass
        passClear = { ...passClear, ...values };
      } else {
        clearDraw(getPassEncoder(), target, values);
      }
    },
    flush: submit,

    readDrawingBuffer: async () => {
      const { color } = getDrawingBuffer();
      const data = await readGPUTexture(color, 4);
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
      if (texture.format !== 'rgba8unorm') notSupported(`reading render targets of format '${texture.format}'`);
      // the rectangle has a bottom-left origin (GL convention of the neutral API)
      return readGPUTexture(texture, 4, { x, y: texture.height - y - height, width, height });
    },
    readTextureAsync: async (texture: Texture, x: number, y: number, width: number, height: number) => {
      const t = asWebGPUTexture(texture);
      const gpuTexture = t.getGPUTexture();
      if (!gpuTexture) throw new Error('texture without storage');
      if (t.format !== 'rgba8unorm') notSupported(`reading textures of format '${t.format}'`);
      // the rectangle has a bottom-left origin (GL convention of the neutral API)
      return readGPUTexture(gpuTexture, 4, { x, y: gpuTexture.height - y - height, width, height });
    },
    readTexture: async (texture: Texture) => {
      const t = asWebGPUTexture(texture);
      const gpuTexture = t.getGPUTexture();
      const info = ReadableFormats[t.format];
      if (!gpuTexture || !info) notSupported(`reading textures of format '${t.format}'`);
      const data = await readGPUTexture(gpuTexture, info.bytes);
      return info.float ? new Float32Array(data.buffer) : data;
    },
    present: () => {
      // rendering goes into the canvas texture, which the browser presents
      submit();
    },
    takeErrors: () => errors.splice(0, errors.length),

    destroy: () => {
      destroyed = true;
      pending.clear();
      submit();
      timer.destroy();
      timerSource?.destroy();
      ring.destroy();
      if (!canvasContext) drawingBuffer?.color.destroy();
      drawingBuffer?.depth.destroy();
      contextLost.complete();
      contextRestored.complete();
    },
  };
}
