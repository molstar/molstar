/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { type GLRenderingContext, isWebGL2 } from './compat.js';
import { checkFramebufferStatus, createNullFramebuffer, getFramebufferHeight } from './framebuffer.js';
import { Scheduler } from '@molstar/core/task';
import { isDebugMode } from '@molstar/core/util/debug';
import { createExtensions, type WebGLExtensions } from './extensions.js';
import { type WebGLState, createState } from './state.js';
import { type WebGLResources, createResources } from './resources.js';
import { createRenderTarget } from './render-target.js';
import type { RenderTarget } from '../render-target.js';
import { Subject } from 'rxjs';
import { now } from '@molstar/core/util/now';
import { PixelData } from '@molstar/core/util/image';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { Viewport } from '@molstar/core/util/viewport';
import { createNullTexture, type Texture, type TextureFilter } from '../texture.js';
import type { Framebuffer } from '../framebuffer.js';
import { createTimer } from './timer.js';
import {
  type GraphicsContext,
  type GraphicsDestroyOptions,
  type GraphicsXRView,
  logGraphicsContextInfo,
  type RenderPassOptions,
  type RenderPassTarget,
  type TargetDescriptor,
} from '../context.js';
import { endRenderPass } from './render-pass.js';
import { createGlobalBindings } from './global-bindings.js';
import { getWebGLCapabilities } from './capabilities.js';
import { createGraphicsStats } from '../stats.js';
import type { ComputeRenderVariant, DrawMode, GraphicsRenderVariant, Transparency } from '../render-item.js';
import type { RenderableSchema, RenderableValues } from '../renderable/schema.js';
import type { ShaderCode } from '../shader-code.js';
import {
  createComputeRenderItem,
  createGraphicsRenderItem,
  createWebGLProgramItem,
  type WebGLRenderItem,
} from './render-item.js';
import { createRenderStateApplier } from './render-state.js';
import { asGLTexture } from './texture.js';
import { getSharedCopyRenderable } from '../compute/util.js';
import { RenderState } from '../render-state.js';
import { registerAtlasIsosurfaceExtraction } from '../compute/marching-cubes/atlas.js';
import { WebGLAtlasIsosurfaceExtraction } from './compute/marching-cubes/isosurface.js';
import { registerGrid3dFragmentImplementation } from '../compute/grid3d.js';
import { WebGLGrid3dFragment } from './compute/grid3d.js';
import { createWebGLTextureReader } from './texture-reader.js';
import { registerRenderPassTargetHeight } from './render-pass-target.js';
import { type DepthRangeConversion, ReversedMinusOneToOne, ReversedZeroToOne } from '../depth-range.js';

export function getGLContext(canvas: HTMLCanvasElement, attribs?: WebGLContextAttributes): GLRenderingContext | null {
  function get(id: 'webgl' | 'experimental-webgl' | 'webgl2') {
    try {
      return canvas.getContext(id, attribs) as GLRenderingContext | null;
    } catch (e) {
      return null;
    }
  }
  const gl = get('webgl2') || get('webgl') || get('experimental-webgl');
  if (isDebugMode) console.log(`isWebgl2: ${isWebGL2(gl)}`);
  return gl;
}

export function getErrorDescription(gl: GLRenderingContext, error: number) {
  switch (error) {
    case gl.NO_ERROR:
      return 'no error';
    case gl.INVALID_ENUM:
      return 'invalid enum';
    case gl.INVALID_VALUE:
      return 'invalid value';
    case gl.INVALID_OPERATION:
      return 'invalid operation';
    case gl.INVALID_FRAMEBUFFER_OPERATION:
      return 'invalid framebuffer operation';
    case gl.OUT_OF_MEMORY:
      return 'out of memory';
    case gl.CONTEXT_LOST_WEBGL:
      return 'context lost';
  }
  return 'unknown error';
}

export function checkError(gl: GLRenderingContext, message?: string) {
  const error = gl.getError();
  if (error !== gl.NO_ERROR) {
    console.log(`WebGL error: '${getErrorDescription(gl, error)}'${message ? ` (${message})` : ''}`);
    // throw new Error(`WebGL error: '${getErrorDescription(gl, error)}'${message ? ` (${message})` : ''}`);
  }
}

export function glEnumToString(gl: GLRenderingContext, value: number) {
  const keys: string[] = [];
  for (const key in gl) {
    if ((gl as any)[key] === value) {
      keys.push(key);
    }
  }
  return keys.length ? keys.join(' | ') : `0x${value.toString(16)}`;
}

function unbindResources(gl: GLRenderingContext) {
  // bind null to all texture units
  const maxTextureImageUnits = gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS);
  for (let i = 0; i < maxTextureImageUnits; ++i) {
    gl.activeTexture(gl.TEXTURE0 + i);
    gl.bindTexture(gl.TEXTURE_2D, null);
    gl.bindTexture(gl.TEXTURE_CUBE_MAP, null);
    if (isWebGL2(gl)) {
      gl.bindTexture(gl.TEXTURE_2D_ARRAY, null);
      gl.bindTexture(gl.TEXTURE_3D, null);
    }
  }

  // assign the smallest possible buffer to all attributes
  const buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  const maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
  for (let i = 0; i < maxVertexAttribs; ++i) {
    gl.vertexAttribPointer(i, 1, gl.FLOAT, false, 0, 0);
  }

  // bind null to all buffers
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
  gl.bindRenderbuffer(gl.RENDERBUFFER, null);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  if (isWebGL2(gl)) {
    gl.bindBuffer(gl.UNIFORM_BUFFER, null);
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
  }
}

/** Programs don't depend on the target on WebGL, all targets are alike */
const WebGLTarget: TargetDescriptor = { key: 'webgl', hasStencil: true };

const tmpPixel = new Uint8Array(1 * 4);

function checkSync(gl: WebGL2RenderingContext, sync: WebGLSync, resolve: () => void) {
  if (gl.getSyncParameter(sync, gl.SYNC_STATUS) === gl.SIGNALED) {
    gl.deleteSync(sync);
    resolve();
  } else {
    Scheduler.setImmediate(checkSync, gl, sync, resolve);
  }
}

function fence(gl: WebGL2RenderingContext, resolve: () => void) {
  const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
  if (!sync) {
    console.warn('Could not create a WebGLSync object');
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel);
    resolve();
  } else {
    Scheduler.setImmediate(checkSync, gl, sync, resolve);
  }
}

let SentWebglSyncObjectNotSupportedInWebglMessage = false;
function waitForGpuCommandsComplete(gl: GLRenderingContext): Promise<void> {
  return new Promise((resolve) => {
    if (isWebGL2(gl)) {
      fence(gl, resolve);
    } else {
      if (!SentWebglSyncObjectNotSupportedInWebglMessage) {
        console.info('Sync object not supported in WebGL');
        SentWebglSyncObjectNotSupportedInWebglMessage = true;
      }
      waitForGpuCommandsCompleteSync(gl);
      resolve();
    }
  });
}

function waitForGpuCommandsCompleteSync(gl: GLRenderingContext): void {
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel);
}

export function readPixels(
  gl: GLRenderingContext,
  x: number,
  y: number,
  width: number,
  height: number,
  buffer: Uint8Array | Float32Array | Int32Array,
) {
  if (isDebugMode) checkFramebufferStatus(gl);
  if (buffer instanceof Uint8Array) {
    gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer);
  } else if (buffer instanceof Float32Array) {
    gl.readPixels(x, y, width, height, gl.RGBA, gl.FLOAT, buffer);
  } else if (buffer instanceof Int32Array && isWebGL2(gl)) {
    gl.readPixels(x, y, width, height, gl.RGBA_INTEGER, gl.INT, buffer);
  } else {
    throw new Error('unsupported readPixels buffer type');
  }
  if (isDebugMode) checkError(gl);
}

function bindDrawingBuffer(gl: GLRenderingContext, xrLayer?: XRWebGLLayer) {
  if (xrLayer) {
    gl.bindFramebuffer(gl.FRAMEBUFFER, xrLayer.framebuffer);
  } else {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  }
}

function getDrawingBufferSize(gl: GLRenderingContext, xrLayer?: XRWebGLLayer, xrInteractionMode?: XRInteractionMode) {
  let width = xrLayer?.framebufferWidth ?? gl.drawingBufferWidth;
  if (xrInteractionMode === 'screen-space') {
    // workaround so XR with a single view behaves simlar to two views
    width *= 2;
  }
  const height = xrLayer?.framebufferHeight ?? gl.drawingBufferHeight;
  return { width, height };
}

function getShaderPrecisionFormat(
  gl: GLRenderingContext,
  shader: 'vertex' | 'fragment',
  precision: 'low' | 'medium' | 'high',
  type: 'float' | 'int',
) {
  const glShader = shader === 'vertex' ? gl.VERTEX_SHADER : gl.FRAGMENT_SHADER;
  const glPrecisionType =
    gl[
      `${precision.toUpperCase()}_${type.toUpperCase()}` as
        | 'LOW_FLOAT'
        | 'MEDIUM_FLOAT'
        | 'HIGH_FLOAT'
        | 'LOW_INT'
        | 'MEDIUM_INT'
        | 'HIGH_INT'
    ];
  return gl.getShaderPrecisionFormat(glShader, glPrecisionType);
}

function getShaderPrecisionFormats(gl: GLRenderingContext, shader: 'vertex' | 'fragment') {
  return {
    lowFloat: getShaderPrecisionFormat(gl, shader, 'low', 'float'),
    mediumFloat: getShaderPrecisionFormat(gl, shader, 'medium', 'float'),
    highFloat: getShaderPrecisionFormat(gl, shader, 'high', 'float'),
    lowInt: getShaderPrecisionFormat(gl, shader, 'low', 'int'),
    mediumInt: getShaderPrecisionFormat(gl, shader, 'medium', 'int'),
    highInt: getShaderPrecisionFormat(gl, shader, 'high', 'int'),
  };
}

type WebGLShaderPrecisionFormats = ReturnType<typeof getShaderPrecisionFormats>;

//

function createParameters(gl: GLRenderingContext, extensions: WebGLExtensions) {
  return {
    maxTextureSize: gl.getParameter(gl.MAX_TEXTURE_SIZE) as number,
    max3dTextureSize: isWebGL2(gl) ? (gl.getParameter(gl.MAX_3D_TEXTURE_SIZE) as number) : 0,
    maxRenderbufferSize: gl.getParameter(gl.MAX_RENDERBUFFER_SIZE) as number,
    maxDrawBuffers: extensions.drawBuffers ? (gl.getParameter(extensions.drawBuffers.MAX_DRAW_BUFFERS) as number) : 0,
    maxTextureImageUnits: gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS) as number,
    maxVertexTextureImageUnits: gl.getParameter(gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS) as number,
    maxVertexAttribs: gl.getParameter(gl.MAX_VERTEX_ATTRIBS) as number,
  };
}

export type WebGLParameters = ReturnType<typeof createParameters>;

//

/** A WebGL context object, including the rendering context, resource caches and counts */
export interface WebGLContext extends GraphicsContext {
  readonly gl: GLRenderingContext;
  readonly isWebGL2: boolean;
  /** Whether uniform blocks are backed by uniform buffers, see `WebGLContextProps` */
  readonly uniformBlocks: boolean;

  readonly extensions: WebGLExtensions;
  readonly state: WebGLState;
  readonly resources: WebGLResources;

  readonly shaderPrecisionFormats: { vertex: WebGLShaderPrecisionFormats; fragment: WebGLShaderPrecisionFormats };
  /**
   * How the projection uniforms are mapped to clip depth (see `DepthRangeTag`): reversed `[0, 1]` with
   * `EXT_clip_control` (as WebGPU, with the precision gain of reversed depth), reversed `[-1, 1]` otherwise; GLSL sees
   * the former as `clipZeroToOne`
   */
  readonly depthRangeConversion: DepthRangeConversion;

  createGraphicsRenderItem: (
    drawMode: DrawMode,
    shaderCode: ShaderCode,
    schema: RenderableSchema,
    values: RenderableValues,
    materialId: number,
    transparency: Transparency,
  ) => WebGLRenderItem<GraphicsRenderVariant>;
  createComputeRenderItem: (
    drawMode: DrawMode,
    shaderCode: ShaderCode,
    schema: RenderableSchema,
    values: RenderableValues,
    materialId?: number,
  ) => WebGLRenderItem<ComputeRenderVariant>;

  /** Binds the drawing buffer outside of render passes, e.g. for readback */
  bindDrawingBuffer: () => void;
  /** Reads from the bound framebuffer, ends the open render pass */
  readPixels: (
    x: number,
    y: number,
    width: number,
    height: number,
    buffer: Uint8Array | Float32Array | Int32Array,
  ) => void;
  waitForGpuCommandsComplete: () => Promise<void>;
  waitForGpuCommandsCompleteSync: () => void;
  getFenceSync: () => WebGLSync | null;
  checkSyncStatus: (sync: WebGLSync) => boolean;
  deleteSync: (sync: WebGLSync) => void;
  checkError: (message?: string) => void;
  checkFramebufferStatus: (message?: string) => void;
}

/** Narrows to the WebGL backend, for code paths that are WebGL-only so far */
export function isWebGLContext(ctx: GraphicsContext): ctx is WebGLContext {
  return (ctx as WebGLContext).gl !== undefined;
}

const ReadRgbaName = 'read-texture-rgba';

/** `texture` (ubyte, single-channel) copied into an rgba texture, the values in `r` */
function copyToRgbaTexture(ctx: WebGLContext, texture: Texture, width: number, height: number) {
  const copy = getSharedCopyRenderable(ctx, texture);
  ctx.namedFramebuffers[ReadRgbaName] ??= ctx.resources.framebuffer();
  const framebuffer = ctx.namedFramebuffers[ReadRgbaName];
  ctx.namedTextures[ReadRgbaName] ??= ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
  const rgba = ctx.namedTextures[ReadRgbaName];
  rgba.define(width, height);
  rgba.attachFramebuffer(framebuffer, 'color0');
  ctx.beginRenderPass(framebuffer, { label: 'read texture (rgba copy)' });
  ctx.setRenderState(RenderState.Default);
  ctx.setViewport(0, 0, width, height);
  ctx.setScissor(0, 0, width, height);
  ctx.clear({ color: [0, 0, 0, 0] });
  copy.render();
  return rgba;
}

export interface WebGLContextProps {
  readonly pixelScale: number;
  /** Uniform blocks backed by uniform buffers (WebGL2 only, on by default); off, all uniforms are set one by one */
  readonly uniformBlocks: boolean;
}

export function createContext(gl: GLRenderingContext, props: Partial<WebGLContextProps> = {}): WebGLContext {
  // GPU compute implemented with fragment shaders, reached through neutral entry points
  registerAtlasIsosurfaceExtraction(WebGLAtlasIsosurfaceExtraction);
  registerGrid3dFragmentImplementation(WebGLGrid3dFragment);
  const extensions = createExtensions(gl);
  const capabilities = getWebGLCapabilities(gl, extensions);
  const state = createState(gl, extensions);
  const stats = createGraphicsStats();
  const parameters = createParameters(gl, extensions);
  const uniformBlocks = isWebGL2(gl) && props.uniformBlocks !== false;
  const resources = createResources(gl, state, stats, extensions, parameters, uniformBlocks, () => isContextLost);
  /** for `readTextureAsync`, created on first use */
  let readFramebuffer: Framebuffer | undefined;
  const timer = createTimer(gl, extensions, stats);
  state.onRenderPassEnd = timer.passEnd;
  const renderState = createRenderStateApplier(gl, state, extensions);
  const drawingBufferHasStencil = !!gl.getContextAttributes()?.stencil;

  if (parameters.maxVertexTextureImageUnits < 8) {
    throw new Error('Need "MAX_VERTEX_TEXTURE_IMAGE_UNITS" >= 8');
  }

  const shaderPrecisionFormats = {
    vertex: getShaderPrecisionFormats(gl, 'vertex'),
    fragment: getShaderPrecisionFormats(gl, 'fragment'),
  };

  // fixed for the context's lifetime, the programs are compiled for it
  const clipZeroToOne = !!extensions.clipControl;
  const depthRangeConversion = clipZeroToOne ? ReversedZeroToOne : ReversedMinusOneToOne;

  // rectangles have a top-left origin in the neutral API, GL's a bottom-left one: viewport and scissor are kept as
  // given and applied flipped with the height of the render pass target, when set in a pass and when a pass begins
  let passTarget: RenderPassTarget = 'drawing-buffer';
  const getRenderPassTargetHeight = () => {
    if (passTarget === 'drawing-buffer') return getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode).height;
    if ('getHeight' in passTarget) return passTarget.getHeight();
    return getFramebufferHeight(passTarget);
  };
  const viewport = { x: 0, y: 0, width: 0, height: 0, set: false };
  const scissor = { x: 0, y: 0, width: 0, height: 0, enabled: false };
  const applyViewport = () => {
    if (!viewport.set) return;
    const { x, y, width, height } = viewport;
    renderState.setViewport(x, getRenderPassTargetHeight() - y - height, width, height);
  };
  const applyScissor = () => {
    if (!scissor.enabled) return;
    const { x, y, width, height } = scissor;
    renderState.setScissor(x, getRenderPassTargetHeight() - y - height, width, height);
  };

  // optimize assuming flats first and last data are same or differences don't matter
  // extension is only available when `FIRST_VERTEX_CONVENTION` is more efficient
  const epv = extensions.provokingVertex;
  epv?.provokingVertex(epv.FIRST_VERTEX_CONVENTION);
  const ecc = extensions.clipControl;
  if (clipZeroToOne && ecc) state.clipControl?.(ecc.LOWER_LEFT, ecc.ZERO_TO_ONE);

  let isContextLost = false;
  const contextLost = new Subject<now.Timestamp>();

  let pixelScale = props.pixelScale || 1;

  const xr = {
    session: undefined as XRSession | undefined,
    layer: undefined as XRWebGLLayer | undefined,
    changed: new Subject<void>(),
    clear: () => {
      xr.layer = undefined;
      xr.session = undefined;
      xr.changed.next();
    },
  };

  // losses and restorations the platform reports for the canvas (none for contexts without one, e.g. headless)
  const canvas = gl.canvas as HTMLCanvasElement | OffscreenCanvas | undefined;
  const hasCanvasEvents = typeof canvas?.addEventListener === 'function';
  /** the pending `recreate`, completed by the restore event */
  let restoring: { resolve: (restored: boolean) => void } | undefined;
  /** lost on purpose (`release`): `recreate` asks for the context back */
  let released: WEBGL_lose_context | undefined;
  const onContextLost = (e: Event) => {
    // allows the platform to restore the context
    e.preventDefault();
    ctx.setContextLost();
  };
  // the restored context of the canvas is used by a new context object (`recreate`), this one stays lost
  const onContextRestored = () => {
    if (!ctx.isContextLost) return;
    const r = restoring;
    restoring = undefined;
    r?.resolve(true);
  };
  if (hasCanvasEvents) {
    canvas!.addEventListener('webglcontextlost', onContextLost, false);
    canvas!.addEventListener('webglcontextrestored', onContextRestored, false);
  }

  const ctx: WebGLContext = {
    gl,
    isWebGL2: isWebGL2(gl),
    uniformBlocks,
    depthRangeConversion,
    get pixelRatio() {
      const dpr = typeof window !== 'undefined' ? window.devicePixelRatio || 1 : 1;
      return dpr * (pixelScale || 1);
    },

    capabilities,
    extensions,
    state,
    stats,
    resources,
    timer,

    get shaderPrecisionFormats() {
      return shaderPrecisionFormats;
    },

    namedComputeRenderables: Object.create(null),
    namedFramebuffers: Object.create(null),
    namedTextures: Object.create(null),

    isProgramSupported: () => true,
    createGraphicsRenderItem: (drawMode, shaderCode, schema, values, materialId, transparency) => {
      return createGraphicsRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, transparency);
    },
    createProgramItem: (_drawMode, shaderCode, schema, defineValues, _linearTextures, transparency) => {
      return createWebGLProgramItem(ctx, shaderCode, schema, defineValues, transparency);
    },
    createComputeRenderItem: (drawMode, shaderCode, schema, values, materialId) => {
      return createComputeRenderItem(ctx, drawMode, shaderCode, schema, values, materialId);
    },
    createGlobalBindings: (schema, uniforms, textures) => {
      return createGlobalBindings(ctx, schema, uniforms, textures);
    },

    get isContextLost() {
      return isContextLost || gl.isContextLost();
    },
    contextLost,
    setContextLost: () => {
      isContextLost = true;
      timer.clear();
      contextLost.next(now());
    },
    recreate: async () => {
      if (!ctx.isContextLost || !hasCanvasEvents) return undefined;
      restoring?.resolve(false);
      const restored = await new Promise<boolean>((resolve) => {
        restoring = { resolve };
        const ext = released;
        released = undefined;
        if (!ext) return;
        try {
          ext.restoreContext();
        } catch {
          // refused, e.g. before the loss was reported
          restoring = undefined;
          resolve(false);
        }
      });
      // the restored context of the canvas, with the options of this one
      return restored ? createContext(gl, props) : undefined;
    },
    simulateContextLoss: () => {
      const ext = gl.getExtension('WEBGL_lose_context');
      if (!ext || ctx.isContextLost) return false;
      ext.loseContext();
      setTimeout(() => {
        if (ctx.isContextLost) ext.restoreContext();
      }, 1000);
      return true;
    },
    release: () => {
      const ext = gl.getExtension('WEBGL_lose_context');
      if (!ext || ctx.isContextLost) return false;
      released = ext;
      ext.loseContext();
      return true;
    },

    xr: {
      get session() {
        return xr.session;
      },
      changed: xr.changed,
      set: async (session: XRSession | undefined, options?: { resolutionScale?: number }) => {
        if (xr.session === session) return;

        await xr.session?.end();
        if (session === undefined) return;

        try {
          await gl.makeXRCompatible();
          xr.session = session;
          xr.layer = new XRWebGLLayer(xr.session, gl, {
            antialias: true,
            alpha: true,
            depth: true,
            stencil: true,
            framebufferScaleFactor: pixelScale * (options?.resolutionScale ?? 1),
          });
          await xr.session.updateRenderState({ baseLayer: xr.layer });

          xr.session.addEventListener('end', xr.clear);
          xr.changed.next();
        } catch (err) {
          if (session) {
            await session.end();
          } else {
            xr.layer = undefined;
            xr.session = undefined;
          }
          throw err;
        }
      },
      end: async () => {
        return xr.session?.end();
      },
      isPresenting: () => !!xr.session && !!xr.layer && xr.session.renderState.baseLayer === xr.layer,
      getView: (view: XRView, out: GraphicsXRView) => {
        const { layer } = xr;
        const vp = layer?.getViewport(view);
        if (!layer || !vp) return false;
        // viewports of XR WebGL layers have GL's bottom-left origin (for `gl.viewport`)
        Viewport.set(out.viewport, vp.x, layer.framebufferHeight - vp.y - vp.height, vp.width, vp.height);
        // WebGL sessions have projections with GL's clip depth
        Mat4.fromArray(out.projection, view.projectionMatrix, 0);
        return true;
      },
    },

    setPixelScale: (value: number) => {
      pixelScale = value;
    },

    createRenderTarget: (
      width: number,
      height: number,
      depthStencil?: 'none' | 'depth' | 'depth-stencil',
      type?: 'uint8' | 'float32' | 'fp16',
      filter?: TextureFilter,
      format?: 'rgba' | 'r',
    ) => {
      return createRenderTarget(gl, resources, width, height, depthStencil, type, filter, format);
    },
    createDrawTarget: () => {
      return {
        id: -1,
        texture: createNullTexture(),
        framebuffer: createNullFramebuffer(),
        depthRenderbuffer: null,

        getByteCount: () => 0,

        getWidth: () => getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode).width,
        getHeight: () => getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode).height,
        bind: () => {
          bindDrawingBuffer(gl, xr.layer);
        },
        setSize: () => {},
        destroy: () => {},
      };
    },
    bindDrawingBuffer: () => bindDrawingBuffer(gl, xr.layer),
    getDrawingBufferSize: () => getDrawingBufferSize(gl, xr.layer, xr.session?.interactionMode),
    readPixels: (
      x: number,
      y: number,
      width: number,
      height: number,
      buffer: Uint8Array | Float32Array | Int32Array,
    ) => {
      endRenderPass(state);
      readPixels(gl, x, y, width, height, buffer);
    },
    waitForGpuCommandsComplete: () => waitForGpuCommandsComplete(gl),
    readRenderTarget: async (target: RenderTarget, x: number, y: number, width: number, height: number) => {
      endRenderPass(state);
      target.bind();
      const array = new Uint8Array(width * height * 4);
      // the rectangle has a top-left origin, GL's a bottom-left one
      readPixels(gl, x, target.getHeight() - y - height, width, height, array);
      // GL reads rows bottom to top
      PixelData.flipY(PixelData.create(array, width, height));
      return array;
    },
    readTextureAsync: async (texture: Texture, x: number, y: number, width: number, height: number) => {
      endRenderPass(state);
      readFramebuffer ??= resources.framebuffer();
      texture.attachFramebuffer(readFramebuffer, 'color0');
      readFramebuffer.bind();
      const array = new Uint8Array(width * height * 4);
      // the rectangle has a top-left origin, GL's a bottom-left one
      readPixels(gl, x, texture.getHeight() - y - height, width, height, array);
      // GL reads rows bottom to top
      PixelData.flipY(PixelData.create(array, width, height));
      return array;
    },
    createTextureReader: (options) => createWebGLTextureReader(ctx, options),
    readTextureData: async (texture: Texture) => {
      endRenderPass(state);
      const { format, type } = asGLTexture(texture);
      const width = texture.getWidth();
      const height = texture.getHeight();
      if (type === gl.INT || (isWebGL2(gl) && type === gl.UNSIGNED_INT)) {
        throw new Error('reading integer textures is not supported');
      }
      // float and half float
      const isFloat = type !== gl.UNSIGNED_BYTE;
      let source = texture;
      if (format !== gl.RGBA) {
        // single-channel formats aren't readable everywhere: copied into rgba, values in `r`
        if (isFloat) throw new Error('reading single-channel float textures is not supported');
        source = copyToRgbaTexture(ctx, texture, width, height);
      }
      readFramebuffer ??= resources.framebuffer();
      source.attachFramebuffer(readFramebuffer, 'color0');
      readFramebuffer.bind();
      const array = isFloat ? new Float32Array(width * height * 4) : new Uint8Array(width * height * 4);
      readPixels(gl, 0, 0, width, height, array);
      return { array, width, height };
    },
    waitForGpuCommandsCompleteSync: () => waitForGpuCommandsCompleteSync(gl),
    getFenceSync: () => {
      return isWebGL2(gl) ? gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0) : null;
    },
    checkSyncStatus: (sync: WebGLSync) => {
      if (!isWebGL2(gl)) return true;

      if (gl.getSyncParameter(sync, gl.SYNC_STATUS) === gl.SIGNALED) {
        gl.deleteSync(sync);
        return true;
      } else {
        return false;
      }
    },
    deleteSync: (sync: WebGLSync) => {
      if (isWebGL2(gl)) gl.deleteSync(sync);
    },
    setRenderState: renderState.setRenderState,
    setBlendConstant: (red: number, green: number, blue: number, alpha: number) => {
      state.blendColor(red, green, blue, alpha);
    },
    setViewport: (x: number, y: number, width: number, height: number) => {
      Object.assign(viewport, { x, y, width, height, set: true });
      if (state.inRenderPass) applyViewport();
    },
    setScissor: (x: number, y: number, width: number, height: number) => {
      Object.assign(scissor, { x, y, width, height, enabled: true });
      if (state.inRenderPass) applyScissor();
    },
    disableScissor: () => {
      scissor.enabled = false;
      renderState.disableScissor();
    },
    describeTarget: () => WebGLTarget,
    beginRenderPass: (target: RenderPassTarget, options?: RenderPassOptions) => {
      timer.passBegin(target, options?.label);
      if (target === 'drawing-buffer') bindDrawingBuffer(gl, xr.layer);
      else target.bind();
      passTarget = target;
      applyViewport();
      applyScissor();
      state.inRenderPass = true;
      // draws of the pass set their render state, which the flag masks
      renderState.setDepthReadOnly(!!options?.depthReadOnly);
      if (options?.clear) renderState.clear(options.clear);
    },
    endRenderPass: () => endRenderPass(state),
    hasRenderPassStencil: () => {
      // offscreen targets used for rendering are created with a stencil, the drawing buffer may lack one
      return drawingBufferHasStencil || gl.getParameter(gl.FRAMEBUFFER_BINDING) !== null;
    },
    clear: renderState.clear,
    flush: () => {
      endRenderPass(state);
      gl.flush();
    },

    checkError: (message?: string) => {
      checkError(gl, message);
    },
    checkFramebufferStatus: (message?: string) => {
      checkFramebufferStatus(gl, message);
    },

    destroy: (options?: GraphicsDestroyOptions) => {
      if (hasCanvasEvents) {
        canvas!.removeEventListener('webglcontextlost', onContextLost, false);
        canvas!.removeEventListener('webglcontextrestored', onContextRestored, false);
      }
      restoring?.resolve(false);
      restoring = undefined;
      // the objects of a lost context are gone with it, and the restored context of the canvas is used by another
      if (!options?.lost) {
        resources.destroy();
        unbindResources(gl);
      }

      xr.session?.removeEventListener('end', xr.clear);
      xr.session?.end();

      contextLost.complete();
      xr.changed.complete();

      // to aid GC
      if (!options?.doNotForceWebGLContextLoss && !options?.lost) {
        gl.getExtension('WEBGL_lose_context')?.loseContext();
        gl.getExtension('STACKGL_destroy_context')?.destroy();
      }
    },
  };
  logGraphicsContextInfo(ctx, {
    parameters,
    extensions: Object.keys(extensions)
      .filter((k) => !!extensions[k as keyof typeof extensions])
      .sort(),
    shaderPrecisionFormats,
  });
  registerRenderPassTargetHeight(ctx, getRenderPassTargetHeight);
  return ctx;
}
