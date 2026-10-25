/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3D, Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { ImagePass } from '@molstar/graphics/canvas3d/passes/image';
import { type AssetManager, DefaultAssetManager } from '@molstar/core/util/assets';
import { isWebGPUContext, type WebGPUFeatures, type WebGPULimits } from '@molstar/graphics/gl/webgpu/context';
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';
import { createWebGPUBackend } from '@molstar/graphics/gl/webgpu/backend';
import { createWebGLBackend } from '@molstar/graphics/gl/webgl/backend';
import { type Capabilities, type CapabilityReport, getCapabilities, getWebGLCapabilityReport } from './capabilities.js';

/** Unported features throw, so that their cases are skipped instead of rendering without them */

export interface BackendProfile {
  readonly id: string;
  readonly kind: 'webgl' | 'webgpu';
  readonly description: string;
  /** Resolves to undefined if the profile cannot be created in this browser */
  create(size: { width: number; height: number }): Promise<BackendInstance | undefined>;
  /**
   * A context of the profile alone on `canvas` (its backend, options and masked extensions), nothing drawn or compiled
   * with it yet (`create` builds a `Canvas3D` and an image pass on it); undefined where the profile can't be created
   */
  createContext(canvas: HTMLCanvasElement, assetManager: AssetManager): Promise<Canvas3DContext | undefined>;
}

/** Renders through `Canvas3D` and `ImagePass` */
export interface BackendInstance {
  readonly kind: 'webgl' | 'webgpu';
  readonly profile: BackendProfile;
  readonly canvas: HTMLCanvasElement;
  readonly context: Canvas3DContext;
  readonly canvas3d: Canvas3D;
  readonly imagePass: ImagePass;
  readonly capabilities: Capabilities;
  readonly report: CapabilityReport;
  /** Returns and clears pending backend errors (e.g. GL errors) */
  takeErrors(): string[];
  dispose(): void;
}

/** Offscreen canvas element in the document (`InputObserver` observes the parent element for resizes) */
function createHostedCanvas(width: number, height: number) {
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;
  const host = document.createElement('div');
  Object.assign(host.style, { position: 'absolute', left: '-10000px', top: '0', overflow: 'hidden' });
  host.appendChild(canvas);
  document.body.appendChild(host);
  return { canvas, host };
}

const ContextAttribs: WebGLContextAttributes = {
  powerPreference: 'high-performance',
  failIfMajorPerformanceCaveat: false,
  // rendering happens offscreen via ImagePass; default framebuffer AA is irrelevant
  antialias: false,
  preserveDrawingBuffer: true,
  alpha: true,
  depth: true,
  premultipliedAlpha: true,
};

/** The context attributes of every profile (`BackendProfile.createContext`) */
const ContextOptions = {
  powerPreference: ContextAttribs.powerPreference,
  failIfMajorPerformanceCaveat: ContextAttribs.failIfMajorPerformanceCaveat,
  antialias: ContextAttribs.antialias,
  preserveDrawingBuffer: ContextAttribs.preserveDrawingBuffer,
};

/** WebGL2 (OpenGL ES 3.0) minimum limits, except `UNIFORM_BUFFER_OFFSET_ALIGNMENT` which is a maximum */
const WebGL2MinimumLimits: { [k: string]: number } = {
  MAX_TEXTURE_SIZE: 2048,
  MAX_3D_TEXTURE_SIZE: 256,
  MAX_CUBE_MAP_TEXTURE_SIZE: 2048,
  MAX_ARRAY_TEXTURE_LAYERS: 256,
  MAX_RENDERBUFFER_SIZE: 2048,
  MAX_DRAW_BUFFERS: 4,
  MAX_COLOR_ATTACHMENTS: 4,
  MAX_SAMPLES: 4,
  MAX_VERTEX_ATTRIBS: 16,
  MAX_VERTEX_UNIFORM_VECTORS: 256,
  MAX_FRAGMENT_UNIFORM_VECTORS: 224,
  MAX_VARYING_VECTORS: 15,
  MAX_TEXTURE_IMAGE_UNITS: 16,
  MAX_VERTEX_TEXTURE_IMAGE_UNITS: 16,
  MAX_COMBINED_TEXTURE_IMAGE_UNITS: 32,
  MAX_VERTEX_UNIFORM_BLOCKS: 12,
  MAX_FRAGMENT_UNIFORM_BLOCKS: 12,
  MAX_COMBINED_UNIFORM_BLOCKS: 24,
  MAX_UNIFORM_BUFFER_BINDINGS: 24,
  MAX_UNIFORM_BLOCK_SIZE: 16384,
  UNIFORM_BUFFER_OFFSET_ALIGNMENT: 256,
};

type ContextMask = {
  /** extension names to hide, '*' hides all */
  readonly hiddenExtensions?: readonly string[] | '*';
  /** reported limits, by GL constant name */
  readonly limits?: { readonly [k: string]: number };
};

/**
 * Creates the context before `Canvas3DContext.create` does so extensions and limits can be masked;
 * later `getContext` calls with the same type return this patched instance.
 */
function getMaskedContext(canvas: HTMLCanvasElement, type: 'webgl' | 'webgl2', mask: ContextMask) {
  const gl = canvas.getContext(type, ContextAttribs) as WebGLRenderingContext | WebGL2RenderingContext | null;
  if (!gl) return gl;

  const { hiddenExtensions = [], limits } = mask;
  if (hiddenExtensions === '*' || hiddenExtensions.length > 0) {
    const isHidden = (name: string) => hiddenExtensions === '*' || hiddenExtensions.includes(name);
    const getExtension = gl.getExtension.bind(gl);
    const getSupportedExtensions = gl.getSupportedExtensions.bind(gl);
    Object.assign(gl, {
      getExtension: (name: string) => (isHidden(name) ? null : getExtension(name)),
      getSupportedExtensions: () => getSupportedExtensions()?.filter((n) => !isHidden(n)) ?? null,
    });
  }
  if (limits) {
    const byValue = new Map<number, number>();
    for (const [name, limit] of Object.entries(limits)) {
      const pname = (gl as any)[name];
      if (typeof pname === 'number') byValue.set(pname, limit);
    }
    const alignment = (gl as any).UNIFORM_BUFFER_OFFSET_ALIGNMENT;
    const getParameter = gl.getParameter.bind(gl);
    Object.assign(gl, {
      getParameter: (pname: number) => {
        const value = getParameter(pname);
        const limit = byValue.get(pname);
        if (limit === undefined || typeof value !== 'number') return value;
        return pname === alignment ? Math.max(value, limit) : Math.min(value, limit);
      },
    });
  }
  return gl;
}

function WebGLProfile(
  id: string,
  type: 'webgl' | 'webgl2',
  description: string,
  mask: ContextMask = {},
  options: { uniformBlocks?: boolean } = {},
): BackendProfile {
  const backend = createWebGLBackend(type, options);
  const profile: BackendProfile = {
    id,
    kind: 'webgl',
    description,
    createContext: async (canvas, assetManager) => {
      if (!getMaskedContext(canvas, type, mask)) return undefined;
      return Canvas3DContext.create(canvas, assetManager, [backend], ContextOptions, {
        pixelScale: 1,
        resolutionMode: 'native',
      });
    },
    create: async ({ width, height }) => {
      const canvas = document.createElement('canvas');
      canvas.width = width;
      canvas.height = height;
      canvas.style.width = `${width}px`;
      canvas.style.height = `${height}px`;

      // InputObserver observes the parent element for resizes
      const host = document.createElement('div');
      Object.assign(host.style, { position: 'absolute', left: '-10000px', top: '0', overflow: 'hidden' });
      host.appendChild(canvas);
      document.body.appendChild(host);

      const assetManager = new DefaultAssetManager();
      const context = await profile.createContext(canvas, assetManager);
      if (!context) {
        assetManager.dispose();
        host.remove();
        return undefined;
      }
      // no idle-time picking warm-up: cases don't hover, and it slows the runs
      const canvas3d = Canvas3D.create(context, {}, { warmUp: false });
      canvas3d.handleResize();
      const imagePass = canvas3d.getImagePass({});

      const webgl = context.graphics;
      if (!isWebGLContext(webgl)) throw new Error('expected a WebGL context');
      const { gl } = webgl;
      const takeErrors = () => {
        const errors: string[] = [];
        for (let i = 0; i < 16; ++i) {
          const e = gl.getError();
          if (e === gl.NO_ERROR) break;
          errors.push(`GL error 0x${e.toString(16)}`);
        }
        // the current context: a restore replaces it (on the same `gl`)
        if (context.graphics.isContextLost) errors.push('context lost');
        return errors;
      };

      return {
        kind: 'webgl',
        profile,
        canvas,
        context,
        canvas3d,
        imagePass,
        capabilities: getCapabilities(webgl),
        report: getWebGLCapabilityReport(webgl),
        takeErrors,
        dispose: () => {
          canvas3d.dispose();
          context.dispose();
          assetManager.dispose();
          host.remove();
        },
      };
    },
  };
  return profile;
}

function WebGPUProfile(
  id: string,
  description: string,
  limits: WebGPULimits,
  features: WebGPUFeatures = 'optional',
): BackendProfile {
  const backend = createWebGPUBackend({ strict: true, limits, features });
  const profile: BackendProfile = {
    id,
    kind: 'webgpu',
    description,
    createContext: async (canvas, assetManager) => {
      if (!backend.isAvailable()) return undefined;
      try {
        return await Canvas3DContext.create(canvas, assetManager, [backend], ContextOptions, {
          pixelScale: 1,
          resolutionMode: 'native',
        });
      } catch {
        return undefined;
      }
    },
    create: async ({ width, height }) => {
      if (!backend.isAvailable()) return undefined;
      const { canvas, host } = createHostedCanvas(width, height);
      const assetManager = new DefaultAssetManager();
      const context = await profile.createContext(canvas, assetManager);
      if (!context) {
        assetManager.dispose();
        host.remove();
        return undefined;
      }
      const { graphics } = context;
      if (!isWebGPUContext(graphics)) throw new Error('expected a WebGPU context');
      const { device } = graphics;
      // no idle-time picking warm-up: cases don't hover, and it slows the runs
      const canvas3d = Canvas3D.create(context, {}, { warmUp: false });
      canvas3d.handleResize();
      const imagePass = canvas3d.getImagePass({});

      // with the strict backend, cases of passes that are not ported throw and are skipped with the reason
      const capabilities = getCapabilities(graphics);

      const info = (device as GPUDevice & { adapterInfo?: GPUAdapterInfo }).adapterInfo;
      const limits: Record<string, number> = {};
      for (const k in device.limits) {
        const v = (device.limits as unknown as Record<string, unknown>)[k];
        if (typeof v === 'number') limits[k] = v;
      }
      const report: CapabilityReport = {
        backend: 'webgpu',
        version: 'WebGPU',
        renderer: [info?.architecture, info?.description].filter(Boolean).join(' '),
        vendor: info?.vendor ?? '',
        userAgent: navigator.userAgent,
        capabilities,
        limits,
        extensions: Array.from(device.features as unknown as Iterable<string>),
      };

      return {
        kind: 'webgpu',
        profile,
        canvas,
        context,
        canvas3d,
        imagePass,
        capabilities,
        report,
        // the current context: a restore replaces it
        takeErrors: () => {
          const current = context.graphics;
          return isWebGPUContext(current) ? current.takeErrors() : [];
        },
        dispose: () => {
          canvas3d.dispose();
          context.dispose();
          assetManager.dispose();
          device.destroy();
          host.remove();
        },
      };
    },
  };
  return profile;
}

export const BackendProfiles: BackendProfile[] = [
  WebGLProfile('webgl2', 'webgl2', 'WebGL2 with all available extensions (reference)'),
  WebGLProfile('webgl1', 'webgl', 'WebGL1 with all available extensions'),
  // not run by default (`scripts/visual-tests.mjs --profiles webgl2-no-ubo`): the `uniformBlocks` option of the render
  // benchmark, uniforms set one by one as on WebGL1
  WebGLProfile('webgl2-no-ubo', 'webgl2', 'WebGL2 without uniform buffers', {}, { uniformBlocks: false }),
  WebGLProfile(
    'webgl2-no-float',
    'webgl2',
    'WebGL2 without float/half-float render targets and linear float filtering',
    {
      hiddenExtensions: [
        'EXT_color_buffer_float',
        'EXT_color_buffer_half_float',
        'OES_texture_float_linear',
        'EXT_float_blend',
      ],
    },
  ),
  WebGLProfile('webgl2-minimum', 'webgl2', 'WebGL2 without extensions and with minimum limits', {
    hiddenExtensions: '*',
    limits: WebGL2MinimumLimits,
  }),
  WebGLProfile('webgl1-minimal', 'webgl', 'WebGL1 with only the extensions required by Mol*', {
    hiddenExtensions: [
      'OES_texture_float_linear',
      'OES_texture_half_float',
      'OES_texture_half_float_linear',
      'WEBGL_color_buffer_float',
      'EXT_color_buffer_half_float',
      'EXT_float_blend',
      'WEBGL_depth_texture',
      'EXT_frag_depth',
      'WEBGL_draw_buffers',
      'EXT_shader_texture_lod',
      'OES_vertex_array_object',
      'WEBGL_multi_draw',
      'EXT_sRGB',
    ],
  }),
  WebGPUProfile('webgpu', 'WebGPU with the limits of the adapter', 'adapter'),
  WebGPUProfile(
    'webgpu-core',
    'WebGPU with the core limits of the spec and no optional features (the minimum spec)',
    'core',
    'none',
  ),
];

export const ReferenceProfileId = 'webgl2';

export function getBackendProfile(id: string) {
  return BackendProfiles.find((p) => p.id === id);
}
