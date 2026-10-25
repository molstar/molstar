/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Jesse Liang <jesse.liang@rcsb.org>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Ke Ma <mark.ma@rcsb.org>
 * @author Adam Midlik <midlik@gmail.com>
 */

import fs from 'fs';
import path from 'path';
import type { PNG } from 'pngjs'; // Only import type here, the actual import must be provided by the caller

import {
  Canvas3D,
  Canvas3DContext,
  type Canvas3DProps,
  DefaultCanvas3DParams,
} from '@molstar/graphics/canvas3d/canvas3d';
import type { ImagePass, ImageProps } from '@molstar/graphics/canvas3d/passes/image';
import { Passes } from '@molstar/graphics/canvas3d/passes/passes';
import { PostprocessingParams, type PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { createContext } from '@molstar/graphics/gl/webgl/context';
import {
  createWebGPUContext,
  requestWebGPUDevice,
  type WebGPUContextOptions,
} from '@molstar/graphics/gl/webgpu/context';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ColorNames } from '@molstar/core/util/color/names';
import { PixelData } from '@molstar/core/util/image';
import { InputObserver } from '@molstar/core/util/input/input-observer';
import { now } from '@molstar/core/util/now';
import { BehaviorSubject, Subject } from 'rxjs';
import { ParamDefinition } from '@molstar/core/util/param-definition';
import type { RuntimeContext } from '@molstar/core/task';

/** Encoded JPEG data returned by an injected codec. */
export interface JpegBufferRet {
  width: number;
  height: number;
  data: Buffer;
}

export interface ExternalModules {
  /** for the WebGL backend (the default) */
  gl?: typeof import('gl');
  /** for the WebGPU backend, e.g. Dawn's Node.js bindings (the `webgpu` package) */
  webgpu?: { create(options: string[]): GPU };
  'jpeg-js'?: { encode(image: RawImageData, quality?: number): JpegBufferRet };
  pngjs?: typeof import('pngjs');
}

export type HeadlessScreenshotHelperOptions = {
  /** `webgpu` needs `ExternalModules.webgpu` and `HeadlessScreenshotHelper.create` (devices are requested asynchronously) */
  backend?: 'webgl' | 'webgpu';
  webgl?: WebGLContextAttributes;
  /** of the WebGPU context, e.g. a lower data binding size (`maxDataBindingSize`) */
  webgpu?: Pick<WebGPUContextOptions, 'maxDataBindingSize'>;
  canvas?: Partial<Canvas3DProps>;
  imagePass?: Partial<ImageProps>;
};

export type RawImageData = {
  data: Uint8ClampedArray;
  width: number;
  height: number;
};

/** The contexts of the headless `Canvas3D`s, see `HeadlessScreenshotHelper.canvas3dContext` */
const HeadlessContexts = new WeakMap<Canvas3D, Canvas3DContext>();

/** `Canvas3D` without canvas and input, rendering with `graphics` */
function createHeadlessCanvas3D(graphics: GraphicsContext, options?: HeadlessScreenshotHelperOptions) {
  const input = InputObserver.create();
  const attribs = { ...Canvas3DContext.DefaultAttribs };
  const props = { ...Canvas3DContext.DefaultProps };
  const assetManager = new DefaultAssetManager();
  let current = graphics;
  let passes = new Passes(current, assetManager, props);
  const pixelScale = 1;
  const syncPixelScale = () => {};
  const changed = new BehaviorSubject<undefined>(undefined);
  /**
   * As the canvas context's: the passes take the pick scale and transparency (`Canvas3D` follows `changed`); the
   * resolution props are kept, a headless canvas draws at its size
   */
  const setProps = (next?: Partial<Canvas3DContext.Props>) => {
    if (!next) return;
    let hasChanged = false;
    if (next.resolutionMode !== undefined) props.resolutionMode = next.resolutionMode;
    if (next.pixelScale !== undefined) props.pixelScale = next.pixelScale;
    if (next.pickScale !== undefined && next.pickScale !== props.pickScale) {
      props.pickScale = next.pickScale;
      passes.setPickScale(next.pickScale);
      hasChanged = true;
    }
    if (next.transparency !== undefined && next.transparency !== props.transparency) {
      props.transparency = next.transparency;
      passes.setTransparency(next.transparency);
      hasChanged = true;
    }
    if (hasChanged) changed.next(undefined);
  };
  // a lost context is replaced by a new one where the backend can create it without a canvas (WebGPU devices), and
  // Canvas3D is built again on it (`contextChanged` with `restore`); a released one (`release`) only on `restore`
  const contextLost = new Subject<now.Timestamp>();
  const contextRestored = new Subject<now.Timestamp>();
  const contextChanged = new Subject<Canvas3DContext.ContextChange>();
  let released = false;
  const replace = (lost: GraphicsContext) =>
    lost.recreate().then((next) => {
      if (!next || lost !== current) {
        next?.destroy({ lost: true });
        return false;
      }
      contextLostSub.unsubscribe();
      current = next;
      passes = new Passes(current, assetManager, props);
      contextLostSub = watch(current);
      contextChanged.next({ previous: lost, reason: 'restore' });
      lost.destroy({ lost: true });
      contextRestored.next(now());
      return true;
    });
  const watch = (lost: GraphicsContext) =>
    lost.contextLost.subscribe((t) => {
      contextLost.next(t);
      if (!released) replace(lost);
    });
  let contextLostSub = watch(current);
  const release = () => {
    if (released || current.isContextLost) return false;
    released = true;
    if (current.release()) return true;
    released = false;
    return false;
  };
  const restore = async () => {
    if (!released || !current.isContextLost) return false;
    released = false;
    const restored = await replace(current);
    if (!restored) released = true;
    return restored;
  };
  const dispose = () => {
    contextLostSub.unsubscribe();
    input.dispose();
    current.destroy();
    contextLost.complete();
    contextRestored.complete();
    contextChanged.complete();
    changed.complete();
  };
  const context: Canvas3DContext = {
    get graphics() {
      return current;
    },
    input,
    get passes() {
      return passes;
    },
    attribs,
    props,
    contextLost,
    contextRestored,
    contextChanged,
    assetManager,
    changed,
    pixelScale,
    syncPixelScale,
    setProps,
    release,
    restore,
    get isReleased() {
      return released;
    },
    dispose,
  };
  const canvas3d = Canvas3D.create(context, options?.canvas ?? defaultCanvas3DParams());
  HeadlessContexts.set(canvas3d, context);
  return canvas3d;
}

/** To render Canvas3D when running in Node.js (without DOM) */
export class HeadlessScreenshotHelper {
  readonly canvas3d: Canvas3D;
  readonly imagePass: ImagePass;

  /** The context of `canvas3d` when the helper made it (losses, `release` and `restore`) */
  get canvas3dContext(): Canvas3DContext | undefined {
    return HeadlessContexts.get(this.canvas3d);
  }

  constructor(
    readonly externalModules: ExternalModules,
    readonly canvasSize: { width: number; height: number },
    canvas3d?: Canvas3D,
    options?: HeadlessScreenshotHelperOptions,
  ) {
    if (canvas3d) {
      this.canvas3d = canvas3d;
    } else {
      if (options?.backend === 'webgpu') {
        throw new Error(`the 'webgpu' backend needs 'HeadlessScreenshotHelper.create'`);
      }
      if (!this.externalModules.gl) throw new Error(`external module 'gl' is required for the 'webgl' backend`);
      const glContext = this.externalModules.gl(
        this.canvasSize.width,
        this.canvasSize.height,
        options?.webgl ?? defaultWebGLAttributes(),
      );
      this.canvas3d = createHeadlessCanvas3D(createContext(glContext), options);
    }

    this.imagePass = this.canvas3d.getImagePass(options?.imagePass ?? defaultImagePassParams());
    this.imagePass.setSize(this.canvasSize.width, this.canvasSize.height);
  }

  /**
   * Creates the helper with a new graphics context of `options.backend` (default `webgl`); needed for `webgpu`, which
   * requests its device asynchronously
   */
  static async create(
    externalModules: ExternalModules,
    canvasSize: { width: number; height: number },
    options?: HeadlessScreenshotHelperOptions,
  ) {
    if (options?.backend !== 'webgpu')
      return new HeadlessScreenshotHelper(externalModules, canvasSize, undefined, options);
    if (!externalModules.webgpu) throw new Error(`external module 'webgpu' is required for the 'webgpu' backend`);
    const gpu = externalModules.webgpu.create([]);
    const device = await requestWebGPUDevice('default', gpu);
    if (!device) throw new Error('no WebGPU adapter available');
    const graphics = createWebGPUContext(device, {
      width: canvasSize.width,
      height: canvasSize.height,
      // also keeps `gpu` alive with the context, collecting it while its device is in use crashes Dawn's Node.js bindings
      requestDevice: () => requestWebGPUDevice('default', gpu),
      ownsDevice: true,
      ...options.webgpu,
    });
    return new HeadlessScreenshotHelper(
      externalModules,
      canvasSize,
      createHeadlessCanvas3D(graphics, options),
      options,
    );
  }

  private async getImageData(runtime: RuntimeContext, width: number, height: number): Promise<RawImageData> {
    this.imagePass.setSize(width, height);
    await this.imagePass.render(runtime);

    // rows top to bottom
    const array = await this.canvas3d.graphics.readRenderTarget(this.imagePass.colorTarget, 0, 0, width, height);
    PixelData.divideByAlpha(PixelData.create(array, width, height));
    // ImageData is not defined in Node.js
    return { data: new Uint8ClampedArray(array), width, height };
  }

  async getImageRaw(
    runtime: RuntimeContext,
    imageSize?: { width: number; height: number },
    postprocessing?: Partial<PostprocessingProps>,
  ): Promise<RawImageData> {
    const width = imageSize?.width ?? this.canvasSize.width;
    const height = imageSize?.height ?? this.canvasSize.height;
    this.canvas3d.commit(true);
    this.imagePass.setProps({
      postprocessing: ParamDefinition.merge(PostprocessingParams, this.canvas3d.props.postprocessing, postprocessing),
    });
    return this.getImageData(runtime, width, height);
  }

  async getImagePng(
    runtime: RuntimeContext,
    imageSize?: { width: number; height: number },
    postprocessing?: Partial<PostprocessingProps>,
  ): Promise<PNG> {
    const imageData = await this.getImageRaw(runtime, imageSize, postprocessing);
    if (!this.externalModules.pngjs) {
      throw new Error(
        "External module 'pngjs' was not provided. If you want to use getImagePng, you must import 'pngjs' and provide it to the HeadlessPluginContext/HeadlessScreenshotHelper constructor.",
      );
    }
    const generatedPng = new this.externalModules.pngjs.PNG({ width: imageData.width, height: imageData.height });
    generatedPng.data = Buffer.from(imageData.data.buffer);
    return generatedPng;
  }

  async getImageJpeg(
    runtime: RuntimeContext,
    imageSize?: { width: number; height: number },
    postprocessing?: Partial<PostprocessingProps>,
    jpegQuality: number = 90,
  ): Promise<JpegBufferRet> {
    const imageData = await this.getImageRaw(runtime, imageSize, postprocessing);
    if (!this.externalModules['jpeg-js']) {
      throw new Error(
        "External module 'jpeg-js' was not provided. If you want to use getImageJpeg, you must import 'jpeg-js' and provide it to the HeadlessPluginContext/HeadlessScreenshotHelper constructor.",
      );
    }
    const generatedJpeg = this.externalModules['jpeg-js'].encode(imageData, jpegQuality);
    return generatedJpeg;
  }

  async saveImage(
    runtime: RuntimeContext,
    outPath: string,
    imageSize?: { width: number; height: number },
    postprocessing?: Partial<PostprocessingProps>,
    format?: 'png' | 'jpeg',
    jpegQuality = 90,
  ) {
    if (!format) {
      const extension = path.extname(outPath).toLowerCase();
      if (extension === '.png') format = 'png';
      else if (extension === '.jpg' || extension === '.jpeg') format = 'jpeg';
      else
        throw new Error(
          `Cannot guess image format from file path '${outPath}'. Specify format explicitly or use path with one of these extensions: .png, .jpg, .jpeg`,
        );
    }
    if (format === 'png') {
      const generatedPng = await this.getImagePng(runtime, imageSize, postprocessing);
      await writePngFile(generatedPng, outPath);
    } else if (format === 'jpeg') {
      const generatedJpeg = await this.getImageJpeg(runtime, imageSize, postprocessing, jpegQuality);
      await writeJpegFile(generatedJpeg, outPath);
    } else {
      throw new Error(`Invalid format: ${format}`);
    }
  }
}

async function writePngFile(png: PNG, outPath: string) {
  await new Promise<void>((resolve) => {
    png.pack().pipe(fs.createWriteStream(outPath)).on('finish', resolve);
  });
}
async function writeJpegFile(jpeg: JpegBufferRet, outPath: string) {
  await new Promise<void>((resolve) => {
    fs.writeFile(outPath, jpeg.data, () => resolve());
  });
}

export function defaultCanvas3DParams(): Partial<Canvas3DProps> {
  return {
    camera: {
      mode: 'orthographic',
      helper: {
        axes: { name: 'off', params: {} },
      },
      stereo: {
        name: 'off',
        params: {},
      },
      fov: 90,
      manualReset: false,
    },
    cameraResetDurationMs: 0,
    cameraFog: {
      name: 'on',
      params: {
        intensity: 50,
      },
    },
    renderer: {
      ...DefaultCanvas3DParams.renderer,
      backgroundColor: ColorNames.white,
    },
    postprocessing: {
      ...DefaultCanvas3DParams.postprocessing,
      occlusion: {
        name: 'off',
        params: {},
      },
      outline: {
        name: 'off',
        params: {},
      },
      antialiasing: {
        name: 'fxaa',
        params: {
          edgeThresholdMin: 0.0312,
          edgeThresholdMax: 0.063,
          iterations: 12,
          subpixelQuality: 0.3,
        },
      },
      background: { variant: { name: 'off', params: {} } },
      shadow: { name: 'off', params: {} },
    },
  };
}

export function defaultWebGLAttributes(): WebGLContextAttributes {
  return {
    antialias: true,
    preserveDrawingBuffer: true,
    alpha: true, // the renderer requires an alpha channel
    depth: true, // the renderer requires a depth buffer
    stencil: true, // the renderer requires a stencil buffer
    premultipliedAlpha: true, // the renderer outputs PMA
  };
}

export function defaultImagePassParams(): Partial<ImageProps> {
  return {
    cameraHelper: {
      axes: { name: 'off', params: {} },
    },
    multiSample: {
      ...DefaultCanvas3DParams.multiSample,
      mode: 'on',
      sampleLevel: 4,
      reuseOcclusion: false,
    },
  };
}

export const STYLIZED_POSTPROCESSING: Partial<PostprocessingProps> = {
  occlusion: {
    name: 'on' as const,
    params: {
      samples: 32,
      multiScale: { name: 'off', params: {} },
      radius: 5,
      bias: 0.8,
      blurKernelSize: 15,
      blurDepthBias: 0.5,
      resolutionScale: 1,
      color: ColorNames.black,
      transparentThreshold: 0.4,
    },
  },
  outline: {
    name: 'on' as const,
    params: {
      scale: 1,
      threshold: 0.95,
      color: ColorNames.black,
      includeTransparent: true,
    },
  },
};
