/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { PickType, type Renderer } from '@molstar/graphics/gl/renderer';
import type { Frame } from '@molstar/graphics/gl/renderable';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { PixelPackBuffer } from '@molstar/graphics/gl/webgl/buffer';
import { isWebGL2 } from '@molstar/graphics/gl/webgl/compat';
import { isWebGLContext, type WebGLContext } from '@molstar/graphics/gl/webgl/context';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Renderbuffer } from '@molstar/graphics/gl/renderbuffer';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import { unpackRGBAToDepth, unpackRGBToInt } from '@molstar/core/util/number-packing';
import type { ICamera } from '../camera.js';
import { Viewport } from '../camera/util.js';
import type { Helper } from '../helper/helper.js';

export type PickData = { id: PickingId; position: Vec3 };

export type AsyncPickData = {
  tryGet: () => 'pending' | PickData | undefined;
};

export const DefaultPickOptions = {
  pickPadding: 1,
  maxAsyncReadLag: 5,
};
export type PickOptions = typeof DefaultPickOptions;

//

export class PickPass {
  private readonly objectPickTarget: RenderTarget;
  private readonly instancePickTarget: RenderTarget;
  private readonly groupPickTarget: RenderTarget;
  private readonly depthPickTarget: RenderTarget;

  private readonly framebuffer: Framebuffer;

  private readonly objectPickTexture: Texture;
  private readonly instancePickTexture: Texture;
  private readonly groupPickTexture: Texture;
  private readonly depthPickTexture: Texture;

  private readonly objectPickFramebuffer: Framebuffer;
  private readonly instancePickFramebuffer: Framebuffer;
  private readonly groupPickFramebuffer: Framebuffer;
  private readonly depthPickFramebuffer: Framebuffer;

  private readonly depthRenderbuffer: Renderbuffer;

  private pickWidth: number;
  private pickHeight: number;

  constructor(
    private graphics: GraphicsContext,
    private width: number,
    private height: number,
    private pickScale: number,
  ) {
    const pickRatio = pickScale / graphics.pixelRatio;
    this.pickWidth = Math.ceil(width * pickRatio);
    this.pickHeight = Math.ceil(height * pickRatio);

    const { resources } = graphics;

    if (graphics.capabilities.multipleRenderTargets) {
      this.objectPickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.objectPickTexture.define(this.pickWidth, this.pickHeight);

      this.instancePickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.instancePickTexture.define(this.pickWidth, this.pickHeight);

      this.groupPickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.groupPickTexture.define(this.pickWidth, this.pickHeight);

      this.depthPickTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      this.depthPickTexture.define(this.pickWidth, this.pickHeight);

      this.framebuffer = resources.framebuffer();

      this.objectPickFramebuffer = resources.framebuffer();
      this.instancePickFramebuffer = resources.framebuffer();
      this.groupPickFramebuffer = resources.framebuffer();
      this.depthPickFramebuffer = resources.framebuffer();

      this.framebuffer.setDrawBuffers(4);

      this.objectPickTexture.attachFramebuffer(this.framebuffer, 'color0');
      this.instancePickTexture.attachFramebuffer(this.framebuffer, 'color1');
      this.groupPickTexture.attachFramebuffer(this.framebuffer, 'color2');
      this.depthPickTexture.attachFramebuffer(this.framebuffer, 'color3');

      this.depthRenderbuffer = resources.renderbuffer(
        'depth32f-stencil8',
        'depth-stencil',
        this.pickWidth,
        this.pickHeight,
      );

      this.depthRenderbuffer.attachFramebuffer(this.framebuffer);

      this.objectPickTexture.attachFramebuffer(this.objectPickFramebuffer, 'color0');
      this.instancePickTexture.attachFramebuffer(this.instancePickFramebuffer, 'color0');
      this.groupPickTexture.attachFramebuffer(this.groupPickFramebuffer, 'color0');
      this.depthPickTexture.attachFramebuffer(this.depthPickFramebuffer, 'color0');
    } else {
      this.objectPickTarget = graphics.createRenderTarget(this.pickWidth, this.pickHeight, 'depth-stencil');
      this.instancePickTarget = graphics.createRenderTarget(this.pickWidth, this.pickHeight, 'depth-stencil');
      this.groupPickTarget = graphics.createRenderTarget(this.pickWidth, this.pickHeight, 'depth-stencil');
      this.depthPickTarget = graphics.createRenderTarget(this.pickWidth, this.pickHeight, 'depth-stencil');
    }
  }

  getByteCount() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      return (
        this.objectPickTexture.getByteCount() +
        this.instancePickTexture.getByteCount() +
        this.groupPickTexture.getByteCount() +
        this.depthPickTexture.getByteCount() +
        this.depthRenderbuffer.getByteCount()
      );
    } else {
      return (
        this.objectPickTarget.getByteCount() +
        this.instancePickTarget.getByteCount() +
        this.groupPickTarget.getByteCount() +
        this.depthPickTarget.getByteCount()
      );
    }
  }

  dispose() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.framebuffer.destroy();

      this.objectPickTexture.destroy();
      this.instancePickTexture.destroy();
      this.groupPickTexture.destroy();
      this.depthPickTexture.destroy();

      this.objectPickFramebuffer.destroy();
      this.instancePickFramebuffer.destroy();
      this.groupPickFramebuffer.destroy();
      this.depthPickFramebuffer.destroy();

      this.depthRenderbuffer.destroy();
    } else {
      this.objectPickTarget.destroy();
      this.instancePickTarget.destroy();
      this.groupPickTarget.destroy();
      this.depthPickTarget.destroy();
    }
  }

  get pickRatio() {
    return this.pickScale / this.graphics.pixelRatio;
  }

  setPickScale(pickScale: number) {
    this.pickScale = pickScale;
    this.setSize(this.width, this.height);
  }

  bindObject() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.objectPickFramebuffer.bind();
    } else {
      this.objectPickTarget.bind();
    }
  }

  bindInstance() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.instancePickFramebuffer.bind();
    } else {
      this.instancePickTarget.bind();
    }
  }

  bindGroup() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.groupPickFramebuffer.bind();
    } else {
      this.groupPickTarget.bind();
    }
  }

  bindDepth() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.depthPickFramebuffer.bind();
    } else {
      this.depthPickTarget.bind();
    }
  }

  get drawingBufferHeight() {
    return this.height;
  }

  /** The object, instance, group and depth pick textures */
  getTextures(): [Texture, Texture, Texture, Texture] {
    if (this.graphics.capabilities.multipleRenderTargets) {
      return [this.objectPickTexture, this.instancePickTexture, this.groupPickTexture, this.depthPickTexture];
    } else {
      return [
        this.objectPickTarget.texture,
        this.instancePickTarget.texture,
        this.groupPickTarget.texture,
        this.depthPickTarget.texture,
      ];
    }
  }

  setSize(width: number, height: number) {
    this.width = width;
    this.height = height;

    const pickRatio = this.pickScale / this.graphics.pixelRatio;
    const pickWidth = Math.ceil(this.width * pickRatio);
    const pickHeight = Math.ceil(this.height * pickRatio);

    if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
      this.pickWidth = pickWidth;
      this.pickHeight = pickHeight;

      if (this.graphics.capabilities.multipleRenderTargets) {
        this.objectPickTexture.define(this.pickWidth, this.pickHeight);
        this.instancePickTexture.define(this.pickWidth, this.pickHeight);
        this.groupPickTexture.define(this.pickWidth, this.pickHeight);
        this.depthPickTexture.define(this.pickWidth, this.pickHeight);

        this.depthRenderbuffer.setSize(this.pickWidth, this.pickHeight);
      } else {
        this.objectPickTarget.setSize(this.pickWidth, this.pickHeight);
        this.instancePickTarget.setSize(this.pickWidth, this.pickHeight);
        this.groupPickTarget.setSize(this.pickWidth, this.pickHeight);
        this.depthPickTarget.setSize(this.pickWidth, this.pickHeight);
      }
    }
  }

  reset() {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.objectPickTexture.attachFramebuffer(this.framebuffer, 'color0');
      this.instancePickTexture.attachFramebuffer(this.framebuffer, 'color1');
      this.groupPickTexture.attachFramebuffer(this.framebuffer, 'color2');
      this.depthPickTexture.attachFramebuffer(this.framebuffer, 'color3');

      this.depthRenderbuffer.attachFramebuffer(this.framebuffer);

      this.objectPickTexture.attachFramebuffer(this.objectPickFramebuffer, 'color0');
      this.instancePickTexture.attachFramebuffer(this.instancePickFramebuffer, 'color0');
      this.groupPickTexture.attachFramebuffer(this.groupPickFramebuffer, 'color0');
      this.depthPickTexture.attachFramebuffer(this.depthPickFramebuffer, 'color0');
    }
  }

  private renderVariant(
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    helper: Helper,
    variant: 'pick' | 'depth',
    pickType: number,
    frame: Frame,
  ) {
    renderer.clear(false);
    renderer.update(camera, scene, frame);
    renderer.renderPick(scene.primitives, camera, variant, pickType);

    if (helper.handle.isEnabled) {
      renderer.renderPick(helper.handle.scene, camera, variant, pickType);
    }

    if (helper.camera.isEnabled) {
      helper.camera.update(camera);
      renderer.update(helper.camera.camera, helper.camera.scene, frame);
      renderer.renderPick(helper.camera.scene, helper.camera.camera, variant, pickType);
    }
  }

  render(renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper, frame: Frame) {
    if (this.graphics.capabilities.multipleRenderTargets) {
      this.graphics.beginRenderPass(this.framebuffer, { label: 'pick' });
      this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.None, frame);
      // if (this.pickWidth < 256) {
      //     printTextureImage(readTexture(this.webgl, this.groupPickTexture, new Uint8Array(this.pickWidth * this.pickHeight * 4)), { scale: 16, id: 'group', pixelated: true, useCanvas: true, flipY: true });
      // }
    } else {
      this.graphics.beginRenderPass(this.objectPickTarget, { label: 'pick object' });
      this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.Object, frame);

      this.graphics.beginRenderPass(this.instancePickTarget, { label: 'pick instance' });
      this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.Instance, frame);

      this.graphics.beginRenderPass(this.groupPickTarget, { label: 'pick group' });
      this.renderVariant(renderer, camera, scene, helper, 'pick', PickType.Group, frame);
      // printTextureImage(readTexture(this.webgl, this.groupPickTarget.texture, new Uint8Array(this.pickWidth * this.pickHeight * 4)), { scale: 16, id: 'group', pixelated: true, useCanvas: true, flipY: true });

      this.graphics.beginRenderPass(this.depthPickTarget, { label: 'pick depth' });
      this.renderVariant(renderer, camera, scene, helper, 'depth', PickType.None, frame);
    }
  }
}

let AsyncPickingWarningShown = false;

export function checkAsyncPickingSupport(graphics: GraphicsContext): boolean {
  // other backends read the pick targets with `readTextureAsync`, see `createPickBuffers`
  if (!isWebGLContext(graphics)) return true;
  if (graphics.isWebGL2) return true;

  if (isDebugMode && !AsyncPickingWarningShown) {
    console.log('WebGL2 required for async picking. Falling back to synchronous picking.');
    AsyncPickingWarningShown = true;
  }
  return false;
}

export enum AsyncPickStatus {
  Pending,
  Resolved,
  Failed,
}

/** CPU-side copies of the pick targets, read synchronously or asynchronously */
export interface PickBuffers {
  /** false where reading the pick targets isn't implemented */
  readonly supported: boolean;
  /** false if `read` only starts an asynchronous read, whose result `check` reports */
  readonly synchronous: boolean;
  setViewport(x: number, y: number, width: number, height: number): void;
  read(): void;
  asyncRead(): void;
  check(): AsyncPickStatus;
  getDepth(x: number, y: number): number;
  getPickingId(x: number, y: number): PickingId | undefined;
  reset(): void;
  dispose(): void;
}

/** Reads the pick textures with `readTextureAsync`; without synchronous reads, `read` starts an asynchronous one */
class AsyncPickBuffers implements PickBuffers {
  readonly supported = true;
  readonly synchronous = false;

  private object: Uint8Array;
  private instance: Uint8Array;
  private group: Uint8Array;
  private depth: Uint8Array;

  private viewport = Viewport.create(0, 0, 0, 0);
  private status = AsyncPickStatus.Failed;
  /** results of earlier reads are dropped */
  private readId = 0;

  setViewport(x: number, y: number, width: number, height: number) {
    Viewport.set(this.viewport, x, y, width, height);
    this.reset();
  }

  read() {
    this.asyncRead();
  }

  asyncRead() {
    const readId = ++this.readId;
    const { x, y, width, height } = this.viewport;
    this.status = AsyncPickStatus.Pending;
    Promise.all(this.pickPass.getTextures().map((t) => this.graphics.readTextureAsync(t, x, y, width, height))).then(
      ([object, instance, group, depth]) => {
        if (readId !== this.readId) return;
        this.object = object;
        this.instance = instance;
        this.group = group;
        this.depth = depth;
        this.status = AsyncPickStatus.Resolved;
      },
      (e) => {
        if (readId !== this.readId) return;
        if (isDebugMode) console.error(e);
        this.status = AsyncPickStatus.Failed;
      },
    );
  }

  check() {
    return this.status;
  }

  /** the rows are read top to bottom */
  private getIdx(x: number, y: number): number {
    const { width, height } = this.viewport;
    return ((height - 1 - y) * width + x) * 4;
  }

  getDepth(x: number, y: number): number {
    if (this.status !== AsyncPickStatus.Resolved) return -1;

    const idx = this.getIdx(x, y);
    const b = this.depth;
    return unpackRGBAToDepth(b[idx], b[idx + 1], b[idx + 2], b[idx + 3]);
  }

  private getId(x: number, y: number, buffer: Uint8Array) {
    const idx = this.getIdx(x, y);
    return unpackRGBToInt(buffer[idx], buffer[idx + 1], buffer[idx + 2]);
  }

  getPickingId(x: number, y: number): PickingId | undefined {
    if (this.status !== AsyncPickStatus.Resolved) return;

    const objectId = this.getId(x, y, this.object);
    if (objectId === -1 || objectId === PickingId.Null) return;

    const instanceId = this.getId(x, y, this.instance);
    if (instanceId === -1 || instanceId === PickingId.Null) return;

    const groupId = this.getId(x, y, this.group);
    if (groupId === -1) return;

    return { objectId, instanceId, groupId };
  }

  reset() {
    this.readId += 1;
    this.status = AsyncPickStatus.Failed;
  }

  dispose() {
    this.reset();
  }

  constructor(
    private graphics: GraphicsContext,
    private pickPass: PickPass,
  ) {}
}

export function createPickBuffers(
  graphics: GraphicsContext,
  pickPass: PickPass,
  maxAsyncReadLag?: number,
): PickBuffers {
  return isWebGLContext(graphics)
    ? new WebGLPickBuffers(graphics, pickPass, maxAsyncReadLag)
    : new AsyncPickBuffers(graphics, pickPass);
}

/** Reads with `readPixels`, asynchronously through pixel pack buffers and fences on WebGL2 */
class WebGLPickBuffers implements PickBuffers {
  readonly supported = true;
  readonly synchronous = true;
  private object: Uint8Array;
  private instance: Uint8Array;
  private group: Uint8Array;
  private depth: Uint8Array;

  private objectBuffer: PixelPackBuffer;
  private instanceBuffer: PixelPackBuffer;
  private groupBuffer: PixelPackBuffer;
  private depthBuffer: PixelPackBuffer;

  private viewport = Viewport.create(0, 0, 0, 0);

  private setup() {
    const size = this.viewport.width * this.viewport.height * 4;
    if (!this.object || this.object.length !== size) {
      this.object = new Uint8Array(size);
      this.instance = new Uint8Array(size);
      this.group = new Uint8Array(size);
      this.depth = new Uint8Array(size);
    }
  }

  setViewport(x: number, y: number, width: number, height: number) {
    Viewport.set(this.viewport, x, y, width, height);
    this.setup();
  }

  read() {
    if (isTimingMode) this.webgl.timer.group('PickBuffers.read');
    const { x, y, width, height } = this.viewport;

    this.pickPass.bindObject();
    this.webgl.readPixels(x, y, width, height, this.object);

    this.pickPass.bindInstance();
    this.webgl.readPixels(x, y, width, height, this.instance);

    this.pickPass.bindGroup();
    this.webgl.readPixels(x, y, width, height, this.group);

    this.pickPass.bindDepth();
    this.webgl.readPixels(x, y, width, height, this.depth);

    this.ready = true;
    if (isTimingMode) this.webgl.timer.groupEnd('PickBuffers.read');
  }

  private fenceSync: WebGLSync | null = null;
  private fenceTimestamp: number = 0;

  private ready = false;
  private lag = 0;

  asyncRead() {
    const { gl } = this.webgl;
    if (!isWebGL2(gl)) return;

    if (isTimingMode) this.webgl.timer.group('PickBuffers.asyncRead');
    if (this.fenceSync !== null) {
      gl.deleteSync(this.fenceSync);
    }
    const { x, y, width, height } = this.viewport;

    this.pickPass.bindObject();
    this.objectBuffer.read(x, y, width, height);

    this.pickPass.bindInstance();
    this.instanceBuffer.read(x, y, width, height);

    this.pickPass.bindGroup();
    this.groupBuffer.read(x, y, width, height);

    this.pickPass.bindDepth();
    this.depthBuffer.read(x, y, width, height);

    this.fenceTimestamp = now();
    this.fenceSync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
    gl.flush();

    this.ready = false;
    if (isTimingMode) this.webgl.timer.groupEnd('PickBuffers.asyncRead');
  }

  check(): AsyncPickStatus {
    if (this.ready) return AsyncPickStatus.Resolved;
    if (this.fenceSync === null) return AsyncPickStatus.Failed;

    const { gl } = this.webgl;
    if (!isWebGL2(gl)) return AsyncPickStatus.Failed;

    const res = gl.clientWaitSync(this.fenceSync, 0, 0);
    if (res === gl.WAIT_FAILED || this.lag >= this.maxAsyncReadLag) {
      // console.log(`failed to get buffer data after ${this.lag + 1} checks`);
      if (res !== gl.WAIT_FAILED && now() - this.fenceTimestamp < 1000 / 60) {
        this.lag += 1;
        return AsyncPickStatus.Pending;
      }
      gl.deleteSync(this.fenceSync);
      this.fenceSync = null;
      this.lag = 0;
      this.ready = false;
      return AsyncPickStatus.Failed;
    } else if (res === gl.TIMEOUT_EXPIRED) {
      this.lag += 1;
      // console.log(`waiting for buffer data for ${this.lag} checks`);
      return AsyncPickStatus.Pending;
    } else {
      this.objectBuffer.getSubData(this.object);
      this.instanceBuffer.getSubData(this.instance);
      this.groupBuffer.getSubData(this.group);
      this.depthBuffer.getSubData(this.depth);

      // console.log(`got buffer data after ${this.lag + 1} checks`);
      gl.deleteSync(this.fenceSync);
      this.fenceSync = null;
      this.lag = 0;
      this.ready = true;

      return AsyncPickStatus.Resolved;
    }
  }

  private getIdx(x: number, y: number): number {
    return (y * this.viewport.width + x) * 4;
  }

  getDepth(x: number, y: number): number {
    if (!this.ready) return -1;

    const idx = this.getIdx(x, y);
    const b = this.depth;
    return unpackRGBAToDepth(b[idx], b[idx + 1], b[idx + 2], b[idx + 3]);
  }

  private getId(x: number, y: number, buffer: Uint8Array) {
    if (!this.ready) return -1;

    const idx = this.getIdx(x, y);
    return unpackRGBToInt(buffer[idx], buffer[idx + 1], buffer[idx + 2]);
  }

  getObjectId(x: number, y: number) {
    return this.getId(x, y, this.object);
  }

  getInstanceId(x: number, y: number) {
    return this.getId(x, y, this.instance);
  }

  getGroupId(x: number, y: number) {
    return this.getId(x, y, this.group);
  }

  getPickingId(x: number, y: number): PickingId | undefined {
    const objectId = this.getObjectId(x, y);
    // console.log('objectId', objectId);
    if (objectId === -1 || objectId === PickingId.Null) return;

    const instanceId = this.getInstanceId(x, y);
    // console.log('instanceId', instanceId);
    if (instanceId === -1 || instanceId === PickingId.Null) return;

    const groupId = this.getGroupId(x, y);
    // console.log('groupId', groupId);
    if (groupId === -1) return;

    return { objectId, instanceId, groupId };
  }

  reset() {
    this.fenceSync = null;
    this.ready = false;
    this.lag = 0;
    this.fenceTimestamp = 0;
  }

  dispose() {
    const { gl } = this.webgl;
    if (!isWebGL2(gl)) return;

    this.objectBuffer.destroy();
    this.instanceBuffer.destroy();
    this.groupBuffer.destroy();
    this.depthBuffer.destroy();

    if (this.fenceSync !== null) {
      gl.deleteSync(this.fenceSync);
      this.fenceSync = null;
    }
  }

  constructor(
    private webgl: WebGLContext,
    private pickPass: PickPass,
    public maxAsyncReadLag = 5,
  ) {
    if (webgl.isWebGL2) {
      this.objectBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
      this.instanceBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
      this.groupBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
      this.depthBuffer = webgl.resources.pixelPack('rgba', 'ubyte');
    }

    this.setup();
  }
}
