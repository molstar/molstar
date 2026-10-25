/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PassBuilder } from './pass-builder.js';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { PickType, type Renderer } from '@molstar/graphics/gl/renderer';
import type { Frame } from '@molstar/graphics/gl/renderable';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Renderbuffer } from '@molstar/graphics/gl/renderbuffer';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { TextureReader } from '@molstar/graphics/gl/texture-reader';
import type { Vec3 } from '@molstar/core/math/linear-algebra';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import { unpackRGBAToDepth, unpackRGBToInt } from '@molstar/core/util/number-packing';
import type { ICamera } from '../camera.js';
import { Viewport } from '../camera/util.js';
import type { Helper } from '../helper/helper.js';

export type PickData = { id: PickingId; position: Vec3 };

/** Pick data at a target: `undefined` if nothing is there, `'pending'` while its asynchronous read hasn't returned */
export type PickResult = PickData | 'pending' | undefined;

/** Pick data being read: poll it with `tryGet` (e.g. per frame, as hover and clicks do) or await `wait` */
export type AsyncPickData = {
  tryGet: () => PickResult;
  /** the pick data once read, polling `tryGet` per animation frame */
  wait: () => Promise<PickData | undefined>;
};

export function createAsyncPickData(tryGet: AsyncPickData['tryGet']): AsyncPickData {
  return {
    tryGet,
    wait: async () => {
      for (;;) {
        const data = tryGet();
        if (data !== 'pending') return data;
        await new Promise<void>((resolve) =>
          typeof requestAnimationFrame !== 'undefined'
            ? requestAnimationFrame(() => resolve())
            : setTimeout(resolve, 16),
        );
      }
    },
  };
}

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

  private renderVariant(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    helper: Helper,
    variant: 'pick' | 'depth',
    pickType: number,
    frame: Frame,
  ) {
    if (b.drawing) {
      renderer.clear(false);
      renderer.update(camera, scene, frame);
    }
    const pass = { name: 'pick', variant } as const;
    b.renderer(pass, scene.primitives, () => renderer.renderPick(scene.primitives, camera, variant, pickType));
    // direct volumes are picked where they look solid (their `surfaceOpacity`), a depth they write
    if (this.graphics.capabilities.fragDepth) {
      b.renderer(pass, scene.volumes, () => renderer.renderPick(scene.volumes, camera, variant, pickType));
    }

    if (helper.handle.isEnabled) {
      b.renderer(pass, helper.handle.scene, () => renderer.renderPick(helper.handle.scene, camera, variant, pickType));
    }

    if (helper.camera.isEnabled) {
      if (b.drawing) {
        helper.camera.update(camera);
        renderer.update(helper.camera.camera, helper.camera.scene, frame);
      }
      b.renderer(pass, helper.camera.scene, () =>
        renderer.renderPick(helper.camera.scene, helper.camera.camera, variant, pickType),
      );
    }
  }

  /** Draws the pick buffers (`b`, see `PassBuilder`) */
  render(b: PassBuilder, renderer: Renderer, camera: ICamera, scene: Scene, helper: Helper, frame: Frame) {
    if (this.graphics.capabilities.multipleRenderTargets) {
      b.begin(this.framebuffer, { label: 'pick' });
      this.renderVariant(b, renderer, camera, scene, helper, 'pick', PickType.None, frame);
      // if (this.pickWidth < 256) {
      //     printTextureImage(readTexture(this.webgl, this.groupPickTexture, new Uint8Array(this.pickWidth * this.pickHeight * 4)), { scale: 16, id: 'group', pixelated: true, useCanvas: true, flipY: true });
      // }
    } else {
      b.begin(this.objectPickTarget, { label: 'pick object' });
      this.renderVariant(b, renderer, camera, scene, helper, 'pick', PickType.Object, frame);

      b.begin(this.instancePickTarget, { label: 'pick instance' });
      this.renderVariant(b, renderer, camera, scene, helper, 'pick', PickType.Instance, frame);

      b.begin(this.groupPickTarget, { label: 'pick group' });
      this.renderVariant(b, renderer, camera, scene, helper, 'pick', PickType.Group, frame);
      // printTextureImage(readTexture(this.webgl, this.groupPickTarget.texture, new Uint8Array(this.pickWidth * this.pickHeight * 4)), { scale: 16, id: 'group', pixelated: true, useCanvas: true, flipY: true });

      b.begin(this.depthPickTarget, { label: 'pick depth' });
      this.renderVariant(b, renderer, camera, scene, helper, 'depth', PickType.None, frame);
    }
  }
}

let AsyncPickingWarningShown = false;

/** Whether `buffers` read asynchronously, otherwise picking falls back to synchronous reads */
export function checkAsyncPickingSupport(buffers: PickBuffers): boolean {
  if (buffers.asynchronous) return true;

  if (isDebugMode && !AsyncPickingWarningShown) {
    console.log('No asynchronous texture reads (WebGL2 required). Falling back to synchronous picking.');
    AsyncPickingWarningShown = true;
  }
  return false;
}

export enum AsyncPickStatus {
  Pending,
  Resolved,
  Failed,
}

const ReadStatus = {
  pending: AsyncPickStatus.Pending,
  resolved: AsyncPickStatus.Resolved,
  failed: AsyncPickStatus.Failed,
} as const;

/** CPU-side copies of the pick targets, read synchronously or asynchronously (see `TextureReader`) */
export class PickBuffers {
  readonly supported = true;

  private readonly arrays: [Uint8Array, Uint8Array, Uint8Array, Uint8Array] = [
    new Uint8Array(0),
    new Uint8Array(0),
    new Uint8Array(0),
    new Uint8Array(0),
  ];
  private viewport = Viewport.create(0, 0, 0, 0);
  private ready = false;
  /** frames the started read has been pending (checks less than a frame apart count once), and the last counted */
  private lag = 0;
  private lagTimestamp = 0;

  /** false if `read` only starts an asynchronous read, whose result `check` reports */
  get synchronous() {
    return this.reader.synchronous;
  }

  /** whether `asyncRead` reads asynchronously */
  get asynchronous() {
    return this.reader.asynchronous;
  }

  setViewport(x: number, y: number, width: number, height: number) {
    Viewport.set(this.viewport, x, y, width, height);
    const size = width * height * 4;
    if (this.arrays[0].length !== size) {
      for (let i = 0; i < 4; ++i) this.arrays[i] = new Uint8Array(size);
    }
    this.reset();
  }

  read() {
    if (!this.reader.synchronous) {
      this.asyncRead();
      return;
    }
    if (isTimingMode) this.graphics.timer.group('PickBuffers.read');
    const { x, y, width, height } = this.viewport;
    this.reader.read(this.pickPass.getTextures(), x, y, width, height, this.arrays);
    this.ready = true;
    if (isTimingMode) this.graphics.timer.groupEnd('PickBuffers.read');
  }

  asyncRead() {
    if (!this.reader.asynchronous) return;
    if (isTimingMode) this.graphics.timer.group('PickBuffers.asyncRead');
    const { x, y, width, height } = this.viewport;
    this.reader.start(this.pickPass.getTextures(), x, y, width, height);
    this.ready = false;
    this.lag = 0;
    this.lagTimestamp = now();
    if (isTimingMode) this.graphics.timer.groupEnd('PickBuffers.asyncRead');
  }

  check(): AsyncPickStatus {
    if (this.ready) return AsyncPickStatus.Resolved;
    const status = ReadStatus[this.reader.check(this.arrays)];
    if (status === AsyncPickStatus.Pending) {
      // given up after `maxAsyncReadLag` frames; a click and the hover polling the same read count a frame once
      const t = now();
      if (t - this.lagTimestamp >= 1000 / 60) {
        this.lag += 1;
        this.lagTimestamp = t;
      }
      if (this.lag > this.maxAsyncReadLag) {
        this.reset();
        return AsyncPickStatus.Failed;
      }
    }
    this.ready = status === AsyncPickStatus.Resolved;
    return status;
  }

  /** the rows are top to bottom */
  private getIdx(x: number, y: number): number {
    return (y * this.viewport.width + x) * 4;
  }

  /** window depth of the camera's GL projections (near 0, far 1), for `cameraUnproject` */
  getDepth(x: number, y: number): number {
    if (!this.ready) return -1;

    const idx = this.getIdx(x, y);
    const b = this.arrays[3];
    // rendered depth is reversed (near 1, far 0)
    return 1 - unpackRGBAToDepth(b[idx], b[idx + 1], b[idx + 2], b[idx + 3]);
  }

  private getId(x: number, y: number, buffer: Uint8Array) {
    const idx = this.getIdx(x, y);
    return unpackRGBToInt(buffer[idx], buffer[idx + 1], buffer[idx + 2]);
  }

  getPickingId(x: number, y: number): PickingId | undefined {
    if (!this.ready) return;
    const [object, instance, group] = this.arrays;

    const objectId = this.getId(x, y, object);
    if (objectId === -1 || objectId === PickingId.Null) return;

    const instanceId = this.getId(x, y, instance);
    if (instanceId === -1 || instanceId === PickingId.Null) return;

    const groupId = this.getId(x, y, group);
    if (groupId === -1) return;

    return { objectId, instanceId, groupId };
  }

  reset() {
    this.reader.reset();
    this.ready = false;
    this.lag = 0;
  }

  dispose() {
    this.reader.dispose();
    this.ready = false;
  }

  constructor(
    private graphics: GraphicsContext,
    private pickPass: PickPass,
    private reader: TextureReader,
    public maxAsyncReadLag = DefaultPickOptions.maxAsyncReadLag,
  ) {}
}

export function createPickBuffers(
  graphics: GraphicsContext,
  pickPass: PickPass,
  maxAsyncReadLag?: number,
): PickBuffers {
  return new PickBuffers(graphics, pickPass, graphics.createTextureReader(), maxAsyncReadLag);
}
