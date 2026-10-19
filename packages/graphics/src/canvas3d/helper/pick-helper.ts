/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderBuilder, type PassBuilder } from '../passes/pass-builder.js';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import { type Frame, createFrame } from '@molstar/graphics/gl/renderable';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { spiral2d } from '@molstar/core/math/misc';
import { isTimingMode } from '@molstar/core/util/debug';
import type { Camera } from '../camera.js';
import { StereoCamera } from '../camera/stereo.js';
import { cameraUnproject, Viewport } from '../camera/util.js';
import type { Helper } from './helper.js';
import {
  type AsyncPickData,
  AsyncPickStatus,
  checkAsyncPickingSupport,
  createAsyncPickData,
  createPickBuffers,
  type PickData,
  type PickResult,
  type PickOptions,
  type PickPass,
} from '../passes/pick.js';

export class PickHelper {
  dirty = true;

  /**
   * The last render skipped draws (pipelines compiling, e.g. right after a restore or load): its read is not reported
   * (`'pending'`), the pick is rendered again when tracked work has settled since
   */
  private incomplete = false;
  private settledCount = 0;
  private renderSettledCount = 0;
  private readonly settledSub = this.graphics.resources.pending.settled.subscribe(() => {
    this.settledCount += 1;
  });

  private pickPadding: number;
  private buffers = createPickBuffers(this.graphics, this.pickPass);
  private viewport = Viewport();

  private pickRatio: number;
  private pickX: number;
  private pickY: number;
  private pickWidth: number;
  private pickHeight: number;
  private halfPickWidth: number;

  private spiral: [number, number][];

  setViewport(x: number, y: number, width: number, height: number) {
    Viewport.set(this.viewport, x, y, width, height);
    this.update();
  }

  setPickPadding(pickPadding: number) {
    if (this.pickPadding !== pickPadding) {
      this.pickPadding = pickPadding;
      this.update();
    }
  }

  private update() {
    const { x, y, width, height } = this.viewport;

    this.pickRatio = this.pickPass.pickRatio;
    this.pickX = Math.ceil(x * this.pickRatio);
    this.pickY = Math.ceil(y * this.pickRatio);

    const pickWidth = Math.floor(width * this.pickRatio);
    const pickHeight = Math.floor(height * this.pickRatio);

    if (pickWidth !== this.pickWidth || pickHeight !== this.pickHeight) {
      this.pickWidth = pickWidth;
      this.pickHeight = pickHeight;
      this.halfPickWidth = Math.floor(this.pickWidth / 2);

      this.buffers.setViewport(this.pickX, this.pickY, this.pickWidth, this.pickHeight);
    }

    this.spiral = spiral2d(Math.ceil(this.pickRatio * this.pickPadding));
    this.dirty = true;
  }

  private render(camera: Camera | StereoCamera, frame?: Frame) {
    this.plan(this.renderBuilder, camera, frame);
  }

  /** What `identify` draws before reading, see `PassBuilder` */
  plan(b: PassBuilder, camera: Camera | StereoCamera, frame?: Frame) {
    const timing = isTimingMode && b.drawing;
    if (timing) this.graphics.timer.group('PickHelper.render', { captureStats: true });
    const { stats } = this.graphics;
    const skippedDraws = stats.skippedDraws;
    const { pickX, pickY, pickWidth, pickHeight, halfPickWidth } = this;
    const { renderer, scene, helper } = this;

    // share one frame across both eyes (and with the caller's frame, e.g. Canvas3D's, when given)
    // so a resting camera can reuse the main draw's cull instead of always recomputing
    const f = frame ?? createFrame();

    renderer.setTransparentBackground(false);
    if (b.drawing) {
      renderer.setDrawingBufferSize(pickWidth, pickHeight);
      renderer.setPixelRatio(this.pickRatio);
    }

    if (StereoCamera.is(camera)) {
      if (b.drawing) renderer.setViewport(pickX, pickY, halfPickWidth, pickHeight);
      this.pickPass.render(b, renderer, camera.left, scene, helper, f);

      // preparing: the right eye draws as the left one
      if (b.drawing) {
        renderer.setViewport(pickX + halfPickWidth, pickY, pickWidth - halfPickWidth, pickHeight);
        this.pickPass.render(b, renderer, camera.right, scene, helper, f);
      }
    } else {
      if (b.drawing) renderer.setViewport(pickX, pickY, pickWidth, pickHeight);
      this.pickPass.render(b, renderer, camera, scene, helper, f);
    }

    if (b.drawing) {
      this.dirty = false;
      this.incomplete = stats.skippedDraws !== skippedDraws;
      this.renderSettledCount = this.settledCount;
    }
    if (timing) this.graphics.timer.groupEnd('PickHelper.render');
  }

  /** Whether `identify` can read pick buffers */
  get supported() {
    return this.buffers.supported;
  }

  private identifyInternal(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
    if (this.graphics.isContextLost) return;

    const { graphics, pickRatio } = this;
    if (graphics.isContextLost) return;

    // window pixels, top-left origin as the pointer coordinates
    x *= graphics.pixelRatio;
    y *= graphics.pixelRatio;

    const { viewport } = this;

    // check if within viewport
    if (x < viewport.x || y < viewport.y || x > viewport.x + viewport.width || y > viewport.y + viewport.height) return;

    const xv = x - viewport.x;
    const yv = y - viewport.y;

    const xp = Math.floor(xv * pickRatio);
    const yp = Math.floor(yv * pickRatio);

    const pickingId = this.buffers.getPickingId(xp, yp);
    if (pickingId === undefined) return;

    const z = this.buffers.getDepth(xp, yp);
    const position = Vec3.create(x, y, z);
    if (StereoCamera.is(camera)) {
      const halfWidth = Math.floor(viewport.width / 2);
      if (x > viewport.x + halfWidth) {
        position[0] = viewport.x + (xv - halfWidth) * 2;
        cameraUnproject(position, position, viewport, camera.right.inverseProjectionView);
      } else {
        position[0] = viewport.x + xv * 2;
        cameraUnproject(position, position, viewport, camera.left.inverseProjectionView);
      }
    } else {
      cameraUnproject(position, position, viewport, camera.inverseProjectionView);
    }

    return { id: pickingId, position };
  }

  private prepare() {
    if (this.pickRatio !== this.pickPass.pickRatio) {
      this.update();
    }
  }

  /** Whether the pick is drawn again: changed since its render, or that render was incomplete and work settled since */
  private get needsRender() {
    return this.dirty || (this.incomplete && this.settledCount !== this.renderSettledCount);
  }

  private getPickData(x: number, y: number, camera: Camera | StereoCamera): PickData | undefined {
    for (const d of this.spiral) {
      const pickData = this.identifyInternal(x + d[0], y + d[1], camera);
      if (pickData) return pickData;
    }
  }

  /**
   * The pick data at (x, y); without synchronous reads `'pending'` until the read started for the current scene and
   * camera has returned (see `asyncIdentify`)
   */
  identify(x: number, y: number, camera: Camera | StereoCamera, frame?: Frame): PickResult {
    if (!this.buffers.supported) return;
    this.prepare();

    if (this.needsRender) {
      if (isTimingMode) this.graphics.timer.group('PickHelper.identify');
      this.graphics.resources.finalizePrograms(['pick'], true);
      this.render(camera, frame);
      this.buffers.read();
      if (isTimingMode) this.graphics.timer.groupEnd('PickHelper.identify');
    }

    if (!this.buffers.synchronous) {
      const status = this.buffers.check();
      if (status === AsyncPickStatus.Pending) return 'pending';
      if (status === AsyncPickStatus.Failed) {
        // read again on the next call
        this.dirty = true;
        return 'pending';
      }
    }
    if (this.incomplete) return 'pending';
    return this.getPickData(x, y, camera);
  }

  asyncIdentify(x: number, y: number, camera: Camera | StereoCamera, frame?: Frame): AsyncPickData | undefined {
    if (!this.buffers.supported) return;
    this.prepare();

    if (this.needsRender) this.asyncRender(camera, frame);

    return createAsyncPickData(() => {
      // nothing once disposed (e.g. the graphics context changed), the scene may hold another context's renderables
      if (this.disposed) return undefined;
      // not when only `dirty` (frames drawn since the call): the read of the call is reported
      if (this.incomplete && this.settledCount !== this.renderSettledCount) this.asyncRender(camera);
      const status = this.buffers.check();
      if (status === AsyncPickStatus.Resolved) {
        return this.incomplete ? 'pending' : this.getPickData(x, y, camera);
      } else if (status === AsyncPickStatus.Failed) {
        // read again (slow under GPU load, e.g. illumination accumulating): a click must not resolve to nothing
        this.asyncRender(camera);
      }
      return 'pending';
    });
  }

  private asyncRender(camera: Camera | StereoCamera, frame?: Frame) {
    if (isTimingMode) this.graphics.timer.group('PickHelper.asyncIdentify');
    this.render(camera, frame);
    this.buffers.asyncRead();
    if (isTimingMode) this.graphics.timer.groupEnd('PickHelper.asyncIdentify');
  }

  /** Whether `identify` reads the pick data synchronously */
  get synchronous() {
    return this.buffers.synchronous;
  }

  private disposed = false;

  dispose() {
    this.disposed = true;
    this.settledSub.unsubscribe();
    this.buffers.dispose();
  }

  private readonly renderBuilder: PassBuilder;

  constructor(
    private graphics: GraphicsContext,
    private renderer: Renderer,
    private scene: Scene,
    private helper: Helper,
    private pickPass: PickPass,
    viewport: Viewport,
    options: PickOptions,
  ) {
    this.renderBuilder = createRenderBuilder(graphics);
    // before the viewport: it sets the spiral of pixels searched around the pick
    this.pickPadding = options.pickPadding;
    this.setViewport(viewport.x, viewport.y, viewport.width, viewport.height);

    if (!checkAsyncPickingSupport(this.buffers)) {
      this.asyncIdentify = (x, y, camera, frame) =>
        createAsyncPickData(() => (this.disposed ? undefined : this.identify(x, y, camera, frame)));
    }
  }
}
