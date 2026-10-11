/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Ray3D } from '@molstar/core/math/geometry/primitives/ray3d';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { degToRad, spiral2d } from '@molstar/core/math/misc';
import { isTimingMode } from '@molstar/core/util/debug';
import { Camera } from '../camera.js';
import { cameraUnproject } from '../camera/util.js';
import { Viewport } from '../camera/util.js';
import type { Helper } from './helper.js';
import {
  type AsyncPickData,
  createPickBuffers,
  type PickBuffers,
  type PickData,
  PickPass,
  type PickOptions,
  checkAsyncPickingSupport,
  AsyncPickStatus,
} from '../passes/pick.js';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';

export class RayHelper {
  private viewport = Viewport();
  private size: number;
  private spiral: [number, number][];

  private pickPadding: number;
  private camera: Camera;
  private pickPass: PickPass;
  private buffers: PickBuffers;

  setPickPadding(pickPadding: number) {
    if (this.pickPadding !== pickPadding) {
      this.pickPadding = pickPadding;
      this.update();
    }
  }

  private update() {
    const size = this.pickPadding * 2 + 1;
    Viewport.set(this.viewport, 0, 0, size, size);
    this.buffers.setViewport(0, 0, size, size);

    this.spiral = spiral2d(this.pickPadding);
    this.size = size;

    this.pickPass.setSize(size, size);
  }

  private render(camera: Camera) {
    if (isTimingMode) this.graphics.timer.group('RayHelper.render', { captureStats: true });
    const { renderer, scene, helper } = this;

    renderer.setTransparentBackground(false);
    renderer.setDrawingBufferSize(this.size, this.size);
    renderer.setPixelRatio(1);

    renderer.setViewport(0, 0, this.size, this.size);
    // this camera is synthetic (ray-aligned), never the main camera - always a fresh,
    // never-reused frame token, so it can never share/collide with the main draw's frame
    this.pickPass.render(renderer, camera, scene, helper, createFrame());

    if (isTimingMode) this.graphics.timer.groupEnd('RayHelper.render');
  }

  private identifyInternal(x: number, y: number): PickData | undefined {
    if (this.graphics.isContextLost) return;

    const { viewport } = this;

    const pickingId = this.buffers.getPickingId(x, y);
    if (pickingId === undefined) return;

    const z = this.buffers.getDepth(x, y);
    const position = Vec3.create(x, y, z);
    cameraUnproject(position, position, viewport, this.camera.inverseProjectionView);

    return { id: pickingId, position };
  }

  private prepare(ray: Ray3D, cam: Camera) {
    this.camera.far = cam.far;
    this.camera.near = cam.near;
    this.camera.fogFar = cam.fogFar;
    this.camera.fogNear = cam.fogNear;
    this.camera.forceFull = cam.forceFull;
    this.camera.scale = cam.scale;
    Viewport.copy(this.camera.viewport, this.viewport);
    Camera.copySnapshot(this.camera.state, { ...cam.state, mode: 'orthographic' });

    updateOrthoRayCamera(this.camera, ray, cam.up);
    Mat4.mul(this.camera.projectionView, this.camera.projection, this.camera.view);
    Mat4.tryInvert(this.camera.inverseProjectionView, this.camera.projectionView);

    Mat4.copy(this.camera.viewEye, cam.view);
  }

  private getPickData(): PickData | undefined {
    const c = this.pickPadding;
    for (const d of this.spiral) {
      const pickData = this.identifyInternal(c + d[0], c + d[1]);
      if (pickData) return pickData;
    }
  }

  sphere = Sphere3D();

  private intersectsScene(ray: Ray3D, scale: number): boolean {
    Sphere3D.scaleNX(this.sphere, this.scene.boundingSphereVisible, scale);
    return Ray3D.isInsideSphere3D(ray, this.sphere) || Ray3D.isIntersectingSphere3D(ray, this.sphere);
  }

  identify(ray: Ray3D, cam: Camera): PickData | undefined {
    // the backend can't read the pick targets (and render the pick variant) yet
    if (!this.buffers.supported) return;
    if (!this.intersectsScene(ray, cam.scale)) return;

    this.prepare(ray, cam);

    if (isTimingMode) this.graphics.timer.group('RayHelper.identify');
    this.render(this.camera);
    this.buffers.read();
    if (isTimingMode) this.graphics.timer.groupEnd('RayHelper.identify');

    return this.getPickData();
  }

  asyncIdentify(ray: Ray3D, cam: Camera): AsyncPickData | undefined {
    if (!this.buffers.supported) return;
    if (!this.intersectsScene(ray, cam.scale)) return;

    this.prepare(ray, cam);

    if (isTimingMode) this.graphics.timer.group('RayHelper.asyncIdentify');
    this.render(this.camera);
    this.buffers.asyncRead();
    if (isTimingMode) this.graphics.timer.groupEnd('RayHelper.asyncIdentify');

    return {
      tryGet: () => {
        const status = this.buffers.check();
        if (status === AsyncPickStatus.Resolved) {
          return this.getPickData();
        } else if (status === AsyncPickStatus.Pending) {
          return 'pending';
        }
      },
    };
  }

  reset() {
    this.buffers.reset();
    this.pickPass.reset();
  }

  dispose() {
    this.buffers.dispose();
    this.pickPass.dispose();
  }

  constructor(
    private graphics: GraphicsContext,
    private renderer: Renderer,
    private scene: Scene,
    private helper: Helper,
    options: PickOptions,
  ) {
    const size = options.pickPadding * 2 + 1;

    this.camera = new Camera();
    this.pickPass = new PickPass(graphics, size, size, 1);
    this.buffers = createPickBuffers(this.graphics, this.pickPass, options.maxAsyncReadLag);
    this.pickPadding = options.pickPadding;

    this.update();

    if (!checkAsyncPickingSupport(graphics)) {
      this.asyncIdentify = (ray, cam) => ({
        tryGet: () => this.identify(ray, cam),
      });
    }
  }
}

//

function updateOrthoRayCamera(camera: Camera, ray: Ray3D, up: Vec3) {
  const { near, far, viewport } = camera;

  const height = 2 * Math.tan(degToRad(0.1) / 2) * Vec3.distance(camera.position, camera.target) * camera.scale;
  const zoom = viewport.height / height;

  const fullLeft = -viewport.width / 2;
  const fullRight = viewport.width / 2;
  const fullTop = viewport.height / 2;
  const fullBottom = -viewport.height / 2;

  const dx = (fullRight - fullLeft) / (2 * zoom);
  const dy = (fullTop - fullBottom) / (2 * zoom);
  const cx = (fullRight + fullLeft) / 2;
  const cy = (fullTop + fullBottom) / 2;

  const left = cx - dx;
  const right = cx + dx;
  const top = cy + dy;
  const bottom = cy - dy;

  // build projection matrix
  Mat4.ortho(camera.projection, left, right, top, bottom, near, far);

  const direction = Vec3.normalize(Vec3(), ray.direction);
  const r = Quat.fromUnitVec3(Quat(), direction, Vec3.negUnitZ);
  Quat.invert(r, r);

  const eye = Vec3.clone(ray.origin);
  const target = Vec3.add(Vec3(), eye, direction);

  // build view matrix
  Mat4.lookAt(camera.view, eye, target, up);
}
