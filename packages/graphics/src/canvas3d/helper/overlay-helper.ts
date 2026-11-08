/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { GraphicsRenderable } from '@molstar/graphics/gl/renderable';
import { Scene } from '@molstar/graphics/gl/scene';
import { getDrawnBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Camera, type ICamera } from '../camera.js';
import { Viewport } from '../camera/util.js';

/** The nearest the near plane of a perspective comes, as a fraction of the far plane (depth precision) */
const MinNearFraction = 1e-3;

const tmpCenter = Vec3();
const tmpOrigin = Vec3();
const tmpTransform = Mat4();

/** What the overlay camera's planes fit: the visible objects' bounds, screen-sized ones (`uScreenSize`) at its scale */
export type OverlayObject = Pick<GraphicsRenderable, 'state' | 'values'>;

/**
 * The world size of one unit of a screen-sized mesh (`uScreenSize`, CSS pixels per unit) with its instance's origin at
 * `origin`, in the view of `camera`: as the mesh vertex stage scales it (`assign_position`)
 */
export function getScreenSizeScale(camera: ICamera, origin: Vec3, screenSize: number, pixelRatio: number) {
  // the view is in scaled units (`Camera.scale`), as the renderer's model matrix
  Vec3.scale(tmpCenter, origin, camera.scale);
  Vec3.transformMat4(tmpCenter, tmpCenter, camera.view);
  const depth = Math.max(-tmpCenter[2], 0);
  // CSS pixels of the camera's viewport, as the renderer's `uCssHeight`
  const cssHeight = camera.viewport.height / pixelRatio;
  const pixelSize = (camera.state.mode === 'orthographic' ? 2 : 2 * depth) / (camera.projection[5] * cssHeight);
  return (screenSize * pixelSize) / camera.scale;
}

/**
 * A camera's view and projection with near and far planes that fit the overlay's objects, without fog: what the
 * overlay is drawn with, so that it is not cut by the scene's planes
 */
export class OverlayCamera implements ICamera {
  readonly viewport = Viewport.create(0, 0, 0, 0);
  readonly view = Mat4();
  readonly projection = Mat4();
  readonly projectionView = Mat4();
  readonly inverseProjectionView = Mat4();
  readonly headRotation = Mat4();
  readonly viewEye = Mat4();
  isAsymmetricProjection = false;

  readonly state = Camera.createDefaultSnapshot();
  readonly viewOffset = Camera.ViewOffset();
  readonly crop = Camera.Crop();
  near = 1;
  far = 10000;
  fogNear = 5000;
  fogFar = 10000;

  forceFull = false;
  scale = 1;
  minTargetDistance = 0;
  disabled = false;

  update(source: ICamera, objects: ArrayLike<OverlayObject>, pixelRatio: number) {
    Viewport.copy(this.viewport, source.viewport);
    Mat4.copy(this.view, source.view);
    Mat4.copy(this.projection, source.projection);
    Mat4.copy(this.headRotation, source.headRotation);
    Mat4.copy(this.viewEye, source.viewEye);
    this.isAsymmetricProjection = source.isAsymmetricProjection;
    Camera.copySnapshot(this.state, source.state);
    this.state.fog = 0;
    Camera.copyViewOffset(this.viewOffset, source.viewOffset);
    if (source.crop) Camera.copyCrop(this.crop, source.crop);
    else this.crop.enabled = false;
    this.forceFull = source.forceFull;
    this.scale = source.scale;
    this.minTargetDistance = source.minTargetDistance;
    this.disabled = source.disabled;

    const [near, far] = getPlanes(source, objects, pixelRatio);
    this.near = near;
    this.far = far;
    this.fogNear = far;
    this.fogFar = far;
    setProjectionPlanes(this.projection, near, far);

    Mat4.mul(this.projectionView, this.projection, this.view);
    Mat4.invert(this.inverseProjectionView, this.projectionView);
    return this;
  }
}

/** The view depths the spheres added (`addSphere`) take */
const depthRange = { min: Infinity, max: -Infinity };

/** Adds a sphere (world units) in the view of `camera` to `depthRange` */
function addSphere(camera: ICamera, center: Vec3, radius: number) {
  Vec3.scale(tmpCenter, center, camera.scale);
  Vec3.transformMat4(tmpCenter, tmpCenter, camera.view);
  const depth = -tmpCenter[2];
  // a little room so that what touches the sphere is not cut
  const r = radius * camera.scale * 1.01;
  depthRange.min = Math.min(depthRange.min, depth - r);
  depthRange.max = Math.max(depthRange.max, depth + r);
}

/**
 * Near and far planes that hold the visible `objects` in the view of `camera`: their bounding spheres, a screen-sized
 * one's per instance at the size the camera gives it (its bounds are of world units); the camera's own when empty
 */
function getPlanes(
  camera: ICamera,
  objects: ArrayLike<OverlayObject>,
  pixelRatio: number,
): [near: number, far: number] {
  depthRange.min = Infinity;
  depthRange.max = -Infinity;
  for (let i = 0, il = objects.length; i < il; ++i) {
    const { state, values } = objects[i];
    if (!state.visible) continue;
    const v = values as { uScreenSize?: { ref: { value: number } } };
    const screenSize = v.uScreenSize?.ref.value ?? 0;
    if (screenSize > 0) {
      const invariant = values.invariantBoundingSphere.ref.value;
      const transforms = values.aTransform.ref.value;
      const objectTransform = values.uObjectTransform?.ref.value;
      for (let j = 0, jl = values.instanceCount.ref.value; j < jl; ++j) {
        Mat4.fromArray(tmpTransform, transforms, j * 16);
        // as drawn: the transient transform after the instance's
        if (objectTransform) Mat4.mul(tmpTransform, objectTransform, tmpTransform);
        Mat4.getTranslation(tmpOrigin, tmpTransform);
        const s = getScreenSizeScale(camera, tmpOrigin, screenSize, pixelRatio);
        const scale = s * Mat4.getMaxScaleOnAxis(tmpTransform);
        const center = Vec3.transformMat4(tmpOrigin, Vec3.scale(tmpOrigin, invariant.center, s), tmpTransform);
        addSphere(camera, center, invariant.radius * scale);
      }
    } else {
      const sphere = getDrawnBoundingSphere(values);
      if (sphere.radius > 0) addSphere(camera, sphere.center, sphere.radius);
    }
  }
  if (depthRange.min > depthRange.max) return [camera.near, camera.far];
  const far = depthRange.max;
  if (camera.state.mode === 'orthographic') return [depthRange.min, far];
  // all behind the camera: nothing to draw
  if (far <= 0) return [camera.near, camera.far];
  return [Math.max(depthRange.min, far * MinNearFraction), far];
}

/**
 * Sets the near and far planes of a GL projection (perspective or orthographic), keeping the rest of it (its frustum
 * at the near plane, offsets, a sample's jitter)
 */
export function setProjectionPlanes(projection: Mat4, near: number, far: number) {
  if (projection[11] === -1) {
    // perspective
    projection[10] = -(far + near) / (far - near);
    projection[14] = (-2 * far * near) / (far - near);
  } else {
    // orthographic
    projection[10] = -2 / (far - near);
    projection[14] = -(far + near) / (far - near);
  }
  return projection;
}

/**
 * The overlay layer (`Representation.Layer`): the representations drawn above the scene, anchored in the world. They
 * have a scene of their own, out of the canvas scene's bounds, drawn after the scene with a depth of their own
 * (cleared) and near and far planes fitting their bounds (`getCamera`), and picked after the scene's.
 */
export class OverlayHelper {
  readonly scene: Scene;
  private cameras = new WeakMap<ICamera, OverlayCamera>();
  private graphics: GraphicsContext;
  /** the scene's opaque depth the overlay was last drawn over (`OverlayPass`), which its picking tests against */
  occluderDepth: Texture | undefined = undefined;
  /** the scene's front transparent depth of that draw, where it wrote one */
  occluderDepthAlpha: Texture | undefined = undefined;

  constructor(graphics: GraphicsContext) {
    this.graphics = graphics;
    this.scene = Scene.create(graphics, 'blended');
  }

  get isEnabled() {
    return this.scene.count > 0;
  }

  /** The overlay camera of `source` (one for each, e.g. a stereo eye, an image pass's own), updated to it */
  getCamera(source: ICamera): OverlayCamera {
    let camera = this.cameras.get(source);
    if (!camera) {
      camera = new OverlayCamera();
      this.cameras.set(source, camera);
    }
    return camera.update(source, this.scene.renderables, this.graphics.pixelRatio);
  }

  /** Whether the render object of `objectId` (a `PickingId`'s) is in the overlay */
  hasObject(objectId: number) {
    let has = false;
    this.scene.forEach((_, o) => {
      if (o.id === objectId) has = true;
    });
    return has;
  }

  /** This overlay, drawn only while `enabled` also says so (e.g. an image pass's own toggle) */
  withEnabled(enabled: () => boolean): OverlayHelper {
    const overlay = this;
    return Object.create(this, {
      isEnabled: { get: () => enabled() && overlay.isEnabled },
    });
  }

  setContext(graphics: GraphicsContext) {
    this.graphics = graphics;
    this.occluderDepth = undefined;
    this.occluderDepthAlpha = undefined;
    return this.scene.setContext(graphics);
  }
}
