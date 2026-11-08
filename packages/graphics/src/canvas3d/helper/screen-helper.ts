/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import { Representation } from '@molstar/graphics/repr/representation';
import { Camera, type ICamera } from '../camera.js';
import { Viewport } from '../camera/util.js';

const tmpRotation = Mat4();
const tmpTranslation = Vec3();
const tmpFrame = Viewport.zero();

/**
 * The camera of one object of the screen layer: orthographic in the CSS pixels of a camera's frame (`Camera.getFrame`:
 * the cropped part of it projected, the jitter's view offset on top), its view the object's place
 * (`Representation.ScreenPlacement`: the anchor's translation, the camera's rotation where followed), so that culling
 * sees where the object is
 */
export class ScreenCamera implements ICamera {
  readonly viewport = Viewport.create(0, 0, 0, 0);
  readonly view = Mat4();
  readonly projection = Mat4();
  readonly projectionView = Mat4();
  readonly inverseProjectionView = Mat4();
  readonly headRotation = Mat4.zero();
  readonly viewEye = Mat4.zero();
  readonly isAsymmetricProjection = false;

  readonly state = Camera.createDefaultSnapshot();
  readonly viewOffset = Camera.ViewOffset();
  readonly crop = Camera.Crop();
  near = 1;
  far = 10000;
  fogNear = 10000;
  fogFar = 10000;

  readonly forceFull = false;
  readonly scale = 1;
  readonly minTargetDistance = 0;
  disabled = false;

  constructor() {
    this.state.mode = 'orthographic';
    this.state.fog = 0;
  }

  /** The camera of an object of `radius` about its origin (CSS pixels) placed by `placement` in the view of `source` */
  update(
    source: ICamera,
    placement: Representation.ScreenPlacement,
    radius: number,
    inset: number,
    pixelRatio: number,
  ) {
    Viewport.copy(this.viewport, source.viewport);
    Camera.copyViewOffset(this.viewOffset, source.viewOffset);
    if (source.crop) Camera.copyCrop(this.crop, source.crop);
    else this.crop.enabled = false;
    this.disabled = source.disabled;

    // anchored in the frame, cropped with the rest where its edge is outside the canvas
    const frame = Camera.getFrame(tmpFrame, source);
    const width = frame.width / pixelRatio;
    const height = frame.height / pixelRatio;
    const { anchor, offsetX, offsetY } = placement;
    // on an edge: inward by the inset and the offset; in the middle of an axis: the offset as is
    const side = Representation.screenAnchorSides(anchor);
    const x = side.x === 0 ? offsetX : side.x * (width / 2 - inset - offsetX);
    const y = side.y === 0 ? offsetY : side.y * (height / 2 - inset - offsetY);

    // the eye in front of the object, as far as it is large
    const r = radius + 1;
    const distance = 2 * r;
    this.near = distance - r;
    this.far = distance + r;
    this.fogNear = this.far;
    this.fogFar = this.far;

    if (placement.followRotation) Mat4.extractRotation(tmpRotation, source.view);
    else Mat4.setIdentity(tmpRotation);
    Mat4.fromTranslation(this.view, Vec3.set(tmpTranslation, x, y, -distance));
    Mat4.mul(this.view, this.view, tmpRotation);

    setScreenProjection(this.projection, this.viewport, this.crop, this.viewOffset, pixelRatio, this.near, this.far);
    Mat4.mul(this.projectionView, this.projection, this.view);
    Mat4.invert(this.inverseProjectionView, this.projectionView);

    // the eye and the object's origin in the object's frame, for the renderer's camera position and direction
    const invView = Mat4.invert(tmpRotation, this.view);
    Mat4.getTranslation(this.state.position, invView);
    Vec3.set(this.state.target, 0, 0, 0);
    return this;
  }
}

/**
 * An orthographic projection in CSS pixels, the frame centered at 0: of the drawn viewport, where it is cropped its part
 * of the frame (window y down, view y up), shifted by the jitter (`viewOffset`, device pixels)
 */
function setScreenProjection(
  out: Mat4,
  viewport: Viewport,
  crop: Camera.Crop,
  viewOffset: Camera.ViewOffset,
  pixelRatio: number,
  near: number,
  far: number,
) {
  let cx = 0;
  let cy = 0;
  if (crop.enabled) {
    cx = crop.offsetX + viewport.width / 2 - crop.fullWidth / 2;
    cy = -(crop.offsetY + viewport.height / 2 - crop.fullHeight / 2);
  }
  if (viewOffset.enabled) {
    cx += viewOffset.offsetX;
    cy -= viewOffset.offsetY;
  }
  const left = (cx - viewport.width / 2) / pixelRatio;
  const right = (cx + viewport.width / 2) / pixelRatio;
  const top = (cy + viewport.height / 2) / pixelRatio;
  const bottom = (cy - viewport.height / 2) / pixelRatio;
  return Mat4.ortho(out, left, right, top, bottom, near, far);
}

interface Entry {
  readonly scene: Scene;
  placement: Representation.ScreenPlacement;
}

/**
 * The screen layer (`Representation.Layer`): objects anchored to the viewport, above the scene and the overlay, their
 * unit a CSS pixel. A scene for each owner (a representation, the corner axes) with its place, drawn with its own
 * camera (`ScreenCamera`) after the overlay, its depth cleared again, and picked last; off in XR (no screen).
 */
export class ScreenHelper {
  private readonly entries = new Map<object, Entry>();
  /** `scenes`, made again where the entries change: read every tick (`Canvas3D` commits the layers) */
  private sceneList: Scene[] | undefined = undefined;
  private readonly cameras = new WeakMap<ICamera, ScreenCamera>();

  constructor(private graphics: GraphicsContext) {}

  get isEnabled() {
    if (this.graphics.xr.session) return false;
    for (const { scene } of this.entries.values()) {
      if (scene.count > 0) return true;
    }
    return false;
  }

  /** The scene of `owner`, created on first use, placed by `placement` from now on */
  getScene(owner: object, placement: Representation.ScreenPlacement): Scene {
    let entry = this.entries.get(owner);
    if (!entry) {
      entry = { scene: Scene.create(this.graphics, 'blended'), placement };
      this.entries.set(owner, entry);
      this.sceneList = undefined;
    }
    entry.placement = placement;
    return entry.scene;
  }

  /** Whether `scene` is one of the layer's */
  has(scene: Scene) {
    for (const entry of this.entries.values()) {
      if (entry.scene === scene) return true;
    }
    return false;
  }

  /** Removes the scene of `owner` and its objects */
  release(owner: object) {
    const entry = this.entries.get(owner);
    if (!entry) return;
    entry.scene.clear();
    this.entries.delete(owner);
    this.sceneList = undefined;
  }

  get scenes(): readonly Scene[] {
    if (!this.sceneList) this.sceneList = Array.from(this.entries.values(), (e) => e.scene);
    return this.sceneList;
  }

  /**
   * Calls `draw` for each scene with objects, with its camera for `source` (one camera, updated before each call:
   * a draw takes its uniforms at once)
   */
  forEach(source: ICamera, draw: (scene: Scene, camera: ScreenCamera) => void) {
    let camera = this.cameras.get(source);
    if (!camera) {
      camera = new ScreenCamera();
      this.cameras.set(source, camera);
    }
    for (const { scene, placement } of this.entries.values()) {
      if (scene.count === 0) continue;
      const { center, radius } = scene.boundingSphere;
      // about its origin: the camera turns it there
      const extent = Vec3.magnitude(center) + radius;
      draw(scene, camera.update(source, placement, extent, placement.inset ?? radius, this.graphics.pixelRatio));
    }
  }

  /** The render objects dropped by the scenes, see `Scene.setContext` */
  setContext(graphics: GraphicsContext): GraphicsRenderObject[] {
    this.graphics = graphics;
    const dropped: GraphicsRenderObject[] = [];
    for (const { scene } of this.entries.values()) dropped.push(...scene.setContext(graphics));
    return dropped;
  }

  clear() {
    for (const { scene } of this.entries.values()) scene.clear();
    this.entries.clear();
    this.sceneList = undefined;
  }
}
