/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The renderer's GPU culling and sphere raster: merged objects culled by compute kernels (two phases where the opaque
// pass is split by occlusion), and the small opaque spheres drawn by the raster and resolved before the impostors

import { Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import type { ValueCell } from '@molstar/core/util';
import type { Viewport } from '@molstar/core/util/viewport';
import type { RenderCamera } from '../camera.js';
import { isDebugMode } from '@molstar/core/util/debug';
import type { OcclusionPyramid } from '../compute/cell-cull.js';
import {
  createSphereRasterPass,
  isSphereRasterEligible,
  type SphereRasterPass,
} from '../compute/sphere-raster-pass.js';
import type { GraphicsContext } from '../context.js';
import { projectSphere } from '../project-sphere.js';
import type { GlobalBindings } from '../render-item.js';
import { type CullBatch, type Frame, type GraphicsRenderable, LodView } from '../renderable.js';
import type { MergedRasterParams } from '../renderable/merged-culler.js';
import type { GlobalUniformValues } from '../renderable/schema.js';
import type { SpheresValues } from '../renderable/spheres.js';
import { getDrawnBoundingSphere } from '../renderable/util.js';
import type { Scene } from '../scene.js';
import { OpaqueState } from './pass-draws.js';

/** The renderer's props the culling reads (`RendererParams`), read when used */
export interface GpuCullingProps {
  readonly sphereRaster: 'off' | 'auto' | 'on' | 'only';
  readonly sphereRasterRadius: number;
  readonly sphereRasterObjects: number;
}

/** What the culling reads of the renderer: its view of the last `update` and the frame (a getter: read when culled) */
export interface GpuCullingHost {
  readonly ctx: GraphicsContext;
  readonly props: GpuCullingProps;
  readonly globals: GlobalBindings;
  readonly globalUniforms: GlobalUniformValues;
  /** the renderer's viewport, the target's region the passes draw */
  readonly viewport: Viewport;
  /** the culling's camera plane and frustum: in the unscaled world (the renderer's `cullFrustum`, any camera scale) */
  readonly cameraPlane: Plane3D;
  readonly frustum: Frustum3D;
  readonly lodView: LodView;
  readonly frame: Frame;
}

export interface GpuCulling {
  /** The objects the last `prepare` culled on the GPU */
  readonly culledCount: number;
  /** the sum of the `recordsVersion`s of the objects of the last `prepare`: the draws a recording holds */
  readonly records: number;
  /** The objects whose second phase is to run after the last `prepare` with `twoPhase` */
  readonly secondPending: number;
  /** The view of `camera` (an `update`): the raster's split for the next `prepare`, the scale culled on the GPU at */
  setView: (camera: RenderCamera) => void;
  setOcclusionPyramid: (get: (() => OcclusionPyramid | undefined) | null) => void;
  /** Culls on the GPU with the planes of the last `update`, outside of passes (`Renderer.prepareCull`) */
  prepare: (scene: Scene.Group, twoPhase?: boolean, clipGhost?: boolean) => void;
  /**
   * Culls on the GPU for a view of its own with the planes of the last `update`, outside of passes, the objects that
   * draw nothing without (`Renderable.prepareCullView`, `Renderer.prepareCullView`)
   */
  prepareView: (scene: Scene.Group) => void;
  /** The second phase against `pyramid`, of the depth so far, and the raster's; none: left out */
  prepareSecond: (scene: Scene.Group, pyramid: OcclusionPyramid | undefined) => void;
  /** The second part drawn: no second phase pending */
  secondDone: () => void;
  /** The raster's spheres of `group` drawn for the next opaque pass, which resolves them (`Renderer.prepareSphereRaster`) */
  prepareRaster: (group: Scene.Group, camera: RenderCamera) => void;
  /** whether the raster is on and the device has it */
  isRasterAvailable: () => boolean;
  /** whether the raster may draw `r`, a precompile list's stand-in too */
  isRasterCandidate: (r: GraphicsRenderable) => boolean;
  /** whether the raster draws `r` instead of the opaque pass in progress (`sphereRaster` 'only') */
  isRasterDrawn: (r: GraphicsRenderable) => boolean;
  /** The raster's spheres resolved into the open pass, before the impostors; the viewport's scissor after */
  resolveRaster: () => void;
  /** An opaque pass drawn after `resolveRaster`: with 'only' the raster's objects not drawn as impostors */
  beginOpaque: () => void;
  endOpaque: () => void;
  /** The raster's values of the objects reset, after the last part of the opaque pass that resolves it */
  finishRaster: () => void;
  /** The raster turned off: its pass released */
  rasterOff: () => void;
  destroy: () => void;
}

export function createGpuCulling(host: GpuCullingHost): GpuCulling {
  const { ctx, props: p, globals, globalUniforms, viewport, cameraPlane, frustum, lodView } = host;

  /**
   * the camera scale of the last `update`: culled in the unscaled world at any (`GpuCullingHost.frustum`), tested
   * against the occlusion pyramid only unscaled (its spheres scaled), as `cull` (`cullObject`)
   */
  let updateScale = 1;
  const getPyramid = () => (updateScale === 1 ? getOcclusionPyramid?.() : undefined);
  let getOcclusionPyramid: (() => OcclusionPyramid | undefined) | null = null;
  let culledCount = 0;
  let records = 0;
  let secondPending = 0;
  /** the objects' culling recorded and run together (`CullBatch`), reused */
  const cullBatch: CullBatch = { clears: [], dispatches: [] };
  const runCullBatch = (label: string) => {
    const { clears, dispatches } = cullBatch;
    for (let i = 0, il = clears.length; i < il; ++i) clears[i]();
    if (dispatches.length > 0) {
      ctx.compute!.beginComputePass({ label });
      for (let i = 0, il = dispatches.length; i < il; ++i) dispatches[i]();
      ctx.compute!.endComputePass();
    }
    clears.length = 0;
    dispatches.length = 0;
  };
  /** the visible objects' `rasterSpheres`, reused */
  const rasterSpheresOf: number[] = [];
  const descending = (a: number, b: number) => b - a;
  /**
   * The spheres from which an object culled on the GPU for the raster only (`Renderable.rasterSpheres`) is culled there
   * this frame: none where the raster does not draw the view (`'off'` too); `'auto'` the `sphereRasterObjects`
   * largest visible ones (each object costs its own dispatches, whatever its size); `'on'` and `'only'` all of them
   */
  const getRasterSpheresMin = (renderables: ReadonlyArray<GraphicsRenderable>) => {
    if (!cullRaster) return Infinity;
    if (p.sphereRaster !== 'auto') return 0;
    const budget = p.sphereRasterObjects;
    if (budget <= 0) return Infinity;
    rasterSpheresOf.length = 0;
    for (let i = 0, il = renderables.length; i < il; ++i) {
      const r = renderables[i];
      const spheres = r.rasterSpheres;
      if (spheres !== undefined && r.state.visible) rasterSpheresOf.push(spheres);
    }
    if (rasterSpheresOf.length <= budget) return 0;
    rasterSpheresOf.sort(descending);
    return rasterSpheresOf[budget - 1];
  };
  const prepare = (scene: Scene.Group, twoPhase = false, clipGhost = false) => {
    culledCount = 0;
    records = 0;
    secondPending = 0;
    const pyramid = getPyramid();
    const { renderables } = scene;
    const rasterSpheresMin = getRasterSpheresMin(renderables);
    const frame = host.frame;
    for (let i = 0, il = renderables.length; i < il; ++i) {
      const r = renderables[i];
      const spheres = r.rasterSpheres;
      // left to the CPU culling: drawn by its lists, as where the kernels are not ready; hidden, not culled at all
      if (spheres !== undefined && (spheres < rasterSpheresMin || !r.state.visible)) continue;
      if (r.prepareCull?.(cameraPlane, frustum, frame, lodView, pyramid, cullRaster, cullBatch, clipGhost, twoPhase)) {
        culledCount += 1;
      }
      records += r.recordsVersion ?? 0;
      if (twoPhase && r.needsCullSecond?.(frame)) secondPending += 1;
    }
    runCullBatch('merged cell culling');
  };
  const prepareView = (scene: Scene.Group) => {
    const { renderables } = scene;
    const frame = host.frame;
    for (let i = 0, il = renderables.length; i < il; ++i) {
      const r = renderables[i];
      if (r.state.visible) r.prepareCullView?.(cameraPlane, frustum, frame, lodView);
    }
  };
  const prepareSecond = (scene: Scene.Group, pyramid: OcclusionPyramid | undefined) => {
    const { renderables } = scene;
    const frame = host.frame;
    for (let i = 0, il = renderables.length; i < il; ++i) {
      renderables[i].prepareCullSecond?.(frame, pyramid, cullBatch);
    }
    runCullBatch('merged cell culling, second phase');
    if (pyramid) sphereRasterPass?.second();
    else sphereRasterPass?.skipSecond();
  };

  let sphereRasterPass: SphereRasterPass | undefined;
  /**
   * the sphere raster's split per instance for the GPU culling of the next `prepare` (merged objects), from the camera
   * of the last `update`; undefined where the raster does not draw the view (`prepareRaster`)
   */
  let cullRaster: MergedRasterParams | undefined;
  const getCullRaster = (camera: RenderCamera): MergedRasterParams | undefined => {
    if (!isSphereRasterView(camera)) return undefined;
    const P = camera.projection;
    const ortho = camera.state.mode === 'orthographic';
    // the culler's sizes unscaled: perspective over the distance (a ratio), orthographic in view units (`camera.scale`)
    const focal = Math.max(P[0] * viewport.width, P[5] * viewport.height) / 2;
    return {
      radius: getSphereRasterRadius(),
      scale: ortho ? focal * camera.scale : focal,
      ortho,
      // the culler's distances unscaled (`GpuCullingHost.cameraPlane`)
      near: camera.near / camera.scale,
    };
  };
  const isRasterAvailable = () => p.sphereRaster !== 'off' && !!ctx.compute && ctx.capabilities.atomicMinMax64;
  /**
   * whether the raster may draw `r` (`isSphereRasterEligible`); a precompile list's stand-in (`collectSignatures`,
   * without the values it reads) where it is of spheres and opaque
   */
  const isRasterCandidate = (r: GraphicsRenderable) => {
    const values = r.values as { [k: string]: ValueCell<any> | undefined };
    if ('dSphereResolve' in values) {
      // `auto`: objects that can be culled on the GPU only
      if (p.sphereRaster === 'auto' && !r.getRasterInput) return false;
      return isSphereRasterEligible(r.values as SpheresValues);
    }
    // `auto`: of objects culled on the GPU (unmerged ones never are, their signatures without the flag)
    if (p.sphereRaster === 'auto' && (r as { gpuCulled?: boolean }).gpuCulled !== true) return false;
    return values.dGeometryType?.ref.value === 'spheres' && values.uAlpha?.ref.value === 1;
  };
  /**
   * whether the raster draws the view: perspective or orthographic, at any scale; a viewport anywhere in its target
   * (stereo's right eye), the visibility buffer the view's, the resolve's pixels offset by its origin; not one starting
   * outside the target (cropped by the backend where the device does not take it)
   */
  const isSphereRasterView = (_camera: RenderCamera) => {
    if (!isRasterAvailable()) return false;
    return viewport.x >= 0 && viewport.y >= 0;
  };
  /** the largest half size in pixels of the target the raster draws */
  const getSphereRasterRadius = () =>
    p.sphereRaster === 'only' ? 1e9 : p.sphereRasterRadius * globalUniforms.uPixelRatio.ref.value;
  /** in the opaque pass with `sphereRaster` 'only': the raster's objects not drawn as impostors */
  let sphereRasterOnly = false;
  const sphereRasterModel = Mat4();
  const sphereRasterModelView = Mat4();
  const sphereRasterProjection = Mat4();
  let sphereRasterNear = 0;
  const sphereRasterCenter = Vec3();
  const sphereRasterRect = Vec4();
  /**
   * the resolve of `r` scissored to its bounding sphere on screen (a pixel's margin), the viewport where the sphere
   * reaches the near plane; false where it is off screen
   */
  const scissorSphereResolve = (r: GraphicsRenderable) => {
    const { center, radius } = getDrawnBoundingSphere(r.values);
    const c = Vec3.transformMat4(sphereRasterCenter, center, sphereRasterModelView);
    const s = radius * Mat4.getMaxScaleOnAxis(sphereRasterModelView);
    const { x, y, width, height } = viewport;
    if (-c[2] - s <= sphereRasterNear) {
      ctx.setScissor(x, y, width, height);
      return true;
    }
    const rect = projectSphere(sphereRasterRect, c, s, sphereRasterProjection);
    const x0 = Math.max(0, Math.floor(rect[0] * width) - 1);
    const y0 = Math.max(0, Math.floor(rect[1] * height) - 1);
    const x1 = Math.min(width, Math.ceil(rect[2] * width) + 1);
    const y1 = Math.min(height, Math.ceil(rect[3] * height) + 1);
    if (x1 <= x0 || y1 <= y0) return false;
    ctx.setScissor(x + x0, y + y0, x1 - x0, y1 - y0);
    return true;
  };
  const prepareRaster = (group: Scene.Group, camera: RenderCamera) => {
    sphereRasterPass?.finish();
    if (!isSphereRasterView(camera)) return;
    sphereRasterPass ??= createSphereRasterPass(ctx.compute!, ctx.resources);
    Mat4.scaleUniformly(sphereRasterModel, group.view, camera.scale);
    Mat4.mul(sphereRasterModelView, camera.view, sphereRasterModel);
    Mat4.copy(sphereRasterProjection, camera.projection);
    sphereRasterNear = camera.near;
    sphereRasterPass.prepare(group.renderables, {
      modelView: sphereRasterModelView,
      projection: camera.projection,
      cameraPlane: globalUniforms.uCameraPlane.ref.value,
      model: sphereRasterModel,
      modelScale: camera.scale,
      near: camera.near,
      width: viewport.width,
      height: viewport.height,
      origin: [viewport.x, viewport.y],
      maxRadius: getSphereRasterRadius(),
      lod: lodView,
      pyramid: getPyramid(),
      twoPhase: secondPending > 0,
      cull: { frame: host.frame, cameraPlane, frustum },
      culledOnly: p.sphereRaster === 'auto',
    });
  };

  let cullRasterLogged = '';
  const setView = (camera: RenderCamera) => {
    cullRaster = getCullRaster(camera);
    updateScale = camera.scale;
    if (isDebugMode && p.sphereRaster !== 'off') {
      const text = `sphere raster culling: ${cullRaster ? JSON.stringify(cullRaster) : 'off for this view'}`;
      if (text !== cullRasterLogged) console.log(text);
      cullRasterLogged = text;
    }
  };

  return {
    get culledCount() {
      return culledCount;
    },
    get records() {
      return records;
    },
    get secondPending() {
      return secondPending;
    },
    setView,
    setOcclusionPyramid: (get) => {
      getOcclusionPyramid = get;
    },
    prepare,
    prepareView,
    prepareSecond,
    secondDone: () => {
      secondPending = 0;
    },
    prepareRaster,
    isRasterAvailable,
    isRasterCandidate,
    isRasterDrawn: (r) => sphereRasterOnly && !!sphereRasterPass?.isDrawn(r),
    resolveRaster: () => {
      if (!sphereRasterPass?.active) return;
      ctx.setRenderState(OpaqueState);
      ctx.stats.sphereRasterResolves += sphereRasterPass.resolve(globals, scissorSphereResolve);
      ctx.setScissor(viewport.x, viewport.y, viewport.width, viewport.height);
    },
    beginOpaque: () => {
      sphereRasterOnly = p.sphereRaster === 'only';
    },
    endOpaque: () => {
      sphereRasterOnly = false;
    },
    finishRaster: () => {
      sphereRasterPass?.finish();
    },
    rasterOff: () => {
      sphereRasterPass?.destroy();
      sphereRasterPass = undefined;
    },
    destroy: () => {
      sphereRasterPass?.destroy();
    },
  };
}
