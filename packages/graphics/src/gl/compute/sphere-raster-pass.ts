/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { isDebugMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import type { ComputeContext, StorageBuffer } from '../compute.js';
import type { GraphicsResources } from '../context.js';
import type { GlobalBindings } from '../render-item.js';
import type { Frame, Renderable } from '../renderable.js';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import type { SpheresValues } from '../renderable/spheres.js';
import { isStorageClipSupported } from '../clip-storage.js';
import {
  canRasterSpheres,
  createSphereRaster,
  SphereRasterKernels,
  type SphereRasterObject,
  SphereRasterVisibilityBytes,
  type SphereRasterView,
} from './sphere-raster.js';

/**
 * Whether the raster can draw the spheres of `values` in the opaque color pass as the impostors draw them: opaque (alpha
 * 1, no transparency layer, not x-ray), clip objects but by instance with a clipping layer (the kernels clip instances,
 * spheres and pixels as the vertex and fragment stages, `isStorageClipSupported`), the exact impostors, no animation,
 * no fade by the object's distance (`fadeLod`'s dither), levels within the kernels'
 */
export function isSphereRasterEligible(values: SpheresValues) {
  if (values.uAlpha.ref.value !== 1 || values.dTransparency.ref.value) return false;
  if (values.dXrayShaded.ref.value !== 'off' || values.dApproximate.ref.value) return false;
  if (!isStorageClipSupported(values)) return false;
  if (values.dWiggleAnimation.ref.value || values.dTumbleAnimation.ref.value) return false;
  const lod = values.uLod.ref.value;
  if (values.dLodFade.ref.value && lod[3] === 0 && (lod[0] !== 0 || lod[1] !== 0)) return false;
  return canRasterSpheres(values);
}

/** What a view draws the raster with, from the renderer's camera and model */
export interface SphereRasterPassView extends Omit<SphereRasterView, 'cameraPlane'> {
  /** the renderer's camera plane (world space); into the objects' space by the model matrix */
  readonly cameraPlane: Vec4;
  readonly model: Mat4;
  /** the view's frame and culling planes: merged objects only in the view culled on the GPU */
  readonly cull: { readonly frame: Frame; readonly cameraPlane: Plane3D; readonly frustum: Frustum3D };
  /** a second phase of GPU culling follows (`SphereRasterPass.second`) */
  readonly twoPhase?: boolean;
  /**
   * objects culled on the GPU only (`sphereRaster` 'auto'): unmerged objects' split per sphere leaves their quads drawn
   * all the same, the raster's work added (`sphere-compute-raster.md` step 6)
   */
  readonly culledOnly?: boolean;
}

/**
 * The sphere raster in the renderer's opaque pass (`sphere-compute-raster.md`): `prepare` before the pass (no pass
 * open) dispatches the raster of the eligible spheres of a view and sets the drawn objects' `uRaster*` values, so their
 * impostors leave out its spheres in the color pass; `resolve` in the opaque pass draws their resolves; `finish`
 * after it, the impostors draw all spheres again (other views and passes of the frame). Merged objects culled on the
 * GPU leave out the raster regions by `Renderable.rasterDrawn` (`resolve` to `finish`) and keep their values until a
 * view no longer draws them: their impostors don't read them, and the opaque pass's recordings replay unchanged
 */
export interface SphereRasterPass {
  prepare: (renderables: readonly Renderable<any>[], view: SphereRasterPassView) => void;
  /**
   * The second phase of GPU culling, after `resolve` (`prepare` announced it, `twoPhase`): what it found drawn too, for a
   * `resolve` again, the first's pixels keeping their depth
   */
  second: () => void;
  /** No second phase after a `prepare` announcing it */
  skipSecond: () => void;
  /**
   * draws the resolves of the drawn objects, each after `scissor` (its screen rectangle: false where it is off screen,
   * not drawn); the caller restores the scissor after. Returns the count drawn
   */
  resolve: (globals: GlobalBindings, scissor?: (r: Renderable<any>) => boolean) => number;
  finish: () => void;
  /** whether the raster draws `r` (its resolve drawn) until `finish` */
  isDrawn: (r: Renderable<any>) => boolean;
  /** whether objects of the last `prepare` are drawn by the raster until `finish` */
  readonly active: boolean;
  destroy: () => void;
}

export function createSphereRasterPass(compute: ComputeContext, resources: GraphicsResources): SphereRasterPass {
  const raster = createSphereRaster(compute, resources);
  let vis: StorageBuffer | undefined;
  /** the renderables drawn by the raster since `prepare` */
  const drawn: Renderable<any>[] = [];
  /** the merged objects (culled on the GPU) of the last views drawn, their values kept (`release` when no longer) */
  const held: Renderable<any>[] = [];
  /** the merged objects whose impostors leave out the raster regions until `finish` (`Renderable.rasterDrawn`) */
  const regionsDrawn: Renderable<any>[] = [];
  /**
   * the objects holding this pass's buffers (the resolve's, in their values): cleared with `destroy`, a value of a
   * destroyed context would make the object foreign to the next (`hasForeignResources`)
   */
  const bound = new Set<Renderable<any>>();
  const objects: SphereRasterObject[] = [];
  const renderablesOf: Renderable<any>[] = [];
  const cameraPlane = Vec4();

  /** the impostors of `r` draw all its spheres */
  function release(r: Renderable<any>) {
    const values = r.values as SpheresValues;
    ValueCell.updateIfChanged(values.uRasterMaxRadius, 0);
    ValueCell.updateIfChanged(values.uRasterChunkCount, 0);
  }
  function finish() {
    for (const r of regionsDrawn) r.rasterDrawn = false;
    regionsDrawn.length = 0;
    for (const r of drawn) {
      if (r.rasterDrawn !== undefined) held.push(r);
      else release(r);
    }
    drawn.length = 0;
  }
  /** the held objects not drawn by the raster now released */
  function releaseHeld() {
    for (const r of held) if (!drawn.includes(r)) release(r);
    held.length = 0;
  }

  // debug mode: why the spheres objects are drawn by the raster or not, logged where it changed; for merged objects
  // the instances the culler put into its raster regions and the impostors' (read back about once a second)
  const report: string[] = [];
  let reported = '';
  let countsAt = 0;
  function logReport(view: SphereRasterPassView) {
    const text = `sphere raster (radius ${view.maxRadius} px, near ${view.near.toFixed(2)}): ${report.join(', ')}`;
    if (text !== reported) console.log(text);
    reported = text;
    if (now() - countsAt < 1000) return;
    countsAt = now();
    for (const o of objects) {
      const input = o.culled;
      if (!input) continue;
      input.records.readAsync().then((data) => {
        const records = new Uint32Array(data);
        const count = (from: number, to: number) => {
          let n = 0;
          for (let k = from; k < to; ++k) n += records[k * input.recordWords + 1];
          return n;
        };
        const raster = count(input.recordBase, 2 * input.recordBase);
        console.log(`sphere raster #${o.item.id}: ${raster} instances raster, ${count(0, input.recordBase)} impostor`);
      });
    }
  }

  // compiled in the background: no raster until they are ready
  let preparing = false;
  function hasKernels() {
    if (compute.hasKernels(SphereRasterKernels)) return true;
    if (!preparing) {
      preparing = true;
      compute.prepareKernels(SphereRasterKernels).finally(() => (preparing = false));
    }
    return false;
  }

  return {
    prepare: (renderables, view) => {
      finish();
      objects.length = 0;
      renderablesOf.length = 0;
      report.length = 0;
      for (const r of bound) if (r.state.disposed) bound.delete(r);
      for (const r of renderables) {
        const values = r.values as SpheresValues;
        if (!('dSphereResolve' in values)) continue;
        const id = `#${r.renderItem.id}`;
        if (!r.state.visible || r.state.disposed || !r.state.opaque || !r.renderItem.renderFullscreen) {
          if (isDebugMode) report.push(`${id} hidden or not opaque`);
          continue;
        }
        if (!isSphereRasterEligible(values)) {
          if (isDebugMode) report.push(`${id} not eligible`);
          continue;
        }
        // merged objects only where culled on the GPU this frame (the CPU culling's lists are not the raster's input)
        const culled = r.getRasterInput?.(view.cull.frame, view.cull.cameraPlane, view.cull.frustum);
        if ((values as { dSegmented?: ValueCell<boolean> }).dSegmented?.ref.value && !culled) {
          if (isDebugMode) report.push(`${id} not culled on the GPU in this view`);
          continue;
        }
        if (view.culledOnly && !culled) {
          if (isDebugMode) report.push(`${id} not culled on the GPU (auto)`);
          continue;
        }
        objects.push({ item: r.renderItem, values, culled });
        renderablesOf.push(r);
      }
      if (objects.length === 0 || !hasKernels()) {
        releaseHeld();
        if (isDebugMode) {
          if (objects.length > 0) report.push('kernels compiling');
          logReport(view);
        }
        return;
      }

      const bytes = view.width * view.height * SphereRasterVisibilityBytes;
      if (!vis || vis.byteLength < bytes) {
        vis?.destroy();
        vis = compute.createBuffer(bytes, 'sphere raster visibility');
      }
      vis.clear(0, bytes);
      // the plane in the objects' space, its value at `model` times a point: the transpose of `model` times it
      const m = view.model;
      const p = view.cameraPlane;
      for (let i = 0; i < 4; ++i) {
        cameraPlane[i] = m[i * 4] * p[0] + m[i * 4 + 1] * p[1] + m[i * 4 + 2] * p[2] + m[i * 4 + 3] * p[3];
      }
      const table = raster.draw(objects, { ...view, cameraPlane }, vis, view.twoPhase);
      for (const range of table.ranges) {
        if (!range.drawn) continue;
        const r = renderablesOf[range.object];
        const values = r.values as SpheresValues;
        // merged objects: their impostors draw the culler's other regions whole (its split per instance)
        const culled = objects[range.object].culled;
        ValueCell.updateIfChanged(values.uRasterMaxRadius, culled ? 0 : view.maxRadius);
        ValueCell.updateIfChanged(values.uRasterChunkBase, range.chunkBase);
        ValueCell.updateIfChanged(values.uRasterChunkCount, range.chunkCapacity);
        ValueCell.updateIfChanged(values.uRasterWidth, view.width);
        if (!Mat4.areEqual(values.uRasterLevels.ref.value, range.levels, 0)) {
          ValueCell.update(values.uRasterLevels, Mat4.copy(values.uRasterLevels.ref.value, range.levels));
        }
        ValueCell.updateIfChanged(values.sRasterVisibility, vis);
        ValueCell.updateIfChanged(values.sRasterChunks, table.chunks);
        ValueCell.updateIfChanged(values.sRasterTransforms, range.transforms!);
        bound.add(r);
        drawn.push(r);
      }
      releaseHeld();
      if (isDebugMode) {
        for (const range of table.ranges) {
          report.push(`#${objects[range.object].item.id} ${range.drawn ? 'drawn' : 'no data drawn'}`);
        }
        logReport(view);
      }
    },
    resolve: (globals, scissor) => {
      let resolved = 0;
      // where the resolve's pipeline is not ready yet (compiling in the background) the impostors draw all spheres
      for (let i = drawn.length - 1; i >= 0; --i) {
        const r = drawn[i];
        if (scissor && !scissor(r)) continue;
        if (!r.renderItem.renderFullscreen!('color', 'dSphereResolve', globals)) {
          release(r);
          drawn.splice(i, 1);
          if (isDebugMode) console.log(`sphere raster #${r.renderItem.id}: resolve not ready`);
        } else {
          resolved += 1;
          if (
            !r.rasterDrawn &&
            r.rasterDrawn !== undefined &&
            (r.values as SpheresValues).uRasterChunkCount.ref.value > 0
          ) {
            r.rasterDrawn = true;
            regionsDrawn.push(r);
          }
        }
      }
      return resolved;
    },
    second: () => {
      if (drawn.length > 0) raster.second();
    },
    skipSecond: () => {
      if (drawn.length > 0) raster.skipSecond();
    },
    finish,
    isDrawn: (r) => drawn.includes(r),
    get active() {
      return drawn.length > 0;
    },
    destroy: () => {
      finish();
      for (const r of held) release(r);
      held.length = 0;
      for (const r of bound) {
        const values = r.values as SpheresValues;
        const none = undefined as unknown as StorageBuffer;
        ValueCell.updateIfChanged(values.sRasterVisibility, none);
        ValueCell.updateIfChanged(values.sRasterChunks, none);
        ValueCell.updateIfChanged(values.sRasterTransforms, none);
      }
      bound.clear();
      raster.destroy();
      vis?.destroy();
    },
  };
}
