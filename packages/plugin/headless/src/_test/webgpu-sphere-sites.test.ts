/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ColorNames } from '@molstar/core/util/color/names';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { createMergedRenderObject, createRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame, LodView } from '@molstar/graphics/gl/renderable';
import type { SpheresValues } from '@molstar/graphics/gl/renderable/spheres';
import type { LodLevels } from '@molstar/graphics/geo/geometry/lod-levels';
import { setLodLevelsByDistance } from './lod-levels-by-distance.js';
import {
  createSiteTypeSize,
  createTextureSize,
  type SizeData,
  sizeDataFactor,
} from '@molstar/graphics/geo/geometry/size-data';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import { packIntToRGBArray } from '@molstar/core/util/number-packing';
import {
  createSphereRaster,
  decodeSphereRasterId,
  SphereRasterVisibilityBytes,
  type SphereRasterView,
} from '@molstar/graphics/gl/compute/sphere-raster';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 64;

/** The camera plane the renderer culls and levels by (`uCameraPlane`): the near plane's normal, through the eye */
function getCameraPlane(camera: Camera) {
  const frustum = Frustum3D.fromProjectionMatrix(Frustum3D(), camera.projectionView);
  const plane = Plane3D.copy(Plane3D(), frustum[Frustum3D.PlaneIndex.Near]);
  plane.constant -= Plane3D.distanceToPoint(plane, camera.position);
  return plane;
}

const Count = 12;
const Origin = Vec3.create(0.5, -0.25, 1);
const Scale = 0.25;

/** Sphere `i`'s center: multiples of 0.25 from the origin, exact as floats and as quantized */
function center(i: number) {
  return Vec3.create(((i % 4) - 1.5) * 2.5, (Math.floor(i / 4) - 1) * 2.5, -i * 1.5);
}

/** Sizes by a sphere's type (`i % 3`, the site's uint16): the table's */
const TypeSizes = [0.6, 0.8, 1];

/**
 * Spheres at several depths as an instance at each of `translations`: their centers as floats or quantized (`sites`,
 * the same centers exactly); `sized` by type: a table by the sites' types (`siteType`), the floats' by group
 */
function createSpheres(translations: Vec3[], sites: boolean, levels?: LodLevels, sized = false) {
  let spheres: Spheres;
  if (sites) {
    const data = new Int16Array(Count * 4);
    for (let i = 0; i < Count; ++i) {
      const c = Vec3.sub(Vec3(), center(i), Origin);
      data.set([c[0] / Scale, c[1] / Scale, c[2] / Scale, i % 3], i * 4);
    }
    spheres = Spheres.fromSites({ data, origin: Origin, scale: Scale }, Count);
  } else {
    const centers = new Float32Array(Count * 3);
    const groups = new Float32Array(Count);
    for (let i = 0; i < Count; ++i) {
      Vec3.toArray(center(i), centers, i * 3);
      groups[i] = i;
    }
    spheres = Spheres.create(centers, groups, Count);
  }
  const n = translations.length;
  const transforms = new Float32Array(n * 16);
  for (let i = 0; i < n; ++i) Mat4.toArray(Mat4.fromTranslation(Mat4(), translations[i]), transforms, i * 16);
  const props = { ...PD.getDefaultValues(Spheres.Params), sizeFactor: 2 };
  const transform = createTransform(transforms, n, Sphere3D.create(Vec3(), 12), props.cellSize, props.batchSize);
  const simple = BaseGeometry.createSimple(ColorNames.orange as never, 1, transform);
  const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, simple.theme, props);
  if (levels) setLodLevelsByDistance(values, levels, props.sizeFactor);
  if (sized) {
    const sizeData = values as unknown as SizeData;
    if (sites) {
      createSiteTypeSize(TypeSizes, sizeData);
    } else {
      const table = createDataArray(Count, 4, Uint8Array);
      for (let i = 0; i < Count; ++i) packIntToRGBArray(TypeSizes[i % 3] * sizeDataFactor, table.array, i * 4);
      createTextureSize(table, 'group', sizeData);
    }
  }
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

// quantized centers (`SphereSites`, `dSites`): the raster draws them as their float twins, unmerged (the object's origin
// and scale) and merged (the segment's, `CulledData`)
describe.skipIf(!WebGPU.module)('sphere sites on WebGPU', () => {
  it('draw quantized centers as their float twins, unmerged and merged, with levels and without, sized by type', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      if (!graphics.capabilities.atomicMinMax64) return;
      const compute = graphics.compute!;
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState({ ...camera.getFocus(Vec3.create(0, 2, -15), 40), radiusMax: 40 }, 0);
      camera.update();
      const plane = getCameraPlane(camera);
      const view: SphereRasterView = {
        modelView: camera.view,
        projection: camera.projection,
        cameraPlane: Plane3D.toArray(plane, Vec4(), 0),
        modelScale: 1,
        near: camera.near,
        width: Size,
        height: Size,
        maxRadius: 1e6,
        lod: LodView(),
      };
      const raster = createSphereRaster(compute, graphics.resources);
      const vis = compute.createBuffer(Size * Size * SphereRasterVisibilityBytes);
      const depths = async (objects: Parameters<typeof raster.draw>[0]) => {
        vis.clear();
        const table = raster.draw(objects, view, vis);
        // both read as drawn: the chunks a scratch buffer, let go once idle (a wait under load past it)
        const [visData, chunkData] = await Promise.all([vis.readAsync(), table.chunks.readAsync()]);
        const data = new Uint32Array(visData);
        const chunks = new Uint32Array(chunkData);
        const out = new Uint32Array(Size * Size);
        const items = new Int32Array(Size * Size).fill(-1);
        for (let o = 0; o < Size * Size; ++o) {
          out[o] = data[o * 2 + 1];
          const sphere = decodeSphereRasterId(table, chunks, data[o * 2]);
          if (sphere) items[o] = sphere.item;
        }
        return { depths: out, items };
      };

      const boundary = Plane3D.distanceToPoint(plane, Vec3.create(0, 6, -15));
      const levels = [
        { minDistance: 0, maxDistance: boundary, stride: 1, scaleBias: 3 },
        { minDistance: boundary, maxDistance: 1e6, stride: 32, scaleBias: 3 },
      ];
      const translations = [[Vec3.create(-6, 0, 10), Vec3.create(6, 0, -40)], [Vec3.create(0, 6, -15)]];
      for (const [lodLevels, sized] of [
        [undefined, false],
        [levels, false],
        [undefined, true],
        [levels, true],
      ] as const) {
        const floats = translations.map((t) => createSpheres(t, false, lodLevels, sized));
        const sites = translations.map((t) => createSpheres(t, true, lodLevels, sized));
        if (sized) expect((sites[0].values as SpheresValues).sizeType.ref.value).toBe('siteType');
        expect((sites[0].values as SpheresValues).dSites.ref.value).toBe(true);
        const scene = (objects: typeof floats) => {
          const s = Scene.create(graphics);
          for (const o of objects) s.add(o);
          s.commit();
          return s;
        };
        const unmerged = (s: Scene) =>
          depths(s.renderables.map((r) => ({ item: r.renderItem, values: r.values as unknown as SpheresValues })));

        // unmerged: the object's origin and scale
        const floatDepths = await unmerged(scene(floats));
        const siteDepths = await unmerged(scene(sites));
        expect(floatDepths.depths.filter((d) => d !== 0).length).toBeGreaterThan(40);
        expect(Array.from(siteDepths.depths)).toEqual(Array.from(floatDepths.depths));
        expect(Array.from(siteDepths.items)).toEqual(Array.from(floatDepths.items));

        // merged and culled on the GPU: each member's origin and scale from the culler's segment data
        const mergedScene = scene([createMergedRenderObject(sites, { gpuCulling: true })]);
        const renderer = Renderer.create(graphics, { sphereRaster: 'only' });
        renderer.setViewport(0, 0, Size, Size);
        const renderable = mergedScene.renderables[0];
        const cullFrustum = Frustum3D.fromProjectionMatrix(Frustum3D(), camera.projectionView);
        let culled;
        for (let i = 0; i < 50 && !culled; ++i) {
          const frame = createFrame();
          renderer.update(camera, mergedScene, frame);
          renderer.prepareCull(mergedScene);
          culled = renderable.getRasterInput?.(frame, plane, cullFrustum);
          await graphics.resources.pending.whenReady();
        }
        expect(culled).toBeDefined();
        const merged = await depths([
          { item: renderable.renderItem, values: renderable.values as unknown as SpheresValues, culled },
        ]);
        expect(Array.from(merged.depths)).toEqual(Array.from(floatDepths.depths));
      }

      raster.destroy();
      vis.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });
});
