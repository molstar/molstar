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
import { ReversedZeroToOne } from '@molstar/graphics/gl/depth-range';
import type { SpheresValues } from '@molstar/graphics/gl/renderable/spheres';
import {
  createSphereRaster,
  decodeSphereRasterId,
  SphereRasterVisibilityBytes,
  type SphereRasterView,
} from '@molstar/graphics/gl/compute/sphere-raster';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
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

/** Spheres at several depths, as an instance at each of `translations` */
function createSpheres(translations: Vec3[], params: Partial<PD.Values<Spheres.Params>> = {}) {
  const count = 12;
  const centers = new Float32Array(count * 3);
  const groups = new Float32Array(count);
  for (let i = 0; i < count; ++i) {
    centers[i * 3] = ((i % 4) - 1.5) * 2.2;
    centers[i * 3 + 1] = (Math.floor(i / 4) - 1) * 2.2;
    centers[i * 3 + 2] = -i * 1.5;
    groups[i] = i;
  }
  const spheres = Spheres.create(centers, groups, count);
  const n = translations.length;
  const transforms = new Float32Array(n * 16);
  for (let i = 0; i < n; ++i) Mat4.toArray(Mat4.fromTranslation(Mat4(), translations[i]), transforms, i * 16);
  const props = { ...PD.getDefaultValues(Spheres.Params), sizeFactor: 2, ...params };
  const transform = createTransform(transforms, n, Sphere3D.create(Vec3(), 12), props.cellSize, props.batchSize);
  const simple = BaseGeometry.createSimple(ColorNames.orange as never, 1, transform);
  const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, simple.theme, props);
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

/**
 * The CPU's raster of `values` as the kernels draw it (no levels, a uniform size): per pixel the nearest sphere's
 * `slot * vertexCount + item + 1` and depth
 */
function referenceRaster(values: SpheresValues, view: SphereRasterView) {
  const ids = new Uint32Array(Size * Size);
  const depths = new Float32Array(Size * Size);
  const { modelView: mv, width, height, near, maxRadius } = view;
  // the depth of the WebGPU programs (reversed, `[0, 1]`)
  const P = ReversedZeroToOne.projection(Mat4(), view.projection);
  const positions = values.tPositionGroup.ref.value.array;
  const transforms = values.aTransform.ref.value;
  const vertexCount = values.uVertexCount.ref.value;
  const radius = values.uSize.ref.value * values.uSizeFactor.ref.value;
  const m = Mat4();
  const p = Vec3();
  let skipped = 0;
  for (let slot = 0; slot < values.instanceCount.ref.value; ++slot) {
    Mat4.mul(m, mv, Mat4.fromArray(Mat4(), transforms, slot * 16));
    for (let i = 0; i < vertexCount; ++i) {
      Vec3.transformMat4(p, Vec3.fromArray(p, positions, i * 4), m);
      const [x, y, z] = p;
      const r = radius;
      // `sphereBounds`, `spherePixelBounds`
      const pzr2 = z * z - r * r;
      const vx = Math.sqrt(x * x + pzr2);
      const minx = ((vx * x - r * z) / (vx * z + r * x)) * P[0];
      const maxx = ((vx * x + r * z) / (vx * z - r * x)) * P[0];
      const vy = Math.sqrt(y * y + pzr2);
      const miny = ((vy * y - r * z) / (vy * z + r * y)) * P[5];
      const maxy = ((vy * y + r * z) / (vy * z - r * y)) * P[5];
      const bx0 = ((-maxx - P[8]) * 0.5 + 0.5) * width;
      const bx1 = ((-minx - P[8]) * 0.5 + 0.5) * width;
      const by0 = ((miny + P[9]) * 0.5 + 0.5) * height;
      const by1 = ((maxy + P[9]) * 0.5 + 0.5) * height;
      if (!(-z - r > near && Math.max(bx1 - bx0, by1 - by0) * 0.5 <= maxRadius)) {
        skipped += 1;
        continue;
      }
      const id = slot * vertexCount + i + 1;
      for (let py = Math.max(0, Math.ceil(by0 - 0.5)); py <= Math.min(height - 1, Math.floor(by1 - 0.5)); ++py) {
        for (let px = Math.max(0, Math.ceil(bx0 - 0.5)); px <= Math.min(width - 1, Math.floor(bx1 - 0.5)); ++px) {
          const nx = ((px + 0.5) / width) * 2 - 1;
          const ny = 1 - ((py + 0.5) / height) * 2;
          const d = Vec3.normalize(Vec3(), Vec3.create((nx + P[8]) / P[0], (ny + P[9]) / P[5], -1));
          const B = Vec3.dot(d, p);
          const det = B * B + r * r - Vec3.dot(p, p);
          if (det < 0) continue;
          const t = B - Math.sqrt(det);
          const h = Vec3.scale(Vec3(), d, t);
          const clip = Vec4.transformMat4(Vec4(), Vec4.create(h[0], h[1], h[2], 1), P);
          const depth = clip[2] / clip[3];
          if (depth <= 0 || depth > 1) continue;
          const o = py * width + px;
          if (depth > depths[o] || (depth === depths[o] && id > ids[o])) {
            depths[o] = depth;
            ids[o] = id;
          }
        }
      }
    }
  }
  return { ids, depths, skipped };
}

// the kernels draw the spheres of a render item from its data pages: each pixel the nearest sphere's id and depth
describe.skipIf(!WebGPU.module)('sphere raster on WebGPU', () => {
  it('writes the ids and depths of the nearest spheres, those within the size limit and unoccluded only', async () => {
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
      const scene = Scene.create(graphics);
      const object = createSpheres([Vec3.create(-5, 0, 0), Vec3.create(5, 0, 0)]);
      scene.add(object);
      scene.commit();
      const values = object.values as SpheresValues;
      const item = scene.renderables[0].renderItem;

      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState(camera.getFocus(Vec3.create(0, 0, -8), 20), 0);
      camera.update();
      const cameraPlane = Plane3D.toArray(getCameraPlane(camera), Vec4(), 0);

      const raster = createSphereRaster(compute, graphics.resources);
      const vis = compute.createBuffer(Size * Size * SphereRasterVisibilityBytes);
      const draw = async (maxRadius: number, pyramid?: OcclusionPyramid) => {
        const view: SphereRasterView = {
          modelView: camera.view,
          projection: camera.projection,
          cameraPlane,
          modelScale: 1,
          near: camera.near,
          width: Size,
          height: Size,
          maxRadius,
          lod: LodView(),
          pyramid,
        };
        vis.clear();
        const table = raster.draw([{ item, values }], view, vis);
        const data = new Uint32Array(await vis.readAsync());
        const chunks = new Uint32Array(await table.chunks.readAsync());
        return { data, table, chunks, reference: referenceRaster(values, view) };
      };

      const check = async (maxRadius: number) => {
        const { data, table, chunks, reference } = await draw(maxRadius);
        const vertexCount = values.uVertexCount.ref.value;
        let covered = 0;
        let differ = 0;
        const drawn = new Set<number>();
        for (let o = 0; o < Size * Size; ++o) {
          const sphere = decodeSphereRasterId(table, chunks, data[o * 2]);
          const id = sphere ? sphere.slot * vertexCount + sphere.item + 1 : 0;
          if (reference.ids[o] !== 0) covered += 1;
          if (id !== reference.ids[o]) differ += 1;
          if (!sphere) continue;
          drawn.add(id);
          // the depth as `frag_depth`, within the rounding of the hit in f32 (the ray's parameter is a difference of
          // squares about the distance's)
          if (id === reference.ids[o]) {
            const depth = new Float32Array(new Uint32Array([data[o * 2 + 1]]).buffer)[0];
            expect(Math.abs(depth - reference.depths[o])).toBeLessThan(2e-4);
          }
          expect(sphere.object).toBe(0);
          expect(sphere.slot).toBeLessThan(2);
          expect(sphere.item).toBeLessThan(vertexCount);
          expect(sphere.level).toBe(0);
          expect(sphere.instance).toBe(values.aInstance.ref.value[sphere.slot]);
        }
        // the same pixels as the CPU's, but where float rounding moves an edge
        expect(covered).toBeGreaterThan(150);
        expect(differ).toBeLessThanOrEqual(covered * 0.02);
        return { drawn, skipped: reference.skipped };
      };

      // the spheres through the near plane left to the impostors
      const all = await check(1e6);
      expect(all.skipped).toBeGreaterThan(0);
      // and those above the limit, none of their pixels written
      const limited = await check(3);
      expect(limited.skipped).toBeGreaterThan(all.skipped);
      expect(limited.skipped).toBeLessThan(24);
      expect(limited.drawn.size).toBeLessThan(all.drawn.size);

      // occlusion of the instances against a pyramid of the view (`HiZPass.gpuPyramid`'s layout, its levels side by side
      // from 64 px to 1), all at one depth: at the near plane (reversed depth 1) every instance is behind it, at the far
      // plane (0) none; from farther away (an instance's sphere off screen in part is not tested, as on the CPU), the
      // pyramid's depths of a near plane at 1 (the test's sphere, enlarged, reaches the camera's)
      camera.setState(camera.getFocus(Vec3.create(0, 0, -8), 60), 0);
      camera.update();
      Plane3D.toArray(getCameraPlane(camera), cameraPlane, 0);
      const levelCount = Math.log2(Size) + 1;
      const texture = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
      let version = 0;
      const pyramid = (depth: number): OcclusionPyramid => {
        const width = 2 * Size - 1;
        texture.load({ array: new Float32Array(width * Size).fill(depth), width, height: Size });
        return {
          texture,
          levelCount,
          minLevel: 0,
          width: Size,
          height: Size,
          view: camera.view,
          projection: camera.projection,
          near: 1,
          far: camera.far,
          version: ++version,
        };
      };
      const visible = (data: Uint32Array) => data.filter((v, i) => i % 2 === 0 && v !== 0).length;
      const unoccluded = await draw(1e6);
      expect(visible(unoccluded.data)).toBeGreaterThan(20);
      expect(visible((await draw(1e6, pyramid(1))).data)).toBe(0);
      expect(Array.from((await draw(1e6, pyramid(0))).data)).toEqual(Array.from(unoccluded.data));
      texture.destroy();

      raster.destroy();
      vis.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });

  it('draws each instance at the levels its distance picks, as the impostors', async () => {
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

      // a near and a far instance side by side, the levels' boundary between them
      const translations = [Vec3.create(-12, 0, 10), Vec3.create(12, 0, -50)];
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      // near and far around the whole scene (`radiusMax`)
      camera.setState({ ...camera.getFocus(Vec3.create(0, 0, -20), 45), radiusMax: 70 }, 0);
      camera.update();
      const plane = getCameraPlane(camera);
      const boundary =
        (Plane3D.distanceToPoint(plane, translations[0]) + Plane3D.distanceToPoint(plane, translations[1])) / 2;
      const object = createSpheres(translations, {
        sizeFactor: 3,
        lodPixelSize: 0,
        lodLevels: [
          { minDistance: 0, maxDistance: boundary, overlap: 0, stride: 1, scaleBias: 3 },
          { minDistance: boundary, maxDistance: 1e6, overlap: 0, stride: 108, scaleBias: 3 },
        ],
      });
      const scene = Scene.create(graphics);
      scene.add(object);
      scene.commit();
      const values = object.values as SpheresValues;
      const item = scene.renderables[0].renderItem;
      const levels = values.lodLevels.ref.value as unknown as number[][];
      expect(levels.length).toBe(2);
      // the coarse level's spheres the first of the level order
      const coarse = levels[1][3] / 6;
      expect(coarse).toBeLessThan(values.uVertexCount.ref.value);

      const raster = createSphereRaster(compute, graphics.resources);
      const vis = compute.createBuffer(Size * Size * SphereRasterVisibilityBytes);
      vis.clear();
      const table = raster.draw(
        [{ item, values }],
        {
          modelView: camera.view,
          projection: camera.projection,
          cameraPlane: Plane3D.toArray(plane, Vec4(), 0),
          modelScale: 1,
          near: camera.near,
          width: Size,
          height: Size,
          maxRadius: 1e6,
          lod: LodView(),
        },
        vis,
      );
      const data = new Uint32Array(await vis.readAsync());
      const chunks = new Uint32Array(await table.chunks.readAsync());
      const slots = [0, 0];
      for (let o = 0; o < Size * Size; ++o) {
        const sphere = decodeSphereRasterId(table, chunks, data[o * 2]);
        if (!sphere) continue;
        slots[sphere.slot] += 1;
        // the instance's spheres at its level: the far one's at the coarse level, the first of the level order
        expect(sphere.level).toBe(sphere.slot);
        if (sphere.level === 1) expect(sphere.item).toBeLessThan(coarse);
      }
      expect(slots[0]).toBeGreaterThan(0);
      expect(slots[1]).toBeGreaterThan(0);

      raster.destroy();
      vis.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });

  it('draws merged objects culled on the GPU as their members unmerged', async () => {
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
      // without levels, and with two by distance faded into each other, the second member at their boundary (drawn at
      // both)
      const boundary = Plane3D.distanceToPoint(plane, Vec3.create(0, 6, -15));
      const levels = [
        { minDistance: 0, maxDistance: boundary + 4, overlap: 4, stride: 1, scaleBias: 3 },
        { minDistance: boundary - 4, maxDistance: 1e6, overlap: 4, stride: 32, scaleBias: 3 },
      ];
      // (the stride divided by the size factor cubed: 4)
      for (const params of [{}, { lodPixelSize: 0, lodLevels: levels }]) {
        const members = [
          createSpheres([Vec3.create(-6, 0, 10), Vec3.create(6, 0, -40)], params),
          createSpheres([Vec3.create(0, 6, -15)], params),
        ];
        const mergedScene = Scene.create(graphics);
        mergedScene.add(createMergedRenderObject(members, { gpuCulling: true }));
        mergedScene.commit();
        const memberScene = Scene.create(graphics);
        for (const m of members) memberScene.add(m);
        memberScene.commit();

        // the culling kernels compile in the background: culled on the CPU until they are ready
        // every instance into the culler's raster regions (its split per instance at any size)
        const renderer = Renderer.create(graphics, { sphereRaster: 'only' });
        renderer.setViewport(0, 0, Size, Size);
        const renderable = mergedScene.renderables[0];
        // the view the renderer culls with
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
        expect(culled!.segmentCount).toBe(2);
        expect(culled!.levelCount).toBe('lodLevels' in params ? 2 : 0);

        const depths = async (objects: Parameters<typeof raster.draw>[0]) => {
          vis.clear();
          const table = raster.draw(objects, view, vis);
          const data = new Uint32Array(await vis.readAsync());
          const chunks = new Uint32Array(await table.chunks.readAsync());
          const out = new Uint32Array(Size * Size);
          for (let o = 0; o < Size * Size; ++o) {
            out[o] = data[o * 2 + 1];
            // the merged object's spheres decode to its compacted instances
            const sphere = decodeSphereRasterId(table, chunks, data[o * 2]);
            if (sphere) expect(sphere.object).toBeLessThan(objects.length);
          }
          return out;
        };
        const merged = await depths([
          { item: renderable.renderItem, values: renderable.values as unknown as SpheresValues, culled },
        ]);
        const unmerged = await depths(
          memberScene.renderables.map((r) => ({ item: r.renderItem, values: r.values as unknown as SpheresValues })),
        );
        // the same spheres at the same depths
        expect(merged.filter((d) => d !== 0).length).toBeGreaterThan(40);
        expect(Array.from(merged)).toEqual(Array.from(unmerged));
      }

      raster.destroy();
      vis.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });
});
