/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
import { Renderer, type RendererProps } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { createMergedRenderObject, createRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 96;

/**
 * Instances of five spheres in a grid of three columns; `single`: of one sphere (the raster's rows of their own) in
 * fifteen, more rows than a workgroup's
 */
function createMember(offset: Vec3, color: number, single: boolean, extra: object = {}, radius = 2.5) {
  const centers = new Float32Array(single ? [0, 0, 0] : [0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, -1.2, 0, 0, 0, -1.2, 0]);
  const count = single ? 1 : 5;
  const spheres = Spheres.create(centers, new Float32Array(Array.from({ length: count }, (_, i) => i)), count);
  const columns = single ? 15 : 3;
  const instances = columns * 12;
  const transforms = new Float32Array(instances * 16);
  for (let i = 0; i < instances; ++i) {
    const t = Vec3.add(
      Vec3(),
      offset,
      Vec3.create((((i % columns) - (columns - 1) / 2) * 80) / columns, 0, -Math.floor(i / columns) * 8),
    );
    Mat4.toArray(Mat4.fromTranslation(Mat4(), t), transforms, i * 16);
  }
  const props = {
    ...PD.getDefaultValues(Spheres.Params),
    cellSize: 8,
    batchSize: 32,
    ...extra,
  } as PD.Values<Spheres.Params>;
  const transform = createTransform(
    transforms,
    instances,
    Sphere3D.create(Vec3(), radius),
    props.cellSize,
    props.batchSize,
  );
  const simple = BaseGeometry.createSimple(color as never, 1, transform);
  const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, simple.theme, props);
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

// the second phase's inputs by hand: a pyramid hiding all stands for a previous frame's depth that hid what is visible
// now, one at the far plane for the frame's depth so far hiding nothing
describe.skipIf(!WebGPU.module)('two-phase GPU culling on WebGPU', () => {
  it('draws in the second phase what the first left out, once, by impostors and by the raster', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const scene = Scene.create(graphics);
      scene.add(
        createMergedRenderObject(
          [
            createMember(Vec3.create(0, 0, 0), ColorNames.orange, false),
            createMember(Vec3.create(0, 4, -3), ColorNames.green, true),
          ],
          { gpuCulling: true },
        ),
      );
      scene.commit();

      // all instances well on screen and behind the near plane: tested for occlusion
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState({ ...camera.getFocus(Vec3.create(0, 2, -44), 90), radius: 120, radiusMax: 120 }, 0);
      camera.update();

      const levelSize = 64;
      const levelCount = Math.log2(levelSize) + 1;
      let version = 0;
      /** `depth` of each level's texel at `x` of its size, the levels side by side */
      const createPyramid = (depth: (x: number) => number) => {
        const texture = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
        const width = 2 * levelSize - 1;
        const array = new Float32Array(width * levelSize);
        for (let i = 0, offset = 0; i < levelCount; offset += levelSize >> i, ++i) {
          const size = levelSize >> i;
          for (let y = 0; y < size; ++y) {
            for (let x = 0; x < size; ++x) array[y * width + offset + x] = depth(x / size);
          }
        }
        texture.load({ array, width, height: levelSize });
        return {
          texture,
          get: (): OcclusionPyramid => ({
            texture,
            levelCount,
            minLevel: 0,
            width: levelSize,
            height: levelSize,
            view: camera.view,
            projection: camera.projection,
            near: camera.near,
            far: camera.far,
            version: ++version,
          }),
        };
      };
      const hiding = createPyramid(() => 1);
      const clear = createPyramid(() => 0);
      const half = createPyramid((x) => (x < 0.5 ? 1 : 0));
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      const read = async () => {
        graphics.flush();
        return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      };

      // the raster's input of its last view
      const renderable = scene.renderables[0];
      const getRasterInput = renderable.getRasterInput!;
      let input: ReturnType<typeof getRasterInput>;
      renderable.getRasterInput = (...args) => {
        input = getRasterInput(...args) ?? input;
        return input;
      };

      for (const props of [{ sphereRaster: 'off' }, { sphereRaster: 'auto', sphereRasterRadius: 8 }] as const) {
        if (props.sphereRaster !== 'off' && !graphics.capabilities.atomicMinMax64) continue;
        const renderer = Renderer.create(graphics, props as Partial<RendererProps>);
        renderer.setViewport(0, 0, Size, Size);
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        const first = (frame = createFrame()) => {
          renderer.update(camera, scene, frame);
          renderer.prepareCull(scene, true);
          renderer.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
          return frame;
        };
        const whole = async () => {
          renderer.setOcclusionPyramid(null);
          first();
          return read();
        };
        // the kernels and the resolve compile in the background
        let reference: number[] = [];
        for (let i = 0; i < 30; ++i) {
          reference = await whole();
          await graphics.resources.pending.whenReady();
          if (i > 5 && (await whole()).every((v, k) => v === reference[k])) break;
        }
        expect(new Set(reference).size).toBeGreaterThan(4);

        // no pyramid: the first phase left out nothing, no second
        expect(renderer.cullSecondPending).toBe(0);

        // all hidden by the previous frame's: nothing drawn yet
        renderer.setOcclusionPyramid(hiding.get);
        first();
        expect(renderer.cullSecondPending).toBe(1);
        const firstPart = await read();
        expect(firstPart).not.toEqual(reference);
        expect(new Set(firstPart).size).toBeLessThanOrEqual(2);
        // nothing hidden by the frame's own: all drawn by the second part, the image the same
        renderer.prepareCullSecond(scene, clear.get());
        graphics.beginRenderPass(target);
        renderer.renderOpaqueSecond(scene, camera);
        expect(await read()).toEqual(reference);

        // the left hidden by the previous frame's: the first part the right, the second the left after it
        await whole();
        renderer.setOcclusionPyramid(half.get);
        const frame = first();
        const rightOnly = await read();
        expect(new Set(rightOnly).size).toBeGreaterThan(2);
        expect(rightOnly).not.toEqual(reference);
        renderer.prepareCullSecond(scene, clear.get());
        graphics.beginRenderPass(target);
        renderer.renderOpaqueSecond(scene, camera);
        expect(await read()).toEqual(reference);
        if (props.sphereRaster !== 'off') {
          // each phase with instances of chunks and of single spheres in the raster's regions, the second's after the
          // first's: else the above does not test drawing the new ones only
          const all = new Uint32Array(await input!.records.readAsync());
          const second = new Uint32Array(await input!.secondRecords.readAsync());
          for (let s = 0; s < 2; ++s) {
            const o = (input!.recordBase + s) * input!.recordWords + 1;
            expect(second[o]).toBeGreaterThan(0);
            expect(all[o]).toBeGreaterThan(second[o]);
          }
        }

        // the frame's later views: all from the lists, no second phase again
        first(frame);
        expect(renderer.cullSecondPending).toBe(0);
        expect(await read()).toEqual(reference);

        // hidden by the frame's own too: stays out
        await whole();
        renderer.setOcclusionPyramid(hiding.get);
        first();
        renderer.prepareCullSecond(scene, hiding.get());
        graphics.beginRenderPass(target);
        renderer.renderOpaqueSecond(scene, camera);
        expect(await read()).toEqual(firstPart);

        // drawn by the first, not again by the second: only the first part draws
        await whole();
        renderer.setOcclusionPyramid(clear.get);
        first();
        const drawnFirst = await read();
        expect(drawnFirst).toEqual(reference);
        const item = scene.renderables[0].renderItem;
        renderer.prepareCullSecond(scene, clear.get());
        graphics.beginRenderPass(target);
        renderer.clear(true);
        const resolve = item.renderFullscreen;
        item.renderFullscreen = () => true;
        renderer.renderOpaqueSecond(scene, camera);
        item.renderFullscreen = resolve;
        expect(new Set(await read()).size).toBeLessThanOrEqual(2);
        // the chunk counts being read back: their buffer goes with the renderer
        await graphics.resources.pending.whenReady();
        renderer.dispose();
      }
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });

  // the second part draws no caps of its own: the first draws those of a mesh's segments whole, whatever is culled
  it('draws clipped objects with their solid interiors as without occlusion', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState({ ...camera.getFocus(Vec3.create(0, 0, -44), 45), radius: 120, radiusMax: 120 }, 0);
      camera.update();

      const levelSize = 64;
      const createPyramid = (depth: (x: number) => number) => {
        const texture = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
        const width = 2 * levelSize - 1;
        const array = new Float32Array(width * levelSize);
        for (let i = 0, offset = 0; levelSize >> i > 0; offset += levelSize >> i, ++i) {
          const size = levelSize >> i;
          for (let y = 0; y < size; ++y) {
            for (let x = 0; x < size; ++x) array[y * width + offset + x] = depth(x / size);
          }
        }
        texture.load({ array, width, height: levelSize });
        let version = 0;
        return (): OcclusionPyramid => ({
          texture,
          levelCount: Math.log2(levelSize) + 1,
          minLevel: 0,
          width: levelSize,
          height: levelSize,
          view: camera.view,
          projection: camera.projection,
          near: camera.near,
          far: camera.far,
          version: ++version,
        });
      };
      const half = createPyramid((x) => (x < 0.5 ? 1 : 0));
      const clear = createPyramid(() => 0);
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      const read = async () => {
        graphics.flush();
        return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      };

      const variants: [string, boolean, Partial<RendererProps>][] = [
        ['spheres', true, { sphereRaster: 'off' }],
        ['spheres', true, { sphereRaster: 'auto', sphereRasterRadius: 8 }],
        ['mesh', true, { sphereRaster: 'off' }],
        ['mesh', false, { sphereRaster: 'off' }],
      ];
      const images = new Map<string, number[]>();
      for (const [geometry, solidInterior, props] of variants) {
        if (props.sphereRaster !== 'off' && !graphics.capabilities.atomicMinMax64) continue;
        const clip = { clip: ClipPlane, solidInterior };
        const scene = Scene.create(graphics);
        scene.add(
          createMergedRenderObject(
            geometry === 'mesh'
              ? [createMeshMember(Vec3.create(0, 0, 0), ColorNames.orange, clip)]
              : [createMember(Vec3.create(0, 0, 0), ColorNames.orange, false, { ...clip, sizeFactor: 6 }, 8)],
            { gpuCulling: true },
          ),
        );
        scene.commit();
        const renderer = Renderer.create(graphics, props);
        renderer.setViewport(0, 0, Size, Size);
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        const first = () => {
          renderer.update(camera, scene, createFrame());
          renderer.prepareCull(scene, true);
          renderer.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
        };
        const whole = async () => {
          renderer.setOcclusionPyramid(null);
          first();
          return read();
        };
        let reference: number[] = [];
        for (let i = 0; i < 30; ++i) {
          reference = await whole();
          await graphics.resources.pending.whenReady();
          if (i > 5 && (await whole()).every((v, k) => v === reference[k])) break;
        }
        expect(new Set(reference).size).toBeGreaterThan(4);

        renderer.setOcclusionPyramid(half);
        first();
        expect(renderer.cullSecondPending).toBe(1);
        expect(await read()).not.toEqual(reference);
        renderer.prepareCullSecond(scene, clear());
        graphics.beginRenderPass(target);
        renderer.renderOpaqueSecond(scene, camera);
        expect(await read()).toEqual(reference);
        images.set(`${geometry} ${solidInterior}`, reference);
        // the chunk counts being read back: their buffer goes with the renderer
        await graphics.resources.pending.whenReady();
        renderer.dispose();
      }
      // the mesh's caps were drawn: else the above does not test them
      expect(images.get('mesh true')).not.toEqual(images.get('mesh false'));
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });
});

const ClipPlane = {
  variant: 'pixel',
  objects: [
    {
      type: 'plane',
      invert: false,
      position: Vec3.create(0, 0, -48),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 70 },
      scale: Vec3.create(1, 1, 1),
      transform: Mat4.identity(),
    },
  ],
};

/** The grid of three columns of `createMember`, of a mesh sphere */
function createMeshMember(offset: Vec3, color: number, extra: object) {
  const state = MeshBuilder.createState();
  addSphere(state, Vec3.origin, 6, 2);
  const mesh = MeshBuilder.getMesh(state);
  const instances = 3 * 12;
  const transforms = new Float32Array(instances * 16);
  for (let i = 0; i < instances; ++i) {
    const t = Vec3.add(Vec3(), offset, Vec3.create(((i % 3) - 1) * 26, 0, -Math.floor(i / 3) * 8));
    Mat4.toArray(Mat4.fromTranslation(Mat4(), t), transforms, i * 16);
  }
  const props = { ...PD.getDefaultValues(Mesh.Params), cellSize: 8, batchSize: 32, ...extra } as PD.Values<Mesh.Params>;
  const transform = createTransform(transforms, instances, Sphere3D.create(Vec3(), 8), props.cellSize, props.batchSize);
  const simple = BaseGeometry.createSimple(color as never, 1, transform);
  const values = Mesh.Utils.createValues(mesh, simple.transform, simple.locationIterator, simple.theme, props);
  return createRenderObject('mesh', values, Mesh.Utils.createRenderableState(props), -1);
}
