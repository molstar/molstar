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
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Renderer, type RendererProps } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { createMergedRenderObject, createRenderObject, hasForeignResources } from '@molstar/graphics/gl/render-object';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 96;

/**
 * Five spheres (`stride`: at a level keeping every so many, 5 one of them, as far instances draw whole units), as
 * instances in three columns receding from the camera, the grid's origin at `offset`
 */
function createMember(offset: Vec3, color: number, stride?: number) {
  const centers = new Float32Array([0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, -1.2, 0, 0, 0, -1.2, 0]);
  const spheres = Spheres.create(centers, new Float32Array([0, 1, 2, 3, 4]), 5);
  const count = 3 * 12;
  const transforms = new Float32Array(count * 16);
  for (let i = 0; i < count; ++i) {
    const t = Vec3.add(Vec3(), offset, Vec3.create(((i % 3) - 1) * 6, 0, -Math.floor(i / 3) * 8));
    Mat4.toArray(Mat4.fromTranslation(Mat4(), t), transforms, i * 16);
  }
  const props = {
    ...PD.getDefaultValues(Spheres.Params),
    cellSize: 8,
    batchSize: 32,
    ...(stride && { lodLevels: [{ minDistance: 0, maxDistance: 1e7, overlap: 0, stride, scaleBias: 1 }] }),
  };
  const transform = createTransform(transforms, count, Sphere3D.create(Vec3(), 2.5), props.cellSize, props.batchSize);
  const simple = BaseGeometry.createSimple(color as never, 1, transform);
  const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, simple.theme, props);
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

// merged objects culled on the GPU: the culler puts the instances the raster draws into regions of their own, the
// impostors draw the others (`sphere-compute-raster.md`, the split per instance)
describe.skipIf(!WebGPU.module)('sphere raster split per instance on WebGPU', () => {
  it('draws the near instances as impostors, the far ones by the raster, as the impostors draw all', async () => {
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
      const members = [
        createMember(Vec3.create(0, 0, 0), ColorNames.orange),
        createMember(Vec3.create(0, 4, -3), ColorNames.green),
      ];
      const scene = Scene.create(graphics);
      const merged = createMergedRenderObject(members, { gpuCulling: true });
      scene.add(merged);
      scene.commit();
      const renderable = scene.renderables[0];

      // along the columns, from just in front of the first instances
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState(
        { ...camera.getFocus(Vec3.create(0, 2, -40), 50), position: Vec3.create(0, 6, 14), radiusMax: 120 },
        0,
      );
      camera.update();

      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      const render = async (props: Partial<RendererProps>) => {
        const renderer = Renderer.create(graphics, props);
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        let resolved = 0;
        const item = renderable.renderItem;
        const renderFullscreen = item.renderFullscreen!;
        item.renderFullscreen = (...args) => {
          const drawn = renderFullscreen(...args);
          if (drawn) resolved += 1;
          return drawn;
        };
        const frame = async () => {
          resolved = 0;
          renderer.setViewport(0, 0, Size, Size);
          renderer.update(camera, scene, createFrame());
          renderer.prepareCull(scene);
          renderer.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
          graphics.flush();
          return Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
        };
        let pixels: number[] = [];
        // the culling and raster kernels and the resolve compile in the background
        for (let i = 0; i < 30; ++i) {
          pixels = await frame();
          await graphics.resources.pending.whenReady();
          if (props.sphereRaster === 'off' ? i > 3 : resolved > 0 && i > 3) break;
        }
        // the same view again: the opaque pass replays its recording (`drawLists`), with the raster on too
        for (let i = 0; i < 4; ++i) await frame();
        const { drawLists } = graphics.stats;
        const replayed = drawLists.replayed;
        const direct = drawLists.chunksDirect;
        const recorded = drawLists.recorded;
        expect(await frame()).toEqual(pixels);
        expect(drawLists.replayed).toBe(replayed + 1);
        expect(drawLists.chunksDirect).toBe(direct);
        expect(drawLists.recorded).toBe(recorded);
        if (props.sphereRaster !== 'off') expect(resolved).toBe(1);
        item.renderFullscreen = renderFullscreen;
        // a context without compute (a switch to WebGL): the raster's buffers leave the object with its renderer, else
        // it counts as holding resources of the previous context and is dropped from the scene for good
        const noCompute = () => ({ resources: { ownsTexture: () => true } }) as unknown as GraphicsContext;
        if (props.sphereRaster !== 'off') expect(hasForeignResources(merged, noCompute())).toBe(true);
        renderer.dispose();
        expect(hasForeignResources(merged, noCompute())).toBe(false);
        return { pixels, resolved };
      };

      // the raster's input of its last view
      const getRasterInput = renderable.getRasterInput!;
      let input: ReturnType<typeof getRasterInput>;
      renderable.getRasterInput = (...args) => {
        input = getRasterInput(...args) ?? input;
        return input;
      };
      const impostors = await render({ sphereRaster: 'off' });
      expect(new Set(impostors.pixels).size).toBeGreaterThan(4);
      const split = await render({ sphereRaster: 'on', sphereRasterRadius: 2 });
      expect(split.resolved).toBeGreaterThan(0);
      // `auto` draws objects culled on the GPU as `on`
      const auto = await render({ sphereRaster: 'auto', sphereRasterRadius: 2 });
      expect(auto.resolved).toBeGreaterThan(0);
      expect(auto.pixels).toEqual(split.pixels);

      // both kinds of region hold instances: the impostors' (near) and the raster's (far)
      renderable.getRasterInput = getRasterInput;
      expect(input).toBeDefined();
      const records = new Uint32Array(await input!.records.readAsync());
      const counts = (from: number, to: number) => {
        let n = 0;
        for (let r = from; r < to; ++r) n += records[r * input!.recordWords + 1];
        return n;
      };
      expect(counts(0, input!.recordBase)).toBeGreaterThan(0);
      expect(counts(input!.recordBase, 2 * input!.recordBase)).toBeGreaterThan(0);

      // the same image, but where the raster's and the fragment stage's arithmetic decide an edge differently
      let differ = 0;
      for (let i = 0; i < impostors.pixels.length; i += 4) {
        const d = Math.max(
          Math.abs(impostors.pixels[i] - split.pixels[i]),
          Math.abs(impostors.pixels[i + 1] - split.pixels[i + 1]),
          Math.abs(impostors.pixels[i + 2] - split.pixels[i + 2]),
        );
        if (d > 8) differ += 1;
      }
      expect(differ).toBeLessThanOrEqual(Size * Size * 0.02);

      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });

  // the chunks of a merged object follow what its views use (`readChunkCounts`): more than it starts with here, 3,000
  // spheres an instance by 4,000 instances, lost in a frame or so, then all drawn
  it('grows its chunks where a view needs more than it started with', async () => {
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
      const n = 3000;
      const centers = new Float32Array(n * 3);
      for (let i = 0; i < n * 3; ++i) centers[i] = (((i * 2654435761) >>> 0) / 4294967296 - 0.5) * 6;
      const spheres = Spheres.create(centers, new Float32Array(n), n);
      const count = 4000;
      const transforms = new Float32Array(count * 16);
      for (let i = 0; i < count; ++i) {
        const t = Vec3.create((i % 63) * 8 - 248, (Math.floor(i / 63) % 64) * 8 - 252, 0);
        Mat4.toArray(Mat4.fromTranslation(Mat4(), t), transforms, i * 16);
      }
      const props = { ...PD.getDefaultValues(Spheres.Params), cellSize: 16, batchSize: 64, sizeFactor: 2 };
      const transform = createTransform(transforms, count, Sphere3D.create(Vec3(), 6), props.cellSize, props.batchSize);
      const simple = BaseGeometry.createSimple(ColorNames.orange as never, 1, transform);
      const values = Spheres.Utils.createValues(
        spheres,
        simple.transform,
        simple.locationIterator,
        simple.theme,
        props,
      );
      const scene = Scene.create(graphics);
      scene.add(
        createMergedRenderObject(
          [createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1)],
          { gpuCulling: true },
        ),
      );
      scene.commit();
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState({ ...camera.getFocus(Vec3(), 300), position: Vec3.create(0, 0, 600), radiusMax: 1000 }, 0);
      camera.update();
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      const render = async (sphereRaster: RendererProps['sphereRaster']) => {
        const renderer = Renderer.create(graphics, { sphereRaster });
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        const item = scene.renderables[0].renderItem;
        const renderFullscreen = item.renderFullscreen!;
        let resolved = 0;
        item.renderFullscreen = (...args) => {
          const drawn = renderFullscreen(...args);
          if (drawn) resolved += 1;
          return drawn;
        };
        const images: number[][] = [];
        for (let i = 0; i < 40; ++i) {
          resolved = 0;
          renderer.setViewport(0, 0, Size, Size);
          renderer.update(camera, scene, createFrame());
          renderer.prepareCull(scene);
          renderer.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
          graphics.flush();
          images.push(Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size)));
          await graphics.resources.pending.whenReady();
          if (sphereRaster === 'off' ? i > 3 : resolved > 0 && i > 8) break;
        }
        if (sphereRaster !== 'off') expect(resolved).toBe(1);
        item.renderFullscreen = renderFullscreen;
        renderer.dispose();
        return images[images.length - 1];
      };
      const impostors = await render('off');
      expect(new Set(impostors).size).toBeGreaterThan(4);
      const raster = await render('on');
      let differ = 0;
      for (let i = 0; i < impostors.length; i += 4) {
        const d = Math.max(
          Math.abs(impostors[i] - raster[i]),
          Math.abs(impostors[i + 1] - raster[i + 1]),
          Math.abs(impostors[i + 2] - raster[i + 2]),
        );
        if (d > 8) differ += 1;
      }
      expect(differ).toBeLessThanOrEqual(Size * Size * 0.02);
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });

  // instances of one sphere take rows of their own, a thread each (`appendChunks`, `uSingles`), beside the chunks of
  // the members with more
  it('draws instances of one sphere as the impostors, beside instances of several', async () => {
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
      const scene = Scene.create(graphics);
      scene.add(
        createMergedRenderObject(
          [
            createMember(Vec3.create(0, 0, 0), ColorNames.orange, 5),
            createMember(Vec3.create(0, 4, -3), ColorNames.green, 1),
            createMember(Vec3.create(0, -4, -3), ColorNames.blue, 5),
          ],
          { gpuCulling: true },
        ),
      );
      scene.commit();
      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      camera.setState(
        // back from all instances, each within the raster's limit (none left to the impostors)
        { ...camera.getFocus(Vec3.create(0, 0, -40), 50), position: Vec3.create(0, 20, 45), radiusMax: 200 },
        0,
      );
      camera.update();
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      const render = async (props: Partial<RendererProps>) => {
        const renderer = Renderer.create(graphics, props);
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        const item = scene.renderables[0].renderItem;
        const renderFullscreen = item.renderFullscreen!;
        let resolved = 0;
        item.renderFullscreen = (...args) => {
          const drawn = renderFullscreen(...args);
          if (drawn) resolved += 1;
          return drawn;
        };
        let pixels: number[] = [];
        // the kernels and the resolve compile in the background
        for (let i = 0; i < 40; ++i) {
          resolved = 0;
          renderer.setViewport(0, 0, Size, Size);
          renderer.update(camera, scene, createFrame());
          renderer.prepareCull(scene);
          renderer.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          renderer.clear(true);
          renderer.renderOpaque(scene, camera);
          graphics.flush();
          pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
          await graphics.resources.pending.whenReady();
          if (i > 3 && (props.sphereRaster === 'off' || resolved > 0)) break;
        }
        if (props.sphereRaster !== 'off') expect(resolved).toBe(1);
        item.renderFullscreen = renderFullscreen;
        renderer.dispose();
        return pixels;
      };
      const impostors = await render({ sphereRaster: 'off' });
      expect(new Set(impostors).size).toBeGreaterThan(4);
      const raster = await render({ sphereRaster: 'auto', sphereRasterRadius: 64 });
      let differ = 0;
      for (let i = 0; i < impostors.length; i += 4) {
        const d = Math.max(
          Math.abs(impostors[i] - raster[i]),
          Math.abs(impostors[i + 1] - raster[i + 1]),
          Math.abs(impostors[i + 2] - raster[i + 2]),
        );
        if (d > 8) differ += 1;
      }
      expect(differ).toBeLessThanOrEqual(Size * Size * 0.02);
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  });
});
