/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterAll, describe, expect, it, vi } from 'vitest';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4, Quat, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { now } from '@molstar/core/util/now';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { Clip } from '@molstar/core/util/clip';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import {
  createPosesTransform,
  createTransform,
  moveTransformData,
  type TransformData,
} from '@molstar/graphics/geo/geometry/transform-data';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import type { MultiDrawBaseData } from '@molstar/graphics/gl/render-item';
import {
  createMergedRenderObject,
  createRenderObject,
  type GraphicsRenderObject,
  isMergedRenderObject,
} from '@molstar/graphics/gl/render-object';
import {
  createMergedCuller,
  MergedCullerCounts,
  MergedCullerLimits,
  MergedCullKernels,
} from '@molstar/graphics/gl/renderable/merged-culler';
import { createMergedValues } from '@molstar/graphics/gl/renderable/merged';
import { moveInstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
// contexts of their own (`own`): two at once compared (what one compiles, the other does not)
import { WebGPU } from './webgpu-module.js';

const Size = 64;
const Copies = 6;

/**
 * Five spheres, as `copies` × `copies` instances on a grid in the xy plane, the grid's origin at `offset`; `params`
 * besides the defaults (clip objects)
 */
function createMember(offset: Vec3, color: number, params: Partial<PD.Values<Spheres.Params>> = {}) {
  const centers = new Float32Array([0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, -1.2, 0, 0, 0, -1.2, 0]);
  const groups = new Float32Array([0, 1, 2, 3, 4]);
  const spheres = Spheres.create(centers, groups, 5);
  const transforms = new Float32Array(Copies * Copies * 16);
  for (let i = 0; i < Copies * Copies; ++i) {
    const t = Vec3.add(Vec3(), offset, Vec3.create((i % Copies) * 6, Math.floor(i / Copies) * 6, 0));
    Mat4.toArray(Mat4.fromTranslation(Mat4(), t), transforms, i * 16);
  }
  // cells of a few instances each (the default 200 Å holds a member's whole grid)
  const props = { ...PD.getDefaultValues(Spheres.Params), cellSize: 8, batchSize: 32, ...params };
  // with the invariant sphere (the five spheres') the instances get a grid
  const sphere = Sphere3D.create(Vec3(), 2.5);
  const transform = createTransform(transforms, Copies * Copies, sphere, props.cellSize, props.batchSize);
  const simple = BaseGeometry.createSimple(color as never, 1, transform);
  const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, simple.theme, props);
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

describe.skipIf(!WebGPU.module)('GPU culling of merged render objects', () => {
  it("writes the identity region in an order set, kept through its regions' layouts", async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const members = [createMember(Vec3.create(0, 0, 0), ColorNames.orange), createMember(Vec3(), ColorNames.blue)];
    const values = members.map((m) => m.values as never);
    const culler = createMergedCuller(
      graphics.compute!,
      graphics.resources,
      createMergedValues('spheres', values[0], values, true),
      values,
      true,
    );
    const n = members.length * Copies * Copies;
    const region = async () => {
      const buffer = culler.values.sCulledIndices.ref.value;
      return Array.from(new Uint32Array(await buffer.readAsync(culler.identityBase * 4, n * 4)));
    };
    const own = Array.from({ length: n }, (_, i) => i);
    expect(await region()).toEqual(own);
    // each member's instances reversed (a grid's order of instances written in another)
    const order = Uint32Array.from(own, (i) => (i < n / 2 ? n / 2 - 1 - i : n - 1 - (i - n / 2)));
    culler.setIdentityOrder(order);
    expect(await region()).toEqual(Array.from(order));
    // laid out again (the members' grids changed: the cells and regions anew): the order kept
    const layouts = MergedCullerCounts.layouts;
    ValueCell.update(members[1].values.instanceGrid, { ...members[1].values.instanceGrid.ref.value });
    culler.sync();
    expect(MergedCullerCounts.layouts).toBe(layouts + 1);
    expect(await region()).toEqual(Array.from(order));
    culler.setIdentityOrder(undefined);
    expect(await region()).toEqual(own);
    culler.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  it('draws the cells the CPU culling draws, from indirect lists', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const members = [
      createMember(Vec3.create(0, 0, 0), ColorNames.orange),
      createMember(Vec3.create(2.5, 2.5, 1), ColorNames.green),
      createMember(Vec3.create(-2.5, 3, -1), ColorNames.blue),
    ];
    const cpu = createMergedRenderObject(members);
    const gpu = createMergedRenderObject(members, { gpuCulling: true });

    const createScene = (o: GraphicsRenderObject) => {
      const scene = Scene.create(graphics);
      scene.add(o);
      scene.commit();
      return scene;
    };
    const cpuScene = createScene(cpu);
    const gpuScene = createScene(gpu);

    // a part of the grid in view: the frustum culls instances
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = cpuScene.boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * 0.3, 0, 0)), radius * 0.35), 0);
    camera.update();

    // the GPU culling's lists as the impostors draw them: without the sphere raster's regions (`sphereRaster` 'auto')
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    for (const scene of [cpuScene, gpuScene]) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
      );
    }

    // the lists each object draws with
    const drawn = new Map<Scene, MultiDrawBaseData[]>();
    for (const scene of [cpuScene, gpuScene]) {
      const item = scene.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, list) => {
        if (list) drawn.set(scene, list);
        render(variant, globals, list);
      };
    }

    const frame = async (scene: Scene, cullOnGpu = true) => {
      // GPU culling dispatches before the passes, with the planes of `update`
      renderer.update(camera, scene, createFrame());
      if (cullOnGpu) renderer.prepareCull(scene);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(scene, camera);
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    // the culling kernels compile in the background: the CPU culling's lists until they are ready
    for (let i = 0; i < 20 && !drawn.get(gpuScene)?.[0].indirect; ++i) await frame(gpuScene);
    const gpuImage = await frame(gpuScene);
    expect(drawn.get(gpuScene)![0].indirect).toBeDefined();
    const cpuImage = await frame(cpuScene);
    const cpuList = drawn.get(cpuScene)![0];
    expect(cpuList.indirect).toBeUndefined();
    let cpuInstances = 0;
    for (let i = 0; i < cpuList.count; ++i) cpuInstances += cpuList.instanceCounts[i];
    // some cells culled, some drawn
    const total = members.length * Copies * Copies;
    expect(cpuInstances).toBeGreaterThan(0);
    expect(cpuInstances).toBeLessThan(total);
    // the GPU at most the CPU's: its cells split those of the grid (finer frustum tests), the same image; the instance
    // counts of the argument records
    const gpuList = drawn.get(gpuScene)![0];
    const { buffer, byteOffset, indexed } = gpuList.indirect!;
    // sphere quads draw indexed by the shared quad indices: indexed records (5 words)
    expect(indexed).toBe(true);
    const words = indexed ? 5 : 4;
    const listInstances = async () => {
      const list = drawn.get(gpuScene)![0];
      const { records } = list.indirect!;
      const args = new Uint32Array(await buffer.readAsync(byteOffset));
      let n = 0;
      for (let i = 0; i < list.count; ++i) n += args[(records ? records[i] : i) * words + 1];
      return n;
    };
    const gpuInstances = await listInstances();
    expect(gpuInstances).toBeGreaterThan(0);
    expect(gpuInstances).toBeLessThan(cpuInstances);
    expect(new Set(cpuImage).size).toBeGreaterThan(2);
    expect(gpuImage).toEqual(cpuImage);

    // the instances of all records of the list, drawn or not: those it leaves out are empty
    const allInstances = async () => {
      const args = new Uint32Array(await buffer.readAsync(byteOffset));
      let n = 0;
      for (let i = 0; i < members.length; ++i) n += args[i * words + 1];
      return n;
    };

    // a still camera: the lists drawing only the records with instances, the same instances and image
    for (let i = 0; i < 10; ++i) await frame(gpuScene);
    const stillImage = await frame(gpuScene);
    const still = drawn.get(gpuScene)![0];
    expect(still.indirect!.records).toBeDefined();
    expect(still.count).toBeLessThanOrEqual(members.length);
    expect(await listInstances()).toBe(gpuInstances);
    expect(await allInstances()).toBe(gpuInstances);
    expect(stillImage).toEqual(gpuImage);

    // moving: the records the bounds do not rule out, none with instances left out; back, the same image at rest
    const home = camera.getSnapshot();
    for (const x of [0.2, 0.6, 1.1, -0.8]) {
      camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * x, 0, 0)), radius * 0.5), 0);
      camera.update();
      await frame(gpuScene);
      expect(await listInstances()).toBe(await allInstances());
    }
    camera.setState(home, 0);
    camera.update();
    for (let i = 0; i < 10; ++i) await frame(gpuScene);
    expect(await frame(gpuScene)).toEqual(gpuImage);

    // another view not culled on the GPU (no `prepareCull`) draws none: culled there alone, no CPU culling
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(-radius * 0.3, 0, 0)), radius * 0.35), 0);
    camera.update();
    const unculled = await frame(gpuScene, false);
    expect(unculled).toEqual(await frame(Scene.create(graphics), false));
    // culled on the GPU for the sphere raster alone, which does not take it here (off): the CPU culling's lists from
    // the instances in their own order after the regions (the identity region), its transforms as the rows the
    // regions hold
    const rasterScene = createScene(createMergedRenderObject(members, { gpuCulling: true, rasterOnly: true }));
    await graphics.resources.preparePipelines(
      renderer.collect({ name: 'opaque' }, rasterScene, graphics.describeTarget(target)),
    );
    const rasterItem = rasterScene.renderables[0].renderItem;
    const rasterRender = rasterItem.render;
    rasterItem.render = (variant, globals, list) => {
      if (list) drawn.set(rasterScene, list);
      rasterRender(variant, globals, list);
    };
    const identityImage = await frame(rasterScene);
    expect(drawn.get(rasterScene)![0].indirect).toBeUndefined();
    const movedImage = await frame(cpuScene);
    expect(new Set(movedImage).size).toBeGreaterThan(2);
    expect(identityImage).toEqual(movedImage);

    // `update` never culls on the GPU: passes call it with their render pass open (pick, marking), a compute pass
    // would end it; `prepareCull` does, before the passes
    const renderable = gpuScene.renderables[0];
    const prepare = renderable.prepareCull!;
    let prepared = 0;
    renderable.prepareCull = (...args) => {
      prepared += 1;
      return prepare(...args);
    };
    renderer.update(camera, gpuScene, createFrame());
    expect(prepared).toBe(0);
    renderer.prepareCull(gpuScene);
    expect(prepared).toBe(1);
  }, 60_000);

  it('culls the same in rows of workgroups past the per-dimension limit', async () => {
    const plugin = await sharedPlugin();
    const { maxWorkgroupsPerDimension, dataBlockInstances, sceneCellInstances } = MergedCullerLimits;
    try {
      // the instances' data written a few instances at a time (`writeCulledData`), each member in several blocks
      MergedCullerLimits.dataBlockInstances = 7;
      // cells of an instance each (`buildSceneCells`): the members in planes apart, small cells
      MergedCullerLimits.sceneCellInstances = 0.1;
      const graphics = plugin.canvas3d!.graphics;
      // more than one workgroup (64) of cells and of instances: split into rows below
      const members = [
        createMember(Vec3.create(0, 0, 0), ColorNames.orange),
        createMember(Vec3.create(2.5, 2.5, 10), ColorNames.green),
        createMember(Vec3.create(-2.5, 3, 20), ColorNames.blue),
      ];
      const cpu = createMergedRenderObject(members);
      const gpu = createMergedRenderObject(members, { gpuCulling: true });
      const createScene = (o: GraphicsRenderObject) => {
        const scene = Scene.create(graphics);
        scene.add(o);
        scene.commit();
        return scene;
      };
      const cpuScene = createScene(cpu);
      const gpuScene = createScene(gpu);

      const camera = new Camera();
      Viewport.set(camera.viewport, 0, 0, Size, Size);
      const { center, radius } = cpuScene.boundingSphere;
      const view = (x: number, distance: number) => {
        camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * x, 0, 0)), radius * distance), 0);
        camera.update();
      };
      view(0.3, 0.35);

      const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      for (const scene of [cpuScene, gpuScene]) {
        await graphics.resources.preparePipelines(
          renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
      }
      let list: MultiDrawBaseData | undefined;
      const item = gpuScene.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, l) => {
        if (l) list = l[0];
        render(variant, globals, l);
      };
      const frame = async (scene: Scene) => {
        renderer.update(camera, scene, createFrame());
        if (scene === gpuScene) renderer.prepareCull(scene);
        graphics.beginRenderPass(target);
        renderer.setViewport(0, 0, Size, Size);
        renderer.clear(true);
        renderer.renderOpaque(scene, camera);
        graphics.flush();
        const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
        await graphics.resources.pending.whenReady();
        return pixels;
      };
      // the instance counts of every record (instances' order within a region is the atomics')
      const records = async () => {
        const { buffer, byteOffset, indexed } = list!.indirect!;
        const words = indexed ? 5 : 4;
        const args = new Uint32Array(await buffer.readAsync(byteOffset));
        return Array.from({ length: args.length / words }, (_, i) => args[i * words + 1]);
      };
      for (let i = 0; i < 20 && !list?.indirect; ++i) await frame(gpuScene);
      expect(list?.indirect).toBeDefined();
      const image = await frame(gpuScene);
      const counts = await records();
      expect(counts.some((n) => n > 0)).toBe(true);

      // a workgroup a row: every dispatch of more than one workgroup in rows; culled again (another view and back)
      MergedCullerLimits.maxWorkgroupsPerDimension = 1;
      const rows = MergedCullerCounts.rows;
      view(-0.5, 0.5);
      await frame(gpuScene);
      view(0.3, 0.35);
      const splitImage = await frame(gpuScene);
      // cells and instances, both phases' dispatches (the first only here)
      expect(MergedCullerCounts.rows - rows).toBeGreaterThanOrEqual(4);
      expect(await records()).toEqual(counts);
      expect(splitImage).toEqual(image);
      expect(new Set(splitImage).size).toBeGreaterThan(2);
      expect(splitImage).toEqual(await frame(cpuScene));
    } finally {
      MergedCullerLimits.maxWorkgroupsPerDimension = maxWorkgroupsPerDimension;
      MergedCullerLimits.dataBlockInstances = dataBlockInstances;
      MergedCullerLimits.sceneCellInstances = sceneCellInstances;
    }
  }, 60_000);

  // marking is redrawn over the frame without culling again: a change of markers keeps the frame's GPU culling, one of
  // what it culled does not
  it('keeps the GPU culling of a frame through an update of what it does not cull by', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const members = [createMember(Vec3.create(0, 0, 0), ColorNames.orange)];
    const scene = Scene.create(graphics);
    scene.add(createMergedRenderObject(members, { gpuCulling: true }));
    scene.commit();
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scene.boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * 0.3, 0, 0)), radius * 0.35), 0);
    camera.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    await graphics.resources.preparePipelines(
      renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
    );

    const renderable = scene.renderables[0];
    let drawn: MultiDrawBaseData[] | undefined;
    const render = renderable.renderItem.render;
    renderable.renderItem.render = (variant, globals, list) => {
      drawn = list;
      render(variant, globals, list);
    };
    const draw = (f: ReturnType<typeof createFrame>, cullOnGpu: boolean) => {
      renderer.update(camera, scene, f);
      if (cullOnGpu) renderer.prepareCull(scene);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(scene, camera);
      graphics.flush();
    };
    // the culling kernels compile in the background: the CPU culling's lists until they are ready
    let frame = createFrame();
    for (let i = 0; i < 20 && !drawn?.[0].indirect; ++i) {
      frame = createFrame();
      draw(frame, true);
      await graphics.resources.pending.whenReady();
    }
    expect(drawn![0].indirect).toBeDefined();

    // markers changed, the same frame drawn again (a redraw): the GPU culling's lists
    const { markerAverage } = members[0].values;
    ValueCell.update(markerAverage, 0.5);
    renderable.update();
    draw(frame, false);
    expect(drawn![0].indirect).toBeDefined();

    // the instances changed: the CPU culling's lists until a frame culls on the GPU again
    const { aTransform } = members[0].values;
    ValueCell.update(aTransform, aTransform.ref.value);
    renderable.update();
    draw(frame, false);
    expect(drawn![0].indirect).toBeUndefined();
  }, 60_000);

  it('culls again on the CPU for a clip moved alone, culled on the GPU for the sphere raster alone (off)', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    // a plane through the grid at `y`, its upper part removed
    const clipAt = (y: number): Clip.Props => ({
      variant: 'pixel',
      objects: [
        {
          type: 'plane',
          invert: false,
          position: Vec3.create(15, y, 0),
          rotation: { axis: Vec3.create(0, 0, 1), angle: 20 },
          scale: Vec3.create(1, 1, 1),
          transform: Mat4.identity(),
        },
      ],
    });
    const member = createMember(Vec3.create(0, 0, 0), ColorNames.orange, { clip: clipAt(1000) });
    const scene = Scene.create(graphics);
    scene.add(createMergedRenderObject([member], { gpuCulling: true, rasterOnly: true }));
    scene.commit();
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scene.boundingSphere;
    camera.setState(camera.getFocus(center, radius), 0);
    camera.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    await graphics.resources.preparePipelines(
      renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
    );
    const renderable = scene.renderables[0];
    let drawn: MultiDrawBaseData[] | undefined;
    const render = renderable.renderItem.render;
    renderable.renderItem.render = (variant, globals, list) => {
      drawn = list;
      render(variant, globals, list);
    };
    const instances = () => {
      renderer.update(camera, scene, createFrame());
      renderer.prepareCull(scene);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(scene, camera);
      graphics.flush();
      expect(drawn![0].indirect).toBeUndefined();
      let n = 0;
      for (const l of drawn!) for (let i = 0; i < l.count; ++i) n += l.instanceCounts[i];
      return n;
    };
    const moveClip = (x: number) => {
      BaseGeometry.updateClip(member.values, clipAt(x));
      renderable.update();
    };
    const all = instances();
    expect(all).toBeGreaterThan(0);
    // the clip moved into the grid: some instances clipped whole
    moveClip(14);
    const clipped = instances();
    expect(clipped).toBeLessThan(all);
    // and away again: all of them back
    moveClip(1000);
    expect(instances()).toBe(all);
  }, 60_000);

  it("tells a frame culled on the GPU stale after a side view's culling (a ray's pick), current culled again", async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const scene = Scene.create(graphics);
    scene.add(createMergedRenderObject([createMember(Vec3.create(0, 0, 0), ColorNames.orange)], { gpuCulling: true }));
    scene.commit();
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scene.boundingSphere;
    camera.setState(camera.getFocus(center, radius), 0);
    camera.update();
    const side = new Camera();
    Viewport.set(side.viewport, 0, 0, Size, Size);
    side.setState(side.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius, 0, 0)), radius * 0.2), 0);
    side.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    await graphics.resources.preparePipelines(
      renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
    );
    const renderable = scene.renderables[0];
    let drawn: MultiDrawBaseData[] | undefined;
    const render = renderable.renderItem.render;
    renderable.renderItem.render = (variant, globals, list) => {
      drawn = list;
      render(variant, globals, list);
    };
    const draw = (f: ReturnType<typeof createFrame>) => {
      renderer.update(camera, scene, f);
      renderer.prepareCull(scene);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(scene, camera);
      graphics.flush();
    };
    let frame = createFrame();
    for (let i = 0; i < 20 && !drawn?.[0].indirect; ++i) {
      frame = createFrame();
      draw(frame);
      await graphics.resources.pending.whenReady();
    }
    expect(drawn![0].indirect).toBeDefined();
    expect(renderer.isGpuCullStale(scene, frame)).toBe(false);

    // a side view culled on the GPU (a ray's pick): the frame's lists that view's now
    renderer.update(side, scene, createFrame());
    renderer.prepareCullView(scene);
    expect(renderer.isGpuCullStale(scene, frame)).toBe(true);

    // the frame's view culled again (as a marking redraw does): current
    renderer.update(camera, scene, frame);
    renderer.prepareCull(scene);
    expect(renderer.isGpuCullStale(scene, frame)).toBe(false);

    // a ray's pick: its view culled on the GPU, then its pass culls with a frame of its own (the object's last cull
    // that frame's, not `frame`'s): stale for `frame` as well
    const ray = createFrame();
    renderer.update(side, scene, ray);
    renderer.prepareCullView(scene);
    graphics.beginRenderPass(target);
    renderer.setViewport(0, 0, Size, Size);
    renderer.renderOpaque(scene, side);
    graphics.flush();
    expect(renderer.isGpuCullStale(scene, frame)).toBe(true);
    renderer.update(camera, scene, frame);
    renderer.prepareCull(scene);
    expect(renderer.isGpuCullStale(scene, frame)).toBe(false);
  }, 60_000);

  it('draws and picks in a canvas with draw lists (passes that update the renderer inside their render pass)', {
    tags: ['medium'],
    timeout: 60_000,
  }, async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: Size, height: Size }, { own: true });
    try {
      await plugin.init();
      const pdb = `${Array.from({ length: 48 }, (_, i) => {
        const f = (v: number) => v.toFixed(3).padStart(8);
        const chain = 'ABCD'[Math.floor(i / 12)];
        // a compact grid filling the view, four chains (units, merged)
        const [x, y, z] = [(i % 4) * 2.5, (Math.floor(i / 4) % 4) * 2.5, Math.floor(i / 16) * 2.5];
        return `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String((i % 12) + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
      }).join('\n')}\nEND\n`;
      const data = await plugin.builders.data.rawData({ data: pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { gpuCulling: true },
      });
      const canvas = plugin.canvas3d!;
      canvas.setProps({ renderer: { ...canvas.props.renderer, drawLists: true } });
      const pixelRatio = canvas.graphics.pixelRatio;
      let time = now();
      for (let i = 0; i < 5; ++i) {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        const pick = await canvas.asyncIdentify(Vec2.create(Size / 2 / pixelRatio, Size / 2 / pixelRatio))?.wait();
        // the first pick's draws compile (skipped), the later ones hit
        if (i > 0) expect(pick?.id.objectId).toBeGreaterThanOrEqual(0);
      }
    } finally {
      plugin.dispose();
    }
  });

  it('occlusion against the depth pyramid draws what is drawn without', {
    tags: ['medium'],
    timeout: 60_000,
  }, async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: Size, height: Size }, { own: true });
    try {
      await plugin.init();
      // two chains, solid walls of atoms (spaced closer than their radii), as 20 copies one behind the other along z
      const f = (v: number) => v.toFixed(3).padStart(8);
      const atoms: string[] = [];
      for (let i = 0; i < 64; ++i) {
        const chain = i < 32 ? 'A' : 'B';
        const j = i % 32;
        const [x, y, z] = [(j % 4) * 1.2 + (i < 32 ? 0 : 4.8), (Math.floor(j / 4) % 4) * 1.2, Math.floor(j / 16) * 1.2];
        atoms.push(
          `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String(j + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`,
        );
      }
      const biomt: string[] = ['REMARK 350 BIOMOLECULE: 1', 'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B'];
      for (let k = 0; k < 20; ++k) {
        const n = String(k + 1).padStart(4);
        biomt.push(`REMARK 350   BIOMT1 ${n}  1.000000  0.000000  0.000000        0.00000`);
        biomt.push(`REMARK 350   BIOMT2 ${n}  0.000000  1.000000  0.000000        0.00000`);
        biomt.push(`REMARK 350   BIOMT3 ${n}  0.000000  0.000000  1.000000 ${f(-4 * k)}`);
      }
      const pdb = `${biomt.join('\n')}\n${atoms.join('\n')}\nEND\n`;
      const data = await plugin.builders.data.rawData({ data: pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model, {
        name: 'assembly',
        params: { id: '1' },
      });
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        // a cell per copy and chain
        typeParams: { gpuCulling: true, cellSize: 3, batchSize: 12 },
      });
      const canvas = plugin.canvas3d!;
      // the GPU's occlusion test alone (by default, `hiZ.gpu`; Hi-Z's CPU test off)
      canvas.setProps({
        multiSample: { ...canvas.props.multiSample, mode: 'off' },
        hiZ: { ...canvas.props.hiZ, enabled: false, gpu: true },
      });
      // the viewports synced (the Hi-Z pass's among them)
      canvas.handleResize();
      const passes = plugin.renderer.canvas3dContext!.passes;
      let time = now();
      const frames = async (count: number) => {
        for (let i = 0; i < count; ++i) {
          canvas.commit(true);
          canvas.requestDraw();
          canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
        }
        return Array.from(await canvas.graphics.readRenderTarget(passes.draw.drawTarget, 0, 0, Size, Size));
      };
      await frames(1);
      // looking along -z at the front copy, the others behind it
      const { center, radius } = canvas.boundingSphere;
      const target = Vec3.clone(center);
      const position = Vec3.add(Vec3(), center, Vec3.create(0, 0, radius * 1.5));
      canvas.camera.setState({ target, position, up: Vec3.create(0, 1, 0) }, 0);
      // a pyramid of the same view first (the culling tests the previous frame's)
      const occluded = await frames(6);
      canvas.setProps({ hiZ: { ...canvas.props.hiZ, enabled: false, gpu: false } });
      const all = await frames(3);
      expect(new Set(all).size).toBeGreaterThan(2);
      expect(occluded).toEqual(all);
    } finally {
      plugin.dispose();
    }
  });
  it('culls a merged object of one unit on the GPU; by default (gpuCulling auto) not with few instances but for the raster', {
    tags: ['medium'],
    timeout: 60_000,
  }, async () => {
    const plugin = await createPlugin();
    try {
      const structure = await copiesStructure(plugin);
      const canvas = plugin.canvas3d!;
      const visible = () => {
        canvas.commit(true);
        return canvas.getRenderObjects().filter((o) => o.state.visible);
      };
      const repr = await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const culledOnGpu = (o: GraphicsRenderObject) =>
        (o.values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value;
      const rasterOnly = (o: GraphicsRenderObject) =>
        (o.values as { gpuCullingRasterOnly?: { ref: { value: boolean } } }).gpuCullingRasterOnly?.ref.value;
      // 40 instances: merged (a units representation's always), fewer than `auto` culls on the GPU by; spheres where the
      // device has the sphere raster culled there for it alone (the renderer culls them on the CPU where it does not draw
      // them), as the raster takes merged objects culled on the GPU only
      const raster = !!canvas.graphics.capabilities.atomicMinMax64;
      expect(visible().map(isMergedRenderObject)).toEqual([true]);
      expect(visible().map(culledOnGpu)).toEqual([raster]);
      if (raster) expect(visible().map(rasterOnly)).toEqual([true]);

      await plugin.build().delete(repr).commit();
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { gpuCulling: true },
      });
      const objects = visible();
      expect(objects.length).toBe(1);
      expect(isMergedRenderObject(objects[0])).toBe(true);
      expect((objects[0].values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value).toBe(true);
    } finally {
      plugin.dispose();
    }
  });

  // an object culled on the GPU draws no CPU lists: the overlay culls its camera's view there for it
  it('draws and picks an object culled on the GPU in the overlay layer', {
    tags: ['medium'],
    timeout: 60_000,
  }, async () => {
    const plugin = await createPlugin();
    try {
      const structure = await copiesStructure(plugin);
      const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { gpuCulling: true },
      });
      const canvas = plugin.canvas3d!;
      const r = repr.data!.repr;
      r.setState({ layer: 'overlay' });
      canvas.update(r);
      canvas.commit(true);
      expect(r.renderObjects.some((o) => isMergedRenderObject(o))).toBe(true);
      const pixelRatio = canvas.graphics.pixelRatio;
      let time = now();
      let picked: number | undefined;
      // its culling's kernels may compile after the frame's
      for (let i = 0; i < 20 && picked === undefined; ++i) {
        canvas.commit(true);
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
        await canvas.graphics.resources.pending.whenReady();
        const pick = await canvas.asyncIdentify(Vec2.create(Size / 2 / pixelRatio, Size / 2 / pixelRatio))?.wait();
        picked = pick?.id.objectId;
      }
      expect(picked).toBeGreaterThanOrEqual(0);
      plugin.renderer.imagePass.setProps({ overlay: true });
      const { data } = await plugin.getImageRaw({ width: Size, height: Size });
      const background = [data[0], data[1], data[2]];
      const center = ((Size / 2) * Size + Size / 2) * 4;
      expect([data[center], data[center + 1], data[center + 2]]).not.toEqual(background);
    } finally {
      plugin.dispose();
    }
  });

  it('prepares the culling kernels by a precompile list naming a GPU-culled object, in a context without one', {
    tags: ['medium'],
    timeout: 60_000,
  }, async () => {
    const a = await createPlugin();
    const b = await createPlugin();
    try {
      const structure = await copiesStructure(a);
      const compute = a.canvas3d!.graphics.compute!;
      expect(compute.hasKernels(MergedCullKernels)).toBe(false);
      await a.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { gpuCulling: true },
      });
      a.canvas3d!.commit(true);
      // ready without a frame drawn here
      await vi.waitFor(() => expect(compute.hasKernels(MergedCullKernels)).toBe(true), { timeout: 20_000 });

      a.canvas3d!.requestDraw();
      a.canvas3d!.tick(now(), { isSynchronous: true });
      const list = a.canvas3d!.getPrecompileList();
      expect(list.objects.some((s) => s.gpuCulled)).toBe(true);

      // another context, nothing in its scene
      const other = b.canvas3d!.graphics.compute!;
      expect(other.hasKernels(MergedCullKernels)).toBe(false);
      await b.canvas3d!.precompile(list);
      expect(other.hasKernels(MergedCullKernels)).toBe(true);
    } finally {
      a.dispose();
      b.dispose();
    }
  });

  // the merged object's transient transform (a drag's frame, `uObjectTransform`): a uniform of the draws and the culling
  // kernels, nothing of the layout redone; drawn as the same move baked into the members' instances
  it('draws a merged object moved by its transient transform as the move baked, culled on the GPU and the CPU', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const move = Vec3.create(4, -14, 0);
    const offsets = [Vec3.create(0, 0, 0), Vec3.create(2.5, 2.5, 1), Vec3.create(-2.5, 3, -1)];
    const colors = [ColorNames.orange, ColorNames.green, ColorNames.blue];
    const members = offsets.map((o, i) => createMember(o, colors[i]));
    const baked = offsets.map((o, i) => createMember(Vec3.add(Vec3(), o, move), colors[i]));
    const scene = (o: GraphicsRenderObject) => {
      const s = Scene.create(graphics);
      s.add(o);
      s.commit();
      return s;
    };
    const gpuObject = createMergedRenderObject(members, { gpuCulling: true });
    const cpuObject = createMergedRenderObject(members);
    const gpuScene = scene(gpuObject);
    const cpuScene = scene(cpuObject);
    const bakedScene = scene(createMergedRenderObject(baked));

    // all of the moved grid in view, the frustum culling some of the members'
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = bakedScene.boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(0, -radius * 0.3, 0)), radius * 0.6), 0);
    camera.update();

    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    for (const s of [gpuScene, cpuScene, bakedScene]) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
      );
    }
    let gpuDrawn: MultiDrawBaseData[] | undefined;
    const item = gpuScene.renderables[0].renderItem;
    const render = item.render;
    item.render = (variant, globals, list) => {
      if (list) gpuDrawn = list;
      render(variant, globals, list);
    };
    const frame = async (s: Scene) => {
      renderer.update(camera, s, createFrame());
      renderer.prepareCull(s);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(s, camera);
      graphics.flush();
      const pixels = await graphics.readRenderTarget(target, 0, 0, Size, Size);
      await graphics.resources.pending.whenReady();
      return Array.from(pixels);
    };
    /** the pixels further than 10 levels apart (any channel): float rounding of the two ways at the edges */
    const differing = (a: number[], b: number[]) => {
      let n = 0;
      for (let i = 0; i < a.length; i += 4) {
        if (Math.max(...[0, 1, 2].map((c) => Math.abs(a[i + c] - b[i + c]))) > 10) n += 1;
      }
      return n;
    };
    for (let i = 0; i < 20 && !gpuDrawn?.[0].indirect; ++i) await frame(gpuScene);
    const still = await frame(gpuScene);

    // the object moved after merging, as a drag moves it
    for (const o of [gpuObject, cpuObject])
      ValueCell.update(o.values.uObjectTransform, Mat4.fromTranslation(Mat4(), move));
    for (const s of [gpuScene, cpuScene]) s.renderables[0].update();
    const expected = await frame(bakedScene);
    const gpuImage = await frame(gpuScene);
    expect(gpuDrawn![0].indirect).toBeDefined();
    const cpuImage = await frame(cpuScene);
    expect(differing(expected, still)).toBeGreaterThan(50);
    expect(differing(gpuImage, expected)).toBeLessThan(6);
    expect(differing(cpuImage, expected)).toBeLessThan(6);
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  // a placement moved as a whole (the release of a drag): the members' grids moved (`moveInstanceGrid`), not built again;
  // the culler moves its cells instead of splitting them again, drawn as members built at the moved place
  it('moves its cells for grids moved from those it split, drawing as members built there', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const move = Vec3.create(4, -14, 0);
    const offsets = [Vec3.create(0, 0, 0), Vec3.create(2.5, 2.5, 1), Vec3.create(-2.5, 3, -1)];
    const colors = [ColorNames.orange, ColorNames.green, ColorNames.blue];
    const members = offsets.map((o, i) => createMember(o, colors[i]));
    const baked = offsets.map((o, i) => createMember(Vec3.add(Vec3(), o, move), colors[i]));
    const scene = (o: GraphicsRenderObject) => {
      const s = Scene.create(graphics);
      s.add(o);
      s.commit();
      return s;
    };
    const movedScene = scene(createMergedRenderObject(members, { gpuCulling: true }));
    const bakedScene = scene(createMergedRenderObject(baked, { gpuCulling: true }));

    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = bakedScene.boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(0, -radius * 0.3, 0)), radius * 0.6), 0);
    camera.update();

    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    for (const s of [movedScene, bakedScene]) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
      );
    }
    let drawn: MultiDrawBaseData[] | undefined;
    const item = movedScene.renderables[0].renderItem;
    const render = item.render;
    item.render = (variant, globals, list) => {
      if (list) drawn = list;
      render(variant, globals, list);
    };
    const frame = async (s: Scene) => {
      renderer.update(camera, s, createFrame());
      renderer.prepareCull(s);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(s, camera);
      graphics.flush();
      const pixels = await graphics.readRenderTarget(target, 0, 0, Size, Size);
      await graphics.resources.pending.whenReady();
      return Array.from(pixels);
    };
    const differing = (a: number[], b: number[]) => {
      let n = 0;
      for (let i = 0; i < a.length; i += 4) {
        if (Math.max(...[0, 1, 2].map((c) => Math.abs(a[i + c] - b[i + c]))) > 10) n += 1;
      }
      return n;
    };
    for (let i = 0; i < 20 && !drawn?.[0].indirect; ++i) await frame(movedScene);
    const still = await frame(movedScene);

    // each member moved as a placement change moves it: its grid, its instances in the grid's order, its bounds
    const m = Mat4.fromTranslation(Mat4(), move);
    for (const { values } of members) {
      const grid = moveInstanceGrid(values.instanceGrid.ref.value, m);
      ValueCell.update(values.instanceGrid, grid);
      ValueCell.update(values.aTransform, grid.cellTransform);
      const b = values.boundingSphere.ref.value;
      ValueCell.update(values.boundingSphere, Sphere3D.create(Vec3.add(Vec3(), b.center, move), b.radius));
    }
    const { layouts, moves } = MergedCullerCounts;
    movedScene.renderables[0].update();
    expect(MergedCullerCounts.layouts).toBe(layouts);
    expect(MergedCullerCounts.moves).toBe(moves + 1);

    const expected = await frame(bakedScene);
    const image = await frame(movedScene);
    expect(drawn![0].indirect).toBeDefined();
    expect(differing(expected, still)).toBeGreaterThan(50);
    expect(differing(image, expected)).toBeLessThan(6);
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  // clip objects (by pixel, without a clipping layer) cull the cells and instances they clip whole on the GPU as on the
  // CPU: fewer instances drawn, the same image
  it('leaves out what clip objects clip whole, drawing what the CPU culling draws', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    // a plane through the grid, its upper part removed
    const clip = {
      variant: 'pixel' as const,
      objects: [
        {
          type: 'plane' as const,
          invert: false,
          position: Vec3.create(15, 14, 0),
          rotation: { axis: Vec3.create(0, 0, 1), angle: 20 },
          scale: Vec3.create(1, 1, 1),
          transform: Mat4.identity(),
        },
      ],
    };
    const members = (params: Partial<PD.Values<Spheres.Params>>) => [
      createMember(Vec3.create(0, 0, 0), ColorNames.orange, params),
      createMember(Vec3.create(2.5, 2.5, 1), ColorNames.green, params),
    ];
    const scene = (o: GraphicsRenderObject) => {
      const s = Scene.create(graphics);
      s.add(o);
      s.commit();
      return s;
    };
    const clipped = members({ clip });
    const cpuScene = scene(createMergedRenderObject(clipped));
    const gpuScene = scene(createMergedRenderObject(clipped, { gpuCulling: true }));
    const wholeScene = scene(createMergedRenderObject(members({}), { gpuCulling: true }));

    // all of the grid in view
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = wholeScene.boundingSphere;
    camera.setState(camera.getFocus(center, radius), 0);
    camera.update();

    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    for (const s of [cpuScene, gpuScene, wholeScene]) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
      );
    }
    const drawn = new Map<Scene, MultiDrawBaseData[]>();
    for (const s of [cpuScene, gpuScene, wholeScene]) {
      const item = s.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, list) => {
        if (list) drawn.set(s, list);
        render(variant, globals, list);
      };
    }
    const frame = async (s: Scene) => {
      renderer.update(camera, s, createFrame());
      renderer.prepareCull(s);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(s, camera);
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    const gpuInstances = async (s: Scene) => {
      const list = drawn.get(s)![0];
      const { buffer, byteOffset, indexed, records } = list.indirect!;
      const words = indexed ? 5 : 4;
      const args = new Uint32Array(await buffer.readAsync(byteOffset));
      let n = 0;
      for (let i = 0; i < list.count; ++i) n += args[(records ? records[i] : i) * words + 1];
      return n;
    };
    for (const s of [gpuScene, wholeScene]) {
      for (let i = 0; i < 20 && !drawn.get(s)?.[0].indirect; ++i) await frame(s);
    }
    const gpuImage = await frame(gpuScene);
    const clippedCount = await gpuInstances(gpuScene);
    await frame(wholeScene);
    const wholeCount = await gpuInstances(wholeScene);
    const cpuImage = await frame(cpuScene);
    expect(clippedCount).toBeGreaterThan(0);
    expect(clippedCount).toBeLessThan(wholeCount * 0.8);
    expect(new Set(cpuImage).size).toBeGreaterThan(2);
    expect(gpuImage).toEqual(cpuImage);
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  });

  // the ghost of what the clip objects remove (`Renderer.renderClipGhost`) culled on the GPU into a region of its own:
  // the instances in view the clip objects reach, drawing what the CPU culling's ghost draws
  it('culls the ghost of what clip objects remove into its own region, drawing what the CPU culling draws', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const clip = {
      variant: 'pixel' as const,
      objects: [
        {
          type: 'plane' as const,
          invert: false,
          position: Vec3.create(15, 14, 0),
          rotation: { axis: Vec3.create(0, 0, 1), angle: 20 },
          scale: Vec3.create(1, 1, 1),
          transform: Mat4.identity(),
        },
      ],
    };
    const members = (params: Partial<PD.Values<Spheres.Params>>) => [
      createMember(Vec3.create(0, 0, 0), ColorNames.orange, params),
      createMember(Vec3.create(2.5, 2.5, 1), ColorNames.green, params),
    ];
    const scene = (o: GraphicsRenderObject) => {
      const s = Scene.create(graphics);
      s.add(o);
      s.commit();
      return s;
    };
    const cpuScene = scene(createMergedRenderObject(members({ clip })));
    const gpuScene = scene(createMergedRenderObject(members({ clip }), { gpuCulling: true }));
    const wholeScene = scene(createMergedRenderObject(members({}), { gpuCulling: true }));

    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = wholeScene.boundingSphere;
    camera.setState(camera.getFocus(center, radius), 0);
    camera.update();

    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    // packed depth, as the ghost pass's target
    const ghostTarget = graphics.createRenderTarget(Size, Size, 'depth-stencil', 'uint8', 'nearest');
    for (const s of [cpuScene, gpuScene, wholeScene]) {
      await graphics.resources.preparePipelines([
        ...renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
        ...renderer.collect({ name: 'clipGhost' }, s, graphics.describeTarget(ghostTarget)),
      ]);
    }
    const drawn = new Map<Scene, MultiDrawBaseData[]>();
    let drawingGhost = false;
    for (const s of [cpuScene, gpuScene, wholeScene]) {
      const item = s.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, list) => {
        if (list && drawingGhost) drawn.set(s, list);
        render(variant, globals, list);
      };
    }
    const frame = async (s: Scene) => {
      renderer.update(camera, s, createFrame());
      renderer.prepareCull(s, false, true);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(s, camera);
      graphics.beginRenderPass(ghostTarget);
      renderer.clearDepth(true);
      drawingGhost = true;
      renderer.renderClipGhost(s, camera);
      drawingGhost = false;
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(ghostTarget, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    const gpuInstances = async (s: Scene) => {
      const list = drawn.get(s)![0];
      const { buffer, byteOffset, indexed, records } = list.indirect!;
      const words = indexed ? 5 : 4;
      const args = new Uint32Array(await buffer.readAsync(byteOffset));
      let n = 0;
      for (let i = 0; i < list.count; ++i) n += args[(records ? records[i] : i) * words + 1];
      return n;
    };
    for (let i = 0; i < 20 && !drawn.get(gpuScene)?.[0].indirect; ++i) await frame(gpuScene);
    const gpuGhost = await frame(gpuScene);
    const ghostCount = await gpuInstances(gpuScene);
    const cpuGhost = await frame(cpuScene);
    const instanceCount = wholeScene.renderables[0].values.instanceCount.ref.value as number;
    expect(ghostCount).toBeGreaterThan(0);
    expect(ghostCount).toBeLessThan(instanceCount);
    // something removed, drawn as packed depth
    expect(new Set(cpuGhost).size).toBeGreaterThan(2);
    expect(gpuGhost).toEqual(cpuGhost);
    target.destroy();
    ghostTarget.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  });

  it('draws members with pending grids as with built ones: their data in their order, the grids built for a CPU cull (raster alone)', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const pending = createPosedMembers(true);
    const built = createPosedMembers();
    for (const t of pending.poseTransforms) {
      expect(t.instanceGrid.ref.value.pending).toBeDefined();
      expect(t.instanceGrid.ref.value.cellCount).toBe(0);
    }
    const scene = (objects: GraphicsRenderObject[]) => {
      const s = Scene.create(graphics);
      for (const o of objects) s.add(o);
      s.commit();
      return s;
    };
    const pendingRaster = createPosedMembers(true);
    // culled on the GPU (alone), and for the sphere raster alone (culled on the CPU where the raster does not take it:
    // here, `sphereRaster` off)
    const scenes = [
      scene([createMergedRenderObject(pending.poses, { gpuCulling: true })]),
      scene([createMergedRenderObject(built.poses, { gpuCulling: true })]),
      scene([createMergedRenderObject(pendingRaster.poses, { gpuCulling: true, rasterOnly: true })]),
      scene([createMergedRenderObject(createPosedMembers().poses, { gpuCulling: true, rasterOnly: true })]),
    ];
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scenes[1].boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * 0.2, 0, 0)), radius * 0.5), 0);
    camera.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    const indirect = new Set<Scene>();
    for (const s of scenes) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
      );
      const item = s.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, list) => {
        if (list?.[0]?.indirect) indirect.add(s);
        render(variant, globals, list);
      };
    }
    const frame = async (s: Scene, cullOnGpu = true) => {
      renderer.update(camera, s, createFrame());
      if (cullOnGpu) renderer.prepareCull(s);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(s, camera);
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    for (let i = 0; i < 20 && !(indirect.has(scenes[0]) && indirect.has(scenes[1])); ++i) {
      await frame(scenes[0]);
      await frame(scenes[1]);
    }
    // culled on the GPU: the same image, no grid built, the members' data in their own order
    const gpu = [await frame(scenes[0]), await frame(scenes[1])];
    expect(indirect.has(scenes[0])).toBe(true);
    expect(new Set(gpu[1]).size).toBeGreaterThan(2);
    expect(countDifferentPixels(gpu[0], gpu[1])).toBeLessThanOrEqual(Size * Size * 0.002);
    for (const t of pending.poseTransforms) {
      expect(t.instanceGrid.ref.value.pending).toBeDefined();
      expect(Array.from(t.aInstance.ref.value)).toEqual(Array.from({ length: t.instanceCount.ref.value }, (_, k) => k));
    }
    // a view not culled on the GPU (no `prepareCull`): none drawn where culled there alone (no CPU culling, no grid
    // built); for the raster alone culled on the CPU, the grids built, the lists their cells' in the identity region,
    // what the built grids draw; the members' data as it was
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(-radius * 0.3, 0, 0)), radius * 0.4), 0);
    camera.update();
    const empty = await frame(scene([]), false);
    expect(await frame(scenes[0], false)).toEqual(empty);
    expect(await frame(scenes[1], false)).toEqual(empty);
    const cpu = [await frame(scenes[2]), await frame(scenes[3])];
    expect(new Set(cpu[1]).size).toBeGreaterThan(2);
    expect(countDifferentPixels(cpu[0], cpu[1])).toBeLessThanOrEqual(Size * Size * 0.002);
    for (const t of pendingRaster.poseTransforms) expect(t.instanceGrid.ref.value.pending).toBeDefined();
    // back on the GPU: as before
    expect(countDifferentPixels(await frame(scenes[0]), await frame(scenes[1]))).toBeLessThanOrEqual(
      Size * Size * 0.002,
    );

    // not merged: nothing draws them by index, the pending grids built and applied (the data in their order)
    const unmerged = createPosedMembers(true);
    const plain = scene(unmerged.poses);
    for (const t of unmerged.poseTransforms) {
      expect(t.instanceGrid.ref.value.pending).toBeUndefined();
      expect(t.instanceGrid.ref.value.cellCount).toBeGreaterThan(0);
    }
    const reference = scene(createPosedMembers().poses);
    for (const s of [plain, reference]) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
      );
    }
    expect(countDifferentPixels(await frame(plain, false), await frame(reference, false))).toBe(0);
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  it("culls on the GPU alone: views of their own (the thickness pass, a ray's pick), a frame's second camera, none meanwhile", async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    // the culling's kernels compiling, as the context tells (an earlier test may have compiled them): from before the
    // object is made, which creates them where compiled
    const hasKernels = vi.spyOn(graphics.compute!, 'hasKernels').mockReturnValue(false);
    try {
      const pending = createPosedMembers(true);
      const scene = Scene.create(graphics);
      scene.add(createMergedRenderObject(pending.poses, { gpuCulling: true }));
      scene.commit();
      const { center, radius } = scene.boundingSphere;
      const cameraAt = (offset: number, distance: number) => {
        const camera = new Camera();
        Viewport.set(camera.viewport, 0, 0, Size, Size);
        camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * offset, 0, 0)), distance), 0);
        camera.update();
        return camera;
      };
      // the main camera at a part of the members, the other at all of them
      const main = cameraAt(0.2, radius * 0.5);
      const other = cameraAt(0, radius * 1.2);
      const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
      const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
      );
      const item = scene.renderables[0].renderItem;
      const render = item.render;
      let indirect = false;
      item.render = (variant, globals, list) => {
        indirect = !!list?.[0]?.indirect;
        render(variant, globals, list);
      };
      /** a frame's view of `camera`: culled by `prepare` (the frame's, a view of its own's, none) */
      const draw = async (
        camera: Camera,
        frame: ReturnType<typeof createFrame>,
        prepare: 'frame' | 'view' | 'none',
      ) => {
        renderer.update(camera, scene, frame);
        if (prepare === 'frame') renderer.prepareCull(scene);
        else if (prepare === 'view') renderer.prepareCullView(scene);
        graphics.beginRenderPass(target);
        renderer.setViewport(0, 0, Size, Size);
        renderer.clear(true);
        renderer.renderOpaque(scene, camera);
        graphics.flush();
        const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
        await graphics.resources.pending.whenReady();
        return pixels;
      };
      const empty = Scene.create(graphics);
      renderer.setViewport(0, 0, Size, Size);
      renderer.update(other, empty, createFrame());
      graphics.beginRenderPass(target);
      renderer.clear(true);
      graphics.flush();
      const background = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      // the kernels waited for: nothing drawn until they are ready, no CPU culling meanwhile
      for (let i = 0; i < 2; ++i) expect(await draw(main, createFrame(), 'frame')).toEqual(background);
      hasKernels.mockRestore();
      for (let i = 0; i < 20 && !indirect; ++i) await draw(main, createFrame(), 'frame');
      const mainImage = await draw(main, createFrame(), 'frame');
      expect(indirect).toBe(true);
      const otherImage = await draw(other, createFrame(), 'frame');
      expect(new Set(mainImage).size).toBeGreaterThan(2);
      expect(countDifferentPixels(mainImage, otherImage)).toBeGreaterThan(20);

      // a frame and camera of its own (a ray's pick): culled on the GPU, the frame's image of that camera
      await draw(main, createFrame(), 'frame');
      expect(countDifferentPixels(await draw(other, createFrame(), 'view'), otherImage)).toBe(0);
      expect(indirect).toBe(true);
      // the main view's next frame culled again
      expect(countDifferentPixels(await draw(main, createFrame(), 'frame'), mainImage)).toBe(0);
      // the same camera's frame of its own (the thickness pass) after the frame's: culled for it, and the frame's
      // passes after it draw that view's lists
      const frame = createFrame();
      await draw(main, frame, 'frame');
      expect(countDifferentPixels(await draw(main, createFrame(), 'view'), mainImage)).toBe(0);
      expect(indirect).toBe(true);
      expect(countDifferentPixels(await draw(main, frame, 'none'), mainImage)).toBe(0);
      expect(indirect).toBe(true);
      // a second camera of a frame (a stereo or XR eye): culled for it, the first's again by the next frame
      const eyes = createFrame();
      expect(countDifferentPixels(await draw(main, eyes, 'frame'), mainImage)).toBe(0);
      expect(countDifferentPixels(await draw(other, eyes, 'frame'), otherImage)).toBe(0);
      expect(indirect).toBe(true);
      expect(countDifferentPixels(await draw(main, createFrame(), 'frame'), mainImage)).toBe(0);
      // each of these views drew the GPU culling's lists: none culled on the CPU (no grid built for it); a frame of
      // another view not culled for it draws none
      expect(await draw(other, createFrame(), 'none')).toEqual(background);

      renderer.dispose();
      scene.clear();
      target.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      hasKernels.mockRestore();
    }
  }, 60_000);

  it('updates pending members whose poses changed alone without a layout: the cells and the poses again', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const shift = Vec3.create(1.5, -2, 0.5);
    const posed = createPosedMembers(true);
    const reference = createPosedMembers(true, shift);
    const scene = (objects: GraphicsRenderObject[]) => {
      const s = Scene.create(graphics);
      for (const o of objects) s.add(o);
      s.commit();
      return s;
    };
    const scenes = [
      scene([createMergedRenderObject(posed.poses, { gpuCulling: true })]),
      scene([createMergedRenderObject(reference.poses, { gpuCulling: true })]),
    ];
    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scenes[1].boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * 0.2, 0, 0)), radius * 0.5), 0);
    camera.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    const indirect = new Set<Scene>();
    for (const s of scenes) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, s, graphics.describeTarget(target)),
      );
      const item = s.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, list) => {
        if (list?.[0]?.indirect) indirect.add(s);
        render(variant, globals, list);
      };
    }
    const frame = async (s: Scene) => {
      renderer.update(camera, s, createFrame());
      renderer.prepareCull(s);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(s, camera);
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    for (let i = 0; i < 20 && !(indirect.has(scenes[0]) && indirect.has(scenes[1])); ++i) {
      await frame(scenes[0]);
      await frame(scenes[1]);
    }
    const before = await frame(scenes[0]);
    const expected = await frame(scenes[1]);
    expect(countDifferentPixels(before, expected)).toBeGreaterThan(20);

    // the members' poses changed alone (the reference's), in place: their pending grids kept
    const { layouts, poses } = MergedCullerCounts;
    posed.poseTransforms.forEach((t, k) => {
      const grid = t.instanceGrid.ref.value;
      const next = reference.poseTransforms[k].poses.ref.value!;
      createPosesTransform(next, Sphere3D.create(Vec3(), 2.5), 8, 32, t, true);
      expect(t.instanceGrid.ref.value).toBe(grid);
      const values = posed.poses[k].values;
      ValueCell.update(values.boundingSphere, reference.poses[k].values.boundingSphere.ref.value);
    });
    scenes[0].renderables[0].update();
    expect(MergedCullerCounts.layouts).toBe(layouts);
    expect(MergedCullerCounts.poses).toBe(poses + 1);
    expect(countDifferentPixels(await frame(scenes[0]), expected)).toBeLessThanOrEqual(Size * Size * 0.002);
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  it('draws members placed by poses, read by index, as members of the same matrices', {
    tags: ['medium'],
  }, async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const { poses, matrices } = createPosedMembers();
    const scenes = [
      createMergedRenderObject(poses, { gpuCulling: true }),
      createMergedRenderObject(matrices, { gpuCulling: true }),
      createMergedRenderObject(matrices),
    ].map((o) => {
      const scene = Scene.create(graphics);
      scene.add(o);
      scene.commit();
      return scene;
    });
    // the culled programs read the poses' data (no rows), the CPU-culled object the matrices' rows
    expect(scenes[0].renderables[0].values.instancePoses?.ref.value?.kind).toBe('segmented');
    expect(scenes[1].renderables[0].values.instancePoses?.ref.value).toBeUndefined();

    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scenes[2].boundingSphere;
    expect(scenes[0].boundingSphere.radius).toBeCloseTo(radius, 3);
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * 0.2, 0, 0)), radius * 0.5), 0);
    camera.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    for (const scene of scenes) {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
      );
    }
    const indirect = new Set<Scene>();
    for (const scene of scenes) {
      const item = scene.renderables[0].renderItem;
      const render = item.render;
      item.render = (variant, globals, list) => {
        if (list?.[0]?.indirect) indirect.add(scene);
        render(variant, globals, list);
      };
    }
    const frame = async (scene: Scene) => {
      renderer.update(camera, scene, createFrame());
      renderer.prepareCull(scene);
      graphics.beginRenderPass(target);
      renderer.setViewport(0, 0, Size, Size);
      renderer.clear(true);
      renderer.renderOpaque(scene, camera);
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    for (let i = 0; i < 20 && !(indirect.has(scenes[0]) && indirect.has(scenes[1])); ++i) {
      await frame(scenes[0]);
      await frame(scenes[1]);
    }
    const [posed, culled, cpu] = [await frame(scenes[0]), await frame(scenes[1]), await frame(scenes[2])];
    expect(indirect.has(scenes[0])).toBe(true);
    expect(indirect.has(scenes[2])).toBe(false);
    expect(new Set(cpu).size).toBeGreaterThan(2);
    expect(culled).toEqual(cpu);
    // the rows built on the GPU from the pose, not packed on the CPU: the same image within float rounding
    expect(countDifferentPixels(posed, cpu)).toBeLessThanOrEqual(Size * Size * 0.002);

    // the sphere raster alone: its kernels and the resolve reading the culler's data, an unmerged object's its own
    if (graphics.capabilities.atomicMinMax64) {
      const unmerged = [poses, matrices].map((members) => {
        const scene = Scene.create(graphics);
        for (const m of members) scene.add(m);
        scene.commit();
        return scene;
      });
      const raster = Renderer.create(graphics, { sphereRaster: 'only' });
      const rasterFrame = async (scene: Scene) => {
        await graphics.resources.preparePipelines(
          raster.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
        );
        let resolved = 0;
        for (const r of scene.renderables) {
          const item = r.renderItem;
          const renderFullscreen = item.renderFullscreen!;
          item.renderFullscreen = (...args) => {
            const drawn = renderFullscreen(...args);
            if (drawn) resolved += 1;
            return drawn;
          };
        }
        let pixels: number[] = [];
        for (let i = 0; i < 30; ++i) {
          resolved = 0;
          raster.setViewport(0, 0, Size, Size);
          raster.update(camera, scene, createFrame());
          raster.prepareCull(scene);
          raster.prepareSphereRaster(scene, camera);
          graphics.beginRenderPass(target);
          raster.clear(true);
          raster.renderOpaque(scene, camera);
          graphics.flush();
          pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
          await graphics.resources.pending.whenReady();
          if (resolved > 0 && i > 3) break;
        }
        expect(resolved).toBeGreaterThan(0);
        return pixels;
      };
      // the raster's regions turned on at its first dispatch: the regions laid out again, not the cells
      const { layouts } = MergedCullerCounts;
      const rasterCulled = await rasterFrame(scenes[1]);
      expect(MergedCullerCounts.layouts).toBe(layouts);
      expect(new Set(rasterCulled).size).toBeGreaterThan(2);
      expect(countDifferentPixels(await rasterFrame(scenes[0]), rasterCulled)).toBeLessThanOrEqual(Size * Size * 0.002);
      const rasterUnmerged = await rasterFrame(unmerged[1]);
      expect(new Set(rasterUnmerged).size).toBeGreaterThan(2);
      expect(countDifferentPixels(await rasterFrame(unmerged[0]), rasterUnmerged)).toBeLessThanOrEqual(
        Size * Size * 0.002,
      );
      raster.dispose();
    }
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  });
});

describe.skipIf(!WebGPU.module)('instances of poses on WebGPU', () => {
  it('draws objects of poses culled on the CPU from their instance data, as objects of the same matrices', {
    tags: ['medium'],
  }, async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const { poses, matrices, poseTransforms, matrixTransforms } = createPosedMembers();
    const createScene = (objects: GraphicsRenderObject[]) => {
      const scene = Scene.create(graphics);
      for (const o of objects) scene.add(o);
      scene.commit();
      return scene;
    };
    // unmerged, and merged culled on the CPU
    const scenes = [
      createScene(poses),
      createScene(matrices),
      createScene([createMergedRenderObject(poses)]),
      createScene([createMergedRenderObject(matrices)]),
    ];
    // the instances read from storage by index (`dInstanceStorage`), not packed as rows into vertex buffers
    for (const scene of [scenes[0], scenes[2]]) {
      for (const r of scene.renderables) expect(r.values.dInstanceStorage?.ref.value).toBe(true);
    }
    for (const scene of [scenes[1], scenes[3]]) {
      for (const r of scene.renderables) expect(r.values.dInstanceStorage?.ref.value).toBeFalsy();
    }

    const camera = new Camera();
    Viewport.set(camera.viewport, 0, 0, Size, Size);
    const { center, radius } = scenes[1].boundingSphere;
    camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(radius * 0.2, 0, 0)), radius * 0.5), 0);
    camera.update();
    const renderer = Renderer.create(graphics, { sphereRaster: 'off' });
    const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
    const draw = async (scene: Scene) => {
      await graphics.resources.preparePipelines(
        renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
      );
      // the viewport first: the clear's scissor and the impostors' rays read it (a new renderer's first frame else not
      // cleared, its spheres unshaded)
      renderer.setViewport(0, 0, Size, Size);
      renderer.update(camera, scene, createFrame());
      graphics.beginRenderPass(target);
      renderer.clear(true);
      renderer.renderOpaque(scene, camera);
      graphics.flush();
      const pixels = Array.from(await graphics.readRenderTarget(target, 0, 0, Size, Size));
      await graphics.resources.pending.whenReady();
      return pixels;
    };
    // drawn again where the frame started work (a pipeline the preparation did not cover compiles at the first draw
    // of this context, for matrices as for poses)
    const frame = async (scene: Scene) => {
      for (let i = 0; i < 5; ++i) {
        const tracked = graphics.resources.pending.trackedCount;
        const pixels = await draw(scene);
        if (graphics.resources.pending.trackedCount === tracked) return pixels;
      }
      throw new Error('the frame keeps starting work');
    };
    const compare = async () => {
      const images = [];
      for (const scene of scenes) images.push(await frame(scene));
      expect(new Set(images[1]).size).toBeGreaterThan(2);
      expect(countDifferentPixels(images[0], images[1])).toBeLessThanOrEqual(Size * Size * 0.002);
      expect(countDifferentPixels(images[2], images[3])).toBeLessThanOrEqual(Size * Size * 0.002);
      return images[1];
    };
    const before = await compare();

    // moved: the poses' matrix (no instance touched), the matrices: the instance data written again on the update
    const move = Mat4.fromTranslation(Mat4(), Vec3.create(4, -3, 0));
    for (const t of [poseTransforms[0], matrixTransforms[0]]) expect(moveTransformData(t, move)).toBe(true);
    for (const scene of scenes) {
      scene.update(undefined, false);
      scene.commit();
    }
    expect(countDifferentPixels(await compare(), before)).toBeGreaterThan(0);
    target.destroy();
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  });
});

/** Pixels (RGBA bytes) of `a` and `b` whose channels differ by more than 2 */
function countDifferentPixels(a: ArrayLike<number>, b: ArrayLike<number>) {
  let n = 0;
  for (let o = 0; o < a.length; o += 4) {
    for (let c = 0; c < 4; ++c) {
      if (Math.abs(a[o + c] - b[o + c]) > 2) {
        n += 1;
        break;
      }
    }
  }
  return n;
}

/**
 * `createMember`'s spheres placed by poses of one particle list (rotated, scaled, about a center): three members, by a
 * range of the list (two, the second with a matrix: rotated and moved) and by indices; and members of the same
 * instances as matrices (`InstanceTransforms.getMatrix`)
 */
/**
 * `pendingGrid`: the poses' grids pending (`InstanceGrid.pending`), built where a CPU cull needs them; `shift`: every
 * other particle moved by it (poses changed alone)
 */
function createPosedMembers(pendingGrid = false, shift?: Vec3) {
  const n = Copies * Copies;
  const total = 3 * n + 2;
  const positions = new Float32Array(total * 3);
  const rotations = new Float32Array(total * 4);
  const scales = new Float32Array(total);
  const q = Quat();
  const axis = Vec3();
  for (let e = 0; e < total; ++e) {
    const i = e % n;
    const m = Math.floor(e / n);
    positions[e * 3] = (i % Copies) * 6 + m * 2.5;
    positions[e * 3 + 1] = Math.floor(i / Copies) * 6 - m;
    positions[e * 3 + 2] = m * 1.5;
    Vec3.normalize(axis, Vec3.create(Math.sin(e), Math.cos(e * 1.3), 0.5));
    Quat.toArray(Quat.setAxisAngle(q, axis, e * 0.7), rotations, e * 4);
    // one not above 0: read as 1
    scales[e] = e === 5 ? 0 : 0.8 + (e % 5) * 0.1;
    if (shift && e % 2 === 0) for (let k = 0; k < 3; ++k) positions[e * 3 + k] += shift[k];
  }
  const center = Vec3.create(0.3, -0.2, 0.4);
  const moved = Mat4.mul(
    Mat4(),
    Mat4.fromTranslation(Mat4(), Vec3.create(-3, 2, -1)),
    Mat4.fromRotation(Mat4(), 0.4, Vec3.create(0, 0, 1)),
  );
  const reversed = Array.from({ length: n }, (_, k) => 2 * n + 1 - k);
  const lists = [
    InstanceTransforms.ofPoses({ count: n, positions, rotations, scales, first: 0, center }),
    InstanceTransforms.ofPoses({ count: n, positions, rotations, scales, first: n, center, matrix: moved }),
    InstanceTransforms.ofPoses({ count: n, positions, rotations, scales, indices: reversed, center }),
  ];
  const colors = [ColorNames.orange, ColorNames.green, ColorNames.blue];
  const centers = new Float32Array([0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, -1.2, 0, 0, 0, -1.2, 0]);
  const groups = new Float32Array([0, 1, 2, 3, 4]);
  const props = { ...PD.getDefaultValues(Spheres.Params), cellSize: 8, batchSize: 32 };
  const create = (transform: TransformData, k: number) => {
    const spheres = Spheres.create(centers, groups, 5);
    const simple = BaseGeometry.createSimple(colors[k] as never, 1, transform);
    const values = Spheres.Utils.createValues(spheres, simple.transform, simple.locationIterator, simple.theme, props);
    return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
  };
  const sphere = Sphere3D.create(Vec3(), 2.5);
  const poseTransforms = lists.map((p) =>
    createPosesTransform(p, sphere, props.cellSize, props.batchSize, undefined, pendingGrid),
  );
  const matrixTransforms = lists.map((p) => {
    const array = new Float32Array(n * 16);
    const m = Mat4();
    for (let i = 0; i < n; ++i) Mat4.toArray(InstanceTransforms.getMatrix(p, i, m), array, i * 16);
    return createTransform(array, n, sphere, props.cellSize, props.batchSize);
  });
  return {
    poses: poseTransforms.map(create),
    matrices: matrixTransforms.map(create),
    poseTransforms,
    matrixTransforms,
  };
}

async function createPlugin() {
  const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: Size, height: Size }, { own: true });
  await plugin.init();
  return plugin;
}

// one context for the tests that draw with its graphics alone (their own scenes, renderers and targets)
let shared: Promise<HeadlessPluginContext> | undefined;
function sharedPlugin() {
  shared ??= createPlugin();
  return shared;
}
afterAll(async () => {
  (await shared)?.dispose();
});

/** One chain of 16 atoms, 40 copies: one render object of 40 instances (fewer than `gpuCulling: auto` culls on the GPU) */
async function copiesStructure(plugin: HeadlessPluginContext) {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const atoms = Array.from({ length: 16 }, (_, i) => {
    const [x, y, z] = [(i % 4) * 1.2, Math.floor(i / 4) * 1.2, 0];
    return `ATOM  ${String(i + 1).padStart(5)}  CA  ALA A${String(i + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
  });
  const biomt: string[] = ['REMARK 350 BIOMOLECULE: 1', 'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A'];
  for (let k = 0; k < 40; ++k) {
    const n = String(k + 1).padStart(4);
    biomt.push(`REMARK 350   BIOMT1 ${n}  1.000000  0.000000  0.000000 ${f((k % 8) * 6)}`);
    biomt.push(`REMARK 350   BIOMT2 ${n}  0.000000  1.000000  0.000000 ${f(Math.floor(k / 8) * 6)}`);
    biomt.push(`REMARK 350   BIOMT3 ${n}  0.000000  0.000000  1.000000        0.00000`);
  }
  const pdb = `${biomt.join('\n')}\n${atoms.join('\n')}\nEND\n`;
  const data = await plugin.builders.data.rawData({ data: pdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  return plugin.builders.structure.createStructure(model, { name: 'assembly', params: { id: '1' } });
}
