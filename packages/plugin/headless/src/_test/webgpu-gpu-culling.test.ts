/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { now } from '@molstar/core/util/now';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
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
  MergedCullerCounts,
  MergedCullerLimits,
  MergedCullKernels,
} from '@molstar/graphics/gl/renderable/merged-culler';
import { moveInstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
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
  it('draws the cells the CPU culling draws, from indirect lists', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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

      // another view not culled on the GPU (no `prepareCull`) draws the CPU culling's lists from the instances in their
      // own order after the regions (the identity region): its transforms as the rows the regions hold
      camera.setState(camera.getFocus(Vec3.add(Vec3(), center, Vec3.create(-radius * 0.3, 0, 0)), radius * 0.35), 0);
      camera.update();
      const identityImage = await frame(gpuScene, false);
      expect(drawn.get(gpuScene)![0].indirect).toBeUndefined();
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
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('culls the same in rows of workgroups past the per-dimension limit', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    const { maxWorkgroupsPerDimension, dataBlockInstances } = MergedCullerLimits;
    try {
      await plugin.init();
      // the instances' data written a few instances at a time (`writeCulledData`), each member in several blocks
      MergedCullerLimits.dataBlockInstances = 7;
      const graphics = plugin.canvas3d!.graphics;
      // more than one workgroup (64) of cells and of instances: split into rows below
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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
      plugin.dispose();
    }
  }, 60_000);

  // marking is redrawn over the frame without culling again: a change of markers keeps the frame's GPU culling, one of
  // what it culled does not
  it('keeps the GPU culling of a frame through an update of what it does not cull by', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('draws and picks in a canvas with draw lists (passes that update the renderer inside their render pass)', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
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
        typeParams: { renderMerged: true, gpuCulling: true },
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
  }, 60_000);

  it('occlusion against the depth pyramid draws what is drawn without', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      // two chains, solid blocks of atoms (spaced closer than their radii), A above B; copies of both one behind the
      // other along z, in front of them a copy of A only: what is behind it above is hidden, below not (a test reading
      // the rows of the pyramid mirrored culls what is visible)
      const f = (v: number) => v.toFixed(3).padStart(8);
      const atoms: string[] = [];
      for (let i = 0; i < 64; ++i) {
        const chain = i < 32 ? 'A' : 'B';
        const j = i % 32;
        const [x, y, z] = [(j % 4) * 1.2, (Math.floor(j / 4) % 4) * 1.2 + (i < 32 ? 1 : -6), Math.floor(j / 16) * 1.2];
        atoms.push(
          `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String(j + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`,
        );
      }
      const operator = (k: number, z: number) => {
        const n = String(k).padStart(4);
        return [
          `REMARK 350   BIOMT1 ${n}  1.000000  0.000000  0.000000        0.00000`,
          `REMARK 350   BIOMT2 ${n}  0.000000  1.000000  0.000000        0.00000`,
          `REMARK 350   BIOMT3 ${n}  0.000000  0.000000  1.000000 ${f(z)}`,
        ];
      };
      const biomt: string[] = ['REMARK 350 BIOMOLECULE: 1', 'REMARK 350 APPLY THE FOLLOWING TO CHAINS: A'];
      biomt.push(...operator(1, 8));
      biomt.push('REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B');
      for (let k = 0; k < 20; ++k) biomt.push(...operator(k + 2, -4 * k));
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
        typeParams: { renderMerged: true, gpuCulling: true },
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
  }, 60_000);

  it('occlusion against the depth pyramid draws what is drawn without', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
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
        typeParams: { renderMerged: true, gpuCulling: true, cellSize: 3, batchSize: 12 },
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
  }, 60_000);
  it('merges one object alone to cull it on the GPU; by default (gpuCulling auto) not with few instances', async () => {
    const plugin = await createPlugin();
    try {
      const structure = await copiesStructure(plugin);
      const canvas = plugin.canvas3d!;
      const visible = () => {
        canvas.commit(true);
        return canvas.getRenderObjects().filter((o) => o.state.visible);
      };
      const repr = await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      // 40 instances: fewer than `auto` culls on the GPU, the object stays as it is
      expect(visible().map(isMergedRenderObject)).toEqual([false]);

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
  }, 60_000);

  it('prepares the culling kernels by a precompile list naming a GPU-culled object, in a context without one', async () => {
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
  }, 60_000);

  // the merged object's transient transform (a drag's frame, `uObjectTransform`): a uniform of the draws and the culling
  // kernels, nothing of the layout redone; drawn as the same move baked into the members' instances
  it('draws a merged object moved by its transient transform as the move baked, culled on the GPU and the CPU', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  // a placement moved as a whole (the release of a drag): the members' grids moved (`moveInstanceGrid`), not built again;
  // the culler moves its cells instead of splitting them again, drawn as members built at the moved place
  it('moves its cells for grids moved from those it split, drawing as members built there', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  // clip objects (by pixel, without a clipping layer) cull the cells and instances they clip whole on the GPU as on the
  // CPU (`clip-culling.md`): fewer instances drawn, the same image
  it('leaves out what clip objects clip whole, drawing what the CPU culling draws', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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
    } finally {
      plugin.dispose();
    }
  });

  // the ghost of what the clip objects remove (`Renderer.renderClipGhost`) culled on the GPU into a region of its own:
  // the instances in view the clip objects reach, drawing what the CPU culling's ghost draws
  it('culls the ghost of what clip objects remove into its own region, drawing what the CPU culling draws', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
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
        renderer.clear(true);
        renderer.setViewport(0, 0, Size, Size);
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
    } finally {
      plugin.dispose();
    }
  });
});

async function createPlugin() {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: Size, height: Size },
    { backend: 'webgpu' },
  );
  await plugin.init();
  return plugin;
}

/** One chain of 16 atoms, 40 copies: one render object of 40 instances (fewer units than `renderMerged: auto` merges) */
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
