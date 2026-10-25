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
import { MergedCullKernels } from '@molstar/graphics/gl/renderable/merged-culler';
import { Scene } from '@molstar/graphics/gl/scene';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 64;
const Copies = 6;

/** Five spheres, as `copies` × `copies` instances on a grid in the xy plane, the grid's origin at `offset` */
function createMember(offset: Vec3, color: number) {
  const centers = new Float32Array([0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, -1.2, 0, 0, 0, -1.2, 0]);
  const groups = new Float32Array([0, 1, 2, 3, 4]);
  const spheres = Spheres.create(centers, groups, 5);
  const transforms = new Float32Array(Copies * Copies * 16);
  for (let i = 0; i < Copies * Copies; ++i) {
    const t = Vec3.add(Vec3(), offset, Vec3.create((i % Copies) * 6, Math.floor(i / Copies) * 6, 0));
    Mat4.toArray(Mat4.fromTranslation(Mat4(), t), transforms, i * 16);
  }
  // cells of a few instances each (the default 200 Å holds a member's whole grid)
  const props = { ...PD.getDefaultValues(Spheres.Params), cellSize: 8, batchSize: 32 };
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

      const renderer = Renderer.create(graphics);
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

      const frame = async (scene: Scene) => {
        // GPU culling dispatches before the passes, with the planes of `update`
        renderer.update(camera, scene, createFrame());
        renderer.prepareCull(scene);
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
      // the GPU the same cells: the instance counts of the argument records (a record per cell)
      const gpuList = drawn.get(gpuScene)![0];
      const { buffer, byteOffset, indexed } = gpuList.indirect!;
      // sphere quads draw indexed by the shared quad indices: indexed records (5 words)
      expect(indexed).toBe(true);
      const words = indexed ? 5 : 4;
      const args = new Uint32Array(await buffer.readAsync(byteOffset, gpuList.count * words * 4));
      let gpuInstances = 0;
      for (let i = 0; i < gpuList.count; ++i) gpuInstances += args[i * words + 1];
      expect(gpuInstances).toBe(cpuInstances);
      expect(new Set(cpuImage).size).toBeGreaterThan(2);
      expect(gpuImage).toEqual(cpuImage);

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
