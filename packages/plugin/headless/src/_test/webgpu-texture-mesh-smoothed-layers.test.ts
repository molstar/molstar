/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// a smoothed layer on a texture mesh (a GPU surface): its grid splatted on the GPU (compute kernels), drawn as the CPU
// mesh's grid, the instances it does not reach left out of the grid
import { afterAll, describe, expect, it } from 'vitest';
import { Color } from '@molstar/core/util/color';
import { Box3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Interval } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { Theme } from '@molstar/graphics/theme/theme';
import { UnitsRepresentation } from '@molstar/graphics/repr/structure/units-representation';
import {
  GaussianSurfaceMeshParams,
  GaussianSurfaceMeshVisual,
  GaussianSurfaceTextureMeshVisual,
} from '@molstar/graphics/repr/structure/visual/gaussian-surface-mesh';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

const Size = 128;
const TwoChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const chain = i < 24 ? 'A' : 'B';
  const k = i % 24;
  const [x, y, z] = [(i < 24 ? 0 : 6) + Math.cos(k * 0.6) * 2, Math.sin(k * 0.6) * 2, k * 0.6];
  return `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String(k + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const Smooth = { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } };
const Off = { name: 'off', params: {} };

async function createPlugin() {
  const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: Size, height: Size });
  await plugin.init();
  return plugin;
}

/** one plugin for the file: its device's pipelines compiled once */
let shared: Promise<HeadlessPluginContext> | undefined;
const sharedPlugin = () => (shared ??= createPlugin());
afterAll(async () => (await shared)?.dispose());

async function setup(plugin: HeadlessPluginContext) {
  const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = (await plugin.builders.structure.createStructure(model)).obj!.data;
  const registries = {
    colorThemeRegistry: plugin.representation.structure.themes.colorThemeRegistry,
    sizeThemeRegistry: plugin.representation.structure.themes.sizeThemeRegistry,
  };
  return { structure, registries };
}

/** an overpaint on the half of `unit` toward the camera (along the helices' axis): the smoothing seen at its edge */
function overpaintOf(structure: Structure, unit: Structure['units'][number]) {
  const n = unit.elements.length;
  const loci = StructureElement.Loci(structure, [{ unit, indices: Interval.ofBounds(n / 2, n) }]);
  return { kind: 'element-loci', layers: [{ loci, color: Color(0xff0000), clear: false }] } as const;
}

/** pixels of `a` and `b` that differ by more than 24 summed over rgb */
function differing(a: Uint8Array, b: Uint8Array) {
  let n = 0;
  for (let i = 0; i < a.length; i += 4) {
    if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 24) n++;
  }
  return n;
}

/** a Gaussian surface of `structure` with `overpaint`, a texture mesh or a CPU mesh per group, its visual's props */
async function createRepr(
  graphics: GraphicsContext,
  registries: any,
  structure: Structure,
  geometry: 'texture-mesh' | 'mesh',
  overpaint: any,
  smoothColors: any,
) {
  // the visuals themselves, not the representation's choice by size (these surfaces are small: on the CPU); the mesh
  // per group (without its merge builder), each grid sampled by the fragments as the texture mesh's
  const visualCtor = (materialId: number) =>
    geometry === 'texture-mesh'
      ? GaussianSurfaceTextureMeshVisual(materialId)
      : new Proxy(GaussianSurfaceMeshVisual(materialId), {
          get: (t, k) => (k === 'builder' ? undefined : Reflect.get(t, k)),
        });
  const repr = UnitsRepresentation(
    'test',
    { graphics, ...registries },
    () => GaussianSurfaceMeshParams as any,
    visualCtor as any,
  );
  repr.setTheme(
    Theme.create(registries, { structure }, { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'physical' } }),
  );
  const props = { ...PD.getDefaultValues(GaussianSurfaceMeshParams), resolution: 0.5, smoothColors };
  await repr.createOrUpdate(props as any, structure).run();
  repr.setState({ overpaint });
  return repr;
}

/** `repr` drawn opaque from the side (the helices' axis across the view) */
async function image(graphics: GraphicsContext, repr: ReturnType<typeof UnitsRepresentation>) {
  const scene = Scene.create(graphics);
  for (const o of repr.renderObjects) scene.add(o);
  scene.commit();
  const camera = new Camera();
  Viewport.set(camera.viewport, 0, 0, Size, Size);
  const { center, radius } = scene.boundingSphere;
  camera.setState(camera.getFocus(center, radius, Vec3.create(0, 0, 1), Vec3.create(0, 1, 0)), 0);
  camera.update();
  const renderer = Renderer.create(graphics);
  const target = graphics.createRenderTarget(Size, Size, 'depth-stencil');
  await graphics.resources.preparePipelines(
    renderer.collect({ name: 'opaque' }, scene, graphics.describeTarget(target)),
  );
  let pixels = new Uint8Array(0);
  for (let i = 0; i < 2; ++i) {
    renderer.update(camera, scene, createFrame());
    graphics.beginRenderPass(target);
    renderer.setViewport(0, 0, Size, Size);
    renderer.clear(true);
    renderer.renderOpaque(scene, camera);
    graphics.flush();
    await graphics.resources.pending.whenReady();
    pixels = new Uint8Array(await graphics.readRenderTarget(target, 0, 0, Size, Size));
  }
  for (const o of repr.renderObjects) scene.remove(o);
  return pixels;
}

describe.skipIf(!WebGPU.module)('smoothed layers on texture meshes (WebGPU)', () => {
  it(
    'draws an overpaint on one chain as the CPU mesh does',
    { tags: ['medium'] },
    async () => {
      const plugin = await sharedPlugin();
      const graphics = plugin.canvas3d!.graphics;
      const { structure, registries } = await setup(plugin);
      const overpaint = overpaintOf(structure, structure.units[0]);
      const draw = async (geometry: 'texture-mesh' | 'mesh', smoothColors: any) => {
        const repr = await createRepr(graphics, registries, structure, geometry, overpaint, smoothColors);
        const types = repr.renderObjects.map((o) => o.type);
        const values = repr.renderObjects[0].values as any;
        const layerType = values.dOverpaintType.ref.value;
        const pixels = await image(graphics, repr);
        repr.destroy();
        return { types, layerType, pixels };
      };
      const smoothed = await draw('texture-mesh', Smooth);
      const cpu = await draw('mesh', Smooth);
      const plain = await draw('texture-mesh', Off);
      const toCpu = differing(smoothed.pixels, cpu.pixels);
      const plainToCpu = differing(plain.pixels, cpu.pixels);
      // a texture mesh a group, its layer a grid (per instance)
      expect(smoothed.types).toEqual(['texture-mesh', 'texture-mesh']);
      expect(smoothed.layerType).toBe('volumeInstance');
      expect(cpu.layerType).toBe('volumeInstance');
      expect(plain.layerType).not.toBe('volumeInstance');
      // the plain layer's hard edge apart from the CPU grid's soft one, the GPU grid nearer it (measured 134 to 171
      // against 283 pixels, from run to run; without a layer the surfaces 32 apart: marching cubes on the GPU's
      // density). Not closer: the grids' lattices differ (the CPU's on the bounding box's cells), and the GPU's samples
      // every ninth vertex of the arena (never indexed), whose order is its kernels' appends
      expect(plainToCpu).toBeGreaterThan(200);
      expect(toCpu).toBeLessThan(plainToCpu * 0.75);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    },
    60_000,
  );

  it('leaves the instances the layer does not reach out of the grid', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const { structure: single, registries } = await setup(plugin);
    // three copies 40 Å apart, a group of three instances a chain
    const transforms = [0, 1, 2].map((i) => Mat4.fromTranslation(Mat4(), Vec3.create(i * 40, 0, 0)));
    const structure = Structure.create([...Structure.instances(single, transforms).units]);
    const reached = structure.units[0];
    const repr = await createRepr(
      graphics,
      registries,
      structure,
      'texture-mesh',
      overpaintOf(structure, reached),
      Smooth,
    );
    const object = repr.renderObjects.find((o) => (o.values as any).dOverpaintType.ref.value === 'volumeInstance')!;
    const values = object.values as any;
    expect(values.instanceCount.ref.value).toBe(3);
    const gridDim: Vec3 = values.uOverpaintGridDim.ref.value;
    const texDim = values.uOverpaintGridTexDim.ref.value;
    const transform = values.uOverpaintGridTransform.ref.value;
    const { array } = await graphics.readTextureData(values.tOverpaintGrid.ref.value);
    // the reached unit's atoms in grid cells, by the surface's reach (radius, offset, smoothness: 5 Å) and the
    // splat's (2 cells) beyond them
    const scale = transform[3];
    // (the unit's boundary is its untransformed elements')
    const box = Box3D.transform(Box3D(), reached.boundary.box, reached.conformation.operator.matrix);
    Box3D.expand(box, box, Vec3.create(5, 5, 5));
    const lo = [0, 1, 2].map((k) => (box.min[k] - transform[k]) * scale - 2);
    const hi = [0, 1, 2].map((k) => (box.max[k] - transform[k]) * scale + 2);
    const columns = Math.floor(texDim[0] / gridDim[0]);
    let inside = 0;
    let outside = 0;
    for (let z = 0; z < gridDim[2]; ++z) {
      const x0 = (z % columns) * gridDim[0];
      const y0 = Math.floor(z / columns) * gridDim[1];
      for (let y = 0; y < gridDim[1]; ++y) {
        for (let x = 0; x < gridDim[0]; ++x) {
          if (array[((y0 + y) * texDim[0] + x0 + x) * 4 + 3] === 0) continue;
          const cell = [x, y, z];
          if ([0, 1, 2].every((k) => cell[k] >= lo[k] && cell[k] <= hi[k])) inside++;
          else outside++;
        }
      }
    }
    expect(inside).toBeGreaterThan(0);
    expect(outside).toBe(0);
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
    repr.destroy();
  }, 60_000);
});
