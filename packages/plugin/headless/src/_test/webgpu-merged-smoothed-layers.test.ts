/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// a smoothed layer on a merged surface: drawn as the per-group visuals' (their grids), not as the plain layer; sampled
// per vertex by a kernel into the storage the merged object draws it from
import { afterAll, describe, expect, it } from 'vitest';
import { Color } from '@molstar/core/util/color';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Interval } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { Renderer } from '@molstar/graphics/gl/renderer';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { getMergedMembers, isMergedRenderObject } from '@molstar/graphics/gl/render-object';
import { type DataArray, readGpuDataArrays } from '@molstar/graphics/gl/renderable/util';
import { SegmentedArray } from '@molstar/graphics/gl/segmented-array';
import type { WebGPUContext } from '@molstar/graphics/gl/contexts/webgpu';
import { KernelCode, type StorageBuffer, StorageBufferSpec } from '@molstar/graphics/gl/compute';
import { LayerSamplingKernels } from '@molstar/graphics/gl/compute/layer-sampling';
import { getGpuMeshVertices } from '@molstar/graphics/geo/geometry/mesh/gpu-vertices';
import { UniformSpec } from '@molstar/graphics/gl/renderable/schema';
import { ValueCell } from '@molstar/core/util';
import { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { Theme } from '@molstar/graphics/theme/theme';
import { UnitsRepresentation } from '@molstar/graphics/repr/structure/units-representation';
import {
  GaussianSurfaceMeshParams,
  GaussianSurfaceMeshVisual,
} from '@molstar/graphics/repr/structure/visual/gaussian-surface-mesh';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
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

/** one plugin for the tests that need no context of their own: its device's pipelines compiled once */
let shared: Promise<HeadlessPluginContext> | undefined;
const sharedPlugin = () => (shared ??= createPlugin());
afterAll(async () => (await shared)?.dispose());

async function setup(plugin: HeadlessPluginContext) {
  const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structureCell = await plugin.builders.structure.createStructure(model);
  const structure = structureCell.obj!.data;
  const registries = {
    colorThemeRegistry: plugin.representation.structure.themes.colorThemeRegistry,
    sizeThemeRegistry: plugin.representation.structure.themes.sizeThemeRegistry,
  };
  const unit = structure.units[0];
  // the half of chain A toward the camera (along the helices' axis): the smoothing seen at the layer's edge
  const n = unit.elements.length;
  const loci = StructureElement.Loci(structure, [{ unit, indices: Interval.ofBounds(n / 2, n) }]);
  const overpaint = { kind: 'element-loci', layers: [{ loci, color: Color(0xff0000), clear: false }] } as const;
  return { structureCell, structure, registries, overpaint };
}

/** the members of `repr`'s merged objects whose overpaint is written on the GPU (`DataArray.gpu`) */
function gpuMembers(repr: { renderObjects: readonly any[] }) {
  const members = repr.renderObjects.filter((o) => isMergedRenderObject(o)).flatMap((o) => getMergedMembers(o));
  return members.filter((m) => !!(m.values.tOverpaint.ref.value as DataArray<Uint8Array>).gpu);
}

describe.skipIf(!WebGPU.module)('smoothed layers on merged surfaces (WebGPU)', () => {
  it('draws an overpaint on one chain as the per-group visuals do', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const { structure, registries, overpaint } = await setup(plugin);
    let written = 0;
    const merged = await image(graphics, registries, structure, 'merged', overpaint, {
      inspect: async (repr) => {
        written = gpuMembers(repr).length;
      },
    });
    const perGroup = await image(graphics, registries, structure, 'per-group', overpaint);
    const plain = await image(graphics, registries, structure, 'merged', overpaint, {
      smoothColors: { name: 'off', params: {} },
    });
    const toPerGroup = differing(merged.pixels, perGroup.pixels);
    const plainToPerGroup = differing(plain.pixels, perGroup.pixels);
    // a part of its own for the smoothed chain beside the plain one, its layer written on the GPU
    expect(merged.objects).toBe(2);
    expect(written).toBe(1);
    // the plain layer's hard edge apart from the grids' soft one; merged, sampled by the same filter per vertex, near
    // them (measured 35 against 288 pixels: the vertices' interpolation of it, not the fragments')
    expect(plainToPerGroup).toBeGreaterThan(150);
    expect(toPerGroup).toBeLessThan(plainToPerGroup / 4);
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  it('reads the layer back as drawn, and samples it again for a geometry made again', async () => {
    const plugin = await sharedPlugin();
    const graphics = plugin.canvas3d!.graphics;
    const { structure, registries, overpaint } = await setup(plugin);
    let readBack: Uint8Array | undefined;
    let drawn: Uint8Array | undefined;
    // another resolution after the layer: the geometry made again, the layer's segments laid out again
    const remade = await image(graphics, registries, structure, 'merged', overpaint, {
      after: (repr, props) => repr.createOrUpdate({ ...props, resolution: 0.6 }, structure).run(),
      inspect: async (repr, scene) => {
        const [member] = gpuMembers(repr);
        const values = await readGpuDataArrays(member.values);
        readBack = (values.tOverpaint.ref.value as DataArray<Uint8Array>).array;
        // the bytes in the storage the object draws its layer from, at the member's part
        const object = repr.renderObjects.find(
          (o: any) => isMergedRenderObject(o) && gpuMembers({ renderObjects: [o] }).length,
        );
        const layer = scene.get(object)!.renderItem.getDataLayer!('tOverpaint')!;
        const array = (object.values.tOverpaint.ref.value as DataArray<any>).array as SegmentedArray;
        const part = array.parts.find((p) => !!p.gpu)!;
        drawn = await readPage(graphics, layer.buffer, layer.offset + part.offset / 4, part.count / 4);
      },
    });
    expect(readBack!.length).toBe(drawn!.length);
    expect(Array.from(readBack!)).toEqual(Array.from(drawn!));
    expect(readBack!.some((b) => b > 0)).toBe(true);
    const perGroup = await image(graphics, registries, structure, 'per-group', overpaint, { resolution: 0.6 });
    const plain = await image(graphics, registries, structure, 'merged', overpaint, {
      resolution: 0.6,
      smoothColors: { name: 'off', params: {} },
    });
    expect(differing(remade.pixels, perGroup.pixels)).toBeLessThan(differing(plain.pixels, perGroup.pixels) / 4);
    expect(graphics.takeErrors?.() ?? []).toEqual([]);
  }, 60_000);

  it(
    'draws the layer again after a device loss, on the new device',
    { tags: ['medium'] },
    async () => {
      const plugin = await createPlugin();
      try {
        const { structureCell, overpaint } = await setup(plugin);
        const cell = await plugin.builders.structure.representation.addRepresentation(structureCell, {
          type: 'gaussian-surface',
          typeParams: { resolution: 0.5, tryUseGpu: false, smoothColors: Smooth } as never,
        });
        const repr = cell!.obj!.data.repr;
        // the sampling's kernels compiled ahead of the layer, the surface smoothing
        await plugin.canvas3d!.graphics.resources.pending.whenReady();
        expect(plugin.canvas3d!.graphics.compute!.hasKernels(LayerSamplingKernels)).toBe(true);
        repr.setState({ overpaint });
        plugin.canvas3d!.update(repr);
        const before = (await plugin.getImageRaw({ width: Size, height: Size })).data;
        expect(gpuMembers(repr)).toHaveLength(1);
        const graphics = plugin.canvas3d!.graphics as WebGPUContext;
        const restored = new Promise<void>((resolve) =>
          plugin.canvas3d!.contextChanged.subscribe(({ reason }) => {
            if (reason === 'restore') resolve();
          }),
        );
        graphics.device.destroy();
        await restored;
        const after = (await plugin.getImageRaw({ width: Size, height: Size })).data;
        expect(Array.from(after)).toEqual(Array.from(before));
        // written by the new device's kernels
        const next = plugin.canvas3d!.graphics;
        const [member] = gpuMembers(repr);
        expect((member.values.tOverpaint.ref.value as DataArray<Uint8Array>).gpu!.compute).toBe(next.compute);
      } finally {
        plugin.dispose();
      }
    },
    60_000,
  );

  it('keeps the vertices on the GPU across layers, released with the representation', async () => {
    const plugin = await sharedPlugin();
    const { structureCell, overpaint } = await setup(plugin);
    const cell = await plugin.builders.structure.representation.addRepresentation(structureCell, {
      type: 'gaussian-surface',
      typeParams: { resolution: 0.5, tryUseGpu: false, smoothColors: Smooth } as never,
    });
    const repr = cell!.obj!.data.repr;
    const compute = plugin.canvas3d!.graphics.compute!;
    const draw = async (color: Color) => {
      const layers = overpaint.layers.map((l) => ({ ...l, color }));
      repr.setState({ overpaint: { ...overpaint, layers } });
      plugin.canvas3d!.update(repr);
      await plugin.getImageRaw({ width: Size, height: Size });
      const [member] = gpuMembers(repr);
      return getGpuMeshVertices(member.values as any, compute);
    };
    const first = await draw(Color(0xff0000));
    const withVertices = compute.getStorageBufferByteCount();
    // another layer: sampled from the same upload of the positions
    const second = await draw(Color(0x0000ff));
    expect(second).toBe(first);
    expect(first.positions.byteLength).toBeGreaterThan(0);
    // the representation removed: its segments' vertices with it
    await plugin.build().delete(cell!).commit();
    expect(withVertices - compute.getStorageBufferByteCount()).toBeGreaterThanOrEqual(first.positions.byteLength);
  }, 60_000);

  it(
    'prepares the sampling kernels by a precompile list naming a merged surface, in another context',
    { tags: ['medium'] },
    async () => {
      const a = await createPlugin();
      const b = await createPlugin();
      try {
        const { structureCell, overpaint } = await setup(a);
        const cell = await a.builders.structure.representation.addRepresentation(structureCell, {
          type: 'gaussian-surface',
          typeParams: { resolution: 0.5, tryUseGpu: false, smoothColors: Smooth } as never,
        });
        const repr = cell!.obj!.data.repr;
        repr.setState({ overpaint });
        a.canvas3d!.update(repr);
        await a.getImageRaw({ width: Size, height: Size });
        const list = a.canvas3d!.getPrecompileList();
        expect(list.objects.some((s) => s.merged && s.defines.dOverpaintType === 'vertexInstance')).toBe(true);

        // another context, nothing in its scene
        const other = b.canvas3d!.graphics.compute!;
        expect(other.hasKernels(LayerSamplingKernels)).toBe(false);
        await b.canvas3d!.precompile(list);
        expect(other.hasKernels(LayerSamplingKernels)).toBe(true);
      } finally {
        a.dispose();
        b.dispose();
      }
    },
    60_000,
  );
});

const CopySchema = {
  uOffset: UniformSpec('i'),
  uCount: UniformSpec('i'),
  page: StorageBufferSpec('u32', 'read'),
  out: StorageBufferSpec('u32', 'read-write'),
};
/** copies `uCount` words of `page` from `uOffset` on into `out` */
const CopyKernel = KernelCode('test-copy-smoothed-layer', [64], {
  wgsl: /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3u) {
    if (gid.x >= u32(uParams.uCount)) { return; }
    out[gid.x] = page[u32(uParams.uOffset) + gid.x];
}
`,
});

/** `words` words of a renderer's storage `page` from word `offset` (a page is not read back directly) */
async function readPage(graphics: GraphicsContext, page: StorageBuffer, offset: number, words: number) {
  const compute = graphics.compute!;
  const out = compute.createBuffer(words * 4);
  const kernel = compute.createKernel(CopyKernel, CopySchema, {
    uOffset: ValueCell.create(offset),
    uCount: ValueCell.create(words),
    page: ValueCell.create(page),
    out: ValueCell.create(out),
  });
  compute.beginComputePass({ label: 'copy' });
  kernel.dispatch(Math.ceil(words / 64));
  compute.endComputePass();
  const bytes = new Uint8Array(await out.readAsync());
  kernel.destroy();
  out.destroy();
  return bytes;
}

/** pixels of `a` and `b` that differ by more than 24 summed over rgb */
function differing(a: Uint8Array, b: Uint8Array) {
  let n = 0;
  for (let i = 0; i < a.length; i += 4) {
    if (Math.abs(a[i] - b[i]) + Math.abs(a[i + 1] - b[i + 1]) + Math.abs(a[i + 2] - b[i + 2]) > 24) n++;
  }
  return n;
}

/**
 * a two-chain Gaussian surface with `overpaint`, merged or per group (a builder-less visual: grids), drawn opaque;
 * `after` the layer is set, `inspect` once drawn
 */
async function image(
  graphics: GraphicsContext,
  registries: any,
  structure: Structure,
  mode: 'merged' | 'per-group',
  overpaint: any,
  options: {
    smoothColors?: any;
    resolution?: number;
    after?: (repr: any, props: any) => Promise<void>;
    inspect?: (repr: any, scene: Scene) => Promise<void>;
  } = {},
) {
  const visualCtor = (materialId: number) => {
    const visual = GaussianSurfaceMeshVisual(materialId);
    return mode === 'merged'
      ? visual
      : new Proxy(visual, { get: (t, k) => (k === 'builder' ? undefined : Reflect.get(t, k)) });
  };
  const repr = UnitsRepresentation(
    'test',
    { graphics, ...registries },
    () => GaussianSurfaceMeshParams as any,
    visualCtor as any,
  );
  repr.setTheme(
    Theme.create(registries, { structure }, { colorTheme: { name: 'uniform' }, sizeTheme: { name: 'physical' } }),
  );
  const props = {
    ...PD.getDefaultValues(GaussianSurfaceMeshParams),
    resolution: options.resolution ?? 0.5,
    smoothColors: options.smoothColors ?? Smooth,
  };
  await repr.createOrUpdate(props as any, structure).run();
  repr.setState({ overpaint });
  await options.after?.(repr, props);
  const scene = Scene.create(graphics);
  for (const o of repr.renderObjects) scene.add(o);
  scene.commit();
  const camera = new Camera();
  Viewport.set(camera.viewport, 0, 0, Size, Size);
  const { center, radius } = scene.boundingSphere;
  // from the side: the helices' axis (z) across the view
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
    renderer.clear(true);
    renderer.setViewport(0, 0, Size, Size);
    renderer.renderOpaque(scene, camera);
    graphics.flush();
    await graphics.resources.pending.whenReady();
    pixels = new Uint8Array(await graphics.readRenderTarget(target, 0, 0, Size, Size));
  }
  const objects = repr.renderObjects.filter((o) => isMergedRenderObject(o)).length;
  await options.inspect?.(repr, scene);
  for (const o of repr.renderObjects) scene.remove(o);
  repr.destroy();
  return { pixels, objects };
}
