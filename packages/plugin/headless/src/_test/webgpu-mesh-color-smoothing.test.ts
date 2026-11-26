/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// a mesh's layer splatted into its grid on the GPU (`computeMeshColorSmoothing`, the mesh's vertices on the GPU) as
// the CPU splats it (`calcMeshColorSmoothing` without them): the same grid within a byte, nothing outside the region
import { describe, expect, it } from 'vitest';
import { Color } from '@molstar/core/util/color';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Interval } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, StructureElement } from '@molstar/model/model/structure';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import { getInstanceMatrices } from '@molstar/graphics/geo/geometry/transform-data';
import { applyMeshColorSmoothing, calcMeshColorSmoothing } from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import { type DataArray, readGpuDataArrays } from '@molstar/graphics/gl/renderable/util';
import { getGpuMeshVertices, releaseGpuMeshVertices } from '@molstar/graphics/geo/geometry/mesh/gpu-vertices';
import { Theme } from '@molstar/graphics/theme/theme';
import { UnitsRepresentation } from '@molstar/graphics/repr/structure/units-representation';
import {
  GaussianSurfaceMeshParams,
  GaussianSurfaceMeshVisual,
} from '@molstar/graphics/repr/structure/visual/gaussian-surface-mesh';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

const TwoChainsPdb = `${Array.from({ length: 48 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const chain = i < 24 ? 'A' : 'B';
  const k = i % 24;
  const [x, y, z] = [(i < 24 ? 0 : 6) + Math.cos(k * 0.6) * 2, Math.sin(k * 0.6) * 2, k * 0.6];
  return `ATOM  ${String(i + 1).padStart(5)}  CA  ALA ${chain}${String(k + 1).padStart(4)}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

/** three copies 30 Å apart, a group of three instances a chain; layers on the half of one chain of the middle copy */
async function setup(plugin: HeadlessPluginContext) {
  const data = await plugin.builders.data.rawData({ data: TwoChainsPdb, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const single = (await plugin.builders.structure.createStructure(model)).obj!.data;
  const transforms = [0, 1, 2].map((i) => Mat4.fromTranslation(Mat4(), Vec3.create(i * 30, 0, 0)));
  const structure = Structure.create([...Structure.instances(single, transforms).units]);
  const registries = {
    colorThemeRegistry: plugin.representation.structure.themes.colorThemeRegistry,
    sizeThemeRegistry: plugin.representation.structure.themes.sizeThemeRegistry,
  };
  // the first chain's middle copy, by its place (the units in the order of their ids: the identity's copies first)
  const x = (u: (typeof structure.units)[number]) => u.conformation.x(u.elements[0]);
  const copies = structure.units
    .filter((u) => u.invariantId === structure.units[0].invariantId)
    .sort((a, b) => x(a) - x(b));
  const unit = copies[1];
  const n = unit.elements.length;
  const loci = StructureElement.Loci(structure, [{ unit, indices: Interval.ofBounds(n / 2, n) }]);
  return { structure, registries, loci };
}

/**
 * the per-group mesh visuals (no merge builder) with the layers, plain (`groupInstance`: their values the input) or
 * smoothed (`smoothColors`)
 */
async function meshValues(
  graphics: GraphicsContext,
  registries: any,
  structure: Structure,
  loci: any,
  smoothColors: object = { name: 'off' },
) {
  const visualCtor = (materialId: number) =>
    new Proxy(GaussianSurfaceMeshVisual(materialId), {
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
  const created = graphics.compute!.getStorageBufferByteCount();
  repr.setState({
    overpaint: { kind: 'element-loci', layers: [{ loci, color: Color(0xff8020), clear: false }] },
    transparency: { kind: 'element-loci', layers: [{ loci, value: 0.6 }] },
  } as any);
  const values = repr.renderObjects
    .map((o) => o.values as unknown as MeshValues)
    .find(
      (v) => v.overpaintType.ref.value === ((smoothColors as any).name === 'off' ? 'groupInstance' : 'volumeInstance'),
    )!;
  return { repr, values, created };
}

/** `calcMeshColorSmoothing` of `layer` on the mesh of `values`, its grid read back (four bytes a texel) */
async function grid(graphics: GraphicsContext, values: MeshValues, layer: 'overpaint' | 'transparency', gpu: boolean) {
  const itemSize = layer === 'overpaint' ? 4 : 1;
  const result = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: getInstanceMatrices(values),
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: layer === 'overpaint' ? values.tOverpaint.ref.value : values.tTransparency.ref.value,
      colorType: 'groupInstance',
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize,
      gpuVertices: gpu ? getGpuMeshVertices(values, graphics.compute!) : undefined,
    },
    { resolution: 1, stride: 3 },
    graphics,
  );
  if (result.kind !== 'volume') throw new Error('expected a grid');
  const { array } = await graphics.readTextureData(result.texture);
  result.texture.destroy();
  return { ...result, array: array as Uint8Array };
}

/**
 * A per-group surface of `structure` colored by `colorTheme`, its colors not smoothed: the input
 * `applyMeshColorSmoothing` takes (by group or by group and instance)
 */
async function colorValues(graphics: GraphicsContext, registries: any, structure: Structure, colorTheme: string) {
  const visualCtor = (materialId: number) =>
    new Proxy(GaussianSurfaceMeshVisual(materialId), {
      get: (t, k) => (k === 'builder' ? undefined : Reflect.get(t, k)),
    });
  const repr = UnitsRepresentation(
    'test',
    { graphics, ...registries },
    () => GaussianSurfaceMeshParams as any,
    visualCtor as any,
  );
  repr.setTheme(
    Theme.create(registries, { structure }, { colorTheme: { name: colorTheme }, sizeTheme: { name: 'uniform' } }),
  );
  const props = { ...PD.getDefaultValues(GaussianSurfaceMeshParams), resolution: 0.5, smoothColors: { name: 'off' } };
  await repr.createOrUpdate(props as any, structure).run();
  return { repr, values: repr.renderObjects[0].values as unknown as MeshValues };
}

describe.skipIf(!WebGPU.module)('mesh color smoothing on the GPU (WebGPU)', () => {
  it('smooths colors per vertex on the GPU as the CPU does, by group and by group and instance', async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: 16, height: 16 });
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const { structure, registries } = await setup(plugin);
      for (const [theme, input, output] of [
        ['sequence-id', 'group', 'vertex'],
        ['element-index', 'groupInstance', 'vertexInstance'],
      ] as const) {
        const cpu = await colorValues(graphics, registries, structure, theme);
        const gpu = await colorValues(graphics, registries, structure, theme);
        expect(cpu.values.colorType.ref.value, theme).toBe(input);
        const options = { resolution: 1, stride: 3 };
        // per vertex on the CPU (no context: a merged surface's before), on the GPU (`perVertexOnGpu`)
        applyMeshColorSmoothing(cpu.values, options);
        applyMeshColorSmoothing(gpu.values, options, graphics, undefined, true);
        expect(cpu.values.colorType.ref.value, theme).toBe(output);
        expect(gpu.values.colorType.ref.value, theme).toBe(output);
        const data = gpu.values.tColor.ref.value as DataArray<Uint8Array>;
        expect(data.gpu, theme).toBeDefined();
        const expected = (cpu.values.tColor.ref.value as DataArray<Uint8Array>).array;
        const count = (cpu.values.tColor.ref.value as DataArray<Uint8Array>).count;
        expect(data.count, theme).toBe(count);
        const read = ((await readGpuDataArrays(gpu.values)).tColor.ref.value as DataArray<Uint8Array>).array;
        // within a byte but where a cell's weight is near zero (the GPU's fixed point a part of it): a few components
        let max = 0;
        let off = 0;
        for (let i = 0; i < count * 4; ++i) {
          const d = Math.abs(read[i] - expected[i]);
          max = Math.max(max, d);
          if (d > 2) off++;
        }
        expect(new Set(expected.subarray(0, count * 4)).size, theme).toBeGreaterThan(4);
        expect(off, theme).toBeLessThanOrEqual(Math.ceil((count * 4) / 1000));
        expect(max, theme).toBeLessThanOrEqual(8);
        cpu.repr.destroy();
        gpu.repr.destroy();
      }
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it('splats a layer into the grid the CPU splats, rgba8 and r8, only the reached region', async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: 16, height: 16 });
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const { structure, registries, loci } = await setup(plugin);
      const { repr, values } = await meshValues(graphics, registries, structure, loci);
      expect(values.uInstanceCount.ref.value).toBe(3);
      for (const layer of ['overpaint', 'transparency'] as const) {
        const cpu = await grid(graphics, values, layer, false);
        const gpu = await grid(graphics, values, layer, true);
        expect(Array.from(gpu.gridDim)).toEqual(Array.from(cpu.gridDim));
        expect(Array.from(gpu.gridTransform)).toEqual(Array.from(cpu.gridTransform));
        expect(gpu.splatCount).toBe(cpu.splatCount);
        // the middle copy and the neighbours its cells overlap, not every instance
        expect(cpu.splatCount).toBeLessThan(3);
        // channels of the layer: red only for single-channel ones (the readers read red)
        const channels = layer === 'overpaint' ? 4 : 1;
        let max = 0;
        let off = 0;
        let nonZero = 0;
        let gpuOnly = 0;
        for (let i = 0; i < cpu.array.length; i += 4) {
          let d = 0;
          for (let c = 0; c < channels; ++c) d = Math.max(d, Math.abs(gpu.array[i + c] - cpu.array[i + c]));
          max = Math.max(max, d);
          if (d > 1) off++;
          if (cpu.array[i] > 0) nonZero++;
          else if (gpu.array[i] > 1) gpuOnly++;
        }
        // within a byte (the GPU's colors unorm, its sums fixed point) but where a cell's weight is near zero (a
        // vertex at the splat's radius): the fixed point's step (1 / 16384) a part of it (measured: 1 cell of 27k off
        // by 4, its color uniform)
        expect(nonZero).toBeGreaterThan(100);
        expect(gpuOnly).toBe(0);
        expect(off).toBeLessThanOrEqual(Math.ceil(nonZero / 1000));
        expect(max).toBeLessThanOrEqual(8);
      }
      releaseGpuMeshVertices(values);
      repr.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it("splats a per-group visual's layer on the GPU, its vertices released with the visual", async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: 16, height: 16 });
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const compute = graphics.compute!;
      const { structure, registries, loci } = await setup(plugin);
      const smooth = { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } };
      const { repr, values, created } = await meshValues(graphics, registries, structure, loci, smooth);
      expect(values).toBeDefined();
      // its positions and groups on the GPU (16 bytes a vertex), kept by its position cell
      const vertexBytes = values.uVertexCount.ref.value * 16;
      expect(compute.getStorageBufferByteCount() - created).toBeGreaterThanOrEqual(vertexBytes);
      const { array } = await graphics.readTextureData(values.tOverpaintGrid.ref.value);
      expect((array as Uint8Array).some((b) => b > 0)).toBe(true);
      const before = compute.getStorageBufferByteCount();
      repr.destroy();
      expect(before - compute.getStorageBufferByteCount()).toBeGreaterThanOrEqual(vertexBytes);
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
