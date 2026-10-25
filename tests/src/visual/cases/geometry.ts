/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { Volume } from '@molstar/model/model/volume';
import { DirectVolumeRepresentationProvider } from '@molstar/graphics/repr/volume/direct-volume';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { SliceRepresentationProvider } from '@molstar/graphics/repr/volume/slice';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type CapabilityKey } from '../framework/capabilities.js';
import { type TestView, type VisualTestCase } from '../framework/test-case.js';
import {
  cylinderStar,
  lineCages,
  meshPrimitives,
  meshTorus,
  pointSphere,
  rowTransforms,
  shapeRepr,
  sphereGrid,
  textLabels,
  textureMeshFromMesh,
} from '../scenes/primitives.js';
import { blobVolume, gpuBlobVolume, negatedVolume, segmentedBlobVolume, volumeRepr } from '../scenes/volume.js';
import { SegmentRepresentationProvider } from '@molstar/graphics/repr/volume/segment';
import { exampleStructure, globuleStructure, helixStructure, structureRepr } from '../scenes/structure.js';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { GaussianVolumeRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-volume';
import { OrderedSet } from '@molstar/core/data/int';
import { ValueCell } from '@molstar/core/util';
import { Material } from '@molstar/core/util/material';
import type { Representation } from '@molstar/graphics/repr/representation';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { StructureElement, type Structure } from '@molstar/model/model/structure';
import type { DrawIndirect } from '@molstar/graphics/gl/compute';
import { FillKernel, FillSchema } from './compute.js';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import {
  getTextureMeshSplatInstances,
  probeZeroWeightNormalization,
} from '@molstar/graphics/geo/geometry/texture-mesh/color-smoothing';
import type { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';

export const ObliqueView: TestView = { direction: [-0.4, -0.5, -1] };

/** Sphere and cylinder impostor shaders require frag depth; there is no fallback at the geometry level */
export const Impostor: CapabilityKey[] = ['fragDepth'];

const textProps: Partial<PD.Values<Text.Params>> = {
  fontQuality: 3,
  background: true,
  backgroundOpacity: 0.8,
  backgroundColor: Color(0xdddddd),
  borderWidth: 0.1,
  attachment: 'middle-center',
};

/**
 * The first half of the triangles of the texture mesh primitives: on WebGPU with all vertices as `drawCount` and half
 * of them as the count of indirect draw arguments written by a kernel, elsewhere with half of them as `drawCount`
 */
async function halfTextureMeshPrimitives(graphics: GraphicsContext) {
  const repr = await shapeRepr('texture-mesh', textureMeshFromMesh(graphics, meshPrimitives()));
  const { values } = repr.renderObjects[0];
  const half = Math.floor(values.drawCount.ref.value / 6) * 3;
  const { compute } = graphics;
  if (!compute) {
    ValueCell.update(values.drawCount, half);
    return repr;
  }
  const buffer = compute.createBuffer(16, 'draw arguments');
  const fill = compute.createKernel(FillKernel, FillSchema, {
    uCount: ValueCell.create(1),
    uOffset: ValueCell.create(0),
    uValue: ValueCell.create(half),
    data: ValueCell.create(buffer),
  });
  compute.beginComputePass({ label: 'draw arguments' });
  fill.dispatch(1);
  compute.endComputePass();
  fill.destroy();
  ValueCell.update(values.drawIndirect as ValueCell<DrawIndirect | undefined>, { buffer, byteOffset: 0 });
  return repr;
}

/** The GPU isosurface of blobs at `first`, then updated to `second` */
async function updatedIsosurface(graphics: GraphicsContext, first: number, second: number) {
  const volume = blobVolume(32);
  const repr = await volumeRepr(graphics, IsosurfaceRepresentationProvider, volume, {
    isoValue: Volume.IsoValue.absolute(first),
    tryUseGpu: true,
  });
  await repr.createOrUpdate({ isoValue: Volume.IsoValue.absolute(second) }, volume).run();
  return repr;
}

/** The atoms of the first unit (chain A of the helix structure) */
function firstUnitLoci(structure: Structure) {
  const unit = structure.units[0];
  return StructureElement.Loci(structure, [
    {
      unit,
      indices: OrderedSet.ofBounds(0 as StructureElement.UnitIndex, unit.elements.length as StructureElement.UnitIndex),
    },
  ]);
}

type SmoothedLayer = 'overpaint' | 'transparency' | 'emissive' | 'substance';

/** Per-group state on chain A, each smoothed into its own grid (`groupInstance` types) */
function smoothedLayerState(structure: Structure, layers: SmoothedLayer[]): Partial<Representation.State> {
  const loci = firstUnitLoci(structure);
  const state: Partial<Representation.State> = {};
  for (const layer of layers) {
    if (layer === 'overpaint')
      state.overpaint = { kind: 'element-loci', layers: [{ loci, color: Color(0xff3311), clear: false }] };
    if (layer === 'transparency') state.transparency = { kind: 'element-loci', layers: [{ loci, value: 0.7 }] };
    if (layer === 'emissive') state.emissive = { kind: 'element-loci', layers: [{ loci, value: 0.8 }] };
    if (layer === 'substance') {
      const material = Material({ metalness: 1, roughness: 0.2 });
      state.substance = { kind: 'element-loci', layers: [{ loci, material, clear: false }] };
    }
  }
  return state;
}

/** A Gaussian surface with smoothed colors and smoothed `layers` on chain A */
async function smoothedLayersRepr(graphics: GraphicsContext, layers: SmoothedLayer[], tryUseGpu: boolean) {
  const structure = await helixStructure(24);
  const repr = await structureRepr(graphics, GaussianSurfaceRepresentationProvider, structure, {
    resolution: 0.5,
    tryUseGpu,
    smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
  });
  repr.setState(smoothedLayerState(structure, layers));
  return repr;
}

const SmoothedLayers: SmoothedLayer[] = ['overpaint', 'transparency', 'emissive', 'substance'];

const LayerGrid = {
  overpaint: { prefix: 'Overpaint', channels: 4 },
  transparency: { prefix: 'Transparency', channels: 1 },
  emissive: { prefix: 'Emissive', channels: 1 },
  substance: { prefix: 'Substance', channels: 4 },
} as const;

/**
 * GPU layer smoothing on an instanced gaussian surface (1cbs over three unit cells, 12 instances) with the layer on one
 * chain: the grids read back with the instances left out (`getTextureMeshSplatInstances`) and with every instance drawn
 * (`meta.vertexBox` cleared) must be the same in the channels the shaders read (rgba, or red of single-channel layers)
 */
async function checkSmoothingInstances(graphics: GraphicsContext): Promise<string[]> {
  const model = await exampleStructure('1cbs_updated.cif');
  const structure = await StructureSymmetry.buildSymmetryRange(model, Vec3.create(0, 0, 0), Vec3.create(2, 0, 0)).run();
  const repr = await structureRepr(graphics, GaussianSurfaceRepresentationProvider, structure, {
    tryUseGpu: true,
    smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
  });
  const metas = repr.renderObjects.map((o) => (o.values as any).meta?.ref.value as TextureMesh['meta'] | undefined);
  for (const m of metas) await m?.extracted;
  // the point sprite path leaves out instances only once the probe has passed
  const leavesOut = graphics.compute ? true : await probeZeroWeightNormalization(graphics);

  const errors: string[] = [];
  if (repr.renderObjects.some((o) => o.type !== 'texture-mesh')) errors.push('expected texture meshes');
  const boxes = metas.map((m) => m?.vertexBox);
  if (boxes.some((b) => !b)) errors.push('expected a vertex box on every surface');

  const states = smoothedLayerState(structure, SmoothedLayers);
  const empty: Partial<Representation.State> = {
    overpaint: { kind: 'element-loci', layers: [] },
    transparency: { kind: 'element-loci', layers: [] },
    emissive: { kind: 'element-loci', layers: [] },
    substance: { kind: 'element-loci', layers: [] },
  };

  const readGrids = async (layer: SmoothedLayer) => {
    const { prefix, channels } = LayerGrid[layer];
    const grids: { id: number; drawn: number; instances: number; read: Uint8Array }[] = [];
    for (const o of repr.renderObjects) {
      const v = o.values as any;
      if (!v[`${layer}Type`].ref.value.startsWith('volume')) continue;
      const data = await graphics.readTextureData(v[`t${prefix}Grid`].ref.value);
      const bytes = new Uint8Array(data.array.buffer, data.array.byteOffset, data.array.byteLength);
      const read = new Uint8Array((bytes.length / 4) * channels);
      for (let t = 0, k = 0; t < bytes.length; t += 4) for (let c = 0; c < channels; ++c) read[k++] = bytes[t + c];
      const transform = v[`u${prefix}GridTransform`].ref.value;
      const { splatCount } = getTextureMeshSplatInstances(
        {
          instanceCount: v.uInstanceCount.ref.value,
          transformBuffer: v.aTransform.ref.value,
          instanceBuffer: v.aInstance.ref.value,
          groupCount: v.uGroupCount.ref.value,
          colorType: 'groupInstance',
          layer: { image: v[`t${prefix}`].ref.value, itemSize: channels },
          vertexBox: v.meta.ref.value.vertexBox,
        } as any,
        v[`u${prefix}GridDim`].ref.value,
        Vec3.create(transform[0], transform[1], transform[2]),
        transform[3],
      );
      grids.push({ id: o.id, drawn: splatCount, instances: v.uInstanceCount.ref.value, read });
    }
    return grids;
  };

  for (const layer of SmoothedLayers) {
    const state = { [layer]: states[layer] };
    repr.setState(state);
    const left = await readGrids(layer);
    repr.setState({ [layer]: empty[layer] });
    metas.forEach((m) => {
      if (m) m.vertexBox = undefined;
    });
    repr.setState(state);
    const all = await readGrids(layer);
    repr.setState({ [layer]: empty[layer] });
    metas.forEach((m, i) => {
      if (m) m.vertexBox = boxes[i];
    });

    if (left.length === 0) errors.push(`${layer}: no smoothed grid`);
    if (leavesOut && !left.some((g) => g.drawn < g.instances)) {
      errors.push(`${layer}: every instance drawn, the case does not leave any out`);
    }
    for (let i = 0; i < left.length; ++i) {
      const a = left[i].read;
      const b = all[i]?.read;
      let diffs = b && a.length === b.length ? 0 : Infinity;
      for (let k = 0; b && k < a.length; ++k) if (a[k] !== b[k]) diffs += 1;
      if (diffs) errors.push(`${layer}: grid ${left[i].id} differs in ${diffs} values from every instance drawn`);
    }
  }
  repr.destroy();
  return errors;
}

export const GeometryCases: VisualTestCase[] = [
  {
    id: 'geometry/mesh/primitives',
    tags: ['smoke'],
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives())],
  },
  {
    id: 'geometry/mesh/flat-shaded',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives(), { props: { flatShaded: true } })],
  },
  {
    id: 'geometry/mesh/instanced',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshTorus(1), { transforms: rowTransforms(3, 3) })],
  },
  {
    id: 'geometry/spheres/grid',
    tags: ['smoke'],
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(3, 3))],
  },
  {
    id: 'geometry/spheres/instanced',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(2, 2.5), { transforms: rowTransforms(3, 6) })],
  },
  {
    id: 'geometry/spheres/approximate',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(3, 3), { props: { approximate: true } })],
  },
  {
    id: 'geometry/spheres/lod-levels',
    description:
      'Per-level `uLod` is set via multi-draw data; everything is in the far level: every 8th sphere at 2x size.',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('spheres', sphereGrid(6, 1.5), {
        props: {
          lodLevels: [
            { minDistance: 0, maxDistance: 1, overlap: 0, stride: 1, scaleBias: 1 },
            { minDistance: 1, maxDistance: 10000000, overlap: 0, stride: 8, scaleBias: 3 },
          ],
        },
      }),
    ],
  },
  {
    id: 'geometry/cylinders/star',
    tags: ['smoke'],
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('cylinders', cylinderStar(4), { size: 0.4 })],
  },
  {
    id: 'geometry/points/square',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'square', sizeFactor: 4 } }),
    ],
  },
  {
    id: 'geometry/points/circle',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'circle', sizeFactor: 6 } }),
    ],
  },
  {
    id: 'geometry/points/fuzzy',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'fuzzy', sizeFactor: 8 } }),
    ],
  },
  {
    id: 'geometry/points/attenuation',
    description: 'Point size attenuated by the view distance.',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), {
        props: { pointStyle: 'circle', sizeFactor: 2, pointSizeAttenuation: true },
      }),
    ],
  },
  {
    id: 'geometry/lines/cages',
    tags: ['smoke'],
    view: ObliqueView,
    scene: async () => [await shapeRepr('lines', lineCages(), { props: { sizeFactor: 3 } })],
  },
  {
    id: 'geometry/text/labels',
    tags: ['smoke'],
    description: 'Font atlas is rasterized with a 2D canvas; expect browser/OS dependent snapshots.',
    scene: async () => [await shapeRepr('text', textLabels(textProps), { props: textProps })],
  },
  {
    id: 'geometry/mesh/isosurface-cpu',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, IsosurfaceRepresentationProvider, blobVolume(32), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/primitives',
    description: 'Texture mesh built on the CPU (vertex data in float textures), as GPU-computed surfaces are.',
    view: ObliqueView,
    scene: async ({ graphics }) => [await shapeRepr('texture-mesh', textureMeshFromMesh(graphics, meshPrimitives()))],
  },
  {
    id: 'geometry/texture-mesh/draw-indirect',
    description:
      'Half of the texture mesh primitives: an indirect draw of a kernel-written vertex count on WebGPU (`drawCount` covers all vertices), `drawCount` elsewhere.',
    view: ObliqueView,
    scene: async ({ graphics }) => [await halfTextureMeshPrimitives(graphics)],
  },
  {
    id: 'geometry/texture-mesh/isosurface-gpu',
    tags: ['smoke'],
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, IsosurfaceRepresentationProvider, blobVolume(32), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: true,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/isosurface-gpu-update',
    description:
      'GPU isosurface updated from iso 0.5 to 0.6: on WebGPU into the capacity of the first extraction, drawn indirectly.',
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ graphics }) => [await updatedIsosurface(graphics, 0.5, 0.6)],
  },
  {
    id: 'geometry/texture-mesh/isosurface-gpu-update-overflow',
    description:
      'GPU isosurface updated from iso 0.9 to 0.5: on WebGPU beyond the capacity of the first extraction, extracted again.',
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ graphics }) => [await updatedIsosurface(graphics, 0.9, 0.5)],
  },
  {
    id: 'geometry/mesh/segment-cpu',
    description: 'Volume segmentation (three labeled blobs), a CPU surface per segment.',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, SegmentRepresentationProvider, segmentedBlobVolume(32), {
        segments: [1, 2, 3],
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/segment-gpu',
    description: 'Volume segmentation (three labeled blobs), a GPU surface per segment (compute kernels on WebGPU).',
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, SegmentRepresentationProvider, segmentedBlobVolume(32), {
        segments: [1, 2, 3],
        tryUseGpu: true,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/isosurface-gpu-negative',
    description: 'GPU isosurface at a negative iso value (winding and normals flipped) of negated blobs.',
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, IsosurfaceRepresentationProvider, negatedVolume(blobVolume(32)), {
        isoValue: Volume.IsoValue.absolute(-0.5),
        tryUseGpu: true,
      }),
    ],
  },
  {
    id: 'geometry/mesh/grid3d-compute',
    description:
      'Volume computed with `createGrid3dComputeRenderable` (loop bound uniform, define on WebGL1; a kernel on WebGPU), CPU isosurface.',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, IsosurfaceRepresentationProvider, await gpuBlobVolume(graphics, 32, false), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/mesh/grid3d-compute-cumulative',
    description: 'As `grid3d-compute`, accumulated per blob (`tCumulativeSum`).',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, IsosurfaceRepresentationProvider, await gpuBlobVolume(graphics, 32, true), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/mesh/gaussian-surface-cpu',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(graphics, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
        tryUseGpu: false,
      }),
    ],
  },
  ...SmoothedLayers.map(
    (layer): VisualTestCase => ({
      id: `geometry/texture-mesh/gaussian-surface-gpu-smooth-${layer}`,
      description: `GPU gaussian surface with ${layer} on chain A, smoothed on the GPU (compute kernels on WebGPU).`,
      requires: ['gaussianSurfaceGpu'],
      view: ObliqueView,
      scene: async ({ graphics }) => [await smoothedLayersRepr(graphics, [layer], true)],
    }),
  ),
  {
    id: 'geometry/texture-mesh/gaussian-surface-gpu-smooth-instances',
    description:
      'GPU layer smoothing on an instanced gaussian surface, a layer on one chain: the grids with the instances the layer does not reach left out are those of every instance drawn.',
    requires: ['gaussianSurfaceGpu'],
    scene: async () => [],
    check: checkSmoothingInstances,
  },
  {
    id: 'geometry/mesh/gaussian-surface-cpu-smooth-layers',
    description:
      'CPU gaussian surface with overpaint, transparency, emissive and substance on chain A, smoothed on the CPU.',
    view: ObliqueView,
    scene: async ({ graphics }) => [await smoothedLayersRepr(graphics, SmoothedLayers, false)],
  },
  {
    id: 'geometry/mesh/gaussian-surface-cpu-smooth-colors',
    description: 'CPU gaussian surface with CPU color smoothing, read from a 2d atlas of the color grid.',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(
        graphics,
        GaussianSurfaceRepresentationProvider,
        await helixStructure(24),
        {
          resolution: 0.5,
          tryUseGpu: false,
          smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
        },
        { name: 'element-index' },
      ),
    ],
  },
  {
    id: 'geometry/direct-volume/gaussian-volume',
    description:
      'Gaussian density computed on the GPU into a 3d texture (2d atlas without 3d textures), drawn as a direct volume.',
    // with 8-bit density targets the volume is blank
    requires: ['halfFloatRenderTarget'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(graphics, GaussianVolumeRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/gaussian-surface-gpu',
    description:
      'Gaussian density computed on the GPU (`tMinDistanceTex`, 2d or 3d by `dGridTexType`), its isosurface too (compute kernels on WebGPU).',
    requires: ['gaussianSurfaceGpu'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(graphics, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
        tryUseGpu: true,
      }),
    ],
  },
  ...[2, 3].map(
    (resolution): VisualTestCase => ({
      id: `geometry/texture-mesh/gaussian-surface-gpu-res${resolution}`,
      description: `A Gaussian surface at a coarse resolution (${resolution} Å) of a ball of 6k atoms: on the GPU with compute kernels (WebGPU), on the CPU with WebGL (resolutions above 1): the GPU surface against the CPU's.`,
      requires: ['gaussianSurfaceGpu'],
      view: ObliqueView,
      scene: async ({ graphics }) => [
        await structureRepr(graphics, GaussianSurfaceRepresentationProvider, await globuleStructure(6000), {
          // the resolution as given ('auto' quality chooses one by the structure's size)
          quality: 'custom',
          resolution,
          tryUseGpu: true,
        }),
      ],
    }),
  ),
  {
    id: 'geometry/texture-mesh/gaussian-surface-gpu-smooth-colors',
    description:
      'GPU gaussian surface with GPU color smoothing (`color-accumulate`, `color-normalize`; compute kernels on WebGPU).',
    requires: ['gaussianSurfaceGpu'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(
        graphics,
        GaussianSurfaceRepresentationProvider,
        await helixStructure(24),
        {
          resolution: 0.5,
          tryUseGpu: true,
          smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
        },
        { name: 'element-index' },
      ),
    ],
  },
  {
    id: 'geometry/direct-volume/blobs',
    tags: ['smoke'],
    view: ObliqueView,
    scene: async ({ graphics }) => [await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(32), {})],
  },
  {
    id: 'geometry/direct-volume/blobs-full-resolution',
    description:
      'As blobs at the full resolution (`volumeResolutionScale` 1): the volume drawn by the scene passes, not upsampled.',
    view: ObliqueView,
    canvasProps: { volumeResolutionScale: 1 },
    scene: async ({ graphics }) => [await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(32), {})],
  },
  {
    id: 'geometry/direct-volume/opaque-edges',
    description:
      'Blobs with an opaque isosurface inside: the rays stop at its depth, and the half-resolution volume is upsampled ' +
      'by the depth at its edges.',
    view: ObliqueView,
    scene: async ({ graphics }) => {
      const volume = blobVolume(32);
      return [
        await volumeRepr(graphics, IsosurfaceRepresentationProvider, volume, {
          isoValue: Volume.IsoValue.absolute(1.15),
        }),
        await volumeRepr(graphics, DirectVolumeRepresentationProvider, volume, {}),
      ];
    },
  },
  {
    id: 'geometry/direct-volume/opaque-transfer',
    description:
      'A transfer function reaching 1: rays opaque after one sample, in the default WBOIT mode (which draws direct ' +
      'volumes also where they are opaque: they have no opaque pass).',
    contextProps: { transparency: 'wboit' },
    requires: ['wboit'],
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(48), {
        controlPoints: [Vec2.create(0.3, 0), Vec2.create(0.35, 1), Vec2.create(1, 1)],
      }),
    ],
  },
  {
    id: 'geometry/direct-volume/bricks',
    description: 'A narrow transfer function at high values: the rays skip most bricks (`createVolumeBrickRange`).',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(48), {
        controlPoints: [Vec2.create(0.7, 0), Vec2.create(0.75, 0.3), Vec2.create(0.85, 0.3), Vec2.create(0.9, 0)],
      }),
    ],
  },
  {
    id: 'geometry/image/slice',
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, SliceRepresentationProvider, blobVolume(32), {
        dimension: { name: 'relativeZ', params: 0.5 },
      }),
    ],
  },
];
