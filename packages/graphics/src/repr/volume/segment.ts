/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Grid, Volume } from '@molstar/model/model/volume';
import type { VisualContext } from '../visual.js';
import type { Theme, ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { VolumeRepresentation, VolumeRepresentationProvider } from './representation.js';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import { logGpuFallback, type VisualUpdateState } from '../util.js';
import { type RepresentationContext, type RepresentationParamsGetter, Representation } from '../representation.js';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import { Mat4, Tensor, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { fillSerial } from '@molstar/core/util/array';
import {
  createSegmentTexture2d,
  eachVolumeLoci,
  getSegmentMaskValue,
  getVolumeGpuFallbackReason,
  getVolumeTexture2dLayout,
} from './util.js';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { ValueCell } from '@molstar/core/util/value-cell';
import { getAtlasIsosurfaceExtraction } from '@molstar/graphics/gl/compute/marching-cubes/atlas';
import { isComputeIsosurfaceSupported } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import type { ComputeContext } from '@molstar/graphics/gl/compute';
import {
  computeTextureMeshIsosurface,
  extractIsosurface,
} from '@molstar/graphics/geo/geometry/texture-mesh/compute-isosurface';
import { Box3D } from '@molstar/core/math/geometry/primitives/box3d';
import { SortedArray } from '@molstar/core/data/int/sorted-array';
import { Interval } from '@molstar/core/data/int/interval';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { type VolumeKey, VolumeVisual } from './visual.js';
import { selectVisual } from '../visual-factory.js';
import { estimateIsosurfaceVertexCount } from '@molstar/graphics/geo/util/marching-cubes/estimate';

export const VolumeSegmentParams = {
  segments: PD.Converted(
    (v: number[]) => v.map((x) => `${x}`),
    (v: string[]) => v.map((x) => parseInt(x)),
    PD.MultiSelect<string>(['0'], PD.arrayToOptions(['0']), {
      isEssential: true,
    }),
  ),
};
export type VolumeSegmentParams = typeof VolumeSegmentParams;
export type VolumeSegmentProps = PD.Values<VolumeSegmentParams>;

const Padding = 1;

const _translate = Mat4();
function getSegmentTransform(grid: Grid, segmentBox: Box3D) {
  const transform = Grid.getGridToCartesianTransform(grid);
  const translate = Mat4.fromTranslation(_translate, segmentBox.min);
  return Mat4.mul(Mat4(), transform, translate);
}

/** Why not on the GPU by the props and the context alone, see `getGpuFallbackReason` */
function getPropsGpuFallbackReason(props: PD.Values<SegmentMeshParams>, graphics?: GraphicsContext) {
  if (!props.tryUseGpu) return 'tryUseGpu is off';
  if (!graphics) return 'no graphics context';
}

/** GPU segments: the WebGL histogram pyramid, or compute kernels where available (WebGPU) */
function getGpuFallbackReason(volume: Volume, props: PD.Values<SegmentMeshParams>, graphics?: GraphicsContext) {
  const reason = getPropsGpuFallbackReason(props, graphics);
  if (reason || !graphics) return reason;
  return getVolumeGpuFallbackReason(volume, graphics, Padding, ['floatRenderTarget', 'multipleRenderTargets']);
}

function useGpu(volume: Volume, props: PD.Values<SegmentMeshParams>, graphics?: GraphicsContext): boolean {
  const reason = getGpuFallbackReason(volume, props, graphics);
  if (reason && props.tryUseGpu) logGpuFallback('Volume segment', reason);
  return !reason;
}

export const SegmentVisual = selectVisual(
  { gpu: SegmentTextureMeshVisual, cpu: SegmentMeshVisual },
  ({ volume }: VolumeKey, props: PD.Values<SegmentMeshParams>, graphics?: GraphicsContext) =>
    useGpu(volume, props, graphics) ? 'gpu' : 'cpu',
  (props, graphics) => (getPropsGpuFallbackReason(props, graphics) ? ['cpu'] : ['gpu', 'cpu']),
);

function getLoci(volume: Volume, props: VolumeSegmentProps) {
  const segments = SortedArray.ofUnsortedArray<Volume.SegmentIndex>(props.segments);
  const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
  return Volume.Segment.Loci(volume, [{ segments, instances }]);
}

function getSegmentLoci(pickingId: PickingId, volume: Volume, key: number, props: VolumeSegmentProps, id: number) {
  const { objectId, groupId, instanceId } = pickingId;

  if (id === objectId) {
    const granularity = Volume.PickingGranularity.get(volume);
    const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
    if (granularity === 'volume') {
      return Volume.Loci(volume, instances);
    } else if (granularity === 'object' || groupId === PickingId.Null) {
      const segments = OrderedSet.ofSingleton(key as Volume.SegmentIndex);
      return Volume.Segment.Loci(volume, [{ segments, instances }]);
    } else {
      const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
      return Volume.Cell.Loci(volume, [{ indices, instances }]);
    }
  }
  return EmptyLoci;
}

export function eachSegment(
  loci: Loci,
  volume: Volume,
  key: number,
  props: VolumeSegmentProps,
  apply: (interval: Interval) => boolean,
) {
  const segments = SortedArray.ofSingleton(key);
  return eachVolumeLoci(loci, volume, { segments }, apply);
}

//

function getSegmentCells(set: number[], bbox: Box3D, cells: Tensor): Tensor {
  const dim = Box3D.size(Vec3(), bbox);
  const [xn, yn, zn] = dim;
  const [minx, miny, minz] = bbox.min;
  const value = getSegmentMaskValue(cells, set, bbox);

  const axisOrder = [...cells.space.axisOrderSlowToFast];
  const segmentSpace = Tensor.Space(dim, axisOrder, Uint8Array);
  const segmentCells = Tensor.create(segmentSpace, segmentSpace.create());

  const segData = segmentCells.data;
  const segSet = segmentSpace.set;

  for (let z = 0; z < zn; ++z) {
    for (let y = 0; y < yn; ++y) {
      for (let x = 0; x < xn; ++x) {
        segSet(segData, x, y, z, value(x + minx, y + miny, z + minz));
      }
    }
  }

  return segmentCells;
}

export async function createVolumeSegmentMesh(
  ctx: VisualContext,
  volume: Volume,
  key: Volume.SegmentIndex,
  theme: Theme,
  props: VolumeSegmentProps,
  mesh?: Mesh,
) {
  const segmentation = Volume.Segmentation.get(volume);
  if (!segmentation) throw new Error('missing volume segmentation');

  ctx.runtime.update({ message: 'Marching cubes...' });

  const bbox = Box3D.clone(segmentation.bounds[key]);
  Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));

  const set = Array.from(segmentation.segments.get(key)!.values());
  const cells = getSegmentCells(set, bbox, volume.grid.cells);
  const ids = fillSerial(new Int32Array(cells.data.length));

  const surface = await computeMarchingCubesMesh(
    {
      isoLevel: 128,
      scalarField: cells,
      idField: Tensor.create(cells.space, Tensor.Data1(ids)),
    },
    mesh,
  ).runAsChild(ctx.runtime);

  const transform = getSegmentTransform(volume.grid, bbox);
  Mesh.transform(surface, transform);
  if (ctx.graphics && !ctx.graphics.capabilities.flatInterpolation) {
    // 2nd arg means not to split triangles based on group id. Splitting triangles
    // is too expensive if each cell has its own group id as is the case here.
    Mesh.uniformTriangleGroup(surface, false);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
  const segments = OrderedSet.ofSingleton(key as Volume.SegmentIndex);
  surface.setBoundingSphere(Volume.Segment.getBoundingSphere(volume, [{ segments, instances }]));

  return surface;
}

export const SegmentMeshParams = {
  ...Mesh.Params,
  ...TextureMesh.Params,
  ...VolumeSegmentParams,
  quality: { ...Mesh.Params.quality, isEssential: false },
  tryUseGpu: PD.Boolean(true),
};
export type SegmentMeshParams = typeof SegmentMeshParams;

export function SegmentMeshVisual(materialId: number): VolumeVisual<SegmentMeshParams> {
  return VolumeVisual<Mesh, SegmentMeshParams>(
    {
      defaultProps: PD.getDefaultValues(SegmentMeshParams),
      createGeometry: createVolumeSegmentMesh,
      createLocationIterator: (volume: Volume, key: number) => {
        const l = Volume.Segment.Location(volume, key);
        return LocationIterator(
          volume.grid.cells.data.length,
          volume.instances ? volume.instances.length : 1,
          1,
          () => l,
        );
      },
      getLoci: getSegmentLoci,
      eachLocation: eachSegment,
      setUpdateState: (
        state: VisualUpdateState,
        newVolume: Volume,
        currentVolume: Volume,
        newProps: PD.Values<SegmentMeshParams>,
        currentProps: PD.Values<SegmentMeshParams>,
      ) => {},
      geometryUtils: Mesh.Utils,
    },
    materialId,
  );
}

//

const SegmentTextureName = 'segment-texture';

function getSegmentTexture(volume: Volume, segment: Volume.SegmentIndex, graphics: GraphicsContext) {
  const segmentation = Volume.Segmentation.get(volume);
  if (!segmentation) throw new Error('missing volume segmentation');

  const { resources } = graphics;

  const bbox = Box3D.clone(segmentation.bounds[segment]);
  Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));

  const transform = getSegmentTransform(volume.grid, bbox);
  const gridDimension = Box3D.size(Vec3(), bbox);
  const { width, height, powerOfTwoSize: texDim } = getVolumeTexture2dLayout(gridDimension, Padding);
  const gridTexDim = Vec3.create(width, height, 0);
  const gridDataDim = Vec3.clone(gridDimension);
  const gridTexScale = Vec2.create(width / texDim, height / texDim);
  // console.log({ texDim, width, height, gridDimension });

  if (texDim > graphics.resources.limits.maxTextureSize / 2) {
    throw new Error('volume too large for gpu segment extraction');
  }

  if (!graphics.namedTextures[SegmentTextureName]) {
    graphics.namedTextures[SegmentTextureName] = resources.texture('image-uint8', 'r', 'ubyte', 'linear');
  }
  const texture = graphics.namedTextures[SegmentTextureName];

  texture.define(texDim, texDim);
  // load volume into sub-section of texture
  const set = Array.from(segmentation.segments.get(segment)!.values());
  texture.load(createSegmentTexture2d(volume, set, bbox, Padding), true);

  gridDimension[0] += Padding;
  gridDimension[1] += Padding;

  return {
    texture,
    transform,
    gridDimension,
    gridTexDim,
    gridDataDim,
    gridTexScale,
  };
}

/**
 * Marching cubes with compute kernels on the smoothed segment mask (as the CPU path, iso 128 of 0–255, groups the cells
 * of the segment box), drawn indirectly
 */
function createVolumeSegmentComputeMesh(
  graphics: GraphicsContext & { compute: ComputeContext },
  volume: Volume,
  segment: Volume.SegmentIndex,
  textureMesh?: TextureMesh,
) {
  const segmentation = Volume.Segmentation.get(volume);
  if (!segmentation) throw new Error('missing volume segmentation');
  const bbox = Box3D.clone(segmentation.bounds[segment]);
  Box3D.expand(bbox, bbox, Vec3.create(2, 2, 2));
  const set = Array.from(segmentation.segments.get(segment)!.values());
  const cells = getSegmentCells(set, bbox, volume.grid.cells);
  const values = Float32Array.from(cells.data as ArrayLike<number>);
  const { space } = cells;
  const strides = Vec3.create(space.dataOffset(1, 0, 0), space.dataOffset(0, 1, 0), space.dataOffset(0, 0, 1));
  const gridDim = Vec3.clone(space.dimensions as Vec3);
  const transform = getSegmentTransform(volume.grid, bbox);

  const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
  const segments = OrderedSet.ofSingleton(segment);
  const boundingSphere = Volume.Segment.getBoundingSphere(volume, [{ segments, instances }]);
  return computeTextureMeshIsosurface(
    graphics,
    async (target) => {
      const buffer = graphics.compute.createBuffer(values.byteLength, 'segment mask');
      buffer.write(values);
      try {
        const input = { kind: 'buffer', volume: buffer, strides } as const;
        return await extractIsosurface(graphics, input, gridDim, 128, transform, false, target);
      } finally {
        buffer.destroy();
      }
    },
    volume.grid.cells.data.length,
    boundingSphere,
    textureMesh,
    () => estimateIsosurfaceVertexCount(values, gridDim, strides, 128),
  );
}

async function createVolumeSegmentTextureMesh(
  ctx: VisualContext,
  volume: Volume,
  segment: Volume.SegmentIndex,
  theme: Theme,
  props: VolumeSegmentProps,
  textureMesh?: TextureMesh,
) {
  const graphics = ctx.graphics;
  if (volume.grid.cells.data.length <= 1) {
    return TextureMesh.createEmpty(textureMesh);
  }
  if (graphics && isComputeIsosurfaceSupported(graphics)) {
    return createVolumeSegmentComputeMesh(graphics, volume, segment, textureMesh);
  }
  const atlas = graphics && getAtlasIsosurfaceExtraction(graphics);
  if (!graphics || !atlas) {
    throw new Error('compute kernels or fragment-shader marching cubes required to create a volume segment');
  }

  const { texture, gridDimension, gridTexDim, gridDataDim, gridTexScale, transform } = getSegmentTexture(
    volume,
    segment,
    graphics,
  );

  const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
  const buffer = textureMesh?.doubleBuffer.get(graphics);
  const gv = atlas.extract(
    graphics,
    {
      texture: texture,
      channel: 'red',
      gridDim: gridDimension,
      gridTexDim: gridTexDim,
      gridDataDim: gridDataDim,
      gridTexScale: gridTexScale,
      transform: transform,
      isoValue: 0.5,
      invert: false,
      packedGroup: false,
      axisOrder: axisOrder,
      constantGroup: true,
    },
    buffer,
  );

  const groupCount = volume.grid.cells.data.length;
  const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
  const segments = OrderedSet.ofSingleton(segment as Volume.SegmentIndex);
  const boundingSphere = Volume.Segment.getBoundingSphere(volume, [{ segments, instances }]);
  const surface = TextureMesh.create(
    gv.vertexCount,
    groupCount,
    gv.vertexTexture,
    gv.groupTexture,
    gv.normalTexture,
    boundingSphere,
    textureMesh,
  );

  return surface;
}

export function SegmentTextureMeshVisual(materialId: number): VolumeVisual<SegmentMeshParams> {
  return VolumeVisual<TextureMesh, SegmentMeshParams>(
    {
      defaultProps: PD.getDefaultValues(SegmentMeshParams),
      createGeometry: createVolumeSegmentTextureMesh,
      createLocationIterator: (volume: Volume, segment: number) => {
        const l = Volume.Segment.Location(volume, segment);
        return LocationIterator(
          volume.grid.cells.data.length,
          volume.instances ? volume.instances.length : 1,
          1,
          () => l,
        );
      },
      getLoci: getSegmentLoci,
      eachLocation: eachSegment,
      setUpdateState: (
        state: VisualUpdateState,
        newVolume: Volume,
        currentVolume: Volume,
        newProps: PD.Values<SegmentMeshParams>,
        currentProps: PD.Values<SegmentMeshParams>,
      ) => {},
      geometryUtils: TextureMesh.Utils,
      dispose: (geometry: TextureMesh) => {
        geometry.vertexTexture.ref.value.destroy();
        geometry.groupTexture.ref.value.destroy();
        geometry.normalTexture.ref.value.destroy();
        geometry.doubleBuffer.destroy();
      },
    },
    materialId,
  );
}

//

function getSegments(props: VolumeSegmentProps, _volume: Volume): SortedArray {
  return SortedArray.ofUnsortedArray(props.segments);
}

const SegmentVisuals = {
  segment: (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, SegmentMeshParams>) =>
    VolumeRepresentation('Segment mesh', ctx, getParams, SegmentVisual, getLoci, getSegments),
};

export const SegmentParams = {
  ...SegmentMeshParams,
  visuals: PD.MultiSelect(['segment'], PD.objectToOptions(SegmentVisuals)),
  bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type SegmentParams = typeof SegmentParams;
export function getSegmentParams(ctx: ThemeRegistryContext, volume: Volume) {
  const p = PD.clone(SegmentParams);

  const segmentation = Volume.Segmentation.get(volume);
  if (segmentation) {
    const segments = Array.from(segmentation.segments.keys());
    p.segments = PD.Converted(
      (v: number[]) => v.map((x) => `${x}`),
      (v: string[]) => v.map((x) => parseInt(x)),
      PD.MultiSelect(
        segments.map((x) => `${x}`),
        PD.arrayToOptions(segments.map((x) => `${x}`)),
        {
          isEssential: true,
        },
      ),
    );
  }
  return p;
}

export type SegmentRepresentation = VolumeRepresentation<SegmentParams>;
export function SegmentRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Volume, SegmentParams>,
): SegmentRepresentation {
  return Representation.createMulti(
    'Segment',
    ctx,
    getParams,
    Representation.StateBuilder,
    SegmentVisuals as unknown as Representation.Def<Volume, SegmentParams>,
  );
}

export const SegmentRepresentationProvider = VolumeRepresentationProvider({
  name: 'segment',
  label: 'Segment',
  description: 'Displays a triangulated segment of volumetric data.',
  factory: SegmentRepresentation,
  getParams: getSegmentParams,
  defaultValues: PD.getDefaultValues(SegmentParams),
  defaultColorTheme: { name: 'volume-segment' },
  defaultSizeTheme: { name: 'uniform' },
  isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !!Volume.Segmentation.get(volume),
});
