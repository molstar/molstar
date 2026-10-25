/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { ComputeContext, StorageBuffer } from '@molstar/graphics/gl/compute';
import { isComputeIsosurfaceSupported } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Grid, Volume } from '@molstar/model/model/volume';
import type { VisualContext } from '../visual.js';
import type { Theme, ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import {
  computeMarchingCubesMesh,
  computeMarchingCubesLines,
} from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { VolumeRepresentation, VolumeRepresentationProvider } from './representation.js';
import { logGpuFallback, type VisualUpdateState } from '../util.js';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { type RepresentationContext, type RepresentationParamsGetter, Representation } from '../representation.js';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { EmptyLoci, type Loci } from '@molstar/model/model/loci';
import { Tensor, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { fillSerial } from '@molstar/core/util/array';
import {
  createVolumeCellLocationIterator,
  createVolumeTexture2d,
  createWrappedTensor,
  createWrappedVolume,
  eachVolumeLoci,
  getVolumeTexture2dLayout,
  getVolumeGpuFallbackReason,
} from './util.js';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import {
  computeTextureMeshIsosurface,
  extractIsosurface,
} from '@molstar/graphics/geo/geometry/texture-mesh/compute-isosurface';
import { getAtlasIsosurfaceExtraction } from '@molstar/graphics/gl/compute/marching-cubes/atlas';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import type { Texture } from '@molstar/graphics/gl/texture';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Interval } from '@molstar/core/data/int/interval';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { type VolumeKey, VolumeVisual } from './visual.js';
import { estimateIsosurfaceVertexCount } from '@molstar/graphics/geo/util/marching-cubes/estimate';
import { selectVisual } from '../visual-factory.js';

export const VolumeIsosurfaceParams = {
  isoValue: Volume.IsoValueParam,
  wrap: PD.Select('auto', PD.arrayToOptions(['off', 'on', 'auto'] as const)),
  floodfill: PD.Select('off', PD.arrayToOptions(['off', 'inside', 'outside']), {
    description: 'If and how to floodfill the volume. Note that this disables GPU support.',
  }),
};
export type VolumeIsosurfaceParams = typeof VolumeIsosurfaceParams;
export type VolumeIsosurfaceProps = PD.Values<VolumeIsosurfaceParams>;

export const VolumeIsosurfaceTextureParams = {
  ...VolumeIsosurfaceParams,
  tryUseGpu: PD.Boolean(true, { hideIf: (p) => p.floodfill !== 'off' }),
  gpuDataType: PD.Select('byte', PD.arrayToOptions(['byte', 'float', 'halfFloat'] as const), {
    hideIf: (p) => !p.tryUseGpu || p.floodfill !== 'off',
  }),
};
export type VolumeIsosurfaceTextureParams = typeof VolumeIsosurfaceTextureParams;
export type VolumeIsosurfaceTextureProps = PD.Values<VolumeIsosurfaceTextureParams>;

function shouldWrap(volume: Volume, wrap: VolumeIsosurfaceProps['wrap']) {
  if (wrap === 'on') return true;
  if (wrap === 'off') return false;
  return Volume.isPeriodic(volume);
}

const Padding = 1;

/** GPU isosurfaces: the WebGL histogram pyramid, or compute kernels where available (WebGPU) */
/** Why not on the GPU by the props and the context alone, see `getGpuFallbackReason` */
/**
 * The fields of the CPU marching cubes: the values, floodfilled and wrapped as the props say, and each cell's own index
 * as its id, wrapped alike. The ids keep the unfilled space: a floodfilled `get` would fill them too.
 */
function getIsosurfaceFields(volume: Volume, isoLevel: number, props: VolumeIsosurfaceProps) {
  const { cells } = volume.grid;
  let scalarField = cells;
  let idField = Tensor.create(cells.space, Tensor.Data1(fillSerial(new Int32Array(cells.data.length))));
  if (props.floodfill !== 'off') {
    scalarField = Tensor.createFloodfilled(scalarField, isoLevel, props.floodfill);
  }
  if (shouldWrap(volume, props.wrap)) {
    scalarField = createWrappedTensor(scalarField);
    idField = createWrappedTensor(idField);
  }
  return { scalarField, idField };
}

function getPropsGpuFallbackReason(props: PD.Values<IsosurfaceMeshParams>, graphics?: GraphicsContext) {
  if (!props.tryUseGpu) return 'tryUseGpu is off';
  if (!graphics) return 'no graphics context';
  if (props.floodfill !== 'off') return `floodfill is '${props.floodfill}'`;
}

function getGpuFallbackReason(volume: Volume, props: PD.Values<IsosurfaceMeshParams>, graphics?: GraphicsContext) {
  const reason = getPropsGpuFallbackReason(props, graphics);
  if (reason || !graphics) return reason;
  return getVolumeGpuFallbackReason(volume, graphics, Padding, ['floatRenderTarget', 'multipleRenderTargets']);
}

function useGpu(volume: Volume, props: PD.Values<IsosurfaceMeshParams>, graphics?: GraphicsContext): boolean {
  const reason = getGpuFallbackReason(volume, props, graphics);
  if (reason && props.tryUseGpu) logGpuFallback('Volume isosurface', reason);
  return !reason;
}

export const IsosurfaceVisual = selectVisual(
  { gpu: IsosurfaceTextureMeshVisual, cpu: IsosurfaceMeshVisual },
  ({ volume }: VolumeKey, props: PD.Values<IsosurfaceMeshParams>, graphics?: GraphicsContext) =>
    useGpu(volume, props, graphics) ? 'gpu' : 'cpu',
  (props, graphics) => (getPropsGpuFallbackReason(props, graphics) ? ['cpu'] : ['gpu', 'cpu']),
);

function getLoci(volume: Volume, props: VolumeIsosurfaceProps) {
  const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
  return Volume.Isosurface.Loci(volume, props.isoValue, instances);
}

function getIsosurfaceLoci(
  pickingId: PickingId,
  volume: Volume,
  key: number,
  props: VolumeIsosurfaceProps,
  id: number,
) {
  const { objectId, groupId, instanceId } = pickingId;
  if (id === objectId) {
    const granularity = Volume.PickingGranularity.get(volume);
    const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
    if (granularity === 'volume') {
      return Volume.Loci(volume, instances);
    } else if (granularity === 'object' || groupId === PickingId.Null) {
      return Volume.Isosurface.Loci(volume, props.isoValue, instances);
    } else {
      const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
      return Volume.Cell.Loci(volume, [{ indices, instances }]);
    }
  }
  return EmptyLoci;
}

export function eachIsosurface(
  loci: Loci,
  volume: Volume,
  key: number,
  props: VolumeIsosurfaceProps,
  apply: (interval: Interval) => boolean,
) {
  return eachVolumeLoci(loci, volume, { isoValue: props.isoValue }, apply);
}

//

export async function createVolumeIsosurfaceMesh(
  ctx: VisualContext,
  volume: Volume,
  key: number,
  theme: Theme,
  props: VolumeIsosurfaceProps,
  mesh?: Mesh,
) {
  ctx.runtime.update({ message: 'Marching cubes...' });

  const isoLevel = Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue;

  const { scalarField, idField } = getIsosurfaceFields(volume, isoLevel, props);
  const surface = await computeMarchingCubesMesh({ isoLevel, scalarField, idField }, mesh).runAsChild(ctx.runtime);

  const transform = Grid.getGridToCartesianTransform(volume.grid);
  Mesh.transform(surface, transform);
  if (ctx.graphics && !ctx.graphics.capabilities.flatInterpolation) {
    // 2nd arg means not to split triangles based on group id. Splitting triangles
    // is too expensive if each cell has its own group id as is the case here.
    Mesh.uniformTriangleGroup(surface, false);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  surface.setBoundingSphere(
    Grid.getIsosurfaceBoundingSphere(
      volume.grid,
      Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue,
    ),
  );

  return surface;
}

export const IsosurfaceMeshParams = {
  ...Mesh.Params,
  ...TextureMesh.Params,
  ...VolumeIsosurfaceParams,
  ...VolumeIsosurfaceTextureParams,
  quality: { ...Mesh.Params.quality, isEssential: false },
};
export type IsosurfaceMeshParams = typeof IsosurfaceMeshParams;

export function IsosurfaceMeshVisual(materialId: number): VolumeVisual<IsosurfaceMeshParams> {
  return VolumeVisual<Mesh, IsosurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
      createGeometry: createVolumeIsosurfaceMesh,
      createLocationIterator: createVolumeCellLocationIterator,
      getLoci: getIsosurfaceLoci,
      eachLocation: eachIsosurface,
      setUpdateState: (
        state: VisualUpdateState,
        newVolume: Volume,
        currentVolume: Volume,
        newProps: PD.Values<IsosurfaceMeshParams>,
        currentProps: PD.Values<IsosurfaceMeshParams>,
      ) => {
        state.createGeometry =
          !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, newVolume.grid.stats) ||
          newProps.wrap !== currentProps.wrap ||
          newProps.floodfill !== currentProps.floodfill;
      },
      geometryUtils: Mesh.Utils,
    },
    materialId,
  );
}

//

namespace VolumeIsosurfaceTexture {
  const name = 'volume-isosurface-texture';
  export const descriptor = CustomPropertyDescriptor({ name });
  export function get(volume: Volume, graphics: GraphicsContext, props: VolumeIsosurfaceTextureProps) {
    const { gpuDataType } = props;
    const wrap = shouldWrap(volume, props.wrap);

    const transform = Grid.getGridToCartesianTransform(volume.grid);
    const gridDimension = Vec3.clone(volume.grid.cells.space.dimensions as Vec3);
    const { width, height, powerOfTwoSize: texDim } = getVolumeTexture2dLayout(gridDimension, Padding);
    const gridTexDim = Vec3.create(width, height, 0);
    const gridDataDim = Vec3.subScalar(Vec3(), gridDimension, wrap ? 1 : 0);
    const gridTexScale = Vec2.create(width / texDim, height / texDim);
    // console.log({ texDim, width, height, gridDimension });

    if (texDim > graphics.resources.limits.maxTextureSize / 2) {
      throw new Error('volume too large for gpu isosurface extraction');
    }

    const dataType = gpuDataType === 'halfFloat' && !graphics.capabilities.halfFloatTexture ? 'float' : gpuDataType;

    const cached = volume._propertyData[name];
    // per context: a texture of a replaced context (backend switch) is not drawn by the new one
    if (cached?.dataType !== dataType || cached?.wrap !== wrap || cached?.graphics !== graphics) {
      const texture =
        dataType === 'byte'
          ? graphics.resources.texture('image-uint8', 'r', 'ubyte', 'linear')
          : dataType === 'halfFloat'
            ? graphics.resources.texture('image-float16', 'r', 'fp16', 'linear')
            : graphics.resources.texture('image-float32', 'r', 'float', 'linear');
      volume._propertyData[name] = { texture, dataType, wrap, graphics };
      texture.define(texDim, texDim);
      // load volume into sub-section of texture
      texture.load(createVolumeTexture2d(volume, 'data', Padding, dataType), true);
      volume.customProperties.add(descriptor);
      volume.customProperties.assets(descriptor, [{ dispose: () => texture.destroy() }]);
    }

    gridDimension[0] += Padding;
    gridDimension[1] += Padding;

    return {
      texture: volume._propertyData[name].texture as Texture,
      transform,
      gridDimension,
      gridTexDim,
      gridDataDim,
      gridTexScale,
    };
  }
}

/** The volume values in a storage buffer for compute kernels, cached on the volume */
namespace VolumeIsosurfaceBuffer {
  const name = 'volume-isosurface-buffer';
  export const descriptor = CustomPropertyDescriptor({ name });
  export function get(volume: Volume, compute: ComputeContext): StorageBuffer {
    const cached = volume._propertyData[name] as { buffer: StorageBuffer; compute: ComputeContext } | undefined;
    if (cached?.compute === compute) return cached.buffer;
    const { data } = volume.grid.cells;
    const values = data instanceof Float32Array ? data : Float32Array.from(data as ArrayLike<number>);
    const buffer = compute.createBuffer(values.byteLength, 'volume isosurface');
    buffer.write(values);
    volume._propertyData[name] = { buffer, compute };
    volume.customProperties.add(descriptor);
    volume.customProperties.assets(descriptor, [{ dispose: () => buffer.destroy() }]);
    return buffer;
  }
}

/** Marching cubes with compute kernels on the volume values (WebGPU) */
async function createVolumeIsosurfaceComputeMesh(
  graphics: GraphicsContext & { compute: ComputeContext },
  volume: Volume,
  props: VolumeIsosurfaceTextureProps,
  textureMesh?: TextureMesh,
) {
  const value = Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue;
  const { space, data } = volume.grid.cells;
  const strides = Vec3.create(space.dataOffset(1, 0, 0), space.dataOffset(0, 1, 0), space.dataOffset(0, 0, 1));
  const dataDim = Vec3.clone(space.dimensions as Vec3);
  const buffer = VolumeIsosurfaceBuffer.get(volume, graphics.compute);
  // periodic: the grid one point larger, wrapped into the data by the kernels
  if (shouldWrap(volume, props.wrap)) {
    volume = createWrappedVolume(volume);
  }
  const gridDim = Vec3.clone(volume.grid.cells.space.dimensions as Vec3);
  const input = { kind: 'buffer', volume: buffer, strides, dataDim } as const;
  const transform = Grid.getGridToCartesianTransform(volume.grid);
  // getting isosurface bounding-sphere is too expensive here
  const boundingSphere = Grid.getBoundingSphere(volume.grid);
  // sized exactly when created, drawn indirectly on updates; the arena grown for the surface beforehand, so that the
  // first extraction gets all its pages
  return computeTextureMeshIsosurface(
    graphics,
    (target) => extractIsosurface(graphics, input, gridDim, value, transform, value < 0, target),
    data.length,
    boundingSphere,
    textureMesh,
    () => estimateIsosurfaceVertexCount(data, gridDim, strides, value, dataDim),
  );
}

function createVolumeIsosurfaceTextureMesh(
  ctx: VisualContext,
  volume: Volume,
  key: number,
  theme: Theme,
  props: VolumeIsosurfaceTextureProps,
  textureMesh?: TextureMesh,
): TextureMesh | Promise<TextureMesh> {
  const graphics = ctx.graphics;
  if (graphics && isComputeIsosurfaceSupported(graphics)) {
    if (volume.grid.cells.data.length <= 1) return TextureMesh.createEmpty(textureMesh);
    return createVolumeIsosurfaceComputeMesh(graphics, volume, props, textureMesh);
  }
  const atlas = graphics && getAtlasIsosurfaceExtraction(graphics);
  if (!graphics || !atlas) {
    throw new Error('compute kernels or fragment-shader marching cubes required to create a volume isosurface');
  }

  if (volume.grid.cells.data.length <= 1) {
    return TextureMesh.createEmpty(textureMesh);
  }

  if (shouldWrap(volume, props.wrap)) {
    volume = createWrappedVolume(volume);
  }

  const { max, min } = volume.grid.stats;
  const diff = max - min;
  const value = Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue;
  const isoLevel = (value - min) / diff;

  const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
  const groupCount = volume.grid.cells.data.length;
  const boundingSphere = Grid.getBoundingSphere(volume.grid); // getting isosurface bounding-sphere is too expensive here

  const create = (textureMesh?: TextureMesh) => {
    const { texture, gridDimension, gridTexDim, gridDataDim, gridTexScale, transform } = VolumeIsosurfaceTexture.get(
      volume,
      graphics,
      props,
    );

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
        isoValue: isoLevel,
        invert: value < 0,
        packedGroup: false,
        axisOrder: axisOrder,
        constantGroup: true,
      },
      buffer,
    );

    return TextureMesh.create(
      gv.vertexCount,
      groupCount,
      gv.vertexTexture,
      gv.groupTexture,
      gv.normalTexture,
      boundingSphere,
      textureMesh,
    );
  };

  return create(textureMesh);
}

export function IsosurfaceTextureMeshVisual(materialId: number): VolumeVisual<IsosurfaceMeshParams> {
  return VolumeVisual<TextureMesh, IsosurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(IsosurfaceMeshParams),
      createGeometry: createVolumeIsosurfaceTextureMesh,
      createLocationIterator: createVolumeCellLocationIterator,
      getLoci: getIsosurfaceLoci,
      eachLocation: eachIsosurface,
      setUpdateState: (
        state: VisualUpdateState,
        newVolume: Volume,
        currentVolume: Volume,
        newProps: PD.Values<IsosurfaceMeshParams>,
        currentProps: PD.Values<IsosurfaceMeshParams>,
      ) => {
        state.createGeometry =
          !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, newVolume.grid.stats) ||
          newProps.gpuDataType !== currentProps.gpuDataType ||
          newProps.wrap !== currentProps.wrap;
      },
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

export async function createVolumeIsosurfaceWireframe(
  ctx: VisualContext,
  volume: Volume,
  key: number,
  theme: Theme,
  props: VolumeIsosurfaceProps,
  lines?: Lines,
) {
  ctx.runtime.update({ message: 'Marching cubes...' });

  const isoLevel = Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue;

  const { scalarField, idField } = getIsosurfaceFields(volume, isoLevel, props);
  const wireframe = await computeMarchingCubesLines({ isoLevel, scalarField, idField }, lines).runAsChild(ctx.runtime);

  const transform = Grid.getGridToCartesianTransform(volume.grid);
  Lines.transform(wireframe, transform);

  wireframe.setBoundingSphere(
    Grid.getIsosurfaceBoundingSphere(
      volume.grid,
      Volume.IsoValue.toAbsolute(props.isoValue, volume.grid.stats).absoluteValue,
    ),
  );

  return wireframe;
}

export const IsosurfaceWireframeParams = {
  ...Lines.Params,
  ...VolumeIsosurfaceParams,
  quality: { ...Lines.Params.quality, isEssential: false },
  sizeFactor: PD.Numeric(3, { min: 0, max: 10, step: 0.1 }),
};
export type IsosurfaceWireframeParams = typeof IsosurfaceWireframeParams;

export function IsosurfaceWireframeVisual(materialId: number): VolumeVisual<IsosurfaceWireframeParams> {
  return VolumeVisual<Lines, IsosurfaceWireframeParams>(
    {
      defaultProps: PD.getDefaultValues(IsosurfaceWireframeParams),
      createGeometry: createVolumeIsosurfaceWireframe,
      createLocationIterator: createVolumeCellLocationIterator,
      getLoci: getIsosurfaceLoci,
      eachLocation: eachIsosurface,
      setUpdateState: (
        state: VisualUpdateState,
        newVolume: Volume,
        currentVolume: Volume,
        newProps: PD.Values<IsosurfaceWireframeParams>,
        currentProps: PD.Values<IsosurfaceWireframeParams>,
      ) => {
        state.createGeometry =
          !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, newVolume.grid.stats) ||
          newProps.wrap !== currentProps.wrap ||
          newProps.floodfill !== currentProps.floodfill;
      },
      geometryUtils: Lines.Utils,
    },
    materialId,
  );
}

//

const IsosurfaceVisuals = {
  solid: (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceMeshParams>) =>
    VolumeRepresentation('Isosurface mesh', ctx, getParams, IsosurfaceVisual, getLoci),
  wireframe: (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, IsosurfaceWireframeParams>) =>
    VolumeRepresentation('Isosurface wireframe', ctx, getParams, IsosurfaceWireframeVisual, getLoci),
};

export const IsosurfaceParams = {
  ...IsosurfaceMeshParams,
  ...IsosurfaceWireframeParams,
  visuals: PD.MultiSelect(['solid'], PD.objectToOptions(IsosurfaceVisuals)),
  bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type IsosurfaceParams = typeof IsosurfaceParams;
export function getIsosurfaceParams(ctx: ThemeRegistryContext, volume: Volume) {
  const p = PD.clone(IsosurfaceParams);
  p.isoValue = Volume.createIsoValueParam(Volume.IsoValue.relative(2), volume.grid.stats);
  return p;
}

export type IsosurfaceRepresentation = VolumeRepresentation<IsosurfaceParams>;
export function IsosurfaceRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Volume, IsosurfaceParams>,
): IsosurfaceRepresentation {
  return Representation.createMulti(
    'Isosurface',
    ctx,
    getParams,
    Representation.StateBuilder,
    IsosurfaceVisuals as unknown as Representation.Def<Volume, IsosurfaceParams>,
  );
}

export const IsosurfaceRepresentationProvider = VolumeRepresentationProvider({
  name: 'isosurface',
  label: 'Isosurface',
  description: 'Displays a triangulated isosurface of volumetric data.',
  factory: IsosurfaceRepresentation,
  getParams: getIsosurfaceParams,
  defaultValues: PD.getDefaultValues(IsosurfaceParams),
  defaultColorTheme: { name: 'uniform' },
  defaultSizeTheme: { name: 'uniform' },
  isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume),
});
