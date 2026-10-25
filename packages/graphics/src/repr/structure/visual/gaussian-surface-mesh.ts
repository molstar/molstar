/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  UnitsMeshParams,
  UnitsTextureMeshParams,
  type UnitsVisual,
  UnitsMeshVisual,
  UnitsTextureMeshVisual,
} from '../units-visual.js';
import {
  GaussianDensityParams,
  computeUnitGaussianDensity,
  computeUnitGaussianDensityTexture,
  computeUnitGaussianDensityTexture2d,
  type GaussianDensityProps,
  computeStructureGaussianDensity,
  computeStructureGaussianDensityTexture,
  computeStructureGaussianDensityTexture2d,
  getComputeDensityMaxCells,
} from './util/gaussian.js';
import type { VisualContext } from '../../visual.js';
import type { Unit, Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import {
  ElementIterator,
  getElementLoci,
  eachElement,
  getSerialElementLoci,
  eachSerialElement,
} from './util/element.js';
import { logGpuFallback, missingCapabilities, type VisualUpdateState } from '../../util.js';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import {
  computeTextureMeshIsosurface,
  extractIsosurface,
} from '@molstar/graphics/geo/geometry/texture-mesh/compute-isosurface';
import { getAtlasIsosurfaceExtraction } from '@molstar/graphics/gl/compute/marching-cubes/atlas';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import {
  type ComplexVisual,
  ComplexMeshParams,
  ComplexMeshVisual,
  ComplexTextureMeshVisual,
  ComplexTextureMeshParams,
} from '../complex-visual.js';
import { getVolumeSliceInfo, type StructureGroup } from './util/common.js';
import type { ComputeContext } from '@molstar/graphics/gl/compute';
import {
  getIsosurfaceKernels,
  isComputeIsosurfaceSupported,
} from '@molstar/graphics/gl/compute/marching-cubes/compute';
import { isGaussianDensityTextureSupported } from '@molstar/graphics/geo/gaussian-density/gpu';
import type { GaussianDensityTextureData } from '@molstar/graphics/geo/gaussian-density';
import { GaussianDensityKernels } from '@molstar/graphics/gl/compute/gaussian-density';
import { getColorSmoothingKernels } from '@molstar/graphics/gl/compute/color-smoothing';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import type { TextureMeshValues } from '@molstar/graphics/gl/renderable/texture-mesh';
import type { Texture } from '@molstar/graphics/gl/texture';
import { applyMeshColorSmoothing } from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import { applyTextureMeshColorSmoothing } from '@molstar/graphics/geo/geometry/texture-mesh/color-smoothing';
import { ColorSmoothingParams, getColorSmoothingProps } from '@molstar/graphics/geo/geometry/base';
import { Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { isTimingMode } from '@molstar/core/util/debug';
import { ValueCell } from '@molstar/core/util/value-cell';
import { getSurfaceSignatureDefines } from '@molstar/graphics/repr/signatures';
import { selectVisual } from '../../visual-factory.js';

const SharedParams = {
  ...GaussianDensityParams,
  ...ColorSmoothingParams,
  tryUseGpu: PD.Boolean(true),
};
type SharedParams = typeof SharedParams;

export const GaussianSurfaceMeshParams = {
  ...UnitsMeshParams,
  ...UnitsTextureMeshParams,
  ...SharedParams,
};
export type GaussianSurfaceMeshParams = typeof GaussianSurfaceMeshParams;

export const StructureGaussianSurfaceMeshParams = {
  ...ComplexMeshParams,
  ...ComplexTextureMeshParams,
  ...SharedParams,
};
export type StructureGaussianSurfaceMeshParams = typeof StructureGaussianSurfaceMeshParams;

/** Of the WebGL path */
const GpuCapabilities = ['floatRenderTarget', 'floatTextureLinear', 'multipleRenderTargets'] as const;

/**
 * What a density is computed over: the unit's box and atoms for the surface of a symmetry group (each unit has its own
 * density), the structure's for the surface of the whole structure
 */
type DensityExtent = { readonly box: Box3D; readonly elementCount: number };

function getDensityExtent(target: Structure | Unit): DensityExtent {
  return 'units' in target
    ? { box: target.boundary.box, elementCount: target.elementCount }
    : { box: target.boundary.box, elementCount: target.elements.length };
}

/** Why the 2d density of the WebGL path doesn't fit, if it doesn't */
function unsuitableForGpu(extent: DensityExtent, props: PD.Values<SharedParams>, graphics: GraphicsContext) {
  // the GPU is much more memory contraint, especially true for integrated GPUs,
  // being conservative here still allows for small and medium sized assemblies
  const d = graphics.resources.limits.maxTextureSize / 3;
  const { areaCells, maxAreaCells } = getVolumeSliceInfo(extent.box, props.resolution, d * d);
  if (areaCells >= maxAreaCells) return `the largest slice has ${areaCells} cells, at most ${maxAreaCells} fit`;
}

/**
 * Why the 3d density of the compute path doesn't fit, if it doesn't: each dimension a 3d texture's, all cells the
 * budget (the WebGL criterion is the 2d atlas of slices, sized by the smaller WebGPU 2d texture limit it would exclude
 * medium structures)
 */
function unsuitableForCompute(
  extent: DensityExtent,
  props: PD.Values<SharedParams>,
  graphics: GraphicsContext & { compute: ComputeContext },
) {
  const { dim, cells } = getComputeDensityGrid(extent.box, props);
  const max3d = graphics.resources.limits.max3dTextureSize;
  const n = Math.max(dim[0], dim[1], dim[2]);
  if (n > max3d) return `the density grid has ${n} cells along an axis, 3d textures at most ${max3d}`;
  const maxCells = getComputeDensityMaxCells(graphics.compute);
  if (cells > maxCells) return `the density grid has ${cells} cells, at most ${maxCells}`;
  // without the counts: logged once, not once per small unit (`logGpuFallback`)
  if (cells < SmallSurfaceCells && extent.elementCount < SmallSurfaceAtoms) {
    return `small (under ${SmallSurfaceCells} density cells and ${SmallSurfaceAtoms} atoms): faster on the CPU`;
  }
}

/** The cells of the compute path's 3d density of `box` */
function getComputeDensityGrid(box: Box3D, props: PD.Values<SharedParams>) {
  // the padding of the density: twice the largest radius (3 Å bounds it) and four cells
  const pad = 2 * 3 + 4 * props.resolution;
  const size = Box3D.size(Vec3(), box);
  const dim = [0, 1, 2].map((i) => Math.ceil((size[i] + 2 * pad) / props.resolution));
  return { dim, cells: dim[0] * dim[1] * dim[2] };
}

/**
 * Below both, the compute path's surfaces are computed on the CPU: it takes 1–6 ms there, less than the compute path's
 * dispatches and readback (about 6–7 ms whatever the size); above, the compute path is faster at any resolution, on a
 * weak integrated GPU too (AMD Radeon of a Ryzen 7900X: 30k atoms at resolution 2 in 11 ms against 30 ms on the CPU).
 * The atoms bound the CPU's density work on coarse grids (10k atoms at resolution 3: 43k cells, 10 ms).
 */
const SmallSurfaceCells = 100_000;
const SmallSurfaceAtoms = 5000;

/** Compute marching cubes on the GPU gaussian density (WebGPU) */
function computeSupport(graphics: GraphicsContext): graphics is GraphicsContext & { compute: ComputeContext } {
  return (
    graphics.capabilities.texture3d &&
    isGaussianDensityTextureSupported(graphics) &&
    isComputeIsosurfaceSupported(graphics)
  );
}

/** Why GPU surfaces are not available, if they are not (compute kernels, or the WebGL path), before what a structure needs */
export function getGaussianSurfaceGpuUnsupportedReason(graphics: GraphicsContext): string | undefined {
  if (computeSupport(graphics)) return;
  if (!getAtlasIsosurfaceExtraction(graphics)) return 'no compute kernels or fragment-shader marching cubes';
  const missing = missingCapabilities(graphics, GpuCapabilities);
  if (missing.length) return `missing ${missing.join(', ')}`;
}

/** Whether GPU surfaces are available (compute kernels, or the WebGL path), before what a structure needs */
export function isGaussianSurfaceGpuSupported(graphics: GraphicsContext) {
  return !getGaussianSurfaceGpuUnsupportedReason(graphics);
}

/** Why the GPU surface is not computed, if it is not */
/** Why not on the GPU by the props and the context alone, see `getGpuFallbackReason` */
function getPropsGpuFallbackReason(props: PD.Values<SharedParams>, graphics?: GraphicsContext) {
  if (!props.tryUseGpu) return 'tryUseGpu is off';
  if (!graphics) return 'no graphics context';
  if (props.includeParent) return 'includeParent is on';
  if (props.floodfill !== 'off') return `floodfill is '${props.floodfill}'`;
  // the compute path decides by the size of the surface instead (`unsuitableForCompute`)
  if (computeSupport(graphics)) return;
  // WebGL: lower resolutions are about as fast on the CPU as with the histogram pyramid on integrated GPUs, and very
  // low resolutions have artifacts there
  if (props.resolution > 1) return `the resolution ${props.resolution} is above 1 (about as fast on the CPU)`;
  if (!getAtlasIsosurfaceExtraction(graphics)) return 'no compute kernels or fragment-shader marching cubes';
  const missing = missingCapabilities(graphics, GpuCapabilities);
  if (missing.length) return `missing ${missing.join(', ')}`;
}

/**
 * Why the Gaussian surface of `target` is computed on the CPU, if it is (the props, the context, the size): a unit's
 * as a surface per symmetry group, a structure's as the surface of the whole structure
 */
export function getGaussianSurfaceGpuFallbackReason(
  target: Structure | Unit,
  props: PD.Values<SharedParams>,
  graphics?: GraphicsContext,
) {
  return getGpuFallbackReason(getDensityExtent(target), props, graphics);
}

function getGpuFallbackReason(extent: DensityExtent, props: PD.Values<SharedParams>, graphics?: GraphicsContext) {
  const reason = getPropsGpuFallbackReason(props, graphics);
  if (reason || !graphics) return reason;
  if (computeSupport(graphics)) return unsuitableForCompute(extent, props, graphics);
  return unsuitableForGpu(extent, props, graphics);
}

/** GPU surfaces: the WebGL histogram pyramid on the 2d density, or compute kernels on the 3d density (WebGPU) */
function useGpu(target: Structure | Unit, props: PD.Values<SharedParams>, graphics?: GraphicsContext): boolean {
  const reason = getGpuFallbackReason(getDensityExtent(target), props, graphics);
  if (reason && props.tryUseGpu) logGpuFallback('Gaussian surface', reason);
  return !reason;
}

/**
 * The 3d density texture of the compute path, shared by the surfaces and only grown: the isosurface reads it by grid
 * coordinates, and defining it to each unit's grid recreated it per unit
 */
function getComputeDensityTexture(graphics: GraphicsContext) {
  const name = 'gaussian-surface-3d';
  if (!graphics.namedTextures[name]) {
    graphics.namedTextures[name] = graphics.resources.texture('volume-float16', 'rgba', 'fp16', 'linear');
  }
  return graphics.namedTextures[name];
}

/** Upper limit of the estimate, about 224 MB of vertex textures; larger surfaces are extracted again */
const MaxVertexEstimate = 8 * 1024 * 1024;

/**
 * Vertex count estimate of a Gaussian surface of `atomCount` atoms (unindexed, three per triangle): polymer chains of
 * PDB entries have 73–86 · atoms / resolution² (cavities inside make it linear in the atoms), with headroom; units of
 * isolated atoms (waters) have about 240 and are extracted again
 */
function estimateVertexCount(atomCount: number, resolution: number) {
  return Math.min(MaxVertexEstimate, (100 * atomCount) / (resolution * resolution));
}

/**
 * The isosurface of a 3d gaussian density with compute kernels, the group packed in the density's rgb; drawn
 * indirectly. `computeDensity` is called by each extraction (the density texture is shared)
 */
/** The visual's props: the density, and the color smoothing applied to the surface once created */
type GaussianSurfaceProps = GaussianDensityProps & Partial<Pick<PD.Values<GaussianSurfaceMeshParams>, 'smoothColors'>>;

/** How the isosurface reads the density texture */
const DensityInput = { kind: 'texture', channel: 'alpha', packedGroup: true } as const;

async function createComputeGaussianSurface(
  graphics: GraphicsContext & { compute: ComputeContext },
  computeDensity: () => GaussianDensityTextureData,
  props: GaussianSurfaceProps,
  sphere: Sphere3D,
  groupCount: number,
  textureMesh?: TextureMesh,
) {
  // compiled in the background, not blocking the queue (and frames) at the first dispatch: all before the density is
  // computed, into a texture shared by the surfaces, which another surface could write while the isosurface waited for
  // its kernels; with color smoothing also its kernels, which the visual's synchronous update (`processValues`)
  // dispatches (the color type not known yet)
  const kernels = [...GaussianDensityKernels, ...getIsosurfaceKernels(DensityInput)];
  if (props.smoothColors && props.smoothColors.name !== 'off') {
    kernels.push(...getColorSmoothingKernels('group'), ...getColorSmoothingKernels('groupInstance'));
  }
  await graphics.compute.prepareKernels(kernels);
  // computed by each extraction, right before its commands: the texture is shared, and the pages of the surface are
  // reserved (awaited) before
  const last: { density?: GaussianDensityTextureData } = {};
  const extraRadius = props.radiusOffset * (1 + Math.exp(-props.smoothness));
  const surface = await computeTextureMeshIsosurface(
    graphics,
    (target) => {
      const d = computeDensity();
      last.density = d;
      const isoLevel = Math.exp(-props.smoothness) / d.radiusFactor;
      const input = { ...DensityInput, texture: d.texture };
      return extractIsosurface(graphics, input, d.gridDim, isoLevel, d.transform, false, target);
    },
    groupCount,
    () => Sphere3D.expand(Sphere3D(), sphere, (last.density?.maxRadius ?? 0) + extraRadius),
    textureMesh,
    // by the resolution asked for (the density's may be coarser for large grids: fewer vertices)
    estimateVertexCount(groupCount, props.resolution),
  );
  const { density } = last;
  if (!density) throw new Error('expected the density of the surface');
  (surface.meta as GaussianSurfaceMeta).resolution = density.resolution;
  // a re-extraction computes the same grid
  surface.meta.vertexBox = TextureMesh.getGridVertexBox(density.gridDim, density.transform);
  return surface;
}

/** The surfaces the props and context allow, before the size decides (`useGpu`) */
function getSurfaceCandidates(props: PD.Values<SharedParams>, graphics?: GraphicsContext) {
  return getPropsGpuFallbackReason(props, graphics) ? (['cpu'] as const) : (['gpu', 'cpu'] as const);
}

export const GaussianSurfaceVisual = selectVisual(
  { gpu: GaussianSurfaceTextureMeshVisual, cpu: GaussianSurfaceMeshVisual },
  // by the unit the geometry is created of (`createGaussianSurfaceMesh`), each unit having its own density
  ({ group }: StructureGroup, props, graphics) => (useGpu(group.units[0], props, graphics) ? 'gpu' : 'cpu'),
  getSurfaceCandidates,
);

export const StructureGaussianSurfaceVisual = selectVisual(
  { gpu: StructureGaussianSurfaceTextureMeshVisual, cpu: StructureGaussianSurfaceMeshVisual },
  (structure: Structure, props, graphics) => (useGpu(structure, props, graphics) ? 'gpu' : 'cpu'),
  getSurfaceCandidates,
);

type GaussianSurfaceMeta = {
  resolution?: number;
  colorTexture?: Texture;
};

//

async function createGaussianSurfaceMesh(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: GaussianDensityProps,
  mesh?: Mesh,
): Promise<Mesh> {
  const { smoothness, floodfill, includeParent, radiusOffset } = props;
  const { transform, field, idField, radiusFactor, resolution, maxRadius } = await computeUnitGaussianDensity(
    structure,
    unit,
    theme.size,
    props,
  ).runInContext(ctx.runtime);

  const isoLevel = Math.exp(-smoothness) / radiusFactor;
  const params = {
    isoLevel,
    scalarField: floodfill !== 'off' ? Tensor.createFloodfilled(field, isoLevel, floodfill) : field,
    idField,
  };
  const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
  (surface.meta.resolution as GaussianSurfaceMeta['resolution']) = resolution;

  if (includeParent) {
    const iterations = Math.ceil(2 / resolution);
    Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
  }

  Mesh.transform(surface, transform);
  if (ctx.graphics && !ctx.graphics.capabilities.flatInterpolation) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
  const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius + extraRadius);
  surface.setBoundingSphere(sphere);

  return surface;
}

export function GaussianSurfaceMeshVisual(materialId: number): UnitsVisual<GaussianSurfaceMeshParams> {
  return UnitsMeshVisual<GaussianSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(GaussianSurfaceMeshParams),
      signatureDefines: getSurfaceSignatureDefines,
      createGeometry: createGaussianSurfaceMesh,
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<GaussianSurfaceMeshParams>,
        currentProps: PD.Values<GaussianSurfaceMeshParams>,
      ) => {
        state.createGeometry =
          newProps.resolution !== currentProps.resolution ||
          newProps.radiusOffset !== currentProps.radiusOffset ||
          newProps.smoothness !== currentProps.smoothness ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.includeParent !== currentProps.includeParent ||
          newProps.floodfill !== currentProps.floodfill;

        if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
          state.updateColor = true;
        } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
          if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor)
            state.updateColor = true;
          if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride)
            state.updateColor = true;
        }
      },
      processValues: (
        values: MeshValues,
        geometry: Mesh,
        props: PD.Values<GaussianSurfaceMeshParams>,
        theme: Theme,
        graphics?: GraphicsContext,
      ) => {
        const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
        if (csp) {
          applyMeshColorSmoothing(values, csp, graphics, colorTexture);
          (geometry.meta.colorTexture as GaussianSurfaceMeta['colorTexture']) = values.tColorGrid.ref.value;
        }
      },
      dispose: (geometry: Mesh) => {
        (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
      },
    },
    materialId,
  );
}

//

async function createStructureGaussianSurfaceMesh(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: GaussianDensityProps,
  mesh?: Mesh,
): Promise<Mesh> {
  const { smoothness, floodfill, includeParent, radiusOffset } = props;
  const { transform, field, idField, radiusFactor, resolution, maxRadius } = await computeStructureGaussianDensity(
    structure,
    theme.size,
    props,
  ).runInContext(ctx.runtime);

  const isoLevel = Math.exp(-smoothness) / radiusFactor;
  const params = {
    isoLevel,
    scalarField: floodfill !== 'off' ? Tensor.createFloodfilled(field, isoLevel, floodfill) : field,
    idField,
  };
  const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
  (surface.meta.resolution as GaussianSurfaceMeta['resolution']) = resolution;

  if (includeParent) {
    const iterations = Math.ceil(2 / resolution);
    Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
  }

  Mesh.transform(surface, transform);
  if (ctx.graphics && !ctx.graphics.capabilities.flatInterpolation) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
  const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius);
  surface.setBoundingSphere(sphere);

  return surface;
}

export function StructureGaussianSurfaceMeshVisual(
  materialId: number,
): ComplexVisual<StructureGaussianSurfaceMeshParams> {
  return ComplexMeshVisual<StructureGaussianSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(StructureGaussianSurfaceMeshParams),
      signatureDefines: getSurfaceSignatureDefines,
      createGeometry: createStructureGaussianSurfaceMesh,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<GaussianSurfaceMeshParams>,
        currentProps: PD.Values<GaussianSurfaceMeshParams>,
      ) => {
        state.createGeometry =
          newProps.resolution !== currentProps.resolution ||
          newProps.radiusOffset !== currentProps.radiusOffset ||
          newProps.smoothness !== currentProps.smoothness ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.includeParent !== currentProps.includeParent ||
          newProps.floodfill !== currentProps.floodfill;

        if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
          state.updateColor = true;
        } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
          if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor)
            state.updateColor = true;
          if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride)
            state.updateColor = true;
        }
      },
      processValues: (
        values: MeshValues,
        geometry: Mesh,
        props: PD.Values<GaussianSurfaceMeshParams>,
        theme: Theme,
        graphics?: GraphicsContext,
      ) => {
        const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
        if (csp) {
          applyMeshColorSmoothing(values, csp, graphics, colorTexture);
          (geometry.meta.colorTexture as GaussianSurfaceMeta['colorTexture']) = values.tColorGrid.ref.value;
        }
      },
      dispose: (geometry: Mesh) => {
        (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
      },
    },
    materialId,
  );
}

//

const GaussianSurfaceName = 'gaussian-surface';

function createGaussianSurfaceTextureMesh(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: GaussianSurfaceProps,
  textureMesh?: TextureMesh,
): TextureMesh | Promise<TextureMesh> {
  const graphics = ctx.graphics;
  if (graphics && computeSupport(graphics)) {
    const computeDensity = () =>
      computeUnitGaussianDensityTexture(
        structure,
        unit,
        theme.size,
        props,
        graphics,
        getComputeDensityTexture(graphics),
        true,
      );
    return createComputeGaussianSurface(
      graphics,
      computeDensity,
      props,
      unit.boundary.sphere,
      unit.elements.length,
      textureMesh,
    );
  }
  const atlas = graphics && getAtlasIsosurfaceExtraction(graphics);
  if (!graphics || !atlas) {
    throw new Error('compute kernels or fragment-shader marching cubes required to create a gaussian surface');
  }

  const {
    namedTextures,
    resources,
    capabilities: { floatBlend, halfFloatRenderTarget },
  } = graphics;
  if (!namedTextures[GaussianSurfaceName]) {
    namedTextures[GaussianSurfaceName] = halfFloatRenderTarget
      ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
      : floatBlend
        ? resources.texture('image-float32', 'rgba', 'float', 'linear')
        : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
  }

  const axisOrder = Vec3.create(0, 1, 2);
  const groupCount = unit.elements.length;

  const create = (textureMesh?: TextureMesh) => {
    if (isTimingMode) graphics.timer.group('createGaussianSurfaceTextureMesh');
    const densityTextureData = computeUnitGaussianDensityTexture2d(
      structure,
      unit,
      theme.size,
      true,
      props,
      graphics,
      namedTextures[GaussianSurfaceName],
    );
    const isoLevel = Math.exp(-props.smoothness) / densityTextureData.radiusFactor;

    const buffer = textureMesh?.doubleBuffer.get(graphics);
    const gv = atlas.extract(
      graphics,
      {
        texture: densityTextureData.texture,
        channel: 'alpha',
        gridDim: densityTextureData.gridDim,
        gridTexDim: densityTextureData.gridTexDim,
        gridDataDim: densityTextureData.gridDataDim,
        gridTexScale: densityTextureData.gridTexScale,
        transform: densityTextureData.transform,
        isoValue: isoLevel,
        invert: false,
        packedGroup: true,
        axisOrder: axisOrder,
        constantGroup: true,
      },
      buffer,
    );
    if (isTimingMode) graphics.timer.groupEnd('createGaussianSurfaceTextureMesh');

    const extraRadius = props.radiusOffset * (1 + Math.exp(-props.smoothness));
    const boundingSphere = Sphere3D.expand(
      Sphere3D(),
      unit.boundary.sphere,
      densityTextureData.maxRadius + extraRadius,
    );
    const surface = TextureMesh.create(
      gv.vertexCount,
      groupCount,
      gv.vertexTexture,
      gv.groupTexture,
      gv.normalTexture,
      boundingSphere,
      textureMesh,
    );
    (surface.meta as GaussianSurfaceMeta).resolution = densityTextureData.resolution;
    surface.meta.vertexBox = TextureMesh.getGridVertexBox(densityTextureData.gridDim, densityTextureData.transform);
    return surface;
  };

  return create(textureMesh);
}

export function GaussianSurfaceTextureMeshVisual(materialId: number): UnitsVisual<GaussianSurfaceMeshParams> {
  return UnitsTextureMeshVisual<GaussianSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(GaussianSurfaceMeshParams),
      createGeometry: createGaussianSurfaceTextureMesh,
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<GaussianSurfaceMeshParams>,
        currentProps: PD.Values<GaussianSurfaceMeshParams>,
      ) => {
        state.createGeometry =
          newProps.resolution !== currentProps.resolution ||
          newProps.radiusOffset !== currentProps.radiusOffset ||
          newProps.smoothness !== currentProps.smoothness ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly;

        if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
          state.updateColor = true;
        } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
          if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor)
            state.updateColor = true;
          if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride)
            state.updateColor = true;
        }
      },
      processValues: (
        values: TextureMeshValues,
        geometry: TextureMesh,
        props: PD.Values<GaussianSurfaceMeshParams>,
        theme: Theme,
        graphics?: GraphicsContext,
      ) => {
        const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
        if (csp && graphics) {
          applyTextureMeshColorSmoothing(values, csp, graphics, colorTexture);
          (geometry.meta as GaussianSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
        }
      },
      dispose: (geometry: TextureMesh) => {
        geometry.vertexTexture.ref.value.destroy();
        geometry.groupTexture.ref.value.destroy();
        geometry.normalTexture.ref.value.destroy();
        geometry.doubleBuffer.destroy();

        (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
      },
    },
    materialId,
  );
}

//

function createStructureGaussianSurfaceTextureMesh(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: GaussianSurfaceProps,
  textureMesh?: TextureMesh,
): TextureMesh | Promise<TextureMesh> {
  const graphics = ctx.graphics;
  if (graphics && computeSupport(graphics)) {
    const computeDensity = () =>
      computeStructureGaussianDensityTexture(
        structure,
        theme.size,
        props,
        graphics,
        getComputeDensityTexture(graphics),
        true,
      );
    return createComputeGaussianSurface(
      graphics,
      computeDensity,
      props,
      structure.boundary.sphere,
      structure.elementCount,
      textureMesh,
    );
  }
  const atlas = graphics && getAtlasIsosurfaceExtraction(graphics);
  if (!graphics || !atlas) {
    throw new Error('compute kernels or fragment-shader marching cubes required to create a gaussian surface');
  }

  const {
    namedTextures,
    resources,
    capabilities: { floatBlend, halfFloatRenderTarget },
  } = graphics;
  if (!namedTextures[GaussianSurfaceName]) {
    namedTextures[GaussianSurfaceName] = halfFloatRenderTarget
      ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
      : floatBlend
        ? resources.texture('image-float32', 'rgba', 'float', 'linear')
        : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
  }

  const axisOrder = Vec3.create(0, 1, 2);
  const groupCount = structure.elementCount;

  const create = (textureMesh?: TextureMesh) => {
    if (isTimingMode) graphics.timer.group('createStructureGaussianSurfaceTextureMesh');
    const densityTextureData = computeStructureGaussianDensityTexture2d(
      structure,
      theme.size,
      true,
      props,
      graphics,
      namedTextures[GaussianSurfaceName],
    );
    const isoLevel = Math.exp(-props.smoothness) / densityTextureData.radiusFactor;

    const buffer = textureMesh?.doubleBuffer.get(graphics);
    const gv = atlas.extract(
      graphics,
      {
        texture: densityTextureData.texture,
        channel: 'alpha',
        gridDim: densityTextureData.gridDim,
        gridTexDim: densityTextureData.gridTexDim,
        gridDataDim: densityTextureData.gridDataDim,
        gridTexScale: densityTextureData.gridTexScale,
        transform: densityTextureData.transform,
        isoValue: isoLevel,
        invert: false,
        packedGroup: true,
        axisOrder: axisOrder,
        constantGroup: true,
      },
      buffer,
    );
    if (isTimingMode) graphics.timer.groupEnd('createStructureGaussianSurfaceTextureMesh');

    const extraRadius = props.radiusOffset * (1 + Math.exp(-props.smoothness));
    const boundingSphere = Sphere3D.expand(
      Sphere3D(),
      structure.boundary.sphere,
      densityTextureData.maxRadius + extraRadius,
    );
    const surface = TextureMesh.create(
      gv.vertexCount,
      groupCount,
      gv.vertexTexture,
      gv.groupTexture,
      gv.normalTexture,
      boundingSphere,
      textureMesh,
    );
    (surface.meta as GaussianSurfaceMeta).resolution = densityTextureData.resolution;
    surface.meta.vertexBox = TextureMesh.getGridVertexBox(densityTextureData.gridDim, densityTextureData.transform);
    return surface;
  };

  return create(textureMesh);
}

export function StructureGaussianSurfaceTextureMeshVisual(
  materialId: number,
): ComplexVisual<StructureGaussianSurfaceMeshParams> {
  return ComplexTextureMeshVisual<StructureGaussianSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(StructureGaussianSurfaceMeshParams),
      createGeometry: createStructureGaussianSurfaceTextureMesh,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<StructureGaussianSurfaceMeshParams>,
        currentProps: PD.Values<StructureGaussianSurfaceMeshParams>,
      ) => {
        state.createGeometry =
          newProps.resolution !== currentProps.resolution ||
          newProps.radiusOffset !== currentProps.radiusOffset ||
          newProps.smoothness !== currentProps.smoothness ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly;

        if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
          state.updateColor = true;
        } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
          if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor)
            state.updateColor = true;
          if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride)
            state.updateColor = true;
        }
      },
      processValues: (
        values: TextureMeshValues,
        geometry: TextureMesh,
        props: PD.Values<GaussianSurfaceMeshParams>,
        theme: Theme,
        graphics?: GraphicsContext,
      ) => {
        const { resolution, colorTexture } = geometry.meta as GaussianSurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
        if (csp && graphics) {
          applyTextureMeshColorSmoothing(values, csp, graphics, colorTexture);
          (geometry.meta as GaussianSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
        }
      },
      dispose: (geometry: TextureMesh) => {
        geometry.vertexTexture.ref.value.destroy();
        geometry.groupTexture.ref.value.destroy();
        geometry.normalTexture.ref.value.destroy();
        geometry.doubleBuffer.destroy();

        (geometry.meta as GaussianSurfaceMeta).colorTexture?.destroy();
      },
    },
    materialId,
  );
}
