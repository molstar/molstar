/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';
import {
  type ParticleList,
  Particle,
  ParticleTarget,
  ParticleValues,
  writeParticleTransform,
} from '@molstar/model/model/particles/particle-list';
import { ParticleTileView } from '@molstar/model/model/particles/particle-tile-view';
import type { Theme } from '@molstar/graphics/theme/theme';
import { type GraphicsRenderObject, createRenderObject } from '@molstar/graphics/gl/render-object';
import { type Loci as ModelLoci, isEveryLoci } from '@molstar/model/model/loci';
import { NullLocation } from '@molstar/model/model/location';
import type { MarkerAction } from '@molstar/core/util/marker-action';
import { VisualUpdateState } from '../../../util.js';
import { Visual, type VisualContext } from '../../../visual.js';
import type { Representation } from '../../../representation.js';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Geometry } from '@molstar/graphics/geo/geometry/geometry';
import { LocationIterator } from '@molstar/graphics/geo/util/location-iterator';
import {
  createTransform,
  createIdentityTransform,
  createPosesTransform,
  moveTransformData,
  type TransformData,
} from '@molstar/graphics/geo/geometry/transform-data';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { createPendingInstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { type ColorData, createColors, createObjectColor } from '@molstar/graphics/geo/geometry/color-data';
import { createSizes, getMaxSize, type SizeData } from '@molstar/graphics/geo/geometry/size-data';
import { createMarkers, type MarkerData } from '@molstar/graphics/geo/geometry/marker-data';
import { Interval, OrderedSet, type SortedArray } from '@molstar/core/data/int';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import {
  calculateTransformDataBoundingSphere,
  setInstanceSphereFast,
  setTransformDataCenters,
} from '@molstar/graphics/gl/renderable/util';
import { isPromiseLike } from '@molstar/core/util/type-helpers';
import type { SizeValues } from '@molstar/graphics/gl/renderable/schema';
import { BaseGeometry, resolveGpuCulling } from '@molstar/graphics/geo/geometry/base';
import { canCullOnGpu } from '@molstar/graphics/gl/renderable/merged';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { Color } from '@molstar/core/util/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { UniformSizeTheme } from '@molstar/graphics/theme/size/uniform';
import {
  createStructureTargetGeometry,
  createStructureTargetSizes,
  type StructureTargetProps,
  structureTargetGeometryPropsChanged,
  structureTargetMustRecreate,
  structureTargetSizeFactorChanged,
  structureTargetUseImpostor,
} from './structure.js';
import {
  blobSurfaceTargetGeometryPropsChanged,
  type BlobSurfaceTargetProps,
  createBlobSurfaceTargetGeometry,
} from './blob-surface.js';
import {
  createSitesBlobSurfaceGeometry,
  createSitesTargetGeometry,
  createSitesTargetSizes,
  setSitesTargetSiteColor,
} from './sites.js';
import { SiteColorParams, type SiteColorProps, setStructureTargetSiteColor } from './site-color.js';
import type { SiteColorValues } from '@molstar/graphics/gl/renderable/spheres';
import { createShapeTargetGeometry } from './shape.js';
import {
  createVolumeTargetGeometry,
  type VolumeTargetProps,
  volumeTargetGeometryPropsChanged,
  volumeTargetMustRecreate,
  volumeTargetSizeFactorChanged,
  volumeTargetUseImpostor,
} from './volume.js';
import type { GeometryUtils } from '@molstar/graphics/geo/geometry/geometry';
import { VolumeIsosurfaceParams } from '../../../volume/isosurface.js';
import { BlobSurfaceCoreParams } from '../../../structure/visual/util/blob-surface.js';
import { getElementCountQuality, getQualityProps } from '../../../util.js';

// ---- Params ----------------------------------------------------------------

/** Removes the params every target kind shares (see `ParticleTargetCommonParams`). */
function withoutCommon<P extends PD.Params>(params: P): P {
  const out: PD.Params = { ...params };
  for (const k of Object.keys(ParticleTargetCommonParams)) delete out[k];
  return out as P;
}

/** Copies `params`, hiding each of them unless the group's `type` is `type`. */
function forType<P extends PD.Params>(type: string, params: P): P {
  const out: PD.Params = {};
  for (const k of Object.keys(params)) {
    const hideIf = params[k].hideIf;
    out[k] = { ...params[k], hideIf: (p: { type: string }) => p.type !== type || !!hideIf?.(p) };
  }
  return out as P;
}

// Quality is resolved per target (see `createParticleTargetPropsProvider`) and is therefore
// part of the per-kind groups, where it also drives the hiding of custom quality params.
const { quality, ...BaseParams } = BaseGeometry.Params;

/** Params shared by all target kinds; the per-kind groups below are merged on top of these. */
export const ParticleTargetCommonParams = {
  ...BaseParams,
  instanceGranularity: PD.Boolean(true, { isHidden: true }), // groupCount is always 1, so granularity is moot
  lodPixelSize: Spheres.getLodPixelSizeParam(Spheres.LodPixelPresets.quality),
};
export type ParticleTargetCommonParams = typeof ParticleTargetCommonParams;

/** Params for targets that are structures, rendered as sphere impostors, a sphere mesh or a blob surface. */
export const ParticleTargetStructureParams = {
  ...withoutCommon(Spheres.Params),
  ...withoutCommon(Mesh.Params),
  solidInterior: PD.Boolean(true, BaseGeometry.ShadingCategory),
  type: PD.Select<string>('spacefill', PD.arrayToOptions(['spacefill', 'blob-surface']), {
    isEssential: true,
    description: 'How the elements of the target are rendered.',
  }),
  ...forType('spacefill', {
    sizeFactor: PD.Numeric(1, { min: 0.01, max: 10, step: 0.01 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true),
    // targets of sites and structures: their atoms by element near the camera
    ...SiteColorParams,
  }),
  ...forType('blob-surface', {
    blobSize: BlobSurfaceCoreParams.blobSize,
    blobMethod: BlobSurfaceCoreParams.blobMethod,
    resolution: BlobSurfaceCoreParams.resolution,
    adjustResolution: BlobSurfaceCoreParams.adjustResolution,
    blobShape: BlobSurfaceCoreParams.blobShape,
    radiusOffset: BlobSurfaceCoreParams.radiusOffset,
    smoothness: BlobSurfaceCoreParams.smoothness,
  }),
};
export type ParticleTargetStructureParams = typeof ParticleTargetStructureParams;

/** Params for targets that are shapes, rendered using the shape's own geometry. */
export const ParticleTargetShapeParams = {
  scaleByRadius: PD.Boolean(true, {
    description: 'Scale each instance of the shape by the radius of the particle it is placed at.',
  }),
  ...withoutCommon(Mesh.Params),
};
export type ParticleTargetShapeParams = typeof ParticleTargetShapeParams;

/** Params for targets that are volumes, rendered as an isosurface mesh or as dots. */
export const ParticleTargetVolumeParams = {
  ...withoutCommon(Spheres.Params),
  ...withoutCommon(Mesh.Params),
  solidInterior: PD.Boolean(true, BaseGeometry.ShadingCategory),
  type: PD.Select<string>('isosurface', PD.arrayToOptions(['isosurface', 'dot']), {
    isEssential: true,
    description: 'How the volume of the target is rendered.',
  }),
  isoValue: VolumeIsosurfaceParams.isoValue,
  ...forType('isosurface', {
    wrap: VolumeIsosurfaceParams.wrap,
    floodfill: VolumeIsosurfaceParams.floodfill,
  }),
  ...forType('dot', {
    perturbPositions: PD.Boolean(false),
    sizeFactor: PD.Numeric(
      1,
      { min: 0, max: 10, step: 0.1 },
      {
        description: 'Radius of a dot. Dots are sized uniformly, the representation size theme does not apply to them.',
      },
    ),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true),
  }),
  scaleByRadius: PD.Boolean(false, {
    description: 'Scale each instance of the volume geometry by the radius of the particle it is placed at.',
  }),
  // dots: levels by size on screen for the dots' own size
  lodRadius: Spheres.getLodRadiusParam(0),
};
export type ParticleTargetVolumeParams = typeof ParticleTargetVolumeParams;

export const ParticleTargetRepresentationParams = {
  ...ParticleTargetCommonParams,
  occluder: PD.Select<'auto' | 'compute' | 'off'>(
    'auto',
    [
      ['auto', 'Given'],
      ['compute', 'Computed where none'],
      ['off', 'Off'],
    ],
    {
      label: 'Occluder Grid',
      description:
        "What the particles hide as optical densities on a grid (the GPU's occlusion primed from it: what lies behind dense content culled through the gaps between molecules). Given: a pack's or its tiles'; computed where none: once a list as it is drawn whole (about a second for 30 million particles); off: none. With clip objects, what they remove hides nothing; none while the object is moved.",
      isEssential: false,
    },
  ),
  structure: PD.Group(ParticleTargetStructureParams, {
    isExpanded: false,
    description: 'Parameters applied to targets that are structures.',
  }),
  shape: PD.Group(ParticleTargetShapeParams, {
    isExpanded: false,
    description: 'Parameters applied to targets that are shapes.',
  }),
  volume: PD.Group(ParticleTargetVolumeParams, {
    isExpanded: false,
    description: 'Parameters applied to targets that are volumes.',
  }),
};
export type ParticleTargetRepresentationParams = typeof ParticleTargetRepresentationParams;
export type ParticleTargetRepresentationProps = PD.Values<ParticleTargetRepresentationParams>;

/**
 * The flat props used for a single target: the common props with the group of the target's
 * kind merged on top. Only the params of the matching kind are present, so kind-specific
 * props are optional here.
 */
export type ParticleTargetProps = PD.Values<ParticleTargetCommonParams> &
  Partial<PD.Values<ParticleTargetStructureParams>> &
  Partial<PD.Values<ParticleTargetShapeParams>> &
  Partial<PD.Values<ParticleTargetVolumeParams>>;

/** Resolves the flat props for a single target, see `createParticleTargetPropsProvider`. */
export type ParticleTargetPropsProvider = (target: ParticleTarget) => ParticleTargetProps;

/**
 * Create a provider that flattens the representation props for a single target: common props
 * with the group for the target's kind merged on top, with quality-derived props resolved
 * against the target's own data (e.g. sphere detail from the number of spheres).
 *
 * The per-kind flattening is done once up front and the quality resolution is memoized per
 * target data object, since targets are commonly shared between many particles and the work
 * would otherwise be repeated for every target on every update.
 *
 * Create one provider per update; it holds the props by value.
 */
export function createParticleTargetPropsProvider(
  props: ParticleTargetRepresentationProps,
): ParticleTargetPropsProvider {
  const { structure, shape, volume, ...common } = props;
  const bases = {
    structure: { ...common, ...structure } as ParticleTargetProps,
    // sites drawn as a structure's elements are
    sites: { ...common, ...structure } as ParticleTargetProps,
    shape: { ...common, ...shape } as ParticleTargetProps,
    volume: { ...common, ...volume } as ParticleTargetProps,
    // a template's: never drawn as one (its leaves are, `templateNotDrawn`), keyed alike
    particles: common as ParticleTargetProps,
  };
  const cache = new Map<unknown, ParticleTargetProps>();
  // the targets' props of one kind and quality one object (thousands of targets, a few qualities): what is keyed by
  // them sees the few (`computeParticleTargetSignatures`)
  const interned = new Map<string, ParticleTargetProps>();

  return (target) => {
    const data = ParticleTarget.data(target);
    let resolved = cache.get(data);
    if (!resolved) {
      const flat = bases[target.kind];
      // sites' auto quality by their count, as a structure's by its elements
      const quality =
        target.kind === 'sites' && flat.quality === 'auto'
          ? { ...flat, quality: getElementCountQuality(target.sites.count) }
          : flat;
      const qualityProps = getQualityProps(quality, data);
      const key = `${target.kind}|${JSON.stringify(qualityProps)}`;
      resolved = interned.get(key);
      if (!resolved) {
        resolved = { ...flat, ...qualityProps };
        interned.set(key, resolved);
      }
      cache.set(data, resolved);
    }
    return resolved;
  };
}

/**
 * Whether the occluder grid (the GPU's occlusion primed from the particles' densities, `occluder`) is off: as asked.
 * With clip objects it is on: the grid's rays skip what they cut away (the object's clip objects, `HiZPass.setOccluder`)
 */
export function isTargetOccluderOff(props: Pick<ParticleTargetRepresentationProps, 'occluder'>) {
  return props.occluder === 'off';
}

// ---- Transform helpers -----------------------------------------------------

/**
 * The instances a target's grid cell holds at least on average where the target is culled on the GPU
 * (`getTargetCellSize`): a target spread over the whole scene (most of a mesoscale model's, a few of each over the cell)
 * has about a cell an instance at the default cell size (the yeast cell: 20 M cells for 31.5 M instances, 551 MB). The
 * GPU culls by the scene's cells (`buildSceneCells`), the grid serving only the frames culled on the CPU (a second view,
 * kernels compiling). Not where the CPU culls by the grid: coarse cells cull less (MG_1184 zoomed: the opaque pass 6.3 ->
 * 18.6 ms), though they draw fewer times (whole: 15.7 -> 5.2 ms a frame)
 */
const MinTargetCellInstances = 32;

/**
 * The cell size of a target's grid: `cellSize`, larger where its instances are sparse, so that its cells hold about
 * `MinTargetCellInstances` on average over the box of their spheres' centers (`min`, `max`)
 */
export function getTargetCellSize(cellSize: number, instanceCount: number, min: Vec3, max: Vec3) {
  if (instanceCount <= MinTargetCellInstances || !(max[0] >= min[0])) return cellSize;
  const volume =
    Math.max(max[0] - min[0], cellSize) * Math.max(max[1] - min[1], cellSize) * Math.max(max[2] - min[2], cellSize);
  const perCell = (instanceCount * cellSize ** 3) / volume;
  return perCell >= MinTargetCellInstances ? cellSize : Math.cbrt((volume * MinTargetCellInstances) / instanceCount);
}

/**
 * Whether a target's visual is drawn in a merged object culled on the GPU by its instances: the merge's rule
 * (`mergeRenderObjectsByType`) applied ahead, the representation's particles over its targets' draws (a target and
 * level each). Not an object merged for the sphere raster alone, which the CPU culls where the raster does not take it
 */
export function isTargetCulledOnGpu(
  particles: ParticleList,
  props: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
) {
  if (props.renderMerged === false || !canCullOnGpu(graphics)) return false;
  const levels = Math.max(1, Spheres.getLodPixelLevelCount(props.lodPixelSize));
  return resolveGpuCulling(props.gpuCulling, particles.count, Math.max(1, particles.targetInfo.size) * levels);
}

/**
 * Particles of a list from which on its targets' bounding spheres are the box sphere at any count
 * (`instanceSphereFast`): the boundary helper's tighter spheres of a few instances cost more than they cull in a list of
 * millions (the yeast cell's 5,639 targets), and with a pack's bounds they are seeded without a pass
 */
export const LargeParticleListCount = 1_000_000;

/**
 * Build per-instance transforms placing the target geometry at each of `indices`.
 *
 * The geometry is built at its actual coordinates, not at the origin. Each particle transform is
 * therefore composed with a pre-translation by `-center` so that the target's `center` (a structure's
 * centroid, a volume's box center, `ParticleTarget.getCenter`; `invariantBoundingSphere.center` when
 * not given) renders at
 * the particle position. When `scaleByRadius`
 * is set (e.g. for shape targets), each instance is additionally scaled by the particle radius;
 * a coarse sample's particles drawn larger too (`ParticleList.scale`, or a run's own: `scales`), as the page store's;
 * the scaling is applied before the center subtraction so the basis already accounts for it.
 * The bounding sphere used for grid/culling purposes is likewise inflated by the largest radius
 * in play, since it is otherwise unaware of the per-instance scale (see `calcInstanceGrid`).
 * `culledOnGpu` (`isTargetCulledOnGpu`): the grid's cells grown where the target is sparse (`getTargetCellSize`).
 * `poses` (`isPosesGeometry`): the instances as the particles' poses, no matrix made (`createPosesTransform`).
 */
export function createTargetParticleTransform(
  particles: ParticleList,
  indices: OrderedSet<number>,
  invariantBoundingSphere: Sphere3D,
  cellSize: number,
  batchSize: number,
  scaleByRadius: boolean,
  transformData?: TransformData,
  center: Vec3 = invariantBoundingSphere.center,
  culledOnGpu = false,
  poses = false,
  /** the same particles posed otherwise (`Particle.withPoses`): a pending grid's cell size kept (what lays them out) */
  posed = false,
): TransformData {
  if (poses) {
    return createTargetParticlePoses(
      particles,
      indices,
      invariantBoundingSphere,
      cellSize,
      batchSize,
      scaleByRadius,
      transformData,
      center,
      culledOnGpu,
      posed,
    );
  }
  const instanceCount = OrderedSet.size(indices);

  // reuse the previous output array when updating in place (createTransform keeps it), a new one handed over; not the
  // grid's one copy (`oneCopy` below: its array the grid's and the merged object's, read until they are rebuilt)
  const prev =
    transformData && transformData.transform.ref.value !== transformData.aTransform.ref.value
      ? transformData.transform.ref.value
      : undefined;
  const transformArray = prev && prev.length >= instanceCount * 16 ? prev : new Float32Array(instanceCount * 16);

  // resolve the OrderedSet union once instead of per-element getAt dispatch in the hot loop
  let sortedIndices: SortedArray | undefined;
  let startIdx = 0;
  if (Interval.is(indices)) startIdx = Interval.start(indices);
  else sortedIndices = indices;

  const radiiArr = scaleByRadius ? particles.radii : undefined;
  const listScales = particles.scales;
  const listScale = getListScale(particles);
  const cx = center[0],
    cy = center[1],
    cz = center[2];
  const hasCenter = cx !== 0 || cy !== 0 || cz !== 0;
  // the box of the instances' spheres' centers, for the grid's cell size (`getTargetCellSize`)
  const [ix, iy, iz] = invariantBoundingSphere.center;
  const boxMin = Vec3.create(Infinity, Infinity, Infinity);
  const boxMax = Vec3.create(-Infinity, -Infinity, -Infinity);

  // each particle's transform written here from its position and rotation (`writeParticleTransform`), not copied from
  // the list's transforms of all particles (`getParticleTransforms`: a copy of them all kept)
  let maxRadius = 1;
  for (let i = 0; i < instanceCount; i++) {
    const pi = sortedIndices ? sortedIndices[i] : startIdx + i;
    const o = i * 16;
    writeParticleTransform(particles, pi, transformArray, o);

    // a run's value read again for the run's next particle (`ParticleValues.get`): no search
    let s = listScales ? ParticleValues.get(listScales, pi) : listScale;
    if (radiiArr) {
      const r = ParticleValues.get(radiiArr, pi);
      if (r > 0) s *= r;
    }
    if (s > maxRadius) maxRadius = s;
    if (s !== 1) for (let j = 0; j < 12; j++) transformArray[o + j] *= s;

    if (hasCenter) {
      transformArray[o + 12] -= transformArray[o + 0] * cx + transformArray[o + 4] * cy + transformArray[o + 8] * cz;
      transformArray[o + 13] -= transformArray[o + 1] * cx + transformArray[o + 5] * cy + transformArray[o + 9] * cz;
      transformArray[o + 14] -= transformArray[o + 2] * cx + transformArray[o + 6] * cy + transformArray[o + 10] * cz;
    }
    const t = transformArray;
    const x = t[o] * ix + t[o + 4] * iy + t[o + 8] * iz + t[o + 12];
    const y = t[o + 1] * ix + t[o + 5] * iy + t[o + 9] * iz + t[o + 13];
    const z = t[o + 2] * ix + t[o + 6] * iy + t[o + 10] * iz + t[o + 14];
    if (x < boxMin[0]) boxMin[0] = x;
    if (y < boxMin[1]) boxMin[1] = y;
    if (z < boxMin[2]) boxMin[2] = z;
    if (x > boxMax[0]) boxMax[0] = x;
    if (y > boxMax[1]) boxMax[1] = y;
    if (z > boxMax[2]) boxMax[2] = z;
  }

  // `calcInstanceGrid` (used inside `createTransform`) builds its grid/culling bounds from
  // `invariantBoundingSphere.radius` alone and does not account for any scale baked into the
  // per-instance transform matrices. Since scaling (by radius, a coarse sample's) applies such a scale above, the
  // sphere passed on for grid/culling purposes must be inflated by the largest scale used,
  // otherwise larger particles get under-sized cells/LOD bounds and can be culled incorrectly.
  const gridBoundingSphere =
    maxRadius !== 1
      ? Sphere3D.create(invariantBoundingSphere.center, invariantBoundingSphere.radius * maxRadius)
      : invariantBoundingSphere;

  const data = createTransform(
    transformArray,
    instanceCount,
    gridBoundingSphere,
    culledOnGpu ? getTargetCellSize(cellSize, instanceCount, boxMin, boxMax) : cellSize,
    batchSize,
    transformData,
    transformArray !== prev,
    // the transforms kept once, in the grid's order (millions of particles)
    true,
  );
  setInstanceSphereFast(data, particles.count >= LargeParticleListCount);
  return data;
}

/** The geometries whose readers take poses: their targets' instances placed by poses */
export function isPosesGeometry(kind: string) {
  return kind === 'spheres' || kind === 'mesh';
}

/**
 * The scales of the particles `indices` by radius (`InstancePoses`): the radii as they are (an array, a scale an
 * element), the one value of their run where they are runs (a pack's target: a radius a target, `ParticleValues`),
 * the runs as an array where `indices` span several
 */
/** a list's scale for all its particles (`ParticleList.scale`), 1 where none or not positive */
function getListScale(particles: ParticleList) {
  return particles.scale !== undefined && particles.scale > 0 ? particles.scale : 1;
}

/**
 * the instances' scales by poses (`InstanceTransforms.ofPoses`: one for all, or per particle): their radii's where
 * `scaleByRadius`, times a coarse sample's (`ParticleList.scale`, or its runs': `scales`); none where all are 1
 */
function getInstanceScales(
  particles: ParticleList,
  indices: OrderedSet<number>,
  scaleByRadius: boolean,
): { scales?: Float32Array; scale?: number } {
  const radius = scaleByRadius ? getRadiusScales(particles, indices) : {};
  const listScales = particles.scales;
  const listScale = getListScale(particles);
  const n = OrderedSet.size(indices);
  const byRadius = scaleByRadius && !!particles.radii;
  if (listScales && n > 0) {
    if (ParticleValues.isRuns(listScales)) {
      const first = ParticleValues.runOf(listScales, OrderedSet.getAt(indices, 0));
      const last = ParticleValues.runOf(listScales, OrderedSet.getAt(indices, n - 1));
      if (first === last && !radius.scales) return { scale: listScales.values[first] * (radius.scale ?? 1) };
    }
    return { scales: getParticleScales(particles, byRadius) };
  }
  if (listScale === 1) return radius;
  if (radius.scales) return { scales: getParticleScales(particles, byRadius) };
  return { scale: (radius.scale ?? 1) * listScale };
}

/** per list its particles' scales (`getParticleScales`), without and with their radii */
const particleScales = [new WeakMap<ParticleList, Float32Array>(), new WeakMap<ParticleList, Float32Array>()];

/**
 * every particle's scale of `particles` (a coarse sample's: `ParticleList.scales` or `scale`), times its radius where
 * `byRadius`: one array a list (and kind), its targets' poses sharing it. Why: poses read a scale by particle (the
 * list's index), an array a target the list's size each (thousands of targets in a tile view)
 */
function getParticleScales(particles: ParticleList, byRadius: boolean) {
  const cache = particleScales[byRadius ? 1 : 0];
  let out = cache.get(particles);
  if (out) return out;
  const count = particles.count;
  const scales = particles.scales;
  const listScale = getListScale(particles);
  out = new Float32Array(count);
  if (scales) out.set(ParticleValues.toArray(scales, count, Float32Array).subarray(0, count));
  else out.fill(listScale);
  if (byRadius) {
    const radii = ParticleValues.toArray(particles.radii!, count, Float32Array);
    for (let i = 0; i < count; ++i) out[i] *= radii[i];
  }
  cache.set(particles, out);
  return out;
}

function getRadiusScales(
  particles: ParticleList,
  indices: OrderedSet<number>,
): { scales?: Float32Array; scale?: number } {
  const radii = particles.radii;
  if (!radii) return {};
  if (!ParticleValues.isRuns(radii)) return { scales: radii };
  const n = OrderedSet.size(indices);
  if (n === 0) return {};
  const first = ParticleValues.runOf(radii, OrderedSet.getAt(indices, 0));
  const last = ParticleValues.runOf(radii, OrderedSet.getAt(indices, n - 1));
  if (first === last) return { scale: radii.values[first] };
  return { scales: ParticleValues.toArray(radii, particles.count, Float32Array) };
}

/**
 * `createTargetParticleTransform` by poses: the particles' positions, rotations and radii (`scaleByRadius`) read where
 * they are, a coarse sample's scales with them (`getInstanceScales`), the instances' particles the range or the indices
 * of `indices`, no matrix an instance made
 */
function createTargetParticlePoses(
  particles: ParticleList,
  indices: OrderedSet<number>,
  invariantBoundingSphere: Sphere3D,
  cellSize: number,
  batchSize: number,
  scaleByRadius: boolean,
  transformData: TransformData | undefined,
  center: Vec3,
  culledOnGpu: boolean,
  posed = false,
): TransformData {
  const instanceCount = OrderedSet.size(indices);
  const range = Interval.is(indices);
  const instanceScales = getInstanceScales(particles, indices, scaleByRadius);
  const scaled = instanceScales.scales !== undefined || instanceScales.scale !== undefined;
  const poses = InstanceTransforms.ofPoses({
    count: instanceCount,
    positions: particles.coordinates,
    rotations: particles.rotations,
    scales: instanceScales.scales,
    scale: instanceScales.scale,
    first: range ? Interval.start(indices) : undefined,
    indices: range ? undefined : (indices as SortedArray<number>),
    center: Vec3.clone(center),
  });
  // the grid's sphere inflated by the largest radius, as the matrices' (`createTargetParticleTransform`)
  const maxScale = scaled ? InstanceTransforms.getMaxScaleOfAll(poses) : 1;
  const maxRadius = Math.max(1, maxScale);
  // the box of the instances' centers from their particles' bounds where they are a target's (a pack's,
  // `Particle.getRangeBounds`): no pass over them. An instance's center is its particle's position moved by its scale
  // times the invariant center's offset from the placed point
  // (a tile-ordered pack's target its runs' united, `getIndicesBounds`)
  const bounds = Particle.getIndicesBounds(particles, indices);
  const offset = maxScale * Vec3.distance(invariantBoundingSphere.center, center);
  const gridBoundingSphere =
    maxRadius !== 1
      ? Sphere3D.create(invariantBoundingSphere.center, invariantBoundingSphere.radius * maxRadius)
      : invariantBoundingSphere;
  let targetCellSize = cellSize;
  if (culledOnGpu && posed && transformData?.instanceGrid.ref.value.pending) {
    // posed otherwise: the cells' size kept, its grid pending as it is (the culler's layout holds)
    targetCellSize = transformData.instanceCellSize.ref.value;
  } else if (culledOnGpu && bounds) {
    const boxMin = Vec3.subScalar(Vec3(), bounds.min, offset);
    const boxMax = Vec3.addScalar(Vec3(), bounds.max, offset);
    targetCellSize = getTargetCellSize(cellSize, instanceCount, boxMin, boxMax);
  } else if (culledOnGpu) {
    const boxMin = Vec3();
    const boxMax = Vec3();
    InstanceTransforms.writePoints(
      poses,
      invariantBoundingSphere.center,
      undefined,
      undefined,
      undefined,
      boxMin,
      boxMax,
    );
    targetCellSize = getTargetCellSize(cellSize, instanceCount, boxMin, boxMax);
  }
  // culled on the GPU: the grid pending, built where a CPU cull needs it
  const data = createPosesTransform(poses, gridBoundingSphere, targetCellSize, batchSize, transformData, culledOnGpu);
  setInstanceSphereFast(data, particles.count >= LargeParticleListCount);
  // no grid's pass over the centers: their sphere from the particles' bounds, the box sphere many instances take
  // (`setTransformDataCenters`): exact where the placed point is the invariant center, else within the offset (an
  // instance's center its particle's position moved by at most its scale times the offset)
  if (culledOnGpu && bounds) {
    const center = Vec3.scale(Vec3(), Vec3.add(Vec3(), bounds.min, bounds.max), 0.5);
    setTransformDataCenters(data, {
      invariantCenter: invariantBoundingSphere.center,
      center,
      distance: bounds.distance + offset,
    });
  }
  return data;
}

// ---- Per-target geometry dispatch ------------------------------------------

/**
 * Whether the theme's color is constant over the target's particles (`ColorTheme.particleRuns`: they lie in one run): a
 * color filled, not evaluated a particle (a pack's 31.5 M particles colored by entity)
 */
function isTargetColorConstant(theme: Theme, indices: OrderedSet<number>) {
  const color = theme.color;
  const runs = color.particleRuns;
  const n = OrderedSet.size(indices);
  if (!runs || color.granularity !== 'instance' || n === 0) return false;
  const first = ParticleValues.runOf(runs, OrderedSet.getAt(indices, 0));
  const last = ParticleValues.runOf(runs, OrderedSet.getAt(indices, n - 1));
  if (last === first) return true;
  // in several runs (a tile-ordered pack's target: a run a tile it is in): of one value where each run it is in is
  // (not the runs between, of other targets)
  const value = runs.values[first];
  let constant = true;
  ParticleValues.forEachRunOf(runs, indices, (r) => {
    if (runs.values[r] === value) return true;
    constant = false;
    return false;
  });
  return constant;
}

/**
 * The theme for the values' first colors (`createValues`): uniform where its color is constant over the target's
 * particles (`isTargetColorConstant`), a color evaluated once, then the object's (`createTargetColors`), else as it is
 */
export function getTargetTheme(theme: Theme, indices: OrderedSet<number>): Theme {
  if (!isTargetColorConstant(theme, indices)) return theme;
  return { ...theme, color: { ...(theme.color as ColorTheme<any>), granularity: 'uniform' } };
}

/**
 * The target's colors: a color constant over its particles the object's (`createObjectColor`: one evaluation, one
 * texel, a merged object keeping each member's; a uniform color would be the merged object's first member's for all),
 * else the theme's (`createColors`)
 */
export function createTargetColors(
  locationIt: LocationIterator,
  positionIt: LocationIterator,
  theme: Theme,
  indices: OrderedSet<number>,
  values: ColorData,
) {
  if (isTargetColorConstant(theme, indices)) {
    locationIt.reset();
    const color = (theme.color as ColorTheme<any>).color(locationIt.getInstanceLocation(0), false);
    return createObjectColor(color, values);
  }
  return createColors(locationIt, positionIt, theme.color, values);
}

/** whether colors or sizes of this granularity are read at positions (per vertex, from a volume): a move changes them */
function isPositional(granularity: string) {
  return granularity === 'vertex' || granularity === 'vertexInstance' || granularity.startsWith('volume');
}

/** The reference data backing a target, used for change detection. */
function targetData(target: ParticleTarget): unknown {
  return ParticleTarget.data(target);
}

/** Whether the target is drawn as spheres of its own, sized each (a structure's elements, sites) */
function isElementTarget(target: ParticleTarget): target is Extract<ParticleTarget, { kind: 'structure' | 'sites' }> {
  return target.kind === 'structure' || target.kind === 'sites';
}

/** Whether the target geometry is scaled by the per-particle radius. */
function targetScaleByRadius(target: ParticleTarget, props: ParticleTargetProps): boolean {
  return isElementTarget(target) ? false : !!props.scaleByRadius;
}

/** Whether a volume target is rendered as dots rather than as an isosurface. */
function isVolumeDot(target: ParticleTarget, props: ParticleTargetProps): boolean {
  return target.kind === 'volume' && (props as VolumeTargetProps).type === 'dot';
}

let _dotSizeTheme: SizeTheme<any> | undefined;
function getDotSizeTheme(): SizeTheme<any> {
  if (!_dotSizeTheme) _dotSizeTheme = UniformSizeTheme({} as any, { value: 1 });
  return _dotSizeTheme;
}

/**
 * The geometry a target is drawn with (`createTargetGeometry` makes it): structures as sphere impostors or a sphere
 * mesh, or a blob surface; volumes as an isosurface or dots (impostors or a mesh); shapes as their own geometry. What
 * the precompile signatures take without building it (`getTargetGeometryUtils`)
 */
export type TargetGeometryChoice =
  | 'structure-spheres'
  | 'structure-mesh'
  | 'blob-surface'
  | 'volume-isosurface'
  | 'volume-dot-spheres'
  | 'volume-dot-mesh'
  | 'shape';

/**
 * A template (`ParticleTarget` kind `particles`, `ParticleNesting`) is drawn by its leaves, never as a geometry of its
 * own: its list placed whole (its particles' targets drawn)
 */
function templateNotDrawn(): never {
  throw new Error('a particles target (a template) is drawn by its leaves, not as a target geometry');
}

export function getTargetGeometryChoice(
  target: ParticleTarget,
  props: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
): TargetGeometryChoice {
  switch (target.kind) {
    case 'structure':
    case 'sites':
      if ((props as StructureTargetProps).type === 'blob-surface') return 'blob-surface';
      return structureTargetUseImpostor(props as StructureTargetProps, graphics)
        ? 'structure-spheres'
        : 'structure-mesh';
    case 'shape':
      return 'shape';
    case 'volume':
      if ((props as VolumeTargetProps).type !== 'dot') return 'volume-isosurface';
      return volumeTargetUseImpostor(props as VolumeTargetProps, graphics) ? 'volume-dot-spheres' : 'volume-dot-mesh';
    case 'particles':
      return templateNotDrawn();
  }
}

/** The utils of the geometry `choice` makes for `target` (a shape's own geometry's) */
export function getTargetGeometryUtils(target: ParticleTarget, choice: TargetGeometryChoice): GeometryUtils<any> {
  switch (choice) {
    case 'structure-spheres':
    case 'volume-dot-spheres':
      return Spheres.Utils;
    case 'shape':
      return Geometry.getUtils((target as Extract<ParticleTarget, { kind: 'shape' }>).shape.geometry as Geometry);
    default:
      return Mesh.Utils;
  }
}

export function createTargetGeometry(
  ctx: VisualContext,
  target: ParticleTarget,
  theme: Theme,
  props: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
  existing?: Geometry,
): Geometry | Promise<Geometry> {
  const spheresOrMesh =
    existing?.kind === 'spheres' || existing?.kind === 'mesh' ? (existing as Spheres | Mesh) : undefined;
  // the choice the precompile signatures take too (`getTargetGeometryChoice`)
  const choice = getTargetGeometryChoice(target, props, graphics);
  switch (target.kind) {
    case 'structure':
      return choice === 'blob-surface'
        ? createBlobSurfaceTargetGeometry(
            ctx,
            target.structure,
            theme,
            props as BlobSurfaceTargetProps,
            existing?.kind === 'mesh' ? (existing as Mesh) : undefined,
          )
        : createStructureTargetGeometry(
            ctx,
            target.structure,
            theme,
            props as StructureTargetProps,
            graphics,
            spheresOrMesh,
          );
    case 'sites':
      return choice === 'blob-surface'
        ? createSitesBlobSurfaceGeometry(
            ctx,
            target.sites,
            theme,
            props as BlobSurfaceTargetProps,
            existing?.kind === 'mesh' ? (existing as Mesh) : undefined,
          )
        : createSitesTargetGeometry(target.sites, theme, props as StructureTargetProps, graphics, spheresOrMesh);
    case 'shape':
      return createShapeTargetGeometry(target.shape, existing);
    case 'volume':
      return createVolumeTargetGeometry(ctx, target.volume, theme, props as VolumeTargetProps, graphics, spheresOrMesh);
    case 'particles':
      return templateNotDrawn();
  }
}

export function targetGeometryPropsChanged(
  target: ParticleTarget,
  oldProps: ParticleTargetProps,
  newProps: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
): boolean {
  switch (target.kind) {
    case 'structure':
    case 'sites':
      return (newProps as StructureTargetProps).type === 'blob-surface'
        ? blobSurfaceTargetGeometryPropsChanged(oldProps as BlobSurfaceTargetProps, newProps as BlobSurfaceTargetProps)
        : structureTargetGeometryPropsChanged(
            oldProps as StructureTargetProps,
            newProps as StructureTargetProps,
            graphics,
          );
    case 'shape':
      return false;
    case 'volume':
      return volumeTargetGeometryPropsChanged(
        target.volume,
        oldProps as VolumeTargetProps,
        newProps as VolumeTargetProps,
        graphics,
      );
    case 'particles':
      return templateNotDrawn();
  }
}

export function targetMustRecreate(
  target: ParticleTarget,
  oldProps: ParticleTargetProps,
  newProps: ParticleTargetProps,
): boolean {
  switch (target.kind) {
    case 'structure':
    case 'sites':
      return structureTargetMustRecreate(oldProps as StructureTargetProps, newProps as StructureTargetProps);
    case 'shape':
      return false;
    case 'volume':
      return volumeTargetMustRecreate(oldProps as VolumeTargetProps, newProps as VolumeTargetProps);
    case 'particles':
      return templateNotDrawn();
  }
}

function targetSizeFactorChanged(
  target: ParticleTarget,
  oldProps: ParticleTargetProps,
  newProps: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
): boolean {
  switch (target.kind) {
    case 'structure':
    case 'sites':
      return structureTargetSizeFactorChanged(
        oldProps as StructureTargetProps,
        newProps as StructureTargetProps,
        graphics,
      );
    case 'shape':
      return false;
    case 'volume':
      return volumeTargetSizeFactorChanged(oldProps as VolumeTargetProps, newProps as VolumeTargetProps, graphics);
    case 'particles':
      return templateNotDrawn();
  }
}

/** The sizes of an element target's spheres, an element or site each */
/** An element target's sizes as `createElementTargetSizes` makes them for impostors `geom` (by the sites' types where
 * drawn from them), on their own */
function createElementSizes(
  target: Extract<ParticleTarget, { kind: 'structure' | 'sites' }>,
  theme: Theme,
  geom: Geometry,
): SizeData {
  if (target.kind === 'structure') return createStructureTargetSizes(target.structure, theme);
  const dSites = (geom as Spheres).shaderData?.dSites?.ref.value;
  return createSitesTargetSizes(target.sites, theme, undefined, !!dSites);
}

function createElementTargetSizes(
  target: Extract<ParticleTarget, { kind: 'structure' | 'sites' }>,
  theme: Theme,
  sizeData: SizeValues,
) {
  if (target.kind === 'structure') createStructureTargetSizes(target.structure, theme, sizeData);
  // by the sites' types where the spheres are drawn from them
  else
    createSitesTargetSizes(
      target.sites,
      theme,
      sizeData,
      !!(sizeData as { dSites?: ValueCell<boolean> }).dSites?.ref.value,
    );
}

/**
 * The second color of a sites target's spheres (`setSitesTargetSiteColor`): by the sites' types where the spheres are
 * drawn from them (`dSites`); nothing for other targets or geometries. The shaders blend it into a texture's color
 * (`dColorType` texture): a uniform one the object's (`createObjectColor`) where it is on
 */
function setTargetSiteColor(target: ParticleTarget, geom: Geometry, props: ParticleTargetProps, values: unknown) {
  if (geom.kind !== 'spheres' || (target.kind !== 'sites' && target.kind !== 'structure')) return;
  const v = values as SiteColorValues & ColorData & { dSites?: ValueCell<boolean> };
  const siteProps = props as unknown as SiteColorProps;
  if (target.kind === 'sites') setSitesTargetSiteColor(target.sites, siteProps, v, !!v.dSites?.ref.value);
  else setStructureTargetSiteColor(target.structure, siteProps, v);
  if (v.dSiteColor.ref.value !== 'off' && v.colorType.ref.value === 'uniform') {
    createObjectColor(Color.fromNormalizedArray(v.uColor.ref.value, 0), v);
  }
}

/** Whether the second color's props changed (`SiteColorParams`) */
function siteColorPropsChanged(a: ParticleTargetProps, b: ParticleTargetProps) {
  const pa = (a as unknown as Partial<SiteColorProps>).siteColor;
  const pb = (b as unknown as Partial<SiteColorProps>).siteColor;
  return !PD.areEqual(SiteColorParams, { siteColor: pa }, { siteColor: pb });
}

/** Whether a size theme change requires rebuilding the geometry (sphere meshes bake in their radii). */
function targetSizeThemeAffectsGeometry(
  target: ParticleTarget,
  props: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
): boolean {
  return isElementTarget(target) && !structureTargetUseImpostor(props as StructureTargetProps, graphics);
}

// ---- Per-target visual -----------------------------------------------------

/**
 * A single visual for one target (structure or shape), managing one renderObject.
 * Geometry = the reference geometry in its invariant frame.
 * Instances = the particles of this target, given as indices into the unfiltered `ParticleList`.
 */
export interface TargetVisual {
  renderObject: GraphicsRenderObject | undefined;
  geometryVersion: number;
  /** Maps local instance index → global particle index in the unfiltered ParticleList. */
  particleIndices: OrderedSet<number> | undefined;
  createOrUpdate(
    ctx: VisualContext,
    theme: Theme,
    props: ParticleTargetProps,
    particles: ParticleList,
    particleIndices: OrderedSet<number>,
    target: ParticleTarget,
  ): Promise<void>;
  /** the representation's transform, on top of the particles' */
  setTransform(matrix: Mat4): void;
  mark(loci: ModelLoci, action: MarkerAction): boolean;
  /** the layers given (each all of its own: cleared, then its layers' loci applied), the others kept */
  setLayers(layers: TargetLayers): void;
  destroy(): void;
}

/** the layers a target's visual takes (`Representation.State`'s), a loci of each by the list's particles */
export type TargetLayers = Partial<
  Pick<Representation.State, 'overpaint' | 'transparency' | 'emissive' | 'substance' | 'clipping' | 'wiggle'>
>;

export const TargetLayerKeys = ['overpaint', 'transparency', 'emissive', 'substance', 'clipping', 'wiggle'] as const;

/** `layers` with each layer's loci `map`ped (a layer whose loci maps to none left out) */
export function mapTargetLayers(layers: TargetLayers, map: (loci: ModelLoci) => ModelLoci | undefined): TargetLayers {
  const mapped: TargetLayers = {};
  for (const k of TargetLayerKeys) {
    const value = layers[k] as { layers: ReadonlyArray<{ loci: ModelLoci }> } | undefined;
    if (!value) continue;
    const out: { loci: ModelLoci }[] = [];
    for (const layer of value.layers) {
      const loci = map(layer.loci);
      if (loci) out.push({ ...layer, loci });
    }
    (mapped as any)[k] = { ...value, layers: out };
  }
  return mapped;
}

/** `layers`' each by `lociApply` on `renderObject` (`Visual.setOverpaint` and the others: each cleared first) */
export function setTargetLayers(renderObject: GraphicsRenderObject, layers: TargetLayers, lociApply: Visual.LociApply) {
  if (layers.overpaint) Visual.setOverpaint(renderObject, layers.overpaint, lociApply, true);
  if (layers.transparency) Visual.setTransparency(renderObject, layers.transparency, lociApply, true);
  if (layers.emissive) Visual.setEmissive(renderObject, layers.emissive, lociApply, true);
  if (layers.substance) Visual.setSubstance(renderObject, layers.substance, lociApply, true);
  if (layers.clipping) Visual.setClipping(renderObject, layers.clipping, lociApply, true);
  if (layers.wiggle) Visual.setWiggle(renderObject, layers.wiggle, lociApply, true);
}

/**
 * A loci of the list on a target's instances `particleIndices` (an instance a particle of them): their runs applied,
 * every loci all of them
 */
export function createTargetLociApply(particleIndices: OrderedSet<number>): Visual.LociApply {
  const instanceCount = OrderedSet.size(particleIndices);
  return (loci, apply) => {
    // Fast path: whole visual selected
    if (isEveryLoci(loci)) {
      return apply(Interval.ofBounds(0, instanceCount));
    }
    if (!Particle.isLoci(loci) || Particle.isLociEmpty(loci)) return false;
    // Intersect loci indices with this visual's particle indices
    const intersection = OrderedSet.intersect(loci.indices, particleIndices);
    const intersectionSize = OrderedSet.size(intersection);
    if (intersectionSize === 0) return false;
    // Fast path: all instances in this target are covered
    if (intersectionSize === instanceCount) {
      return apply(Interval.ofBounds(0, instanceCount));
    }
    // Iterate intersection in order; indexOf gives monotonically increasing instance indices → merge runs
    let changed = false;
    let start = -1,
      end = -1;
    OrderedSet.forEach(intersection, (particleIdx) => {
      const instanceIdx = OrderedSet.indexOf(particleIndices, particleIdx);
      if (start === -1) {
        start = end = instanceIdx;
      } else if (instanceIdx === end + 1) {
        end = instanceIdx;
      } else {
        if (apply(Interval.ofBounds(start, end + 1))) changed = true;
        start = end = instanceIdx;
      }
    });
    if (start !== -1 && apply(Interval.ofBounds(start, end + 1))) changed = true;
    return changed;
  };
}

/**
 * Whether `target`'s instances in `b` are those of `a`: both views of one tile set at one frame
 * (`ParticleTileView`), its part of them the same tiles with the same counts. The same particles in the same order,
 * at other indices of the list where the tiles before them changed: what is built of them kept, their indices taken
 */
function isSameTilePart(a: ParticleList, b: ParticleList, target: number) {
  const va = ParticleTileView.of(a);
  const vb = ParticleTileView.of(b);
  return !!va && !!vb && va.set === vb.set && va.frame === vb.frame && ParticleTileView.isSamePart(va, vb, target);
}

export function createTargetVisual(targetId: number, materialId: number, graphics?: GraphicsContext): TargetVisual {
  const updateState = VisualUpdateState.create();

  let renderObject: GraphicsRenderObject | undefined;
  let geometry: Geometry | undefined;
  let geometryVersion = -1;
  let locationIt: LocationIterator;
  let positionIt: LocationIterator;

  let currentTargetData: unknown;
  let currentTargetKind: ParticleTarget['kind'] | undefined;
  let currentProps: ParticleTargetProps | undefined;
  let currentTheme: Theme | undefined;
  let currentParticles: ParticleList | undefined;
  let currentParticleIndices: OrderedSet<number> | undefined;
  let currentTarget: ParticleTarget | undefined;

  /** the particles' transforms into the render object in place (its matrix kept), and the bounds they give */
  function updateTransform(
    particles: ParticleList,
    particleIndices: OrderedSet<number>,
    target: ParticleTarget,
    props: ParticleTargetProps,
    scaleByRadius: boolean,
    posed = false,
  ) {
    const values = renderObject!.values;
    const gridVersion = values.instanceGrid.ref.version;
    createTargetParticleTransform(
      particles,
      particleIndices,
      geometry!.boundingSphere,
      props.cellSize,
      props.batchSize,
      scaleByRadius,
      values as unknown as TransformData,
      ParticleTarget.getCenter(target),
      isTargetCulledOnGpu(particles, props, graphics),
      isPosesGeometry(geometry!.kind),
      posed,
    );
    // from the transforms as drawn, as `Visual.setTransform` does: stale bounds cull a moved list where it was
    const boundingSphere = calculateTransformDataBoundingSphere(
      values.invariantBoundingSphere.ref.value,
      values as unknown as TransformData,
    );
    ValueCell.update(values.boundingSphere, boundingSphere);
    if ('lodLevels' in values && values.instanceGrid.ref.version !== gridVersion) {
      // instanceGrid (and its cellCount) was just rebuilt above; bump the `lodLevels`
      // version to force `renderable.cull` to resize/refresh its per-level multidraw
      // data (`mdbDataList`) instead of reusing stale, potentially undersized arrays -
      // otherwise cells beyond the old cellCount are silently dropped from culling,
      // which can leave the wrong (smaller) LOD level's spheres rendered.
      ValueCell.update(values.lodLevels, values.lodLevels.ref.value);
    }
  }

  /**
   * The instances moved as a whole by `move` (the list moved by it, `Particle.getMove`): the transforms and their grid
   * moved, not built again (`moveTransformData`), so the GPU culling keeps its split of them; the bounds moved. False
   * where the transforms can't be moved (a representation transform around them)
   */
  function moveTransform(move: Mat4) {
    const values = renderObject!.values;
    if (!moveTransformData(values as unknown as TransformData, move)) return false;
    const b = values.boundingSphere.ref.value;
    const moved = Sphere3D.create(Vec3.transformMat4(Vec3(), b.center, move), b.radius * Mat4.getMaxScaleOnAxis(move));
    ValueCell.update(values.boundingSphere, moved);
    return true;
  }

  /** Location iterator: groupCount = 1 (particle-granularity), instanceCount = particle count for this target.
   *  Uses the original (unfiltered) `ParticleList` so that themes receive the correct particle index. */
  function createLocIt(particles: ParticleList, indices: OrderedSet<number>): LocationIterator {
    const loc = Particle.Location(particles, 0);
    // the set's kind resolved once, not per instance (`OrderedSet.getAt`): a theme reads every particle's location
    let getLocation: (groupIndex: number, instanceIndex: number) => Particle.Location;
    if (OrderedSet.isInterval(indices)) {
      const start = Interval.start(indices);
      getLocation = (_groupIndex, instanceIndex) => {
        loc.index = start + instanceIndex;
        return loc;
      };
    } else {
      const array = indices as SortedArray<number>;
      getLocation = (_groupIndex, instanceIndex) => {
        loc.index = array[instanceIndex];
        return loc;
      };
    }
    return LocationIterator(1, OrderedSet.size(indices), 1, getLocation, false);
  }

  /** Create a renderObject from scratch for a geometry + the target's particles. */
  function createRO(
    geom: Geometry,
    particles: ParticleList,
    particleIndices: OrderedSet<number>,
    theme: Theme,
    props: ParticleTargetProps,
    target: ParticleTarget,
    scaleByRadius: boolean,
  ): GraphicsRenderObject {
    const { createValues, createRenderableState } = Geometry.getUtils(geom);
    const instanceCount = OrderedSet.size(particleIndices);
    const transform =
      instanceCount > 0
        ? createTargetParticleTransform(
            particles,
            particleIndices,
            geom.boundingSphere,
            props.cellSize,
            props.batchSize,
            scaleByRadius,
            undefined,
            ParticleTarget.getCenter(target),
            isTargetCulledOnGpu(particles, props, graphics),
            isPosesGeometry(geom.kind),
          )
        : createIdentityTransform();
    locationIt = createLocIt(particles, particleIndices);
    // Sizes of a structure target are per-element, unlike the per-particle colors and picking: an impostors' target's
    // made first, the values with one uniform size, their largest (the bounds' padding as the final one, no bounds
    // computed twice), not the representation's size a particle (31.5 M of a pack) made and overwritten
    const elementSizes =
      isElementTarget(target) && geom.kind === 'spheres' ? createElementSizes(target, theme, geom) : undefined;
    let valuesTheme = getTargetTheme(theme, particleIndices);
    if (elementSizes)
      valuesTheme = { ...valuesTheme, size: UniformSizeTheme({} as any, { value: getMaxSize(elementSizes) }) };
    const values = createValues(geom, transform, locationIt, valuesTheme, props as any);
    const state = createRenderableState(props as any);
    positionIt = Geometry.getUtils(geom).createPositionIterator(geom, values as any);
    // a color constant over the particles: filled per instance (the values' first one uniform, `getTargetTheme`)
    if (isTargetColorConstant(theme, particleIndices)) {
      createTargetColors(locationIt, positionIt, theme, particleIndices, values as any);
    }
    if (elementSizes && 'uSize' in values) {
      const v = values as unknown as { [k: string]: ValueCell<unknown> };
      for (const [k, cell] of Object.entries(elementSizes))
        if (v[k]) ValueCell.update(v[k], (cell as ValueCell<unknown>).ref.value);
      // the padding by the sizes as they are (the same largest: the bounds kept, `calculateTransformDataBoundingSphere`)
      Geometry.getUtils(geom).updateBoundingSphere(values as any, geom as any);
    } else if (isElementTarget(target) && 'uSize' in values) {
      createElementTargetSizes(target, theme, values as unknown as SizeValues);
      // the padding and bounds `createValues` took from the representation's (per-particle) size theme
      Geometry.getUtils(geom).updateBoundingSphere(values as any, geom as any);
    }
    setTargetSiteColor(target, geom, props, values);
    return createRenderObject(geom.kind, values, state, materialId);
  }

  async function createOrUpdate(
    ctx: VisualContext,
    representationTheme: Theme,
    inputProps: ParticleTargetProps,
    particles: ParticleList,
    particleIndices: OrderedSet<number>,
    target: ParticleTarget,
  ): Promise<void> {
    VisualUpdateState.reset(updateState);

    // Marking is per particle, and this visual's location iterator is per particle too
    // (groupCount 1). A structure target's geometry groups are its elements though, so
    // `groupInstance` granularity would index marker data by particle id + element id.
    const props = { ...inputProps, instanceGranularity: true };

    // Dots are sized uniformly and scaled by `sizeFactor`; the representation's per-particle
    // size theme would make every dot as large as the particle it belongs to.
    const theme = isVolumeDot(target, props)
      ? { color: representationTheme.color, size: getDotSizeTheme() }
      : representationTheme;

    const data = targetData(target);
    const targetKindChanged = currentTargetKind !== undefined && target.kind !== currentTargetKind;
    const targetChanged = data !== currentTargetData;
    // `particleIndices` is usually an `Interval`, which is packed into a plain number, so two
    // particle lists with the same target layout (e.g. two frames of a trajectory) compare `===`.
    // The instance data therefore has to be keyed on the `ParticleList` itself as well.
    // the list moved as a whole (a placement's release): its instances moved, not built again, where its particles are
    // the same (the groups computed again give equal index sets)
    const move =
      currentParticles && particles !== currentParticles ? Particle.getMove(currentParticles, particles) : undefined;
    const moved =
      !!move &&
      !!currentParticleIndices &&
      (particleIndices === currentParticleIndices || OrderedSet.areEqual(particleIndices, currentParticleIndices));
    // the same particles posed otherwise (`Particle.withPoses`: a trajectory's frame, a drag): their transforms built
    // again, what their identity gives kept (colors and sizes by particle, as a move keeps them)
    const posed =
      !moved &&
      !!currentParticles &&
      particles !== currentParticles &&
      Particle.isSamePoses(currentParticles, particles) &&
      !!currentParticleIndices &&
      (particleIndices === currentParticleIndices || OrderedSet.areEqual(particleIndices, currentParticleIndices));
    // a view of other tiles where this target's are the same (a tile of other targets loaded or let go): its instances
    // the same particles, only their indices in the list others
    const kept =
      !moved &&
      !posed &&
      !!currentParticles &&
      particles !== currentParticles &&
      !!currentParticleIndices &&
      OrderedSet.size(particleIndices) === OrderedSet.size(currentParticleIndices) &&
      isSameTilePart(currentParticles, particles, targetId);
    const indicesChanged = !moved && !posed && !kept && particleIndices !== currentParticleIndices;
    const particlesChanged = indicesChanged || (!moved && !posed && !kept && particles !== currentParticles);
    const scaleByRadius = targetScaleByRadius(target, props);

    const geometryPropsChanged =
      currentProps && !targetKindChanged ? targetGeometryPropsChanged(target, currentProps, props, graphics) : false;
    const colorThemeChanged = !currentTheme || !ColorTheme.areEqual(theme.color, currentTheme.color);
    const sizeThemeChanged = !currentTheme || !SizeTheme.areEqual(theme.size, currentTheme.size);
    const sizeFactorChanged =
      !!currentProps && !targetKindChanged && targetSizeFactorChanged(target, currentProps, props, graphics);

    if (
      !renderObject ||
      targetKindChanged ||
      targetChanged ||
      (currentProps && targetMustRecreate(target, currentProps, props))
    ) {
      updateState.createNew = true;
    } else {
      updateState.createGeometry =
        targetChanged ||
        geometryPropsChanged ||
        (sizeThemeChanged && targetSizeThemeAffectsGeometry(target, props, graphics));
      updateState.updateMatrix =
        particlesChanged || posed || targetChanged || targetScaleByRadius(target, currentProps!) !== scaleByRadius;
      // colors/sizes are keyed by instance (via `locationIt`/`positionIt`, built from `particleIndices`),
      // so a changed particle set requires rebuilding them even if the theme itself is unchanged; a moved one only
      // where they read positions
      updateState.updateColor =
        colorThemeChanged ||
        particlesChanged ||
        updateState.createGeometry ||
        ((moved || posed) && isPositional(theme.color.granularity));
      updateState.updateSize =
        sizeThemeChanged ||
        sizeFactorChanged ||
        updateState.createGeometry ||
        (particlesChanged && !isElementTarget(target)) ||
        ((moved || posed) && !isElementTarget(target) && isPositional(theme.size.granularity));
    }

    if (updateState.createNew) {
      updateState.createGeometry = true;
      updateState.updateColor = true;
      updateState.updateSize = true;
      updateState.updateMatrix = true;
    }

    let newGeom: Geometry | undefined;
    if (updateState.createGeometry) {
      const built = createTargetGeometry(ctx, target, theme, props, graphics, geometry);
      newGeom = isPromiseLike(built) ? await built : built;
    }

    if (updateState.createNew || newGeom) {
      const g = (newGeom ?? geometry)!;
      geometry = g;
      geometryVersion++;
      renderObject = createRO(g, particles, particleIndices, theme, props, target, scaleByRadius);
    } else if (renderObject) {
      if (indicesChanged) {
        // Markers are per-instance; a new particle set changes the instance count and the
        // meaning of each index. Stale markers make `getMarkersAverage` read out of bounds
        // (NaN `markerAverage`), which silently drops the visual from the marking pass.
        const instanceCount = OrderedSet.size(particleIndices);
        const type = renderObject.values.markerType.ref.value as 'instance' | 'groupInstance';
        createMarkers(instanceCount, type, renderObject.values as unknown as MarkerData);
      }

      if (moved && !updateState.updateMatrix && OrderedSet.size(particleIndices) > 0) {
        if (!moveTransform(move!)) updateTransform(particles, particleIndices, target, props, scaleByRadius);
        positionIt = Geometry.getUtils(geometry!).createPositionIterator(geometry!, renderObject.values as any);
      }
      if (updateState.updateMatrix && OrderedSet.size(particleIndices) > 0) {
        updateTransform(particles, particleIndices, target, props, scaleByRadius, posed);
        positionIt = Geometry.getUtils(geometry!).createPositionIterator(geometry!, renderObject.values as any);
      }

      // kept: the location iterator of the new indices (marking, a later color update), nothing built again
      if (kept) locationIt = createLocIt(particles, particleIndices);

      if (updateState.updateColor) {
        locationIt = createLocIt(particles, particleIndices);
        createTargetColors(locationIt, positionIt, theme, particleIndices, renderObject.values as any);
      }

      if (updateState.updateSize && 'uSize' in renderObject.values) {
        if (isElementTarget(target)) {
          createElementTargetSizes(target, theme, renderObject.values as SizeValues);
          Geometry.getUtils(geometry!).updateBoundingSphere(renderObject.values as any, geometry! as any);
        } else {
          locationIt = createLocIt(particles, particleIndices);
          createSizes(locationIt, positionIt, theme.size, renderObject.values as SizeValues);
        }
      }
      const geomUtils = Geometry.getUtils(geometry!);
      geomUtils.updateValues(renderObject.values as any, props as any);
      geomUtils.updateRenderableState(renderObject.state, props as any);
      // a color update may have made the colors uniform again
      if (!currentProps || siteColorPropsChanged(currentProps, props) || updateState.updateColor) {
        setTargetSiteColor(target, geometry!, props, renderObject.values);
      }
    }

    currentTargetData = data;
    currentTargetKind = target.kind;
    // the instances' order and markers are the old set's: an equal one moved keeps it
    if (moved) particleIndices = currentParticleIndices!;
    currentProps = { ...props };
    currentTheme = theme;
    currentParticles = particles;
    currentParticleIndices = particleIndices;
    currentTarget = target;
  }

  function setTransform(matrix: Mat4) {
    if (!renderObject || !currentParticles || !currentParticleIndices || !currentTarget || !currentProps) return;
    if (OrderedSet.size(currentParticleIndices) === 0) return;
    const m = renderObject.values.matrix;
    if (Mat4.areEqual(m.ref.value, matrix, 0)) return;
    ValueCell.update(m, Mat4.copy(m.ref.value, matrix));
    // rebuilt from the particles: the grid needs the sphere `createTargetParticleTransform` inflates by the radii
    updateTransform(
      currentParticles,
      currentParticleIndices,
      currentTarget,
      currentProps,
      targetScaleByRadius(currentTarget, currentProps),
    );
  }

  function mark(loci: ModelLoci, action: MarkerAction): boolean {
    if (!renderObject || !currentParticleIndices) return false;
    return Visual.mark(renderObject, loci, action, createTargetLociApply(currentParticleIndices));
  }

  function setLayers(layers: TargetLayers) {
    if (!renderObject || !currentParticleIndices) return;
    setTargetLayers(renderObject, layers, createTargetLociApply(currentParticleIndices));
  }

  function destroy() {
    renderObject = undefined;
    geometry = undefined;
    currentTargetData = undefined;
    currentTargetKind = undefined;
    currentProps = undefined;
    currentTheme = undefined;
    currentParticles = undefined;
    currentParticleIndices = undefined;
    currentTarget = undefined;
  }

  return {
    get renderObject() {
      return renderObject;
    },
    get geometryVersion() {
      return geometryVersion;
    },
    get particleIndices() {
      return currentParticleIndices;
    },
    createOrUpdate,
    setTransform,
    mark,
    setLayers,
    destroy,
  };
}

// ---- A member of a store-backed object -------------------------------------

/** A target's member of a store-backed object (`slots.ts`): its geometry and values, its instances the store's */
export interface TargetMember {
  readonly renderObject: GraphicsRenderObject;
  readonly geometry: Geometry;
}

/**
 * `target`'s geometry and values as a member of an object whose instances are kept elsewhere (a page store,
 * `MergedInstanceSource`): one instance of its own, placed at the identity (its transforms, poses and grid are not read),
 * its colors and markers of instance granularity (the source's layers; `theme`'s palette as it is), its sizes as a
 * visual's (`createTargetVisual`: an element target's by element). Its instance count, bounds and markers are set by the
 * store's owner
 */
export async function createTargetMember(
  ctx: VisualContext,
  target: ParticleTarget,
  theme: Theme,
  inputProps: ParticleTargetProps,
  graphics: GraphicsContext | undefined,
  materialId: number,
): Promise<TargetMember> {
  const props = { ...inputProps, instanceGranularity: true };
  const built = createTargetGeometry(ctx, target, theme, props, graphics);
  const geometry = isPromiseLike(built) ? await built : built;
  const { createValues, createRenderableState } = Geometry.getUtils(geometry);
  const locationIt = LocationIterator(1, 1, 1, () => NullLocation, false);
  // the layer's colors the source's: a placeholder of the granularity they are read by
  const color = {
    ...ColorTheme.Empty,
    granularity: 'instance',
    color: () => Color(0xffffff),
    palette: theme.color.palette,
  } as ColorTheme<any>;
  const elementSizes =
    isElementTarget(target) && geometry.kind === 'spheres' ? createElementSizes(target, theme, geometry) : undefined;
  const size = elementSizes ? UniformSizeTheme({} as any, { value: getMaxSize(elementSizes) }) : theme.size;
  const values = createValues(geometry, createIdentityTransform(), locationIt, { color, size }, props as any);
  if (elementSizes && 'uSize' in values) {
    const v = values as unknown as { [k: string]: ValueCell<unknown> };
    for (const [k, cell] of Object.entries(elementSizes))
      if (v[k]) ValueCell.update(v[k], (cell as ValueCell<unknown>).ref.value);
    Geometry.getUtils(geometry).updateBoundingSphere(values as any, geometry as any);
  } else if (isElementTarget(target) && 'uSize' in values) {
    createElementTargetSizes(target, theme, values as unknown as SizeValues);
    Geometry.getUtils(geometry).updateBoundingSphere(values as any, geometry as any);
  }
  setTargetSiteColor(target, geometry, props, values);
  // its grid pending: built from the source's poses where a CPU cull needs its cells (`MergedInstanceSource.getMember`)
  ValueCell.update(
    values.instanceGrid,
    createPendingInstanceGrid(props.cellSize, props.batchSize, values.invariantBoundingSphere.ref.value),
  );
  // a grid of cells for what culls by it (`hasInstanceCells`: a CPU cull by its cells and levels, not whole)
  ValueCell.updateIfChanged(values.instanceCellSize, props.cellSize);
  const renderObject = createRenderObject(geometry.kind, values, createRenderableState(props as any), materialId);
  return { renderObject, geometry };
}
