/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import { type Boundary, getFastBoundary } from '@molstar/core/math/geometry/boundary';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { fillSerial } from '@molstar/core/util/array';
import { packIntToRGBArray } from '@molstar/core/util/number-packing';
import { NullLocation } from '@molstar/model/model/location';
import type { ParticleSites } from '@molstar/model/model/particles/particle-list';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import { sphereVertexCount } from '@molstar/graphics/geo/primitive/sphere';
import type { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import {
  type SizeData,
  createTextureSize,
  createValueSize,
  sizeDataFactor,
} from '@molstar/graphics/geo/geometry/size-data';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Theme } from '@molstar/graphics/theme/theme';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import type { VisualContext } from '../../../visual.js';
import { computePositionBlobSurface, DefaultBlobDensityProps } from '../../../structure/visual/util/blob-surface.js';
import { getSizeScale, type StructureTargetProps, structureTargetUseImpostor } from './structure.js';
import type { BlobSurfaceTargetProps } from './blob-surface.js';

// A target of sites (`ParticleTarget` kind `sites`): drawn as a structure target is (`structure.ts`), from the sites'
// positions and radii, not a structure's elements; a site's index its group.

/** A site's radius as the size theme gives it: the physical radius by the theme's scale, or the uniform size */
function getSiteRadius(sites: ParticleSites, theme: Theme): (index: number) => number {
  if (theme.size.granularity === 'uniform') {
    const size = theme.size.size(NullLocation);
    return () => size;
  }
  const scale = getSizeScale(theme);
  const { radii } = sites;
  return (i) => radii[i] * scale;
}

/** The sites' positions as their coordinates' arrays (`PositionData`) */
function getSitesPosition(sites: ParticleSites) {
  const { positions, count } = sites;
  const x = new Float32Array(count);
  const y = new Float32Array(count);
  const z = new Float32Array(count);
  for (let i = 0, o = 0; i < count; ++i, o += 3) {
    x[i] = positions[o];
    y[i] = positions[o + 1];
    z[i] = positions[o + 2];
  }
  return { x, y, z, indices: OrderedSet.ofBounds(0, count), id: fillSerial(new Uint32Array(count)) };
}

const boundaries = new WeakMap<ParticleSites, Boundary>();
const boundaryHelperCoarse = new BoundaryHelper('14');
const boundaryHelperFine = new BoundaryHelper('98');
/**
 * The boundary of the sites, as a target structure's (`computeStructureBoundary`: over its units' fast boundaries, their
 * boxes' spheres, a chain a unit): the bounds the impostors' level of detail is chosen by (their size on screen) and the
 * scene's sphere merges, and the blob surface's grid
 */
function getSitesBoundary(sites: ParticleSites): Boundary {
  let boundary = boundaries.get(sites);
  if (!boundary) {
    const position = getSitesPosition(sites);
    const offsets = sites.chains?.offsets;
    const spheres: Sphere3D[] = [];
    if (offsets && offsets.length > 2) {
      for (let c = 0; c + 1 < offsets.length; ++c) {
        if (offsets[c + 1] === offsets[c]) continue;
        spheres.push(getFastBoundary({ ...position, indices: OrderedSet.ofBounds(offsets[c], offsets[c + 1]) }).sphere);
      }
    } else {
      spheres.push(getFastBoundary(position).sphere);
    }
    // the coarse helper above 10k extrema, as a structure's (8 a box's sphere)
    const helper = spheres.length * 8 > 10_000 ? boundaryHelperCoarse : boundaryHelperFine;
    helper.reset();
    for (const s of spheres) helper.includeSphere(s);
    helper.finishedIncludeStep();
    for (const s of spheres) helper.radiusSphere(s);
    boundary = { box: helper.getBox(), sphere: helper.getSphere() };
    boundaries.set(sites, boundary);
  }
  return boundary;
}

/** Build the geometry for a sites target: impostor spheres or a sphere mesh, as `createStructureTargetGeometry` */
export function createSitesTargetGeometry(
  sites: ParticleSites,
  theme: Theme,
  props: StructureTargetProps,
  graphics?: GraphicsContext,
  existing?: Spheres | Mesh,
): Spheres | Mesh {
  const { positions, count } = sites;
  const radius = getSiteRadius(sites, theme);
  let maxRadius = 0;
  for (let i = 0; i < count; ++i) maxRadius = Math.max(maxRadius, radius(i));
  const boundingSphere = Sphere3D.expand(
    Sphere3D(),
    getSitesBoundary(sites).sphere,
    maxRadius * props.sizeFactor + 0.05,
  );

  if (structureTargetUseImpostor(props, graphics)) {
    // radii are supplied as size data (`createSitesTargetSizes`), the geometry only holds centers
    const builder = SpheresBuilder.create(count, count / 2, existing?.kind === 'spheres' ? existing : undefined);
    for (let i = 0, o = 0; i < count; ++i, o += 3) builder.add(positions[o], positions[o + 1], positions[o + 2], i);
    const s = builder.getSpheres();
    if (count > 0) s.setBoundingSphere(boundingSphere);
    return s;
  }
  // meshes carry no size data: the radius baked into the vertices
  const { detail, sizeFactor } = props;
  const vertexCount = count * sphereVertexCount(detail);
  const state = MeshBuilder.createState(vertexCount, vertexCount / 2, existing?.kind === 'mesh' ? existing : undefined);
  const v = Vec3();
  for (let i = 0, o = 0; i < count; ++i, o += 3) {
    state.currentGroup = i;
    addSphere(state, Vec3.set(v, positions[o], positions[o + 1], positions[o + 2]), radius(i) * sizeFactor, detail);
  }
  const m = MeshBuilder.getMesh(state);
  if (count > 0) m.setBoundingSphere(boundingSphere);
  return m;
}

/** The sizes of a sites target a site each, as `createStructureTargetSizes`: its physical radius, or a uniform size */
export function createSitesTargetSizes(sites: ParticleSites, theme: Theme, sizeData?: SizeData): SizeData {
  if (theme.size.granularity === 'uniform') {
    return createValueSize(theme.size.size(NullLocation), sizeData);
  }
  const scale = getSizeScale(theme);
  const { radii, count } = sites;
  const sizes = createDataArray(Math.max(1, count), 4, Uint8Array, sizeData?.tSize.ref.value.array);
  for (let i = 0; i < count; ++i) packIntToRGBArray(radii[i] * scale * sizeDataFactor, sizes.array, i * 4);
  return createTextureSize(sizes, 'group', sizeData);
}

/** The blob surface of a sites target, as `createBlobSurfaceTargetGeometry` of a structure's */
export async function createSitesBlobSurfaceGeometry(
  ctx: VisualContext,
  sites: ParticleSites,
  theme: Theme,
  props: BlobSurfaceTargetProps,
  existing?: Mesh,
): Promise<Mesh> {
  const { blobSize, blobMethod, resolution, adjustResolution, blobShape, radiusOffset, smoothness } = props;
  // the positions as their coordinates' arrays for the density only
  const position = getSitesPosition(sites);
  const boundary = getSitesBoundary(sites);
  const { transform, field, idField, radiusFactor, maxRadius } = await computePositionBlobSurface(
    position,
    boundary,
    getSiteRadius(sites, theme),
    {
      ...DefaultBlobDensityProps,
      blobSize,
      blobMethod,
      resolution,
      adjustResolution,
      blobShape,
      radiusOffset,
      smoothness,
    },
  ).runInContext(ctx.runtime);

  const isoLevel = Math.exp(-smoothness) / radiusFactor;
  const surface = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }, existing).runAsChild(
    ctx.runtime,
  );

  Mesh.transform(surface, transform);
  if (ctx.graphics && !ctx.graphics.capabilities.flatInterpolation) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }
  const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
  surface.setBoundingSphere(Sphere3D.expand(Sphere3D(), boundary.sphere, maxRadius + extraRadius));
  return surface;
}
