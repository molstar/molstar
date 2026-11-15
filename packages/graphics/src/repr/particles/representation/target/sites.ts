/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { fillSerial } from '@molstar/core/util/array';
import { packIntToRGBArray } from '@molstar/core/util/number-packing';
import { NullLocation } from '@molstar/model/model/location';
import { ParticleSites } from '@molstar/model/model/particles/particle-list';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import { sphereVertexCount } from '@molstar/graphics/geo/primitive/sphere';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import {
  type SizeData,
  createTextureSize,
  createValueSize,
  createSiteTypeSize,
  sizeDataFactor,
} from '@molstar/graphics/geo/geometry/size-data';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import type { SiteColorValues } from '@molstar/graphics/gl/renderable/spheres';
import { ElementSymbolColors, elementNumberColor } from '@molstar/graphics/theme/color/element-symbol';
import { Color } from '@molstar/core/util/color';
import { type SiteColorProps, setSiteColorTable } from './site-color.js';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Theme } from '@molstar/graphics/theme/theme';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import type { VisualContext } from '../../../visual.js';
import { computePositionBlobSurface, DefaultBlobDensityProps } from '../../../structure/visual/util/blob-surface.js';
import { getSizeScale, type StructureTargetProps, structureTargetUseImpostor } from './structure.js';
import type { BlobSurfaceTargetProps } from './blob-surface.js';

// A target of sites (`ParticleTarget` kind `sites`): drawn as a structure target is (`structure.ts`), from the sites'
// positions and radii, not a structure's elements; a site's index its group. Read a site at a time
// (`ParticleSites.getPosition`, `getRadius`): a pack's quantized sites read in place, no decoded arrays made; the impostors
// drawn from them as they are (`Spheres.fromSites`, WebGL2 and WebGPU), no copy of the centers.

/** A site's radius as the size theme gives it: the physical radius by the theme's scale, or the uniform size */
function getSiteRadius(sites: ParticleSites, theme: Theme): (index: number) => number {
  if (theme.size.granularity === 'uniform') {
    const size = theme.size.size(NullLocation);
    return () => size;
  }
  const scale = getSizeScale(theme);
  return (i) => ParticleSites.getRadius(sites, i) * scale;
}

/**
 * The largest of the sites' radii as the size theme gives them (`getSiteRadius`): a uniform size as it is, a pack's
 * quantized sites by their types read in place (the largest type present, by the theme's scale), not a radius read a
 * site (20.3 M of the yeast cell)
 */
function getSitesMaxRadius(sites: ParticleSites, theme: Theme, radius: (index: number) => number): number {
  const { count } = sites;
  if (count === 0) return 0;
  if (theme.size.granularity === 'uniform') return radius(0);
  const q = sites.quantized;
  if (!q) {
    let max = 0;
    for (let i = 0; i < count; ++i) max = Math.max(max, radius(i));
    return max;
  }
  const { data, typeRadii } = q;
  let max = 0;
  for (let o = 3, ol = count * 4; o < ol; o += 4) {
    const r = typeRadii[data[o] & 0xffff];
    if (r > max) max = r;
  }
  return max * getSizeScale(theme);
}

/** The sites' positions as their coordinates' arrays (`PositionData`) */
function getSitesPosition(sites: ParticleSites) {
  const { count } = sites;
  const x = new Float32Array(count);
  const y = new Float32Array(count);
  const z = new Float32Array(count);
  const p = Vec3();
  for (let i = 0; i < count; ++i) {
    ParticleSites.getPosition(sites, i, p);
    x[i] = p[0];
    y[i] = p[1];
    z[i] = p[2];
  }
  return { x, y, z, indices: OrderedSet.ofBounds(0, count), id: fillSerial(new Uint32Array(count)) };
}

/**
 * The box of sites `[start, end)`, as `Box3D.computeBounding` of their float32 positions (`getFastBoundary`'s): a
 * pack's quantized ones by their integers' extremes, decoded once (float32's rounding is monotonic: the same box), no
 * decoded array made (the yeast cell's 20.3 M sites)
 */
function getSitesBox(sites: ParticleSites, start: number, end: number, chain: number): Box3D {
  const q = sites.quantized;
  const min = [Infinity, Infinity, Infinity];
  const max = [-Infinity, -Infinity, -Infinity];
  if (q) {
    const { data, origin, scale } = q;
    // the chain's extremes as a pack has them, else from its sites
    const e = q.extremes ?? ParticleSites.getQuantizedExtremes(data, start, end, extremesScratch);
    const o = q.extremes ? chain * 6 : 0;
    for (let a = 0; a < 3; ++a) {
      min[a] = e[o + a];
      max[a] = e[o + 3 + a];
    }
    // the decoded float32s of the extremes (`getSitesPosition`'s values); a scale not above 0 flips nothing: `ofQuantized`'s
    // scale is a step, positive
    for (let a = 0; a < 3; ++a) {
      min[a] = Math.fround(origin[a] + min[a] * scale);
      max[a] = Math.fround(origin[a] + max[a] * scale);
    }
  } else {
    const p = sites.positions;
    for (let o = start * 3, ol = end * 3; o < ol; o += 3) {
      for (let a = 0; a < 3; ++a) {
        const v = p[o + a];
        if (v < min[a]) min[a] = v;
        if (v > max[a]) max[a] = v;
      }
    }
  }
  return Box3D.create(Vec3.create(min[0], min[1], min[2]), Vec3.create(max[0], max[1], max[2]));
}

const extremesScratch = new Int16Array(6);
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
    const offsets = sites.chains?.offsets;
    const spheres: Sphere3D[] = [];
    if (offsets && offsets.length > 2) {
      for (let c = 0; c + 1 < offsets.length; ++c) {
        if (offsets[c + 1] === offsets[c]) continue;
        spheres.push(Sphere3D.fromBox3D(Sphere3D(), getSitesBox(sites, offsets[c], offsets[c + 1], c)));
      }
    } else {
      spheres.push(Sphere3D.fromBox3D(Sphere3D(), getSitesBox(sites, 0, sites.count, 0)));
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
  const { count } = sites;
  const radius = getSiteRadius(sites, theme);
  const p = Vec3();
  const maxRadius = getSitesMaxRadius(sites, theme, radius);
  const boundingSphere = Sphere3D.expand(
    Sphere3D(),
    getSitesBoundary(sites).sphere,
    maxRadius * props.sizeFactor + 0.05,
  );

  if (structureTargetUseImpostor(props, graphics)) {
    // radii are supplied as size data (`createSitesTargetSizes`), the geometry only holds centers: a pack's quantized
    // ones as they are (8 bytes a site, a view of the file), else a copy as floats
    const q = sites.quantized;
    if (q && Spheres.isSitesSupported(graphics)) {
      const reuse = existing?.kind === 'spheres' ? existing : undefined;
      // coarse sites (`ParticleSites.stride`) drawn as their stride's level of detail
      const s = Spheres.fromSites(
        { data: q.data, origin: q.origin, scale: q.scale, stride: sites.stride },
        count,
        reuse,
      );
      if (count > 0) s.setBoundingSphere(boundingSphere);
      return s;
    }
    const builder = SpheresBuilder.create(count, count / 2, existing?.kind === 'spheres' ? existing : undefined);
    for (let i = 0; i < count; ++i) {
      ParticleSites.getPosition(sites, i, p);
      builder.add(p[0], p[1], p[2], i);
    }
    const s = builder.getSpheres();
    if (count > 0) s.setBoundingSphere(boundingSphere);
    return s;
  }
  // meshes carry no size data: the radius baked into the vertices
  const { detail, sizeFactor } = props;
  const vertexCount = count * sphereVertexCount(detail);
  const state = MeshBuilder.createState(vertexCount, vertexCount / 2, existing?.kind === 'mesh' ? existing : undefined);
  for (let i = 0; i < count; ++i) {
    state.currentGroup = i;
    addSphere(state, ParticleSites.getPosition(sites, i, p), radius(i) * sizeFactor, detail);
  }
  const m = MeshBuilder.getMesh(state);
  if (count > 0) m.setBoundingSphere(boundingSphere);
  return m;
}

/**
 * The sizes of a sites target a site each, as `createStructureTargetSizes`: its physical radius, or a uniform size;
 * `bySiteType` (its spheres drawn from a pack's quantized sites, `dSites`) a table of the radii by type, read by the
 * site's type (`siteType` sizes: a few entries, not one a site)
 */
export function createSitesTargetSizes(
  sites: ParticleSites,
  theme: Theme,
  sizeData?: SizeData,
  bySiteType = false,
): SizeData {
  if (theme.size.granularity === 'uniform') {
    return createValueSize(theme.size.size(NullLocation), sizeData);
  }
  const scale = getSizeScale(theme);
  const q = sites.quantized;
  if (bySiteType && q) {
    return createSiteTypeSize(
      Array.from(q.typeRadii, (r) => r * scale),
      sizeData,
    );
  }
  const { count } = sites;
  const sizes = createDataArray(Math.max(1, count), 4, Uint8Array, sizeData?.tSize.ref.value.array);
  for (let i = 0; i < count; ++i) {
    packIntToRGBArray(ParticleSites.getRadius(sites, i) * scale * sizeDataFactor, sizes.array, i * 4);
  }
  return createTextureSize(sizes, 'group', sizeData);
}

const CarbonNumber = 6;

/**
 * The second color of a sites target's spheres (`SiteColorValues`): off where asked or where the sites' elements are not
 * known; else by the sites' types where the spheres are drawn from them (`bySiteType`: a texel a type) or by group (a
 * texel a site); a carbon's alpha 0 where it keeps the particle's color
 */
export function setSitesTargetSiteColor(
  sites: ParticleSites,
  props: SiteColorProps,
  values: SiteColorValues,
  bySiteType: boolean,
) {
  const p = props.siteColor;
  if (!p || p.name === 'off' || !ParticleSites.hasElements(sites)) {
    ValueCell.updateIfChanged(values.dSiteColor, 'off');
    return;
  }
  const keepCarbon = p.params.carbonColor === 'particle';
  const typeElements = bySiteType ? sites.quantized?.typeElements : undefined;
  const count = typeElements ? typeElements.length : sites.count;
  const table = createDataArray(Math.max(1, count), 4, Uint8Array, values.tSiteColor.ref.value.array);
  const { array } = table;
  for (let i = 0; i < count; ++i) {
    const element = typeElements ? typeElements[i] : ParticleSites.getElement(sites, i);
    Color.toArray(elementNumberColor(ElementSymbolColors, element), array, i * 4);
    array[i * 4 + 3] = keepCarbon && element === CarbonNumber ? 0 : 255;
  }
  setSiteColorTable(values, table, typeElements ? 'siteType' : 'group', p.params);
}

/**
 * The defines a sites target's impostors take from its data (the precompile's signatures, which see none:
 * `getGeometrySignatures` of an empty geometry): drawn from the quantized sites where supported (`dSites`,
 * `createSitesTargetGeometry`), their second color by type there, else by group, off without elements
 * (`setSitesTargetSiteColor`)
 */
export function getSitesSignatureDefines(
  sites: ParticleSites,
  props: SiteColorProps,
  graphics: GraphicsContext | undefined,
): { readonly dSites: readonly boolean[]; readonly dSiteColor: readonly string[] } {
  const quantized = !!sites.quantized && Spheres.isSitesSupported(graphics);
  const second = props.siteColor?.name === 'element' && ParticleSites.hasElements(sites);
  return { dSites: [quantized], dSiteColor: [second ? (quantized ? 'siteType' : 'group') : 'off'] };
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
