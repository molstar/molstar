/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import type { Column } from '@molstar/core/data/db';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import { GridLookup3D } from '@molstar/core/math/geometry/lookup3d/grid';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { CentroidHelper } from '@molstar/core/math/geometry/centroid-helper';
import { Grid } from '../volume/grid.js';
import type { ModelFormat } from '@molstar/model/formats/format';
import { CustomProperties, CustomPropertyDescriptor } from '../custom-property.js';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import { fillIdentityTransform } from '@molstar/core/math/linear-algebra/3d/transform-array';
import type { Cell } from '@molstar/core/math/geometry/spacegroup/cell';

/**
 * A single named per-particle scalar attribute (e.g. `'cc'`, `'class'`, `'score'`).
 *
 * `column` provides the per-particle values (indexed by particle position, length = `count`).
 * Formats should reuse an existing `Column` (e.g. via `Column.view`) or provide a lambda-based
 * accessor instead of copying data whenever possible.
 */
export interface ParticleAttribute {
  readonly column: Column<any>;
  readonly label: string;
  readonly min: number;
  readonly max: number;
}

export interface ParticleCompartmentInfo {
  /** Compartment name/path (e.g. `"root.mge.surface.proteins"`). */
  readonly name: string;
}

export interface ParticleEntityInfo {
  /** Entity name (e.g. `"MG_191_192_NAP"`). */
  readonly name: string;
  /** Optional functional annotation of the entity (e.g. `"Lipid"`). */
  readonly function?: string;
}

/** Lightweight metadata for a single target id, always available (unlike `targetMapping`). */
export interface ParticleTargetInfo {
  /** Optional display name for this target (e.g. chain id, tomogram id, agent type name). */
  readonly name?: string;
  /** Index into `entityInfo`; states/enforces that every particle of this target shares this entity. */
  readonly entity?: number;
}

export interface ParticleList {
  readonly entryId?: string;
  readonly label?: string;

  readonly count: number;

  /**
   * Optional periodic/bounding cell for the particle system (e.g. the simulation box
   * size of a Simularium trajectory).
   */
  readonly cell?: Cell;

  /** Unique keys for each particle for mapping to source data. */
  readonly keys: Int32Array;

  /**
   * Per-particle target index (length = `count`). Each value identifies which
   * target structure/volume/shape this particle belongs to. Use 0 for
   * single-target data. The distinct values in this array correspond to the
   * keys of `targetMapping` when present.
   */
  readonly targets: Int32Array;

  /**
   * Metadata for each unique target ID in `targets`. Always present (unlike `targetMapping`),
   * so target ids can be enumerated cheaply without scanning `targets`.
   */
  readonly targetInfo: ReadonlyMap<number, ParticleTargetInfo>;

  /**
   * Optional mapping from each unique target ID in `targets` to the reference object
   * instanced at every particle with that target ID. The objects are owned by the
   * particle list itself and are not part of the plugin state tree.
   */
  readonly targetMapping?: ReadonlyMap<number, ParticleTarget>;

  /**
   * Optional per-particle compartment index (length = `count`). Each value identifies
   * which compartment this particle belongs to. A value of -1 means "no compartment".
   * The distinct non-negative values correspond to keys of `compartmentInfo`.
   */
  readonly compartments?: Int32Array;

  /**
   * Optional mapping from each unique compartment index in `compartments` to the
   * compartment info.
   */
  readonly compartmentInfo?: ReadonlyMap<number, ParticleCompartmentInfo>;

  /**
   * Optional per-particle entity index (length = `count`). Each value identifies
   * which entity (molecule type) this particle belongs to. A value of -1 means "no entity".
   * The distinct non-negative values correspond to keys of `entityInfo`.
   */
  readonly entities?: Int32Array;

  /**
   * Optional mapping from each unique entity index in `entities` to the entity info.
   */
  readonly entityInfo?: ReadonlyMap<number, ParticleEntityInfo>;

  /** Particle positions in angstrom, packed as `[x0, y0, z0, x1, y1, z1, ...]`. */
  readonly coordinates: Float32Array;
  /** Optional per-particle orientations as unit quaternions, packed as `[x0, y0, z0, w0, ...]`. */
  readonly rotations?: Float32Array;
  /** Optional per-particle bounding sphere radii in angstrom (length = `count`). */
  readonly radii?: Float32Array;

  /**
   * Optional polyline (fiber) connectivity over particles, stored in compressed-sparse-row
   * form. Fiber `f` (for `0 <= f < count`) is the ordered polyline through the particles
   * `indices[offsets[f]]` .. `indices[offsets[f + 1]) - 1]`. Used by formats such as
   * Simularium where an agent expands into a chain of particles.
   */
  readonly fibers?: {
    readonly count: number;
    readonly offsets: Int32Array;
    readonly indices: Int32Array;
  };

  /**
   * Named per-particle scalar attributes, indexed by particle position (length = count).
   * Keys are short identifiers (e.g. 'cc', 'class', 'score').
   */
  readonly attributes?: ReadonlyMap<string, ParticleAttribute>;

  /** A particle in words, its index from 1 first: "Particle 165 · tomo 1 · class 2 · cc 0.4213" */
  readonly getParticleLabel: (index: number) => string;

  readonly sourceData: ModelFormat;

  /**
   * The list this one is a copy of moved as a whole (`Particle.transform`) and the move from it: the first of a line
   * of moves, the moves composed, so two lists of a line give the move between them (`Particle.getMove`) and what was
   * built for one can be moved for the other
   */
  readonly movedFrom?: { readonly particles: ParticleList; readonly transform: Mat4 };

  customProperties: CustomProperties;
  _propertyData: { [name: string]: any };
}

/**
 * A reference object instanced at each particle of a given target id. Each particle has
 * exactly one target id (see `ParticleList.targets`); the distinct target ids map to these
 * targets via `ParticleList.targetMapping`.
 */
export type ParticleTarget =
  | {
      readonly kind: 'structure';
      readonly structure: import('../structure/structure.js').Structure;
      /**
       * The point of the structure placed at a particle's coordinates: the centroid of its atoms (the default,
       * `ParticleTarget.getCenter`), given when the format has it already (the centroid it placed its particles by)
       */
      readonly center?: Vec3;
    }
  | {
      /** typed spheres (atoms, coarse beads) without a structure: what a target's geometry needs, read from its file */
      readonly kind: 'sites';
      readonly sites: ParticleSites;
      /** as a structure target's: the point placed at a particle's coordinates, the sites' centroid by default */
      readonly center?: Vec3;
    }
  | { readonly kind: 'shape'; readonly shape: import('../shape/shape.js').Shape }
  | { readonly kind: 'volume'; readonly volume: import('../volume/volume.js').Volume };

/**
 * The sites of a target (`ParticleTarget` kind `sites`): positions and radii, a site's index its group in the target's
 * geometry. A few bytes a site, not a model and a structure of them (millions of atoms in a mesoscale model's targets)
 */
export interface ParticleSites {
  readonly count: number;
  /** x, y, z a site */
  readonly positions: Float32Array;
  /** the physical radius a site (an atom's van der Waals radius, a bead's own) */
  readonly radii: Float32Array;
  readonly label?: string;
  /** the sites' chains, their sites consecutive (as a structure's units): runs, a few bytes a chain, not an index a site */
  readonly chains?: ParticleSiteChains;
}

/** Chain `c` of `ParticleSites`: its sites `[offsets[c], offsets[c + 1])`, its label (`label_asym_id`) */
export interface ParticleSiteChains {
  /** a chain's first site, the site count last (chain count + 1 values) */
  readonly offsets: Int32Array;
  readonly labels: readonly string[];
}

export namespace ParticleTarget {
  export function data(target: ParticleTarget) {
    switch (target.kind) {
      case 'structure':
        return target.structure;
      case 'sites':
        return target.sites;
      case 'shape':
        return target.shape;
      case 'volume':
        return target.volume;
    }
  }

  const centroids = new WeakMap<object, Vec3>();
  const position = Vec3();

  /**
   * The point of the target placed at a particle's coordinates: a structure's centroid (`center` when given, cached
   * per structure otherwise); a volume's box center (voxel N/2, the subtomogram formats' convention), not its
   * isosurface's, which moves with the iso value; none for shapes, placed by their geometry's bounding sphere center
   */
  export function getCenter(target: ParticleTarget): Vec3 | undefined {
    if (target.kind === 'shape') return undefined;
    if (target.kind === 'volume') return Grid.getBoundingSphere(target.volume.grid).center;
    if (target.center) return target.center;
    if (target.kind === 'sites') {
      const { sites } = target;
      let center = centroids.get(sites);
      if (!center) {
        const helper = new CentroidHelper();
        const p = sites.positions;
        for (let i = 0, o = 0; i < sites.count; ++i, o += 3) {
          helper.includeStep(Vec3.set(position, p[o], p[o + 1], p[o + 2]));
        }
        helper.finishedIncludeStep();
        center = Vec3.clone(helper.center);
        centroids.set(sites, center);
      }
      return center;
    }
    const { structure } = target;
    let center = centroids.get(structure);
    if (!center) {
      const helper = new CentroidHelper();
      for (const unit of structure.units) {
        const { elements, conformation } = unit;
        for (let i = 0, il = elements.length; i < il; ++i) {
          helper.includeStep(conformation.position(elements[i], position));
        }
      }
      helper.finishedIncludeStep();
      center = Vec3.clone(helper.center);
      centroids.set(structure, center);
    }
    return center;
  }
}

const ParticleSpacingDescriptor = CustomPropertyDescriptor({ name: 'particle-spacing-radius' });

/** How many particles the spacing is measured over (evenly spread): the median settles well before */
const SpacingSamples = 2000;

/**
 * A radius for particles that carry none (picks, oriented points): half the median distance from a particle to its
 * nearest neighbor, so neighbors about touch, as cryo-ET picks of one kind do. Measured on up to 2,000 particles,
 * cached on the list; 1 for fewer than two particles or all at one point.
 */
export function getParticleSpacingRadius(data: ParticleList): number {
  const cached = data._propertyData[ParticleSpacingDescriptor.name] as number | undefined;
  if (cached !== undefined) return cached;
  const { count, coordinates } = data;
  let radius = 1;
  if (count > 1) {
    const x = new Float32Array(count);
    const y = new Float32Array(count);
    const z = new Float32Array(count);
    for (let i = 0; i < count; i++) {
      x[i] = coordinates[3 * i];
      y[i] = coordinates[3 * i + 1];
      z[i] = coordinates[3 * i + 2];
    }
    const lookup = GridLookup3D({ x, y, z, indices: OrderedSet.ofBounds(0, count) }, Particle.getBoundary(data));
    const step = Math.max(1, Math.floor(count / SpacingSamples));
    const nearest: number[] = [];
    for (let i = 0; i < count; i += step) {
      const r = lookup.nearest(x[i], y[i], z[i], 2);
      let best = Infinity;
      for (let k = 0; k < r.count; k++) {
        if (r.indices[k] !== i && r.squaredDistances[k] > 0) best = Math.min(best, r.squaredDistances[k]);
      }
      if (best < Infinity) nearest.push(Math.sqrt(best));
    }
    if (nearest.length) {
      nearest.sort((a, b) => a - b);
      radius = nearest[nearest.length >> 1] / 2;
    }
  }
  data.customProperties.add(ParticleSpacingDescriptor);
  data._propertyData[ParticleSpacingDescriptor.name] = radius;
  return radius;
}

/**
 * Particle `i`'s transform (column-major mat4, its rotation then its position) into `out` at `offset`, all 16 values:
 * as `getParticleTransforms` holds it, without the list's copy of all (a target's own transforms, millions of
 * particles)
 */
export function writeParticleTransform(data: ParticleList, i: number, out: Float32Array, offset: number) {
  const { rotations, coordinates } = data;
  const o = offset;
  const c = i * 3;
  if (rotations) {
    // quat → mat4 (see Mat4.fromQuat); the identity quaternion gives the identity
    const q = i * 4;
    const x = rotations[q],
      y = rotations[q + 1],
      z = rotations[q + 2],
      w = rotations[q + 3];
    const x2 = x + x,
      y2 = y + y,
      z2 = z + z;
    const xx = x * x2,
      yx = y * x2,
      yy = y * y2;
    const zx = z * x2,
      zy = z * y2,
      zz = z * z2;
    const wx = w * x2,
      wy = w * y2,
      wz = w * z2;
    out[o] = 1 - yy - zz;
    out[o + 1] = yx + wz;
    out[o + 2] = zx - wy;
    out[o + 4] = yx - wz;
    out[o + 5] = 1 - xx - zz;
    out[o + 6] = zy + wx;
    out[o + 8] = zx + wy;
    out[o + 9] = zy - wx;
    out[o + 10] = 1 - xx - yy;
  } else {
    out[o] = 1;
    out[o + 1] = 0;
    out[o + 2] = 0;
    out[o + 4] = 0;
    out[o + 5] = 1;
    out[o + 6] = 0;
    out[o + 8] = 0;
    out[o + 9] = 0;
    out[o + 10] = 1;
  }
  out[o + 3] = 0;
  out[o + 7] = 0;
  out[o + 11] = 0;
  out[o + 12] = coordinates[c];
  out[o + 13] = coordinates[c + 1];
  out[o + 14] = coordinates[c + 2];
  out[o + 15] = 1;
}

const ParticleTransformsDescriptor = CustomPropertyDescriptor({ name: 'particle-transforms' });

/**
 * Per-particle transforms as a flat array, 16 consecutive floats (column-major mat4) per
 * particle. Computed once and cached on the `ParticleList`.
 */
export function getParticleTransforms(data: ParticleList): Float32Array {
  let cached = data._propertyData[ParticleTransformsDescriptor.name] as Float32Array | undefined;
  if (!cached) {
    const particleCount = data.count;
    const transformArray = new Float32Array(particleCount * 16);
    const { rotations, coordinates } = data;

    let hasRotations = false;
    if (rotations) {
      for (let i = 0; i < particleCount; ++i) {
        const o = i * 4;
        if (rotations[o + 0] !== 0 || rotations[o + 1] !== 0 || rotations[o + 2] !== 0 || rotations[o + 3] !== 1) {
          hasRotations = true;
          break;
        }
      }
    }

    if (rotations && hasRotations) {
      for (let i = 0; i < particleCount; ++i) writeParticleTransform(data, i, transformArray, i * 16);
    } else {
      fillIdentityTransform(transformArray, particleCount);
      for (let i = 0; i < particleCount; ++i) {
        const cOffset = i * 3;
        transformArray[i * 16 + 12] = coordinates[cOffset + 0];
        transformArray[i * 16 + 13] = coordinates[cOffset + 1];
        transformArray[i * 16 + 14] = coordinates[cOffset + 2];
      }
    }

    data.customProperties.add(ParticleTransformsDescriptor);
    data._propertyData[ParticleTransformsDescriptor.name] = transformArray;
    cached = transformArray;
  }
  return cached;
}

const ParticleTargetGroupsDescriptor = CustomPropertyDescriptor({ name: 'particle-target-groups' });

/** The particles of a `ParticleList` grouped by their target id. */
export interface ParticleTargetGroups {
  /** The distinct target ids, ascending. */
  readonly targetIds: Int32Array;
  /** The particles of group `g`, as indices into the unfiltered `ParticleList`. */
  readonly sets: ReadonlyArray<OrderedSet<number>>;
  readonly groupOfTarget: ReadonlyMap<number, number>;
}

/**
 * Group the particles by target id. Computed once and cached on the `ParticleList`.
 *
 * When `targets` is already non-decreasing - the usual case, since formats emit particles
 * grouped by target - each group is a contiguous range and is returned as an `Interval`,
 * which needs no per-group storage and makes `OrderedSet.getAt`/`indexOf` arithmetic rather
 * than a binary search. Otherwise the groups are built with a counting sort.
 */
export function getParticleTargetGroups(data: ParticleList): ParticleTargetGroups {
  if (!data._propertyData[ParticleTargetGroupsDescriptor.name]) {
    data.customProperties.add(ParticleTargetGroupsDescriptor);
    data._propertyData[ParticleTargetGroupsDescriptor.name] = computeParticleTargetGroups(data);
  }
  return data._propertyData[ParticleTargetGroupsDescriptor.name];
}

function computeParticleTargetGroups(data: ParticleList): ParticleTargetGroups {
  const { targets, count } = data;

  const groupOfTarget = new Map<number, number>();

  if (count === 0) {
    return { targetIds: new Int32Array(0), sets: [], groupOfTarget };
  }

  let sorted = true;
  let minId = targets[0];
  let maxId = targets[0];
  let runCount = 1;
  for (let i = 1; i < count; ++i) {
    const t = targets[i];
    const p = targets[i - 1];
    if (t < p) sorted = false;
    if (t !== p) runCount += 1;
    if (t < minId) minId = t;
    if (t > maxId) maxId = t;
  }

  if (sorted) {
    // runs of equal ids are the distinct ids and are already ascending
    const targetIds = new Int32Array(runCount);
    const sets: OrderedSet<number>[] = new Array(runCount);
    let g = 0;
    let start = 0;
    for (let i = 1; i <= count; ++i) {
      if (i === count || targets[i] !== targets[start]) {
        targetIds[g] = targets[start];
        sets[g] = OrderedSet.ofBounds(start, i);
        groupOfTarget.set(targets[start], g);
        g += 1;
        start = i;
      }
    }
    return { targetIds, sets, groupOfTarget };
  }

  // Counting sort keyed on the target id, offset so that negative ids are handled too.
  let targetIds: Int32Array;
  const span = maxId - minId + 1;
  const useHistogram = span <= 4 * count + 1024;
  if (useHistogram) {
    const histogram = new Int32Array(span);
    for (let i = 0; i < count; ++i) histogram[targets[i] - minId] += 1;
    let distinctCount = 0;
    for (let s = 0; s < span; ++s) {
      if (histogram[s] !== 0) distinctCount += 1;
    }
    targetIds = new Int32Array(distinctCount);
    let g = 0;
    for (let s = 0; s < span; ++s) {
      if (histogram[s] !== 0) {
        targetIds[g] = s + minId;
        groupOfTarget.set(s + minId, g);
        g += 1;
      }
    }
  } else {
    const distinct = new Set<number>();
    for (let i = 0; i < count; ++i) distinct.add(targets[i]);
    const ids = Array.from(distinct).sort((a, b) => a - b);
    targetIds = new Int32Array(ids.length);
    for (let g = 0; g < ids.length; ++g) {
      targetIds[g] = ids[g];
      groupOfTarget.set(ids[g], g);
    }
  }

  const groupCount = targetIds.length;
  const sets: OrderedSet<number>[] = new Array(groupCount);
  const offsets = new Int32Array(groupCount + 1);
  for (let i = 0; i < count; ++i) offsets[groupOfTarget.get(targets[i])! + 1] += 1;
  for (let g = 0; g < groupCount; ++g) offsets[g + 1] += offsets[g];

  const indices = new Int32Array(count);
  const cursor = Int32Array.from(offsets.subarray(0, groupCount));
  for (let i = 0; i < count; ++i) {
    const g = groupOfTarget.get(targets[i])!;
    indices[cursor[g]] = i;
    cursor[g] += 1;
  }
  for (let g = 0; g < groupCount; ++g) {
    sets[g] = OrderedSet.ofSortedArray(indices.subarray(offsets[g], offsets[g + 1]));
  }

  return { targetIds, sets, groupOfTarget };
}

const FiberParticleMaskDescriptor = CustomPropertyDescriptor({ name: 'particle-fiber-mask' });

/**
 * Per-particle mask (length = `count`) marking fiber particles, i.e. particles referenced
 * by `fibers.indices`. Returns `undefined` if the particle list has no `fibers` data.
 * Computed once and cached on the `ParticleList`.
 */
export function getFiberParticleMask(data: ParticleList): Uint8Array | undefined {
  const { fibers } = data;
  if (!fibers) return undefined;
  if (!data._propertyData[FiberParticleMaskDescriptor.name]) {
    const mask = new Uint8Array(data.count);
    const { offsets, indices } = fibers;
    for (let f = 0; f < fibers.count; ++f) {
      const start = offsets[f];
      const end = offsets[f + 1];
      for (let i = start + 0; i < end; ++i) {
        mask[indices[i]] = 1;
      }
    }
    data.customProperties.add(FiberParticleMaskDescriptor);
    data._propertyData[FiberParticleMaskDescriptor.name] = mask;
  }
  return data._propertyData[FiberParticleMaskDescriptor.name];
}

export namespace Particle {
  export function is(x: any): x is ParticleList {
    return (
      typeof x?.count === 'number' &&
      x?.keys instanceof Int32Array &&
      x?.targets instanceof Int32Array &&
      x?.coordinates instanceof Float32Array &&
      !!x?.sourceData &&
      !!x?.customProperties &&
      !!x?._propertyData
    );
  }

  /** A single particle within a `ParticleList`. */
  export interface Location {
    readonly kind: 'particle-location';
    particles: ParticleList;
    /** Particle index in the list. */
    index: number;
  }

  export function Location(particles?: ParticleList, index = 0): Location {
    return { kind: 'particle-location', particles: particles!, index };
  }
  export function isLocation(x: any): x is Location {
    return !!x && x.kind === 'particle-location';
  }
  /** Write the particle's position into `out`. */
  export function position(out: Vec3, location: Location): Vec3 {
    const i = location.index * 3;
    const { coordinates } = location.particles;
    return Vec3.set(out, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
  }

  /** A loci over one or more particles in a `ParticleList`. */
  export interface Loci {
    readonly kind: 'particle-loci';
    readonly particles: ParticleList;
    readonly indices: OrderedSet<number>;
  }
  export function Loci(particles: ParticleList, indices: OrderedSet<number>): Loci {
    return { kind: 'particle-loci', particles, indices };
  }
  export function isLoci(x: any): x is Loci {
    return !!x && x.kind === 'particle-loci';
  }
  export function areLociEqual(a: Loci, b: Loci) {
    return a.particles === b.particles && OrderedSet.areEqual(a.indices, b.indices);
  }
  export function isLociEmpty(loci: Loci) {
    return OrderedSet.isEmpty(loci.indices);
  }
  export function lociSize(loci: Loci) {
    return OrderedSet.size(loci.indices);
  }
  /** Remap a loci to a new `ParticleList`; indices outside the new range are dropped. */
  export function remapLoci(loci: Loci, particles: ParticleList): Loci {
    if (loci.particles === particles) return loci;
    const { count } = particles;
    if (count === 0) return Loci(particles, OrderedSet.Empty);
    const filtered: number[] = [];
    OrderedSet.forEach(loci.indices, (v) => {
      if (v < count) filtered.push(v);
    });
    return Loci(particles, OrderedSet.ofSortedArray(filtered));
  }

  const boundaryHelperCoarse = new BoundaryHelper('14');
  const boundaryHelperFine = new BoundaryHelper('98');
  function getBoundaryHelper(count: number) {
    return count > 10_000 ? boundaryHelperCoarse : boundaryHelperFine;
  }
  const _tmpPos = Vec3();

  export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D): Sphere3D {
    const { particles, indices } = loci;
    const { coordinates, radii } = particles;
    if (isLociEmpty(loci)) {
      return boundingSphere ? Sphere3D.setZero(boundingSphere) : Sphere3D();
    } else if (lociSize(loci) === particles.count) {
      const sphere = getBoundary(particles).sphere;
      return boundingSphere ? Sphere3D.copy(sphere, boundingSphere) : sphere;
    }
    const boundaryHelper = getBoundaryHelper(OrderedSet.size(indices));
    boundaryHelper.reset();
    if (radii) {
      OrderedSet.forEach(indices, (v) => {
        const i = v * 3;
        Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
        boundaryHelper.includePositionRadius(_tmpPos, radii[v]);
      });
      boundaryHelper.finishedIncludeStep();
      OrderedSet.forEach(indices, (v) => {
        const i = v * 3;
        Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
        boundaryHelper.radiusPositionRadius(_tmpPos, radii[v]);
      });
    } else {
      OrderedSet.forEach(indices, (v) => {
        const i = v * 3;
        Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
        boundaryHelper.includePosition(_tmpPos);
      });
      boundaryHelper.finishedIncludeStep();
      OrderedSet.forEach(indices, (v) => {
        const i = v * 3;
        Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
        boundaryHelper.radiusPosition(_tmpPos);
      });
    }
    return boundaryHelper.getSphere(boundingSphere);
  }

  export function getLabel(loci: Loci): string {
    const size = OrderedSet.size(loci.indices);
    if (size === 0) return 'None';
    if (size === 1) {
      const index = OrderedSet.start(loci.indices);
      return loci.particles.getParticleLabel(index);
    }
    return `${size.toLocaleString('en')} particles`;
  }

  /**
   * A copy of `particles` moved by `transform`: the coordinates transformed, the orientations turned by its rotation
   * (given where it turns and the list has none); the radii, the cell and the targets as they were. What is cached on
   * the list by position (the boundary, the transforms) is computed again, what is not kept (`keptByMove`). The copy
   * records the move (`movedFrom`). For rigid moves: the coordinates take any affine `transform`, but an orientation is
   * a unit quaternion, turned by the rotation part alone (`Mat4.extractRotation`), so a mirroring transform (negative
   * determinant) or a shear leaves the orientations turned, not mirrored or sheared, and a scale leaves the radii.
   */
  export function transform(particles: ParticleList, transform: Mat4): ParticleList {
    const { count, coordinates: src } = particles;
    const coordinates = new Float32Array(src.length);
    // `Vec3.transformMat4` and `Quat.multiply` inlined (millions of particles); the transform affine, w 1
    const [m0, m1, m2, , m4, m5, m6, , m8, m9, m10, , m12, m13, m14, m15] = transform;
    const w = m15 || 1;
    for (let o = 0, ol = count * 3; o < ol; o += 3) {
      const x = src[o],
        y = src[o + 1],
        z = src[o + 2];
      coordinates[o] = (m0 * x + m4 * y + m8 * z + m12) / w;
      coordinates[o + 1] = (m1 * x + m5 * y + m9 * z + m13) / w;
      coordinates[o + 2] = (m2 * x + m6 * y + m10 * z + m14) / w;
    }
    const turn = Quat.normalize(Quat(), Quat.fromMat4(Quat(), Mat4.extractRotation(Mat4(), transform)));
    let rotations = particles.rotations;
    // q and −q are the same turn
    if (rotations || 1 - Math.abs(turn[3]) > 1e-9) {
      const turned = new Float32Array(count * 4);
      const [ax, ay, az, aw] = turn;
      for (let o = 0, ol = count * 4; o < ol; o += 4) {
        // the identity where the list has none
        const bx = rotations ? rotations[o] : 0,
          by = rotations ? rotations[o + 1] : 0,
          bz = rotations ? rotations[o + 2] : 0,
          bw = rotations ? rotations[o + 3] : 1;
        turned[o] = ax * bw + aw * bx + ay * bz - az * by;
        turned[o + 1] = ay * bw + aw * by + az * bx - ax * bz;
        turned[o + 2] = az * bw + aw * bz + ax * by - ay * bx;
        turned[o + 3] = aw * bw - ax * bx - ay * by - az * bz;
      }
      rotations = turned;
    }
    const from = particles.movedFrom;
    const customProperties = new CustomProperties();
    const propertyData = Object.create(null);
    const rigid = Math.abs(Mat4.getMaxScaleOnAxis(transform) - 1) < 1e-6;
    for (const descriptor of keptByMove) {
      const value = particles._propertyData[descriptor.name];
      if (value === undefined || (descriptor === ParticleSpacingDescriptor && !rigid)) continue;
      customProperties.add(descriptor);
      propertyData[descriptor.name] = value;
    }
    return {
      ...particles,
      coordinates,
      rotations,
      movedFrom: {
        particles: from?.particles ?? particles,
        transform: from ? Mat4.mul(Mat4(), transform, from.transform) : Mat4.copy(Mat4(), transform),
      },
      customProperties,
      _propertyData: propertyData,
    };
  }

  /**
   * What `transform` keeps of the list's cached data: what does not depend on positions (the target groups, the fiber
   * mask: the arrays they are of are shared), the spacing of a move without scale (distances kept); the release of a
   * move handle moves a list of millions, the groups alone cost 60 ms again
   */
  const keptByMove = [ParticleTargetGroupsDescriptor, FiberParticleMaskDescriptor, ParticleSpacingDescriptor];

  /**
   * The move from `from` to `to` where both are of one line of moves (`movedFrom`, either the first of it): `to` is
   * `from` moved by it, its particles the same; undefined otherwise
   */
  export function getMove(from: ParticleList, to: ParticleList): Mat4 | undefined {
    if ((from.movedFrom?.particles ?? from) !== (to.movedFrom?.particles ?? to)) return undefined;
    const t = to.movedFrom?.transform ?? Mat4.identity();
    const f = from.movedFrom?.transform;
    return f ? Mat4.mul(Mat4(), t, Mat4.invert(Mat4(), f)) : Mat4.copy(Mat4(), t);
  }

  /**
   * Return a copy of `particles` with the given per-target reference objects attached; of the same line of moves
   * (`movedFrom`: the coordinates the same), so a list moved under it gives its move
   */
  export function withTargets(
    particles: ParticleList,
    targetMapping: ReadonlyMap<number, ParticleTarget>,
  ): ParticleList {
    return { ...particles, targetMapping, movedFrom: particles.movedFrom ?? { particles, transform: Mat4.identity() } };
  }

  export const BoundaryDescriptor: CustomPropertyDescriptor<Boundary> = CustomPropertyDescriptor({
    name: 'particle-boundary',
  });
  export function setBoundary(particles: ParticleList, boundary: Boundary): void {
    particles.customProperties.add(BoundaryDescriptor);
    particles._propertyData[BoundaryDescriptor.name] = boundary;
  }
  export function getBoundary(particles: ParticleList): Boundary {
    if (!particles._propertyData[BoundaryDescriptor.name]) {
      // Compute boundary from particle positions and radii, and store it in the particle list for later retrieval.
      // loop over positions and radii to compute the boundary
      const { count, coordinates, radii } = particles;
      const boundaryHelper = getBoundaryHelper(count);
      boundaryHelper.reset();
      if (radii) {
        for (let i = 0; i < count; i++) {
          const cOffset = i * 3;
          Vec3.set(_tmpPos, coordinates[cOffset], coordinates[cOffset + 1], coordinates[cOffset + 2]);
          boundaryHelper.includePositionRadius(_tmpPos, radii[i]);
        }
        boundaryHelper.finishedIncludeStep();
        for (let i = 0; i < count; i++) {
          const cOffset = i * 3;
          Vec3.set(_tmpPos, coordinates[cOffset], coordinates[cOffset + 1], coordinates[cOffset + 2]);
          boundaryHelper.radiusPositionRadius(_tmpPos, radii[i]);
        }
      } else {
        for (let i = 0; i < count; i++) {
          const cOffset = i * 3;
          Vec3.set(_tmpPos, coordinates[cOffset], coordinates[cOffset + 1], coordinates[cOffset + 2]);
          boundaryHelper.includePosition(_tmpPos);
        }
        boundaryHelper.finishedIncludeStep();
        for (let i = 0; i < count; i++) {
          const cOffset = i * 3;
          Vec3.set(_tmpPos, coordinates[cOffset], coordinates[cOffset + 1], coordinates[cOffset + 2]);
          boundaryHelper.radiusPosition(_tmpPos);
        }
      }
      setBoundary(particles, { box: boundaryHelper.getBox(), sphere: boundaryHelper.getSphere() });
    }
    return particles._propertyData[BoundaryDescriptor.name];
  }
}
