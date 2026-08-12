/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../mol-data/int';
import { Column } from '../../mol-data/db';
import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { Sphere3D } from '../../mol-math/geometry';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { ModelFormat } from '../../mol-model-formats/format';
import { CustomProperties, CustomPropertyDescriptor } from '../custom-property';
import { Boundary } from '../../mol-math/geometry/boundary';
import { fillIdentityTransform } from '../../mol-geo/geometry/transform-data';
import { Cell } from '../../mol-math/geometry/spacegroup/cell';

/**
 * A single named per-particle scalar attribute (e.g. `'cc'`, `'class'`, `'score'`).
 *
 * `column` provides the per-particle values (indexed by particle position, length = `count`).
 * Formats should reuse an existing `Column` (e.g. via `Column.view`) or provide a lambda-based
 * accessor instead of copying data whenever possible.
 */
export interface ParticleAttribute {
    readonly column: Column<any>
    readonly label: string
    readonly min: number
    readonly max: number
}

export interface ParticleCompartmentInfo {
    /** Compartment name/path (e.g. `"root.mge.surface.proteins"`). */
    readonly name: string
}

export interface ParticleEntityInfo {
    /** Entity name (e.g. `"MG_191_192_NAP"`). */
    readonly name: string
    /** Optional functional annotation of the entity (e.g. `"Lipid"`). */
    readonly function?: string
}

export interface ParticleList {
    readonly entryId?: string
    readonly label?: string

    readonly count: number

    /**
     * Optional periodic/bounding cell for the particle system (e.g. the simulation box
     * size of a Simularium trajectory).
     */
    readonly cell?: Cell

    /** Unique keys for each particle for mapping to source data. */
    readonly keys: Int32Array

    /**
     * Per-particle target index (length = `count`). Each value identifies which
     * target structure (or later volume) this particle belongs to.  Use 0 for
     * single-target data.  The distinct values in this array correspond to the
     * keys of `targetMapping` when present.
     */
    readonly targets: Int32Array

    /**
     * Optional mapping from each unique target ID in `targets` to the reference object
     * instanced at every particle with that target ID. The objects are owned by the
     * particle list itself and are not part of the plugin state tree.
     */
    readonly targetMapping?: ReadonlyMap<number, ParticleTarget>

    /**
     * Optional per-particle compartment index (length = `count`). Each value identifies
     * which compartment this particle belongs to. A value of -1 means "no compartment".
     * The distinct non-negative values correspond to keys of `compartmentInfo`.
     */
    readonly compartments?: Int32Array

    /**
     * Optional mapping from each unique compartment index in `compartments` to the
     * compartment info.
     */
    readonly compartmentInfo?: ReadonlyMap<number, ParticleCompartmentInfo>

    /**
     * Optional per-particle entity index (length = `count`). Each value identifies
     * which entity (molecule type) this particle belongs to. A value of -1 means "no entity".
     * The distinct non-negative values correspond to keys of `entityInfo`.
     */
    readonly entities?: Int32Array

    /**
     * Optional mapping from each unique entity index in `entities` to the entity info.
     */
    readonly entityInfo?: ReadonlyMap<number, ParticleEntityInfo>

    /** Particle positions in angstrom, packed as `[x0, y0, z0, x1, y1, z1, ...]`. */
    readonly coordinates: Float32Array
    /** Optional per-particle orientations as unit quaternions, packed as `[x0, y0, z0, w0, ...]`. */
    readonly rotations?: Float32Array
    /** Optional per-particle bounding sphere radii in angstrom (length = `count`). */
    readonly radii?: Float32Array

    /**
     * Optional polyline (fiber) connectivity over particles, stored in compressed-sparse-row
     * form. Fiber `f` (for `0 <= f < count`) is the ordered polyline through the particles
     * `indices[offsets[f]]` .. `indices[offsets[f + 1]) - 1]`. Used by formats such as
     * Simularium where an agent expands into a chain of particles.
     */
    readonly fibers?: {
        readonly count: number
        readonly offsets: Int32Array
        readonly indices: Int32Array
    }

    /**
     * Named per-particle scalar attributes, indexed by particle position (length = count).
     * Keys are short identifiers (e.g. 'cc', 'class', 'score').
     */
    readonly attributes?: ReadonlyMap<string, ParticleAttribute>

    readonly getParticleLabel: (index: number) => string

    readonly sourceData: ModelFormat

    customProperties: CustomProperties
    _propertyData: { [name: string]: any }
}

/**
 * A reference object instanced at each particle of a given target id. Each particle has
 * exactly one target id (see `ParticleList.targets`); the distinct target ids map to these
 * targets via `ParticleList.targetMapping`.
 */
export type ParticleTarget =
    | { readonly kind: 'structure', readonly structure: import('../structure/structure').Structure }
    | { readonly kind: 'shape', readonly shape: import('../shape/shape').Shape }
    | { readonly kind: 'volume', readonly volume: import('../volume/volume').Volume }

export namespace ParticleTarget {
    export function data(target: ParticleTarget) {
        switch (target.kind) {
            case 'structure': return target.structure;
            case 'shape': return target.shape;
            case 'volume': return target.volume;
        }
    }
}

const ParticleTransformsDescriptor = CustomPropertyDescriptor({ name: 'particle-transforms' });

/**
 * Per-particle transforms as a flat array, 16 consecutive floats (column-major mat4) per
 * particle. Computed once and cached on the `ParticleList`.
 */
export function getParticleTransforms(data: ParticleList): Float32Array {
    if (!data._propertyData[ParticleTransformsDescriptor.name]) {
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
            const m = Mat4.identity();
            const q = Quat();
            for (let i = 0; i < particleCount; ++i) {
                const cOffset = i * 3;
                const qOffset = i * 4;
                Quat.set(q,
                    rotations[qOffset + 0],
                    rotations[qOffset + 1],
                    rotations[qOffset + 2],
                    rotations[qOffset + 3],
                );
                Mat4.fromQuat(m, q);
                m[12] = coordinates[cOffset + 0];
                m[13] = coordinates[cOffset + 1];
                m[14] = coordinates[cOffset + 2];
                for (let j = 0; j < 16; j++) {
                    transformArray[i * 16 + j] = m[j];
                }
                // transformArray.set(m, i * 16);
            }
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
    }
    return data._propertyData[ParticleTransformsDescriptor.name];
}

const ParticleTargetGroupsDescriptor = CustomPropertyDescriptor({ name: 'particle-target-groups' });

/** The particles of a `ParticleList` grouped by their target id. */
export interface ParticleTargetGroups {
    /** The distinct target ids, ascending. */
    readonly targetIds: Int32Array
    /** The particles of group `g`, as indices into the unfiltered `ParticleList`. */
    readonly sets: ReadonlyArray<OrderedSet<number>>
    readonly groupOfTarget: ReadonlyMap<number, number>
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

    let sorted = true;
    let maxId = 0;
    let minId = 0;
    let runCount = 0;
    for (let i = 0; i < count; ++i) {
        const t = targets[i];
        if (i === 0) {
            minId = maxId = t;
            runCount = 1;
        } else {
            const p = targets[i - 1];
            if (t < p) sorted = false;
            if (t !== p) runCount += 1;
            if (t < minId) minId = t;
            if (t > maxId) maxId = t;
        }
    }

    const targetIds = new Int32Array(runCount);
    const sets: OrderedSet<number>[] = new Array(runCount);
    const groupOfTarget = new Map<number, number>();

    if (count === 0) {
        return { targetIds, sets, groupOfTarget };
    }

    if (sorted) {
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
    const span = maxId - minId + 1;
    const useHistogram = span <= 4 * count + 1024;
    const histogram = useHistogram ? new Int32Array(span) : undefined;
    if (histogram) {
        for (let i = 0; i < count; ++i) histogram[targets[i] - minId] += 1;
        let g = 0;
        for (let s = 0; s < span; ++s) {
            if (histogram[s] !== 0) {
                targetIds[g] = s + minId;
                groupOfTarget.set(s + minId, g);
                g += 1;
            }
        }
    } else {
        const counts = new Map<number, number>();
        for (let i = 0; i < count; ++i) {
            const t = targets[i];
            counts.set(t, (counts.get(t) ?? 0) + 1);
        }
        const ids = Array.from(counts.keys()).sort((a, b) => a - b);
        for (let g = 0; g < ids.length; ++g) {
            targetIds[g] = ids[g];
            groupOfTarget.set(ids[g], g);
        }
    }

    const groupCount = targetIds.length;
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
        readonly kind: 'particle-location'
        particles: ParticleList
        /** Particle index in the list. */
        index: number
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
        readonly kind: 'particle-loci'
        readonly particles: ParticleList
        readonly indices: OrderedSet<number>
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
        OrderedSet.forEach(loci.indices, v => { if (v < count) filtered.push(v); });
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
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                boundaryHelper.includePositionRadius(_tmpPos, radii[v]);
            });
            boundaryHelper.finishedIncludeStep();
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                boundaryHelper.radiusPositionRadius(_tmpPos, radii[v]);
            });
        } else {
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                boundaryHelper.includePosition(_tmpPos);
            });
            boundaryHelper.finishedIncludeStep();
            OrderedSet.forEach(indices, v => {
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
        return `${size} Particles`;
    }

    /** Return a copy of `particles` with the given per-target reference objects attached. */
    export function withTargets(particles: ParticleList, targetMapping: ReadonlyMap<number, ParticleTarget>): ParticleList {
        return { ...particles, targetMapping };
    }

    export const BoundaryDescriptor: CustomPropertyDescriptor<Boundary> = CustomPropertyDescriptor({ name: 'particle-boundary' });
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
        };
        return particles._propertyData[BoundaryDescriptor.name];
    }
}
