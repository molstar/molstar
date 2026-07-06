/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../mol-data/int';
import { Mat4, Quat, Vec3 } from '../../mol-math/linear-algebra';
import { Sphere3D } from '../../mol-math/geometry';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { ModelFormat } from '../../mol-model-formats/format';
import { CustomProperties, CustomPropertyDescriptor } from '../custom-property';
import { Boundary } from '../../mol-math/geometry/boundary';
import { fillIdentityTransform } from '../../mol-geo/geometry/transform-data';

export interface ParticleList {
    readonly entryId?: string
    readonly label?: string

    readonly count: number

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
     * Optional mapping from each unique target ID in `targets` to the list of
     * canonical chain IDs (label_asym_id values) that make up the reference
     * structure for that target.  Used by `buildTargetStructuresFromMapping`
     * in `src/mol-model/particles/particle-structure-registry.ts` to automatically
     * split a parent structure into per-target sub-structures.
     */
    readonly targetMapping?: ReadonlyMap<number, ReadonlyArray<string>>

    /**
     * Optional mapping from each unique target ID in `targets` to a trajectory model index.
     * When present, each target's reference structure is the full structure built from that
     * trajectory model, rather than a chain-split sub-structure of a single parent structure.
     * Used by the petworld mmCIF variant where each molecule type is stored as a separate model
     * (`pdbx_PDB_model_num`) and chain IDs are reused across models. Takes precedence over
     * `targetMapping`.
     */
    readonly targetModels?: ReadonlyMap<number, number>

    /**
     * Optional per-particle compartment index (length = `count`). Each value identifies
     * which compartment this particle belongs to. A value of -1 means "no compartment".
     * The distinct non-negative values correspond to keys of `compartmentInfo`.
     */
    readonly compartments?: Int32Array

    /**
     * Optional mapping from each unique compartment index in `compartments` to the
     * compartment name/path string (e.g. `"root.mge.surface.proteins"`).
     */
    readonly compartmentInfo?: ReadonlyMap<number, string>

    /**
     * Optional per-particle entity index (length = `count`). Each value identifies
     * which entity (molecule type) this particle belongs to. A value of -1 means "no entity".
     * The distinct non-negative values correspond to keys of `entityInfo`.
     */
    readonly entities?: Int32Array

    /**
     * Optional mapping from each unique entity index in `entities` to the entity
     * name string (e.g. `"MG_191_192_NAP"`).
     */
    readonly entityInfo?: ReadonlyMap<number, string>

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
     * Simularium where an agent expands into a chain of particles. Currently informational
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
    readonly attributes?: ReadonlyMap<string, Float32Array>

    /** Metadata for each key in `attributes`. */
    readonly attributeInfo?: ReadonlyMap<string, {
        readonly label: string
        readonly min: number
        readonly max: number
    }>

    readonly getParticleLabel: (index: number) => string

    readonly sourceData: ModelFormat

    customProperties: CustomProperties
    _propertyData: { [name: string]: any }
}

/**
 * A reference object instanced at each particle of a given target id. Each particle has
 * exactly one target id (see `ParticleList.targets`); the distinct target ids map to these
 * targets via `Particle.setParticleTargets` / `Particle.getParticleTargets`.
 */
export type ParticleTarget =
    | { readonly kind: 'structure', readonly structure: import('../structure/structure').Structure }
    | { readonly kind: 'shape', readonly shape: import('../shape/shape').Shape }

const ParticleTransformsDescriptor = CustomPropertyDescriptor({ name: 'particle-transforms' });
const ParticleTransformsAsMat4Descriptor = CustomPropertyDescriptor({ name: 'particle-transforms-as-mat4' });

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

/**
 * Per-particle transforms as `Mat4` instances. Computed from `getParticleTransforms` and
 * cached on the `ParticleList`.
 *
 * Note: the returned matrices are shared/cached, do not mutate them in place; clone before
 * mutating (e.g. `Mat4.mul(Mat4(), transform, offset)`).
 */
export function getParticleTransformsAsMat4(data: ParticleList): Mat4[] {
    if (!data._propertyData[ParticleTransformsAsMat4Descriptor.name]) {
        const transformArray = getParticleTransforms(data);
        const particleCount = data.count;
        const transforms: Mat4[] = new Array(particleCount);
        for (let i = 0; i < particleCount; ++i) {
            transforms[i] = Mat4.fromArray(Mat4(), transformArray, i * 16);
        }

        data.customProperties.add(ParticleTransformsAsMat4Descriptor);
        data._propertyData[ParticleTransformsAsMat4Descriptor.name] = transforms;
    }
    return data._propertyData[ParticleTransformsAsMat4Descriptor.name];
}

export namespace Particle {
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

    const _boundaryHelper = new BoundaryHelper('98');
    const _tmpPos = Vec3();
    export function getBoundingSphere(loci: Loci, boundingSphere?: Sphere3D): Sphere3D {
        if (!boundingSphere) boundingSphere = Sphere3D();
        const { particles, indices } = loci;
        const { coordinates, radii } = particles;
        if (OrderedSet.isEmpty(indices)) {
            boundingSphere.center[0] = boundingSphere.center[1] = boundingSphere.center[2] = 0;
            boundingSphere.radius = 0;
            return boundingSphere;
        }
        _boundaryHelper.reset();
        if (radii) {
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                _boundaryHelper.includePositionRadius(_tmpPos, radii[v]);
            });
            _boundaryHelper.finishedIncludeStep();
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                _boundaryHelper.radiusPositionRadius(_tmpPos, radii[v]);
            });
        } else {
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                _boundaryHelper.includePosition(_tmpPos);
            });
            _boundaryHelper.finishedIncludeStep();
            OrderedSet.forEach(indices, v => {
                const i = v * 3;
                Vec3.set(_tmpPos, coordinates[i], coordinates[i + 1], coordinates[i + 2]);
                _boundaryHelper.radiusPosition(_tmpPos);
            });
        }
        const sphere = _boundaryHelper.getSphere();
        Sphere3D.copy(boundingSphere, sphere);
        return boundingSphere;
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

    const ParticleTargetsDescriptor = CustomPropertyDescriptor({ name: 'particle-targets' });
    export function setParticleTargets(
        particles: ParticleList,
        map: ReadonlyMap<number, ParticleTarget>
    ): void {
        particles.customProperties.add(ParticleTargetsDescriptor);
        particles._propertyData[ParticleTargetsDescriptor.name] = map;
    }
    export function getParticleTargets(
        particles: ParticleList
    ): ReadonlyMap<number, ParticleTarget> | undefined {
        return particles._propertyData[ParticleTargetsDescriptor.name];
    }

    export const BoundaryDescriptor: CustomPropertyDescriptor<Boundary> = CustomPropertyDescriptor({ name: 'particle-boundary' });
    export function setBoundary(particles: ParticleList, boundary: Boundary): void {
        particles.customProperties.add(BoundaryDescriptor);
        particles._propertyData[BoundaryDescriptor.name] = boundary;
    }
    const boundaryHelperCoarse = new BoundaryHelper('14');
    const boundaryHelperFine = new BoundaryHelper('98');
    function getBoundaryHelper(count: number) {
        return count > 10_000 ? boundaryHelperCoarse : boundaryHelperFine;
    }
    export function getBoundary(particles: ParticleList): Boundary {
        if (!particles._propertyData[BoundaryDescriptor.name]) {
            // Compute boundary from particle positions and radii, and store it in the particle list for later retrieval.
            // loop over positions and radii to compute the boundary
            const { count, coordinates, radii } = particles;
            const boundaryHelper = getBoundaryHelper(count);
            const _tmpPos = Vec3();
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
            const sphere = boundaryHelper.getSphere();
            particles._propertyData[BoundaryDescriptor.name] = { box: boundaryHelper.getBox(), sphere };
        };
        return particles._propertyData[BoundaryDescriptor.name];
    }
}
