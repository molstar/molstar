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

const TargetStructuresDescriptor = CustomPropertyDescriptor({ name: 'particle-target-structures' });

export function getParticleTransforms(data: ParticleList) {
    const particleCount = data.count;
    const transforms: Mat4[] = [];
    const { rotations } = data;

    for (let i = 0; i < particleCount; ++i) {
        const cOffset = i * 3;

        let transform: Mat4;
        if (rotations) {
            const qOffset = i * 4;
            const q = Quat.create(
                rotations[qOffset + 0],
                rotations[qOffset + 1],
                rotations[qOffset + 2],
                rotations[qOffset + 3],
            );
            transform = Mat4.fromQuat(Mat4(), q);
        } else {
            transform = Mat4.identity();
        }
        transform[12] = data.coordinates[cOffset + 0];
        transform[13] = data.coordinates[cOffset + 1];
        transform[14] = data.coordinates[cOffset + 2];
        transforms.push(transform);
    }

    return transforms;
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

    export function setTargetStructures(
        particles: ParticleList,
        map: ReadonlyMap<number, import('../structure').Structure>
    ): void {
        particles.customProperties.add(TargetStructuresDescriptor);
        particles._propertyData[TargetStructuresDescriptor.name] = map;
    }

    export function getTargetStructures(
        particles: ParticleList
    ): ReadonlyMap<number, import('../structure').Structure> | undefined {
        return particles._propertyData[TargetStructuresDescriptor.name];
    }
}
