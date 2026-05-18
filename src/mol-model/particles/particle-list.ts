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

export interface ParticleList {
    readonly entryId?: string
    readonly label?: string

    readonly count: number

    /** Unique keys for each particle for mapping to source data. */
    readonly keys: Int32Array

    /** Particle positions in angstrom, packed as `[x0, y0, z0, x1, y1, z1, ...]`. */
    readonly coordinates: Float32Array
    /** Optional per-particle orientations as unit quaternions, packed as `[x0, y0, z0, w0, ...]`. */
    readonly rotations?: Float32Array

    readonly getParticleLabel: (index: number) => string

    // TODO: add common data fields, anything format-specific should be accessed via `sourceData`

    readonly sourceData: ModelFormat
}

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
        const { coordinates } = particles;
        if (OrderedSet.isEmpty(indices)) {
            boundingSphere.center[0] = boundingSphere.center[1] = boundingSphere.center[2] = 0;
            boundingSphere.radius = 0;
            return boundingSphere;
        }
        _boundaryHelper.reset();
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
}
