/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet, SortedArray } from '../../../../mol-data/int';
import { CustomProperties } from '../../../../mol-model/custom-property';
import { ParticleList } from '../../../../mol-model/particles/particle-list';
import { Sphere3D } from '../../../../mol-math/geometry';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { createTargetParticleTransform } from '../target/common';

function createTestParticles(): ParticleList {
    return {
        count: 3,
        keys: new Int32Array([0, 1, 2]),
        targets: new Int32Array(3),
        targetInfo: new Map([[0, {}]]),
        coordinates: new Float32Array([
            1, 2, 3,
            4, 5, 6,
            7, 8, 9,
        ]),
        radii: new Float32Array([2, 3, 4]),
        getParticleLabel: index => `${index}`,
        sourceData: { kind: 'test', name: 'test', data: {} },
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

const sphereAtOrigin = Sphere3D.create(Vec3.create(0, 0, 0), 1);
const sphereOffCenter = Sphere3D.create(Vec3.create(10, 20, 30), 1);

describe('createTargetParticleTransform', () => {
    it('places each particle (contiguous indices, origin-centered geometry)', () => {
        const particles = createTestParticles();
        const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereAtOrigin, 100, 1000, false);
        const t = data.transform.ref.value;

        expect(data.instanceCount.ref.value).toBe(3);
        for (let i = 0; i < 3; ++i) {
            // identity basis
            expect(t[i * 16 + 0]).toBe(1);
            expect(t[i * 16 + 5]).toBe(1);
            expect(t[i * 16 + 10]).toBe(1);
            expect(t[i * 16 + 15]).toBe(1);
            // translation = particle position
            expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0]);
            expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1]);
            expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2]);
        }
    });

    it('subtracts the geometry center from the translation', () => {
        const particles = createTestParticles();
        const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, false);
        const t = data.transform.ref.value;

        for (let i = 0; i < 3; ++i) {
            expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0] - 10);
            expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1] - 20);
            expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2] - 30);
        }
    });

    it('gathers non-contiguous indices', () => {
        const particles = createTestParticles();
        const data = createTargetParticleTransform(particles, OrderedSet.ofSortedArray([0, 2]), sphereOffCenter, 100, 1000, false);
        const t = data.transform.ref.value;

        expect(data.instanceCount.ref.value).toBe(2);
        expect(t[12]).toBe(1 - 10);
        expect(t[16 + 12]).toBe(7 - 10);
    });

    it('scales by the particle radius', () => {
        const particles = createTestParticles();
        const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, true);
        const t = data.transform.ref.value;

        for (let i = 0; i < 3; ++i) {
            const s = particles.radii![i];
            expect(t[i * 16 + 0]).toBe(s);
            expect(t[i * 16 + 5]).toBe(s);
            expect(t[i * 16 + 10]).toBe(s);
            // translation accounts for the scaled basis
            expect(t[i * 16 + 12]).toBe(particles.coordinates[i * 3 + 0] - s * 10);
            expect(t[i * 16 + 13]).toBe(particles.coordinates[i * 3 + 1] - s * 20);
            expect(t[i * 16 + 14]).toBe(particles.coordinates[i * 3 + 2] - s * 30);
        }
    });

    it('matches between the bulk and gather paths', () => {
        const particles = createTestParticles();
        const bulk = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, false);
        // a raw SortedArray is not normalized to an Interval, forcing the per-instance gather path
        const gather = createTargetParticleTransform(particles, SortedArray.ofSortedArray([0, 1, 2]), sphereOffCenter, 100, 1000, false);

        expect(Array.from(gather.transform.ref.value)).toEqual(Array.from(bulk.transform.ref.value));
    });

    it('reuses the output array on in-place updates', () => {
        const particles = createTestParticles();
        const data = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereOffCenter, 100, 1000, false);
        const before = data.transform.ref.value;

        const updated = createTargetParticleTransform(particles, OrderedSet.ofBounds(0, 3), sphereAtOrigin, 100, 1000, false, data);
        expect(updated).toBe(data);
        expect(updated.transform.ref.value).toBe(before);
        expect(updated.transform.ref.value[12]).toBe(1);
    });
});
