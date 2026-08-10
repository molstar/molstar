/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../../../mol-data/int';
import { CustomProperties } from '../../../../mol-model/custom-property';
import { Particle, ParticleList } from '../../../../mol-model/particles/particle-list';
import { ColorTheme } from '../../../../mol-theme/color';
import { SizeTheme } from '../../../../mol-theme/size';
import { FibersRepresentation, getFibersParams } from '../fibers';

function createTestParticles(): ParticleList {
    return {
        count: 3,
        keys: new Int32Array([0, 1, 2]),
        targets: new Int32Array(3),
        coordinates: new Float32Array([
            0, 0, 0,
            1, 0, 0,
            2, 0, 0,
        ]),
        fibers: {
            count: 1,
            offsets: new Int32Array([0, 2]),
            indices: new Int32Array([1, 2]),
        },
        getParticleLabel: index => `${index}`,
        sourceData: { kind: 'test', name: 'test', data: {} },
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };
}

describe('fibers representation', () => {
    it('includes only fiber particles in its all-loci', async () => {
        const particles = createTestParticles();
        const representation = FibersRepresentation({
            colorThemeRegistry: ColorTheme.createRegistry(),
            sizeThemeRegistry: SizeTheme.createRegistry(),
        }, getFibersParams);

        await representation.createOrUpdate({ visuals: ['lines'] }, particles).run();

        const [loci] = representation.getAllLoci();
        expect(Particle.isLoci(loci)).toBe(true);
        if (!Particle.isLoci(loci)) return;

        const indices: number[] = [];
        OrderedSet.forEach(loci.indices, index => indices.push(index));
        expect(indices).toEqual([1, 2]);

        representation.destroy();
    });
});