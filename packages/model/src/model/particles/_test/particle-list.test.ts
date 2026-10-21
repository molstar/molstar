/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '../../custom-property.js';
import { OrderedSet } from '@molstar/core/data/int';
import { lociLabel } from '../../label.js';
import { getParticleSpacingRadius, getParticleTransforms, Particle, ParticleList } from '../particle-list.js';

function createParticles(coordinates: Float32Array, rotations?: Float32Array): ParticleList {
  const count = coordinates.length / 3;
  return {
    count,
    keys: Int32Array.from({ length: count }, (_, i) => i),
    targets: new Int32Array(count),
    targetInfo: new Map([[0, {}]]),
    coordinates,
    rotations,
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

describe('getParticleTransforms', () => {
  it('uses identity bases without rotations', () => {
    const particles = createParticles(new Float32Array([1, 2, 3, 4, 5, 6]));
    const t = getParticleTransforms(particles);

    const expected = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
    for (let j = 0; j < 16; ++j) expect(t[j]).toBeCloseTo(expected[j], 6);
    expect(t[16 + 12]).toBe(4);
    expect(t[16 + 13]).toBe(5);
    expect(t[16 + 14]).toBe(6);
  });

  it('matches Mat4.fromQuat for rotated particles', () => {
    // deterministic LCG so the test is reproducible
    let seed = 3;
    const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;

    const count = 25;
    const coordinates = new Float32Array(count * 3);
    const rotations = new Float32Array(count * 4);
    const q = Quat();
    for (let i = 0; i < count; ++i) {
      coordinates[i * 3] = rand() * 100;
      coordinates[i * 3 + 1] = rand() * 100;
      coordinates[i * 3 + 2] = rand() * 100;
      Quat.set(q, rand() - 0.5, rand() - 0.5, rand() - 0.5, rand() - 0.5);
      Quat.normalize(q, q);
      rotations.set(q, i * 4);
    }

    const particles = createParticles(coordinates, rotations);
    const t = getParticleTransforms(particles);

    const m = Mat4();
    for (let i = 0; i < count; ++i) {
      Quat.set(q, rotations[i * 4], rotations[i * 4 + 1], rotations[i * 4 + 2], rotations[i * 4 + 3]);
      Mat4.fromQuat(m, q);
      m[12] = coordinates[i * 3];
      m[13] = coordinates[i * 3 + 1];
      m[14] = coordinates[i * 3 + 2];
      for (let j = 0; j < 16; ++j) {
        expect(t[i * 16 + j]).toBeCloseTo(m[j], 5);
      }
    }
  });

  it('caches the result on the particle list', () => {
    const particles = createParticles(new Float32Array([1, 2, 3]));
    expect(getParticleTransforms(particles)).toBe(getParticleTransforms(particles));
  });
});

describe('getParticleSpacingRadius', () => {
  it('is half the distance to the nearest particle, for particles on a lattice', () => {
    const coords: number[] = [];
    for (let i = 0; i < 6; i++)
      for (let j = 0; j < 6; j++) for (let k = 0; k < 6; k++) coords.push(i * 10, j * 10, k * 10);
    const list = createParticles(new Float32Array(coords));
    expect(getParticleSpacingRadius(list)).toBeCloseTo(5, 4);
    // cached
    expect(getParticleSpacingRadius(list)).toBe(getParticleSpacingRadius(list));
  });

  it('is 1 for a single particle or particles all at one point', () => {
    expect(getParticleSpacingRadius(createParticles(new Float32Array([1, 2, 3])))).toBe(1);
    expect(getParticleSpacingRadius(createParticles(new Float32Array([1, 2, 3, 1, 2, 3])))).toBe(1);
  });
});

describe('a particle loci label', () => {
  it('says the particle in bold, its list beneath, the list name escaped', () => {
    const list = { ...createParticles(new Float32Array([0, 0, 0, 5, 0, 0])), label: 'picks <1>.em' };
    const one = lociLabel(Particle.Loci(list, OrderedSet.ofSingleton(1)));
    expect(one).toBe('<b>1</b></br><small>picks &lt;1&gt;.em</small>');
    expect(lociLabel(Particle.Loci(list, OrderedSet.ofBounds(0, 2)))).toBe(
      '<b>2 particles</b></br><small>picks &lt;1&gt;.em</small>',
    );
  });
});
