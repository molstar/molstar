/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';
import {
  ParticleEntityColorTheme,
  ParticleEntityColorThemeParams,
} from '@molstar/graphics/theme/color/particle-entity';
import {
  ParticleHierarchyColorTheme,
  ParticleHierarchyColorThemeParams,
} from '@molstar/graphics/theme/color/particle-hierarchy';

/** Counts the reads of an index property (the per-particle passes) */
function counted(array: Int32Array, reads: { count: number }): Int32Array {
  return new Proxy(array, {
    get(t, k) {
      if (typeof k === 'string' && k.length > 0 && k.charCodeAt(0) >= 48 && k.charCodeAt(0) <= 57) reads.count++;
      const v = Reflect.get(t, k, t);
      return typeof v === 'function' ? v.bind(t) : v;
    },
  });
}

/** Two compartments; compartment 0 has entities 0 and 1, compartment 1 has entity 2 and a particle without one */
function particleList(reads: { count: number }): ParticleList {
  const compartments = Int32Array.from([0, 0, 0, 1, 1, -1]);
  const entities = Int32Array.from([0, 1, 1, 2, -1, 0]);
  return {
    count: 6,
    keys: Int32Array.from([0, 1, 2, 3, 4, 5]),
    targets: new Int32Array(6),
    targetInfo: new Map(),
    compartments: counted(compartments, reads),
    compartmentInfo: new Map([
      [0, { name: 'a' }],
      [1, { name: 'b' }],
    ]),
    entities: counted(entities, reads),
    entityInfo: new Map([
      [0, { name: 'e0', function: 'f' }],
      [1, { name: 'e1', function: 'f' }],
      [2, { name: 'e2' }],
    ]),
  } as unknown as ParticleList;
}

function colorsOf(theme: { color: (l: any, isSecondary: boolean) => number }, particles: ParticleList) {
  const ret: number[] = [];
  for (let i = 0; i < particles.count; i++) ret.push(theme.color(Particle.Location(particles, i), false));
  return ret;
}

describe('particle color themes', () => {
  it('particle-hierarchy indexes the particles once per list and grouping (perf-audit-2 09/P6)', () => {
    const reads = { count: 0 };
    const particles = particleList(reads);
    const ctx = { particles };
    const props = PD.getDefaultValues(ParticleHierarchyColorThemeParams);

    const first = colorsOf(ParticleHierarchyColorTheme(ctx, props), particles);
    expect(reads.count).toBeGreaterThan(0);
    // entities 0 and 1 of compartment 0 differ; compartment 1's particle without an entity has its base color
    expect(first[0]).not.toBe(first[1]);
    expect(first[1]).toBe(first[2]);
    expect(first[4]).not.toBe(first[3]);
    expect(first[5]).toBe(0xcccccc);

    // another palette: the colors change, the particles are not read again
    reads.count = 0;
    const shifted = ParticleHierarchyColorTheme(ctx, { ...props, shift: 50 });
    expect(reads.count).toBe(0);
    expect(colorsOf(shifted, particles)).not.toEqual(first);
    // the same props: the same colors
    expect(colorsOf(ParticleHierarchyColorTheme(ctx, props), particles)).toEqual(first);

    // another grouping is its own index
    const byFunction = colorsOf(ParticleHierarchyColorTheme(ctx, { ...props, by: 'function' }), particles);
    expect(reads.count).toBeGreaterThan(0);
    expect(byFunction[0]).toBe(byFunction[1]);
  });

  it('particle-entity indexes the particles once per list and grouping (perf-audit-2 09/P6)', () => {
    const reads = { count: 0 };
    const particles = particleList(reads);
    const ctx = { particles };
    const props = PD.getDefaultValues(ParticleEntityColorThemeParams);

    const first = colorsOf(ParticleEntityColorTheme(ctx, props), particles);
    expect(reads.count).toBeGreaterThan(0);
    expect(first[0]).not.toBe(first[1]);
    expect(first[0]).toBe(first[5]);

    reads.count = 0;
    expect(colorsOf(ParticleEntityColorTheme(ctx, props), particles)).toEqual(first);
    expect(reads.count).toBe(0);

    // an entity added to the shared info (Simularium frames do): indexed again, the palette spans it
    (particles.entityInfo as Map<number, { name: string }>).set(3, { name: 'e3' });
    ParticleEntityColorTheme(ctx, props);
    expect(reads.count).toBeGreaterThan(0);
  });
});
