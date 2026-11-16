/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OrderedSet, SortedArray } from '@molstar/core/data/int';
import type { Theme } from '@molstar/graphics/theme/theme';
import { createTargetColors, getTargetTheme } from '../../repr/particles/representation/target/common.js';
import { LocationIterator } from '../../geo/util/location-iterator.js';
import { createValueColor } from '../../geo/geometry/color-data.js';
import { Color } from '@molstar/core/util/color';
import { Particle, type ParticleList, ParticleValues } from '@molstar/model/model/particles/particle-list';
import {
  ParticleCompartmentColorTheme,
  ParticleCompartmentColorThemeParams,
} from '@molstar/graphics/theme/color/particle-compartment';
import { ParticleSizeTheme, ParticleSizeThemeParams } from '@molstar/graphics/theme/size/particle-size';
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

  it('color and size particles of runs (a pack list) as of arrays', () => {
    // three targets as runs, their compartments, entities and radii a value each
    const offsets = Int32Array.from([0, 2, 5, 6]);
    const list = (runs: boolean): ParticleList => {
      const of = <A extends Int32Array | Float32Array>(values: A, ctor: { new (n: number): A }) =>
        runs
          ? ParticleValues.ofRuns(offsets, values)
          : ParticleValues.toArray(ParticleValues.ofRuns(offsets, values), 6, ctor);
      return {
        ...particleList({ count: 0 }),
        targets: of(Int32Array.from([0, 1, 2]), Int32Array),
        compartments: of(Int32Array.from([1, 0, -1]), Int32Array),
        entities: of(Int32Array.from([2, 0, 1]), Int32Array),
        radii: of(Float32Array.from([3, 5, 7]), Float32Array),
      } as ParticleList;
    };
    const arrays = list(false);
    const runs = list(true);
    const themes = [
      (p: ParticleList) =>
        ParticleEntityColorTheme({ particles: p }, PD.getDefaultValues(ParticleEntityColorThemeParams)),
      (p: ParticleList) =>
        ParticleCompartmentColorTheme({ particles: p }, PD.getDefaultValues(ParticleCompartmentColorThemeParams)),
      (p: ParticleList) =>
        ParticleHierarchyColorTheme({ particles: p }, PD.getDefaultValues(ParticleHierarchyColorThemeParams)),
    ];
    for (const theme of themes) expect(colorsOf(theme(runs), runs)).toEqual(colorsOf(theme(arrays), arrays));
    const sizes = (p: ParticleList) => {
      const theme = ParticleSizeTheme({ particles: p }, PD.getDefaultValues(ParticleSizeThemeParams));
      return Array.from({ length: 6 }, (_, i) => theme.size(Particle.Location(p, i)));
    };
    expect(sizes(runs)).toEqual(sizes(arrays));
    expect(sizes(runs)).toEqual([3, 3, 5, 5, 5, 7]);
  });

  it('colors a target of one run of entities with one color, uniform; one of several or of arrays a particle each', () => {
    const offsets = Int32Array.from([0, 2, 5, 6]);
    const entities = Int32Array.from([2, 0, 1]);
    const runs = {
      ...particleList({ count: 0 }),
      count: 6,
      entities: ParticleValues.ofRuns(offsets, entities),
    } as ParticleList;
    const arrays = {
      ...runs,
      entities: ParticleValues.toArray(ParticleValues.ofRuns(offsets, entities), 6, Int32Array),
    };
    const props = PD.getDefaultValues(ParticleEntityColorThemeParams);
    const color = ParticleEntityColorTheme({ particles: runs }, props);
    expect(Array.from(color.particleRuns!.offsets)).toEqual(Array.from(offsets));
    expect(ParticleEntityColorTheme({ particles: arrays }, props).particleRuns).toBeUndefined();
    const theme = { color, size: {} } as unknown as Theme;
    // the second run's particles: uniform, the color of its first
    const one = getTargetTheme(theme, OrderedSet.ofBounds(2, 5)).color;
    expect(one.granularity).toBe('uniform');
    expect(one.color(Particle.Location(runs, 2), false)).toBe(color.color(Particle.Location(runs, 4), false));
    // across two runs, or no runs: a color a particle
    expect(getTargetTheme(theme, OrderedSet.ofBounds(1, 5)).color.granularity).toBe('instance');
    expect(getTargetTheme(theme, SortedArray.ofSortedArray([2, 5])).color.granularity).toBe('instance');
  });

  it("colors a target of one run as the object (a texel), as a merged object keeps its members' colors", () => {
    const offsets = Int32Array.from([0, 2, 5, 6]);
    const runs = {
      ...particleList({ count: 0 }),
      count: 6,
      entities: ParticleValues.ofRuns(offsets, Int32Array.from([2, 0, 1])),
    } as ParticleList;
    const color = ParticleEntityColorTheme({ particles: runs }, PD.getDefaultValues(ParticleEntityColorThemeParams));
    const theme = { color, size: {} } as unknown as Theme;
    const indices = OrderedSet.ofBounds(2, 5);
    const it = LocationIterator(1, 3, 1, (_g, i) => Particle.Location(runs, 2 + i));
    const data = createTargetColors(it, it, theme, indices, createValueColor(Color(0)));
    expect(data.colorType.ref.value).toBe('object');
    expect(data.dColorType.ref.value).toBe('texture');
    expect(data.tColor.ref.value.count).toBe(1);
    const expected = Color.toArray(color.color(Particle.Location(runs, 3), false), [0, 0, 0], 0);
    const texels = data.tColor.ref.value.array;
    expect(Array.from(texels.subarray(0, 4))).toEqual([...expected, 255]);
  });
});
