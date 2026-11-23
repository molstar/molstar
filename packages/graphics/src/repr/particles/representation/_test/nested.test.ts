/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { ParticleHierarchyColorThemeProvider } from '@molstar/graphics/theme/color/particle-hierarchy';
import { ParticleCompartmentColorThemeProvider } from '@molstar/graphics/theme/color/particle-compartment';
import { SpacefillParticlesParamsGet, SpacefillParticlesRepresentation } from '../spacefill.js';

function createList(fields: Partial<ParticleList> & Pick<ParticleList, 'count' | 'coordinates'>): ParticleList {
  return {
    targets: new Int32Array(fields.count),
    targetInfo: new Map([[0, {}]]),
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    ...fields,
  };
}

describe('particle representation of a nested list', () => {
  it('draws a list with templates by its leaves (the flattened list), its loci the list’s by path', async () => {
    const template = createList({ count: 3, coordinates: Float32Array.of(0, 0, 0, 1, 0, 0, 2, 0, 0) });
    const nested = createList({
      count: 2,
      coordinates: Float32Array.of(0, 0, 0, 0, 10, 0),
      targetMapping: new Map<number, ParticleTarget>([
        [0, { kind: 'particles', particles: template, center: Vec3.create(0, 0, 0) }],
      ]),
    });
    const representation = SpacefillParticlesRepresentation(
      { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() },
      SpacefillParticlesParamsGet,
    );
    await representation.createOrUpdate({}, nested).run();
    const [loci] = representation.getAllLoci();
    expect(Particle.isLoci(loci)).toBe(true);
    if (!Particle.isLoci(loci)) return;
    // every instance: both copies whole
    expect(loci.particles).toBe(nested);
    expect(OrderedSet.toArray(loci.indices)).toEqual([0, 1]);
    expect(Particle.lociSize(loci)).toBe(6);
    // the flattened list's loci of the same instances
    expect(ParticleNesting.toFlatLoci(loci).particles).toBe(ParticleNesting.flatten(nested));
    expect(OrderedSet.size(ParticleNesting.toFlatLoci(loci).indices)).toBe(6);
    representation.destroy();
  });

  it('offers the hierarchy and compartment themes where its templates’ leaves have compartments and entities', () => {
    const template = createList({
      count: 2,
      coordinates: Float32Array.of(0, 0, 0, 1, 0, 0),
      compartments: Int32Array.of(0, 1),
      entities: Int32Array.of(0, 1),
      compartmentInfo: new Map([
        [0, { name: 'inner' }],
        [1, { name: 'outer' }],
      ]),
      entityInfo: new Map([
        [0, { name: 'a' }],
        [1, { name: 'b' }],
      ]),
    });
    const nested = createList({
      count: 2,
      coordinates: Float32Array.of(0, 0, 0, 0, 10, 0),
      targetMapping: new Map<number, ParticleTarget>([
        [0, { kind: 'particles', particles: template, center: Vec3.create(0, 0, 0) }],
      ]),
    });
    // the theme list's first page its top particles (no compartment of their own)
    const ctx = { particles: ParticleNesting.getThemeList(nested) };
    expect(ctx.particles.pages!.pages[0].list.compartments).toBeUndefined();
    expect(ParticleHierarchyColorThemeProvider.isApplicable!(ctx)).toBe(true);
    expect(ParticleCompartmentColorThemeProvider.isApplicable!(ctx)).toBe(true);
    // the leaves colored by their compartments: the two apart
    const theme = ParticleHierarchyColorThemeProvider.factory(ctx, ParticleHierarchyColorThemeProvider.defaultValues);
    const at = (index: number) => theme.color(Particle.Location(ctx.particles, index), false);
    expect(at(2)).not.toBe(at(3));
  });
});
