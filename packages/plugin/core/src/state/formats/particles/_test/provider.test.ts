/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticlePageCache } from '@molstar/model/model/particles/particle-page-cache';
import { getComplexVisualsPlan } from '../provider.js';

// The default visuals of a particle list: by its leaves where it has templates (their targets and infos), not by the
// top list's alone (its particles placing copies)

const leaf = (name: string): ParticleTarget => ({ kind: 'shape', shape: { name } as never });

function createList(fields: Partial<ParticleList> & Pick<ParticleList, 'count' | 'targets'>) {
  return {
    targetInfo: new Map(),
    coordinates: new Float32Array(fields.count * 3),
    getParticleLabel: (i: number) => `particle ${i + 1}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    ...fields,
  } as ParticleList;
}

/** two targeted particles, one without (target 1 unmapped), entities on them */
function createTemplate() {
  return createList({
    label: 'template',
    count: 3,
    targets: Int32Array.of(0, 0, 1),
    targetMapping: new Map([[0, leaf('a')]]),
    entities: Int32Array.of(0, 0, 0),
    entityInfo: new Map([[0, { name: 'A' }]]),
  });
}

/** two particles placing the template (target 2), no infos of its own */
function createTop(template: ParticleList) {
  return createList({
    label: 'top',
    count: 2,
    targets: Int32Array.of(2, 2),
    targetMapping: new Map<number, ParticleTarget>([[2, { kind: 'particles', particles: template }]]),
    fibers: { count: 1 } as never,
  });
}

describe('the default particle visuals', () => {
  it('draws a list with templates by its leaves: their untargeted, colored by their infos, no fibers (not kept)', () => {
    const plan = getComplexVisualsPlan(createTop(createTemplate()));
    expect(plan).toEqual({ hasTargets: true, hasUntargeted: true, hasFibers: false, colorTheme: 'particle-entity' });
  });

  it('leaves out the untargeted of a template read by range (flattening refuses it)', () => {
    const template = createTemplate();
    Particle.setDeferred(
      template,
      new ParticlePageCache().deferred(template, () => Promise.reject(new Error('not read'))),
    );
    expect(getComplexVisualsPlan(createTop(template)).hasUntargeted).toBe(false);
  });

  it('draws a flat list as it is', () => {
    const plan = getComplexVisualsPlan({ ...createTemplate(), fibers: { count: 2 } as never });
    expect(plan).toEqual({ hasTargets: true, hasUntargeted: true, hasFibers: true, colorTheme: 'particle-entity' });
  });
});
