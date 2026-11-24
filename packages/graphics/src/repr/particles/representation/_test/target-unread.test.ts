/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticlePageCache } from '@molstar/model/model/particles/particle-page-cache';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { ParticleTargetRepresentationProvider } from '../target/representation.js';

// A list with templates whose particles are not read (a pack opened by range) where the GPU's cull is not there (a
// backend switched to WebGL): drawn as nothing, not refused (the plugin reads it whole then); read, drawn.

const site = (): ParticleTarget => ({
  kind: 'sites',
  sites: { count: 1, positions: new Float32Array([0, 0, 0]), radii: new Float32Array([1]) },
});

const listOf = (count: number, targetMapping: Map<number, ParticleTarget>): ParticleList => ({
  count,
  targets: new Int32Array(count),
  targetInfo: new Map([[0, {}]]),
  targetMapping,
  coordinates: Float32Array.from({ length: count * 3 }, (_, i) => i),
  radii: new Float32Array(count).fill(1),
  getParticleLabel: (index) => `${index}`,
  sourceData: { kind: 'test', name: 'test', data: {} },
  customProperties: new CustomProperties(),
  _propertyData: Object.create(null),
});

/** two copies of a template of three particles, its particles not read where `deferred` */
function nested(deferred: boolean) {
  const template = listOf(3, new Map([[0, site()]]));
  if (deferred) {
    Particle.setDeferred(
      template,
      new ParticlePageCache().deferred(template, () => Promise.reject(new Error('not read'))),
    );
  }
  return listOf(2, new Map<number, ParticleTarget>([[0, { kind: 'particles', particles: template }]]));
}

function create() {
  const provider = ParticleTargetRepresentationProvider;
  return provider.factory(
    { colorThemeRegistry: ColorTheme.createRegistry(), sizeThemeRegistry: SizeTheme.createRegistry() },
    provider.getParams,
  );
}

describe('particle targets of templates not read, without the GPU cull', () => {
  it('draws nothing, not refused; the list read, drawn', async () => {
    const repr = create();
    await repr.createOrUpdate({ gpuCulling: false }, nested(true)).run();
    expect(repr.renderObjects).toHaveLength(0);
    await repr.createOrUpdate({}, nested(false)).run();
    expect(repr.renderObjects.length).toBeGreaterThan(0);
    repr.destroy();
  });
});
