/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Asset } from '@molstar/core/util/assets';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Particle, type ParticleList, type ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import { writeParticlePack } from '@molstar/model/formats/particles/pack';
import { PluginContext } from '@molstar/plugin/context';
import { needsPackStreamReread, ParticleListFromPackStream } from '../pack.js';

// A nested pack's by-range transform where no GPU cull draws it (a session restored on WebGL, no canvas): read whole

function createNested(): ParticleList {
  const n = 300;
  const base = {
    getParticleLabel: (i: number) => `${i}`,
    sourceData: { kind: 'test', name: 'test', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
  const template: ParticleList = {
    ...base,
    label: 'template',
    count: n,
    coordinates: Float32Array.from({ length: n * 3 }, (_, k) => (k * 7) % 50),
    targets: Int32Array.from({ length: n }, (_, i) => i % 2),
    targetInfo: new Map([
      [0, { name: 'a' }],
      [1, { name: 'b' }],
    ]),
  };
  return {
    ...base,
    _propertyData: Object.create(null),
    customProperties: new CustomProperties(),
    label: 'nested',
    count: 2,
    coordinates: Float32Array.of(0, 0, 0, 100, 0, 0),
    targets: new Int32Array(2),
    targetInfo: new Map([[0, { name: 'copy' }]]),
    targetMapping: new Map<number, ParticleTarget>([[0, { kind: 'particles', particles: template }]]),
  };
}

describe('a nested particle pack by range', () => {
  it('is read whole where no GPU cull draws it by range', async () => {
    const plugin = new PluginContext({ behaviors: [], registry: [] });
    await plugin.init();
    try {
      const bytes = writeParticlePack(createNested());
      const file = Asset.File(new File([bytes as BlobPart], 'nested.mspk'));
      const list = await plugin
        .build()
        .toRoot()
        .apply(ParticleListFromPackStream, { source: { name: 'file', params: file } })
        .commit();
      const data = list.data!;
      expect(data.count).toBe(2);
      const template = ParticleNesting.getTemplate(data, 0)!.particles;
      expect(template.count).toBe(300);
      expect(Particle.getDeferred(template)).toBeUndefined();
      expect(list.obj!.description).toBe('Particle Pack');
      // read again (`reread`, a backend change): as the backend draws it, whatever was saved before
      await plugin
        .build()
        .to(list.ref)
        .update({ source: { name: 'file', params: file }, reread: 1 })
        .commit();
      expect(Particle.getDeferred(ParticleNesting.getTemplate(list.data!, 0)!.particles)).toBeUndefined();
    } finally {
      plugin.dispose();
    }
  });

  it('is read again after a backend change where its list is not as the new backend draws it', () => {
    // the GPU's cull gone: its templates not read read again (whole), a list read whole as it is
    expect(needsPackStreamReread(false, true)).toBe(true);
    expect(needsPackStreamReread(false, false)).toBe(false);
    // the cull there: a list read whole read again (by range), one read by range as it is
    expect(needsPackStreamReread(true, false)).toBe(true);
    expect(needsPackStreamReread(true, true)).toBe(false);
  });
});
