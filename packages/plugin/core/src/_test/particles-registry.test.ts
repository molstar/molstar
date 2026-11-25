/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultRegistry } from '@molstar/plugin/default-registry';
import {
  ParticleActions,
  ParticleAnimations,
  ParticleColorThemes,
  ParticleFormats,
  ParticleRepresentations,
  ParticleSizeThemes,
  ParticleThemes,
  ParticlesRegistry,
} from '@molstar/plugin/particles-registry';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { writeParticlePack } from '@molstar/model/formats/particles/pack';

/** Five oriented particles in a CryoET Data Portal NDJSON file */
const Ndjson = [0, 1, 2, 3, 4]
  .map((i) =>
    JSON.stringify({
      type: 'orientedPoint',
      location: { x: i * 10, y: i * 5, z: 1 },
      xyz_rotation_matrix: [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ],
    }),
  )
  .join('\n');

async function createPlugin(withParticles: boolean) {
  const plugin = new PluginContext({
    behaviors: [],
    registry: withParticles ? [...DefaultRegistry, ...ParticlesRegistry] : DefaultRegistry,
  });
  await plugin.init();
  return plugin;
}

describe('ParticlesRegistry', () => {
  it('lists its entries in order', () => {
    expect(ParticlesRegistry).toEqual([
      ParticleThemes,
      ParticleRepresentations,
      ParticleActions,
      ParticleFormats,
      ParticleAnimations,
    ]);
    expect(ParticleThemes.particles!.themes!.color).toEqual(expect.arrayContaining([...ParticleColorThemes]));
    expect(ParticleThemes.particles!.themes!.size).toEqual(expect.arrayContaining([...ParticleSizeThemes]));
  });

  it('registers the particles scope, its formats and its animation in a plugin that composes it', async () => {
    const plugin = await createPlugin(true);
    const { particles } = plugin.representation;
    expect(particles.registry.default?.name).toBe('spacefill');
    expect(particles.registry.list.map((e) => e.name)).toEqual(['spacefill', 'orientation', 'fibers', 'target']);
    for (const name of ['particle-index', 'particle-entity', 'uniform', 'external-structure']) {
      expect(particles.themes.colorThemeRegistry.has(name)).toBe(true);
    }
    expect(particles.themes.sizeThemeRegistry.has('particle-size')).toBe(true);
    // the particle themes are the particles scope's only
    expect(plugin.representation.structure.themes.colorThemeRegistry.has('particle-index')).toBe(false);
    for (const f of ParticleFormats.formats!) expect(plugin.dataFormats.has(f.name)).toBe(true);
    expect(plugin.managers.animation.animations.map((a) => a.name)).toContain('built-in.animate-particle-trajectory');
    plugin.dispose();
  });

  it('loads a particle file, which a plugin without it does not recognize', async () => {
    const without = await createPlugin(false);
    const file = { name: 'picks.ndjson', ext: 'ndjson' } as any;
    expect(without.dataFormats.auto(file, { data: Ndjson } as any)).toBeUndefined();
    without.dispose();

    const plugin = await createPlugin(true);
    const data = await plugin.builders.data.rawData({ data: Ndjson, label: 'picks.ndjson' });
    const provider = plugin.dataFormats.get('cryoet_ndjson_particles')!;
    const parsed = await provider.parse(plugin, data);
    await provider.visuals?.(plugin, parsed);

    const { lists } = plugin.managers.particles.hierarchy.current;
    expect(lists).toHaveLength(1);
    expect(lists[0].cell.obj!.data.count).toBe(5);
    expect(lists[0].representations.map((r) => r.cell.transform.params?.type.name)).toEqual(['spacefill']);
    plugin.dispose();
  });

  it('loads a particle pack (its extension, its magic checked), its list a view of the bytes', async () => {
    const plugin = await createPlugin(true);
    const bytes = writeParticlePack({
      count: 3,
      keys: Int32Array.from([0, 1, 2]),
      targets: Int32Array.from([0, 0, 1]),
      targetInfo: new Map([
        [0, { name: 'a' }],
        [1, { name: 'b' }],
      ]),
      coordinates: Float32Array.from([0, 0, 0, 10, 0, 0, 0, 10, 0]),
      rotations: Float32Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]),
      getParticleLabel: (i) => String(i),
      sourceData: { kind: 'test', name: 'test', data: undefined },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    });
    expect(plugin.dataFormats.auto({ name: 'cell.mspk', ext: 'mspk' } as any, { data: bytes } as any)?.name).toBe(
      'particle_pack',
    );
    // by its extension, its magic checked: other bytes named .mspk are not one
    expect(
      plugin.dataFormats.auto({ name: 'cell.mspk', ext: 'mspk' } as any, { data: new Uint8Array(128) } as any),
    ).toBeUndefined();
    const data = await plugin.builders.data.rawData({ data: bytes, label: 'cell.mspk' });
    const provider = plugin.dataFormats.get('particle_pack')!;
    const parsed = await provider.parse(plugin, data);
    await provider.visuals?.(plugin, parsed);

    const { lists } = plugin.managers.particles.hierarchy.current;
    expect(lists).toHaveLength(1);
    const list = lists[0].cell.obj!.data;
    expect(list.count).toBe(3);
    expect(GridPositions.is(list.coordinates) && list.coordinates.data.buffer).toBe(bytes.buffer);
    expect(lists[0].representations.length).toBeGreaterThan(0);
    plugin.dispose();
  });
});
