/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList, ParticleTarget } from '@molstar/model/model/particles/particle-list';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultRegistry } from '@molstar/plugin/default-registry';
import { ParticlesRegistry } from '@molstar/plugin/particles-registry';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParticlesRepresentation3D } from '../representation.js';

// A list with templates' representation types: those of its theme list (a slot list of its leaves: the target's) and
// of the list itself (drawn flattened: spacefill's), not the theme list's alone (a spacefill normalized to a target)

function createNested(): ParticleList {
  const base = () => ({
    getParticleLabel: (i: number) => `${i}`,
    sourceData: { kind: 'test', name: 'test', data: undefined },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  });
  const template: ParticleList = {
    ...base(),
    label: 'template',
    count: 4,
    coordinates: Float32Array.of(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1),
    targets: Int32Array.of(0, 0, 1, 1),
    targetInfo: new Map([[0, { name: 'a' }]]),
    targetMapping: new Map<number, ParticleTarget>([[0, { kind: 'shape', shape: { name: 'a' } as never }]]),
  };
  return {
    ...base(),
    label: 'nested',
    count: 2,
    coordinates: Float32Array.of(0, 0, 0, 10, 0, 0),
    targets: new Int32Array(2),
    targetInfo: new Map([[0, { name: 'copy' }]]),
    targetMapping: new Map<number, ParticleTarget>([[0, { kind: 'particles', particles: template }]]),
  };
}

describe("a nested list's particle representation", () => {
  it('offers the types of its leaves and of the list itself, a spacefill kept as one', async () => {
    const plugin = new PluginContext({ behaviors: [], registry: [...DefaultRegistry, ...ParticlesRegistry] });
    await plugin.init();
    try {
      const a = new SO.Particle.List(createNested(), { label: 'nested' });
      const params = ParticlesRepresentation3D.definition.params!(a, plugin) as PD.Params;
      const type = params.type as PD.Mapped<any>;
      const names = type.select.options.map(([name]: [string, string]) => name);
      expect(names).toContain('target');
      expect(names).toContain('spacefill');
      const value = PD.normalizeParams(params, { type: { name: 'spacefill', params: {} } }, 'children');
      expect(value.type.name).toBe('spacefill');
    } finally {
      plugin.dispose();
    }
  });
});
