/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { parseCryoEtDataPortalNdjson } from '@molstar/io/reader/cryoet/ndjson';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { createParticleListFromCryoEtDataPortalNdjson } from '@molstar/model/formats/particles/ndjson';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { ParticlesFormatCategory } from './category.js';
import { simpleVisuals, SimpleParticleVisuals } from './provider.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { ParseCryoEtDataPortalNdjson };
type ParseCryoEtDataPortalNdjson = typeof ParseCryoEtDataPortalNdjson;
const ParseCryoEtDataPortalNdjson = PluginStateTransform.BuiltIn({
  name: 'parse-cryoet-data-portal-ndjson',
  display: { name: 'Parse CryoET NDJSON', description: 'Parse CryoET Data Portal NDJSON from String data' },
  from: [SO.Data.String],
  to: SO.Format.CryoEtDataPortalNdjson,
})({
  apply({ a }) {
    return Task.create('Parse CryoET NDJSON', async (ctx) => {
      const parsed = await parseCryoEtDataPortalNdjson(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.CryoEtDataPortalNdjson(parsed.result, { label: a.label || 'CryoET NDJSON' });
    });
  },
});

export { ParticleListFromCryoEtDataPortalNdjson };
type ParticleListFromCryoEtDataPortalNdjson = typeof ParticleListFromCryoEtDataPortalNdjson;
const ParticleListFromCryoEtDataPortalNdjson = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-cryoet-data-portal-ndjson',
  display: {
    name: 'Particle List from CryoET NDJSON',
    description: 'Create ParticleList from CryoET Data Portal NDJSON data.',
  },
  from: SO.Format.CryoEtDataPortalNdjson,
  to: SO.Particle.List,
  params: {
    pixelSize: PD.Numeric(
      1,
      { min: 0, step: 0.001 },
      {
        description:
          'Pixel size in Å/pixel used to convert pixel-space NDJSON coordinates to angstrom. Required because CryoET Data Portal NDJSON does not encode distance units.',
      },
    ),
    type: PD.Optional(PD.Text('')),
    particleRadius: PD.Numeric(
      0,
      { min: 0, step: 0.1 },
      { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' },
    ),
  },
})({
  apply({ a, params }) {
    return Task.create('Create Particle List from CryoET NDJSON', async (ctx) => {
      const list = createParticleListFromCryoEtDataPortalNdjson(a.data, {
        pixelSize: params.pixelSize,
        type: params.type || void 0,
        particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
      });
      return new SO.Particle.List(list, {
        label: list.label || 'Particles',
        description: 'CryoET NDJSON Particle List',
      });
    });
  },
});

export const CryoEtDataPortalNdjsonParticlesProvider = DataFormatProvider({
  name: 'cryoet_ndjson_particles',
  label: 'CryoET NDJSON Particles',
  description: 'CryoET NDJSON Particles',
  category: ParticlesFormatCategory,
  stringExtensions: ['ndjson'],
  parse: async (plugin, data, params?: { label?: string; type?: string }) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseCryoEtDataPortalNdjson, void 0, { state: { isGhost: true } });

    const list = format.apply(ParticleListFromCryoEtDataPortalNdjson, {
      type: params?.type,
    });

    await format.commit({ revertOnError: true });

    return { format: format.selector, list: list.selector };
  },
  visuals: simpleVisuals,
});

/** The CryoEtDataPortalNdjsonParticles data format with its actions and the representations its visuals use. */
export const CryoEtDataPortalNdjsonParticles: PluginRegistryEntry = {
  formats: [CryoEtDataPortalNdjsonParticlesProvider],
  actions: [ParticleListFromCryoEtDataPortalNdjson],
  ...SimpleParticleVisuals,
};
