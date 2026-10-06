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
import { parseDynamoTbl } from '@molstar/io/reader/dynamo/tbl';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { getDynamoTblTomogramIds, createParticleListFromDynamoTbl } from '@molstar/model/formats/particles/tbl';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { ParticlesFormatCategory } from './category.js';
import { simpleVisuals, SimpleParticleVisuals } from './provider.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { ParseDynamoTbl };
type ParseDynamoTbl = typeof ParseDynamoTbl;
const ParseDynamoTbl = PluginStateTransform.BuiltIn({
  name: 'parse-dynamo-tbl',
  display: { name: 'Parse Dynamo TBL', description: 'Parse Dynamo TBL from String data' },
  from: [SO.Data.String],
  to: SO.Format.DynamoTbl,
})({
  apply({ a }) {
    return Task.create('Parse Dynamo TBL', async (ctx) => {
      const parsed = await parseDynamoTbl(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.DynamoTbl(parsed.result, { label: a.label || 'Dynamo TBL' });
    });
  },
});

export { ParticleListFromDynamoTbl };
type ParticleListFromDynamoTbl = typeof ParticleListFromDynamoTbl;
const ParticleListFromDynamoTbl = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-dynamo-tbl',
  display: { name: 'Particle List from Dynamo TBL', description: 'Create ParticleList from Dynamo TBL data.' },
  from: SO.Format.DynamoTbl,
  to: SO.Particle.List,
  params: (a) => {
    if (!a) {
      return {
        tomos: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all tomograms.' }),
        pixelSize: PD.Optional(
          PD.Numeric(
            0,
            { min: 0, step: 0.001 },
            {
              description:
                'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from the table’s `apix` field.',
            },
          ),
        ),
        particleRadius: PD.Numeric(
          0,
          { min: 0, step: 0.1 },
          { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' },
        ),
      };
    }
    const ids = getDynamoTblTomogramIds(a.data);
    const options = ids.map((id) => [String(id), String(id)] as [string, string]);
    const defaultValue = ids.length > 0 ? [String(ids[0])] : [];
    return {
      tomos: PD.MultiSelect<string>(defaultValue, options, { description: 'Empty selection includes all tomograms.' }),
      pixelSize: PD.Optional(
        PD.Numeric(
          0,
          { min: 0, step: 0.001 },
          {
            description:
              'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from the table’s `apix` field.',
          },
        ),
      ),
      particleRadius: PD.Numeric(
        0,
        { min: 0, step: 0.1 },
        { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' },
      ),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create Particle List from Dynamo TBL', async (ctx) => {
      const list = createParticleListFromDynamoTbl(a.data, {
        tomos: params.tomos.map((v) => Number(v)),
        pixelSize: params.pixelSize && params.pixelSize > 0 ? params.pixelSize : void 0,
        particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
      });
      return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'Dynamo Particle List' });
    });
  },
});

export const DynamoTblParticlesProvider = DataFormatProvider({
  name: 'dynamo_tbl_particles',
  label: 'Dynamo TBL Particles',
  description: 'Dynamo TBL Particles',
  category: ParticlesFormatCategory,
  stringExtensions: ['tbl'],
  parse: async (plugin, data, params?: { label?: string; tomo?: number }) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseDynamoTbl, void 0, { state: { isGhost: true } });

    const list = format.apply(ParticleListFromDynamoTbl, {
      tomos: params?.tomo !== void 0 ? [String(params.tomo)] : [],
    });

    await format.commit({ revertOnError: true });

    return { format: format.selector, list: list.selector };
  },
  visuals: simpleVisuals,
});

/** The DynamoTblParticles data format with its actions and the representations its visuals use. */
export const DynamoTblParticles: PluginRegistryEntry = {
  formats: [DynamoTblParticlesProvider],
  actions: [ParticleListFromDynamoTbl],
  ...SimpleParticleVisuals,
};
