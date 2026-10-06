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
import { parseArtiatomiEm } from '@molstar/io/reader/artiatomi/em';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  getArtiatomiMotivelistTomogramIds,
  createParticleListFromArtiatomiEm,
} from '@molstar/model/formats/particles/em';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { ParticlesFormatCategory } from './category.js';
import { simpleVisuals } from './provider.js';

export { ParseArtiatomiEm };
type ParseArtiatomiEm = typeof ParseArtiatomiEm;
const ParseArtiatomiEm = PluginStateTransform.BuiltIn({
  name: 'parse-artiatomi-em',
  display: { name: 'Parse Artiatomi EM', description: 'Parse Artiatomi EM motivelist from Binary data' },
  from: [SO.Data.Binary],
  to: SO.Format.ArtiatomiEm,
})({
  apply({ a }) {
    return Task.create('Parse Artiatomi EM', async () => {
      const parsed = await parseArtiatomiEm(a.data).run();
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.ArtiatomiEm(parsed.result, { label: a.label || 'Artiatomi EM' });
    });
  },
});

export { ParticleListFromArtiatomiEm };
type ParticleListFromArtiatomiEm = typeof ParticleListFromArtiatomiEm;
const ParticleListFromArtiatomiEm = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-artiatomi-em',
  display: {
    name: 'Particle List from Artiatomi EM',
    description: 'Create ParticleList from Artiatomi EM motivelist data.',
  },
  from: SO.Format.ArtiatomiEm,
  to: SO.Particle.List,
  params: (a) => {
    if (!a) {
      return {
        tomos: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all tomograms.' }),
        pixelSize: PD.Numeric(
          1,
          { min: 0, step: 0.001 },
          {
            description:
              'Pixel size in Å/pixel used to convert voxel-space coordinates to angstrom. Required because Artiatomi EM files do not encode distance units.',
          },
        ),
        particleRadius: PD.Numeric(
          0,
          { min: 0, step: 0.1 },
          { description: 'Uniform particle radius in angstrom. Leave 0 to omit.' },
        ),
      };
    }
    const ids = getArtiatomiMotivelistTomogramIds(a.data);
    const options = ids.map((id) => [String(id), String(id)] as [string, string]);
    const defaultValue = ids.length > 0 ? [String(ids[0])] : [];
    return {
      tomos: PD.MultiSelect<string>(defaultValue, options, { description: 'Empty selection includes all tomograms.' }),
      pixelSize: PD.Numeric(
        1,
        { min: 0, step: 0.001 },
        {
          description:
            'Pixel size in Å/pixel used to convert voxel-space coordinates to angstrom. Required because Artiatomi EM files do not encode distance units.',
        },
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
    return Task.create('Create Particle List from Artiatomi EM', async () => {
      const list = createParticleListFromArtiatomiEm(a.data, {
        tomos: params.tomos.map((v) => Number(v)),
        pixelSize: params.pixelSize,
        label: a.label,
        particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
      });
      return new SO.Particle.List(list, {
        label: list.label || 'Particles',
        description: 'Artiatomi EM Particle List',
      });
    });
  },
});

export const ArtiatomiEmParticlesProvider = DataFormatProvider({
  name: 'artiatomi_em_particles',
  label: 'Artiatomi EM Particles',
  description: 'Artiatomi EM Particles',
  category: ParticlesFormatCategory,
  binaryExtensions: ['em'],
  parse: async (plugin, data, params?: { label?: string; tomo?: number; pixelSize?: number }) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseArtiatomiEm, void 0, { state: { isGhost: true } });

    const list = format.apply(ParticleListFromArtiatomiEm, {
      tomos: params?.tomo !== void 0 ? [String(params.tomo)] : [],
      pixelSize: params?.pixelSize ?? 1,
    });

    await format.commit({ revertOnError: true });

    return { format: format.selector, list: list.selector };
  },
  visuals: simpleVisuals,
});
