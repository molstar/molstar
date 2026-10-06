/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  getRelionStarTomogramNames,
  getRelionStarMicrographNames,
  parseRelionStar,
} from '@molstar/io/reader/relion/star';
import { Task } from '@molstar/core/task';
import { createParticleListFromRelionStar } from '@molstar/model/formats/particles/star';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { ParticlesFormatCategory } from './category.js';
import { simpleVisuals } from './provider.js';
import { ParseCif } from '@molstar/plugin/state/formats/cif';

export { ParticleListFromRelionStar };
type ParticleListFromRelionStar = typeof ParticleListFromRelionStar;
const ParticleListFromRelionStar = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-relion-star',
  display: { name: 'Particle List from RELION STAR', description: 'Create ParticleList from RELION STAR data.' },
  from: SO.Format.Cif,
  to: SO.Particle.List,
  params: (a) => {
    if (!a) {
      return {
        tomograms: PD.MultiSelect<string>([], [], { description: 'Empty selection includes all tomograms.' }),
        micrographs: PD.MultiSelect<string>([], [], {
          description: 'Empty selection includes all micrographs. Combined with the tomogram filter using AND.',
        }),
        pixelSize: PD.Optional(
          PD.Numeric(
            0,
            { min: 0, step: 0.001 },
            {
              description:
                'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from STAR optics/particle metadata.',
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
    let tomoNames: string[] = [];
    let micrographNames: string[] = [];
    try {
      tomoNames = getRelionStarTomogramNames(a.data);
      micrographNames = getRelionStarMicrographNames(a.data);
    } catch {
      // ignore; apply will surface parse errors
    }
    const tomoOptions = tomoNames.map((n) => [n, n] as [string, string]);
    const micrographOptions = micrographNames.map((n) => [n, n] as [string, string]);
    const tomoDefault = tomoNames.length > 0 ? [tomoNames[0]] : [];
    const micrographDefault = micrographNames.length > 0 ? [micrographNames[0]] : [];
    return {
      tomograms: PD.MultiSelect<string>(tomoDefault, tomoOptions, {
        description: 'Empty selection includes all tomograms.',
      }),
      micrographs: PD.MultiSelect<string>(micrographDefault, micrographOptions, {
        description: 'Empty selection includes all micrographs. Combined with the tomogram filter using AND.',
      }),
      pixelSize: PD.Optional(
        PD.Numeric(
          0,
          { min: 0, step: 0.001 },
          {
            description:
              'Override pixel size in Å/pixel for converting pixel-space coordinates to angstrom. Leave 0 to auto-detect from STAR optics/particle metadata.',
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
    return Task.create('Create Particle List from RELION STAR', async (ctx) => {
      const relion = parseRelionStar(a.data);
      if (relion.isError) throw new Error(relion.message);

      const list = createParticleListFromRelionStar(relion.result, {
        tomograms: params.tomograms,
        micrographs: params.micrographs,
        pixelSize: params.pixelSize && params.pixelSize > 0 ? params.pixelSize : void 0,
        particleRadius: params.particleRadius > 0 ? params.particleRadius : void 0,
      });

      return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'RELION Particle List' });
    });
  },
});

export const RelionStarParticlesProvider = DataFormatProvider({
  name: 'relion_star_particles',
  label: 'RELION STAR Particles',
  description: 'RELION STAR Particles',
  category: ParticlesFormatCategory,
  stringExtensions: ['star'],
  parse: async (plugin, data, params?: { label?: string; tomogram?: string }) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseCif, void 0, { state: { isGhost: true } });

    const list = format.apply(ParticleListFromRelionStar, {
      tomograms: params?.tomogram ? [params.tomogram] : [],
    });

    await format.commit({ revertOnError: true });

    return { format: format.selector, list: list.selector };
  },
  visuals: simpleVisuals,
});
