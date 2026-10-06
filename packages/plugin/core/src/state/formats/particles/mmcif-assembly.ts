/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  type MmcifVariant,
  getAssemblyIdsFromMmcif,
  getAsymIdsFromMmcif,
  createParticleListFromMmcifAssembly,
  looksLikeMmcifParticles,
} from '@molstar/model/formats/particles/mmcif';
import { Task } from '@molstar/core/task';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { ParticlesFormatCategory } from './category.js';
import type { PluginContext } from '@molstar/plugin/context';
import { type ParticleFormatData, complexVisuals } from './provider.js';
import { ParseCif } from '@molstar/plugin/state/formats/cif';

export { ParticleListFromMmcifAssembly };
type ParticleListFromMmcifAssembly = typeof ParticleListFromMmcifAssembly;
const ParticleListFromMmcifAssembly = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-mmcif-assembly',
  display: {
    name: 'Particle List from mmCIF Assembly',
    description: 'Create ParticleList from a mmCIF (CellPack/PetWorld) assembly.',
  },
  from: SO.Format.Cif,
  to: SO.Particle.List,
  params: (a) => {
    const variant = PD.Select<MmcifVariant>(
      'auto',
      [
        ['auto', 'Auto'],
        ['cellpack', 'CellPack'],
        ['standard', 'Standard'],
        ['petworld', 'PetWorld'],
      ],
      {
        description:
          'mmCIF variant used to interpret entity names/compartments. Auto detects from the file header and categories.',
      },
    );
    const label = PD.Optional(PD.Text(''));
    const resolveTargets = PD.Boolean(true, {
      description:
        'Build the per-target reference structures from the same mmCIF data and instance them at each particle.',
    });

    if (!a) {
      return {
        assemblyId: PD.Text('', { description: 'Assembly identifier from _pdbx_struct_assembly.id.' }),
        asymIds: PD.MultiSelect<string>([], [], {
          description: 'Empty selection includes all chains for the assembly.',
        }),
        variant,
        label,
        resolveTargets,
      };
    }

    let assemblyIds: string[] = [];
    try {
      assemblyIds = getAssemblyIdsFromMmcif(a.data);
    } catch {
      // ignore; apply will surface parse errors
    }
    const assemblyOptions = assemblyIds.map((id) => [id, id] as [string, string]);
    const defaultAssemblyId = assemblyIds.length > 0 ? assemblyIds[0] : '';

    let asymIds: string[] = [];
    try {
      asymIds = defaultAssemblyId ? getAsymIdsFromMmcif(a.data, defaultAssemblyId) : [];
    } catch {
      // ignore; apply will surface parse errors
    }
    const asymOptions = asymIds.map((id) => [id, id] as [string, string]);

    return {
      assemblyId: PD.Select(defaultAssemblyId, assemblyOptions, {
        description: 'Assembly identifier from _pdbx_struct_assembly.id.',
      }),
      asymIds: PD.MultiSelect<string>([], asymOptions, {
        description: 'Empty selection includes all chains for the assembly.',
      }),
      variant,
      label,
      resolveTargets,
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Create Particle List from mmCIF Assembly', async (ctx) => {
      // Params default to the first assembly ID only when interactively edited; fall back
      // here too so auto-loading a file (which calls parse() with no params) still works.
      const assemblyId = params.assemblyId || getAssemblyIdsFromMmcif(a.data)[0];
      if (!assemblyId) {
        throw new Error('mmCIF file contains no _pdbx_struct_assembly_gen assemblies; cannot create a particle list.');
      }
      const list = await createParticleListFromMmcifAssembly(ctx, a.data, {
        assemblyId,
        asymIds: params.asymIds.length > 0 ? params.asymIds : void 0,
        variant: params.variant,
        label: params.label || a.label,
        resolveTargets: params.resolveTargets,
      });
      return new SO.Particle.List(list, { label: list.label || 'Particles', description: 'mmCIF Particle List' });
    });
  },
});

export const MmcifParticlesProvider = DataFormatProvider({
  label: 'mmCIF Particles',
  description: 'mmCIF Particles (CellPack / PetWorld assemblies)',
  category: ParticlesFormatCategory,
  stringExtensions: ['cif'],
  binaryExtensions: ['bcif'],
  /**
   * Higher than the default (0) mmCIF/CifCore trajectory providers so that CellPack/PetWorld
   * assemblies are recognized ahead of them during auto-detection, regardless of registration order.
   */
  priority: 10,
  isApplicable: (info, data) => looksLikeMmcifParticles(info, data),
  parse: async (
    plugin,
    data,
    params?: { label?: string; assemblyId?: string; asymIds?: string[]; variant?: MmcifVariant },
  ) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseCif, void 0, { state: { isGhost: true } });

    const list = format.apply(ParticleListFromMmcifAssembly, {
      assemblyId: params?.assemblyId,
      asymIds: params?.asymIds,
      variant: params?.variant,
      label: params?.label,
    });

    await format.commit({ revertOnError: true });

    return { format: format.selector, list: list.selector };
  },
  visuals: (plugin: PluginContext, data: ParticleFormatData) => complexVisuals(plugin, data),
});
