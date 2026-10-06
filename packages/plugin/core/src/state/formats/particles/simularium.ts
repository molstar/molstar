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
import { parseSimularium } from '@molstar/io/reader/simularium/parser';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  getSimulariumAgentTypeNames,
  createSimulariumParticleTrajectory,
  getSimulariumFrameCount,
} from '@molstar/model/formats/particles/simularium';
import type { PluginContext } from '@molstar/plugin/context';
import type { Asset } from '@molstar/core/util/assets';
import { createSimulariumGeometryResolver } from '@molstar/plugin/state/helpers/particle-targets';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { ParticlesFormatCategory } from './category.js';
import { type ParticleFormatData, complexVisuals } from './provider.js';
import { ParticleListFromTrajectory } from '@molstar/plugin/state/transforms/particles/ops';

export { ParseSimularium };
type ParseSimularium = typeof ParseSimularium;
const ParseSimularium = PluginStateTransform.BuiltIn({
  name: 'parse-simularium',
  display: { name: 'Parse Simularium', description: 'Parse a Simularium trajectory (JSON or binary) from Binary data' },
  from: [SO.Data.Binary],
  to: SO.Format.Simularium,
})({
  apply({ a }) {
    return Task.create('Parse Simularium', async (ctx) => {
      const parsed = await parseSimularium(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Simularium(parsed.result, { label: a.label || 'Simularium' });
    });
  },
});

export { ParticleTrajectoryFromSimularium };
type ParticleTrajectoryFromSimularium = typeof ParticleTrajectoryFromSimularium;
const ParticleTrajectoryFromSimularium = PluginStateTransform.BuiltIn({
  name: 'particle-trajectory-from-simularium',
  display: {
    name: 'Particle Trajectory from Simularium',
    description: 'Create a ParticleTrajectory wrapping all frames of a Simularium file.',
  },
  from: SO.Format.Simularium,
  to: SO.Particle.Trajectory,
  params: (a) => {
    if (!a) {
      return {
        types: PD.MultiSelect<string>([], [], {
          description: 'Agent types to include. Empty selection includes all types.',
        }),
        scale: PD.Numeric(
          0,
          { min: 0, step: 0.001 },
          { description: 'Spatial scale to angstrom. Leave 0 to auto-detect from the file spatial units.' },
        ),
        loadGeometries: PD.Boolean(true, {
          description:
            'Load the PDB structures and OBJ meshes referenced by the agent types and instance them at each particle of the matching type.',
        }),
      };
    }
    const typeOptions = getSimulariumAgentTypeNames(a.data).map((t) => [String(t.id), t.name] as [string, string]);
    return {
      types: PD.MultiSelect<string>([], typeOptions, {
        description: 'Agent types to include. Empty selection includes all types.',
      }),
      scale: PD.Numeric(
        0,
        { min: 0, step: 0.001 },
        { description: 'Spatial scale to angstrom. Leave 0 to auto-detect from the file spatial units.' },
      ),
      loadGeometries: PD.Boolean(true, {
        description:
          'Load the PDB structures and OBJ meshes referenced by the agent types and instance them at each particle of the matching type.',
      }),
    };
  },
})({
  apply({ a, params, cache }, plugin: PluginContext) {
    return Task.create('Particle Trajectory from Simularium', async (ctx) => {
      const assets: Asset.Wrapper[] = [];
      (cache as any).assets = assets;
      const traj = await createSimulariumParticleTrajectory(a.data, {
        scale: params.scale && params.scale > 0 ? params.scale : void 0,
        typeFilter: params.types.length > 0 ? params.types.map((v) => Number(v)) : void 0,
        resolveGeometry: params.loadGeometries ? createSimulariumGeometryResolver(plugin, ctx, assets) : void 0,
      });
      const frameCount = getSimulariumFrameCount(a.data);
      return new SO.Particle.Trajectory(traj, {
        label: a.label,
        description: `${frameCount} frame${frameCount !== 1 ? 's' : ''}`,
      });
    });
  },
  dispose({ cache }) {
    for (const asset of ((cache as any)?.assets as Asset.Wrapper[] | undefined) ?? []) asset.dispose();
  },
});

export const SimulariumParticlesProvider = DataFormatProvider({
  name: 'simularium_particles',
  label: 'Simularium Particles',
  description: 'Simularium Particles',
  category: ParticlesFormatCategory,
  binaryExtensions: ['simularium'],
  parse: async (plugin, data, params?: { label?: string; frameIndex?: number; loadGeometries?: boolean }) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseSimularium, void 0, { state: { isGhost: false } });

    const trajectory = format.apply(ParticleTrajectoryFromSimularium, {
      loadGeometries: params?.loadGeometries ?? true,
    });

    const list = trajectory.apply(ParticleListFromTrajectory, {
      frameIndex: params?.frameIndex ?? 0,
    });

    await format.commit({ revertOnError: true });

    return { format: format.selector, trajectory: trajectory.selector, list: list.selector };
  },
  visuals: (plugin: PluginContext, data: ParticleFormatData) => complexVisuals(plugin, data),
});
