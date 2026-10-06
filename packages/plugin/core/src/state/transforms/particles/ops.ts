/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Task } from '@molstar/core/task';
import type { PluginContext } from '@molstar/plugin/context';
import {
  particleTargetFormatOptions,
  particleTargetFileExtensions,
  type ParticleTargetSource,
  matchParticleTargetFiles,
  loadParticleTarget,
  loadParticleTargets,
} from '@molstar/plugin/state/helpers/particle-targets';
import type { Asset } from '@molstar/core/util/assets';
import { type ParticleList, type ParticleTarget, Particle } from '@molstar/model/model/particles/particle-list';
import { StateTransformer } from '@molstar/core/state';

const plus1 = (v: number) => v + 1,
  minus1 = (v: number) => v - 1;

export { ParticleListFromTrajectory };
type ParticleListFromTrajectory = typeof ParticleListFromTrajectory;
const ParticleListFromTrajectory = PluginStateTransform.BuiltIn({
  name: 'particle-list-from-trajectory',
  display: { name: 'Particle List from Trajectory', description: 'Extract a single frame from a ParticleTrajectory.' },
  from: SO.Particle.Trajectory,
  to: SO.Particle.List,
  params: (a) => {
    if (!a) {
      return { frameIndex: PD.Numeric(0, {}, { description: 'Zero-based index of the frame', immediateUpdate: true }) };
    }
    return {
      frameIndex: PD.Converted(
        plus1,
        minus1,
        PD.Numeric(
          1,
          { min: 1, max: a.data.frameCount, step: 1 },
          { description: 'Frame Index', immediateUpdate: true },
        ),
      ),
    };
  },
})({
  isApplicable: (a) => a.data.frameCount > 0,
  apply({ a, params }) {
    return Task.create('Particle List from Trajectory', async (ctx) => {
      const idx = Math.max(0, Math.min(params.frameIndex, a.data.frameCount - 1));
      const list = await Task.resolveInContext(a.data.getFrameAtIndex(idx), ctx);
      return new SO.Particle.List(list, {
        label: list.label || 'Particles',
        description: `Frame ${params.frameIndex + 1} of ${a.data.frameCount}`,
      });
    });
  },
});

function ParticleListWithTargetsParams(a: SO.Particle.List | undefined, plugin: PluginContext) {
  const targetIds = a ? Array.from(a.data.targetInfo.keys()).sort((x, y) => x - y) : [];
  const targetIdDescription =
    targetIds.length > 0
      ? `Target ID in the particle list this object maps to. Available IDs: ${targetIds.slice(0, 32).join(', ')}${targetIds.length > 32 ? ', …' : ''}.`
      : 'Target ID in the particle list this object maps to (matches ParticleList.targets).';
  const formatOptions = particleTargetFormatOptions(plugin);

  return {
    files: PD.FileList({
      accept: particleTargetFileExtensions(plugin),
      multiple: true,
      description:
        'Reference objects matched to particle targets by exact filename (without extension) and entity name.',
    }),
    targets: PD.ObjectList(
      {
        targetId: PD.Numeric(0, { min: 0, step: 1 }, { description: targetIdDescription }),
        source: PD.MappedStatic(
          'file',
          {
            file: PD.Group(
              {
                file: PD.File({ accept: particleTargetFileExtensions(plugin) }),
                format: PD.Select('auto', formatOptions),
              },
              { isFlat: true },
            ),
            url: PD.Group(
              {
                url: PD.Url('', { label: 'URL' }),
                format: PD.Select('auto', formatOptions),
                isBinary: PD.Optional(
                  PD.Boolean(false, { description: 'Leave unset to derive from the file extension.' }),
                ),
              },
              { isFlat: true },
            ),
          },
          {
            options: [
              ['file', 'File'],
              ['url', 'URL'],
            ] as ['url' | 'file', string][],
          },
        ),
      },
      (e) => `Target ${e.targetId}`,
      {
        description: 'Reference objects instanced at each particle, mapped per particle target ID.',
      },
    ),
    includeParent: PD.Boolean(true, {
      description: 'Keep the reference objects the particle list already provides (e.g. built by its format).',
    }),
  };
}

interface ParticleListWithTargetsCache {
  assets?: Asset.Wrapper[];
  sourceTargetMapping?: ParticleList['targetMapping'];
  sourceEntityInfo?: ParticleList['entityInfo'];
}

export { ParticleListWithTargets };
type ParticleListWithTargets = typeof ParticleListWithTargets;
const ParticleListWithTargets = PluginStateTransform.BuiltIn({
  name: 'particle-list-with-targets',
  display: {
    name: 'Particle List with Targets',
    description:
      'Attach reference structures, shapes, or volumes that are instanced at each particle of the matching target ID.',
  },
  isDecorator: true,
  from: SO.Particle.List,
  to: SO.Particle.List,
  params: ParticleListWithTargetsParams,
})({
  apply({ a, params, cache }, plugin: PluginContext) {
    return Task.create('Particle List with Targets', async (ctx) => {
      const transformCache = cache as ParticleListWithTargetsCache;
      const sources: { targetId: number; source: ParticleTargetSource }[] = [];
      const explicitTargetIds = new Set<number>();
      for (const { targetId, source } of params.targets) {
        explicitTargetIds.add(targetId);
        if (source.name === 'url') {
          sources.push({
            targetId,
            source: {
              kind: 'url',
              url: source.params.url,
              format: source.params.format,
              isBinary: source.params.isBinary,
            },
          });
        } else if (source.params.file) {
          sources.push({ targetId, source: { kind: 'file', file: source.params.file, format: source.params.format } });
        }
      }

      const targetMapping = new Map<number, ParticleTarget>();
      if (params.includeParent && a.data.targetMapping) {
        for (const [targetId, target] of a.data.targetMapping) targetMapping.set(targetId, target);
      }
      const assets: Asset.Wrapper[] = [];
      transformCache.assets = assets;

      const fileMatches = matchParticleTargetFiles(a.data, params.files ?? [], explicitTargetIds);
      for (const warning of fileMatches.warnings) plugin.log.warn(warning);
      for (const { file, targetIds } of fileMatches.matches) {
        try {
          const target = await loadParticleTarget(plugin, ctx, { kind: 'file', file }, assets);
          for (const targetId of targetIds) targetMapping.set(targetId, target);
        } catch (e) {
          console.error(e);
          plugin.log.warn(`Could not load particle target file '${file.name}'.`);
        }
      }

      for (const { targetId, target } of await loadParticleTargets(plugin, ctx, sources, assets)) {
        targetMapping.set(targetId, target);
      }

      transformCache.sourceTargetMapping = a.data.targetMapping;
      transformCache.sourceEntityInfo = a.data.entityInfo;
      return new SO.Particle.List(Particle.withTargets(a.data, targetMapping), {
        label: a.label,
        description: a.description,
      });
    });
  },
  update({ a, b, oldParams, newParams, cache }, plugin: PluginContext) {
    const transformCache = cache as ParticleListWithTargetsCache;
    if (
      !PD.areEqual(ParticleListWithTargetsParams(a, plugin), oldParams, newParams) ||
      a.data.targetMapping !== transformCache.sourceTargetMapping ||
      a.data.entityInfo !== transformCache.sourceEntityInfo ||
      !b.data.targetMapping
    ) {
      return StateTransformer.UpdateResult.Recreate;
    }

    b.data = Particle.withTargets(a.data, b.data.targetMapping);
    b.label = a.label;
    b.description = a.description;
    return StateTransformer.UpdateResult.Updated;
  },
  dispose({ cache }) {
    for (const asset of (cache as ParticleListWithTargetsCache).assets ?? []) asset.dispose();
  },
});
