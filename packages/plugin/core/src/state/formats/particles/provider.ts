/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateObjectRef } from '@molstar/core/state';
import type { PluginStateObject } from '@molstar/plugin/state/objects';
import { type ParticleList, getParticleTargetGroups } from '@molstar/model/model/particles/particle-list';
import type { PluginContext } from '@molstar/plugin/context';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticleListUnitcell3D } from '@molstar/plugin/state/transforms/particles/unitcell';

export interface ParticleFormatData {
  format: StateObjectRef;
  list: StateObjectRef<PluginStateObject.Particle.List>;
}

/** Whether the particle list has target objects and whether some particles are without one. */
function getParticleTargetCoverage(particles: ParticleList | undefined) {
  const mapping = particles?.targetMapping;
  if (!particles || !mapping || mapping.size === 0) {
    return { hasTargets: false, hasUntargeted: !!particles && particles.count > 0 };
  }

  const { targetIds } = getParticleTargetGroups(particles);
  let hasTargets = false;
  let hasUntargeted = false;
  for (let i = 0; i < targetIds.length; ++i) {
    if (mapping.has(targetIds[i])) hasTargets = true;
    else hasUntargeted = true;
  }
  return { hasTargets, hasUntargeted };
}

export function complexVisuals(
  plugin: PluginContext,
  data: ParticleFormatData,
  targetProps?: { params?: {}; colorTheme?: { name: string; params: {} }; sizeTheme?: { name: string; params: {} } },
) {
  const builder = plugin.state.data.build();
  const particleList = StateObjectRef.resolveAndCheck(plugin.state.data, data.list)?.obj?.data;

  const { hasTargets, hasUntargeted } = getParticleTargetCoverage(particleList);
  const hasEntities = !!particleList?.entityInfo;
  const hasCompartments = !!particleList?.compartmentInfo;

  const colorTheme =
    hasEntities && hasCompartments
      ? 'particle-hierarchy'
      : hasEntities
        ? 'particle-entity'
        : hasCompartments
          ? 'particle-compartment'
          : 'particle-index';

  if (hasTargets) {
    builder.to(data.list).apply(ParticlesRepresentation3D, {
      type: { name: 'target', params: targetProps?.params ?? {} },
      colorTheme: targetProps?.colorTheme ?? { name: colorTheme, params: {} },
    });
  }

  if (hasUntargeted) {
    builder.to(data.list).apply(ParticlesRepresentation3D, {
      type: { name: 'spacefill', params: { excludeTargets: hasTargets } },
      colorTheme: targetProps?.colorTheme ?? { name: colorTheme, params: {} },
      ...(targetProps?.sizeTheme && { sizeTheme: targetProps.sizeTheme }),
    });

    if (particleList?.fibers && particleList.fibers.count > 0) {
      builder.to(data.list).apply(ParticlesRepresentation3D, {
        type: { name: 'fibers', params: {} },
        colorTheme: targetProps?.colorTheme ?? { name: colorTheme, params: {} },
        ...(targetProps?.sizeTheme && { sizeTheme: targetProps.sizeTheme }),
      });
    }
  }

  builder.to(data.list).apply(ParticleListUnitcell3D, { attachment: 'center' });

  return builder.commit();
}

export function simpleVisuals(plugin: PluginContext, data: ParticleFormatData) {
  const builder = plugin.state.data.build();

  builder.to(data.list).apply(ParticlesRepresentation3D, { type: { name: 'spacefill', params: {} } });

  builder.to(data.list).apply(ParticleListUnitcell3D, { attachment: 'center' });

  return builder.commit();
}
