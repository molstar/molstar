/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateObjectRef } from '@molstar/core/state';
import type { PluginStateObject } from '@molstar/plugin/state/objects';
import { Particle, type ParticleList, getParticleTargetGroups } from '@molstar/model/model/particles/particle-list';
import { ParticleNesting } from '@molstar/model/model/particles/particle-nesting';
import type { PluginContext } from '@molstar/plugin/context';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { ParticleListUnitcell3D } from '@molstar/plugin/state/transforms/particles/unitcell';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { ParticleSpacefill } from '@molstar/plugin/registry/particles/spacefill';
import { ParticleFibers } from '@molstar/plugin/registry/particles/fibers';
import { ParticleTarget } from '@molstar/plugin/registry/particles/target';
import { SpacefillParticlesRepresentationProvider } from '@molstar/graphics/repr/particles/representation/spacefill';
import { FibersRepresentationProvider } from '@molstar/graphics/repr/particles/representation/fibers';
import { ParticleTargetRepresentationProvider } from '@molstar/graphics/repr/particles/representation/target/representation';
import { ParticleIndexColorThemeProvider } from '@molstar/graphics/theme/color/particle-index';
import { ParticleEntityColorThemeProvider } from '@molstar/graphics/theme/color/particle-entity';
import { ParticleCompartmentColorThemeProvider } from '@molstar/graphics/theme/color/particle-compartment';
import { ParticleHierarchyColorThemeProvider } from '@molstar/graphics/theme/color/particle-hierarchy';

export interface ParticleFormatData {
  format: StateObjectRef;
  list: StateObjectRef<PluginStateObject.Particle.List>;
}

/**
 * Whether the particle list has target objects and whether some particles are without one. A list with templates by
 * its leaves: its own particles (a template's copy drawn by the target representation) and its templates' (the
 * untargeted drawn flattened; not those of a template read by range, which flattening refuses)
 */
function getParticleTargetCoverage(particles: ParticleList | undefined) {
  if (!particles) return { hasTargets: false, hasUntargeted: false };
  let hasTargets = false;
  let hasUntargeted = false;
  const visit = (list: ParticleList, drawsUntargeted = true) => {
    const mapping = list.targetMapping;
    if (!mapping || mapping.size === 0) {
      if (list.count > 0 && drawsUntargeted) hasUntargeted = true;
      return;
    }
    const { targetIds } = getParticleTargetGroups(list);
    for (let i = 0; i < targetIds.length; ++i) {
      if (mapping.has(targetIds[i])) hasTargets = true;
      else if (drawsUntargeted) hasUntargeted = true;
    }
  };
  visit(particles);
  if (ParticleNesting.hasTemplates(particles)) {
    for (const template of ParticleNesting.getIdSpace(particles).moves.keys()) {
      visit(template, !Particle.getDeferred(template));
    }
  }
  return { hasTargets, hasUntargeted };
}

/**
 * What `complexVisuals` draws of `particleList`: the target representation where it has targets, spacefill where some
 * particles have none, fibers with them where it has fibers (not kept for a list with templates:
 * `ParticleNesting.flatten`); colored by what infos it has (a list with templates by its leaves',
 * `ParticleNesting.getThemeList`)
 */
export function getComplexVisualsPlan(particleList: ParticleList | undefined) {
  const { hasTargets, hasUntargeted } = getParticleTargetCoverage(particleList);
  const nested = !!particleList && ParticleNesting.hasTemplates(particleList);
  const themed = particleList && nested ? ParticleNesting.getThemeList(particleList) : particleList;
  const hasEntities = !!themed?.entityInfo;
  const hasCompartments = !!themed?.compartmentInfo;
  const colorTheme =
    hasEntities && hasCompartments
      ? 'particle-hierarchy'
      : hasEntities
        ? 'particle-entity'
        : hasCompartments
          ? 'particle-compartment'
          : 'particle-index';
  const hasFibers = hasUntargeted && !nested && !!particleList?.fibers && particleList.fibers.count > 0;
  return { hasTargets, hasUntargeted, hasFibers, colorTheme };
}

export function complexVisuals(
  plugin: PluginContext,
  data: ParticleFormatData,
  targetProps?: { params?: {}; colorTheme?: { name: string; params: {} }; sizeTheme?: { name: string; params: {} } },
) {
  const builder = plugin.state.data.build();
  const particleList = StateObjectRef.resolveAndCheck(plugin.state.data, data.list)?.obj?.data;
  const { hasTargets, hasUntargeted, hasFibers, colorTheme } = getComplexVisualsPlan(particleList);

  if (hasTargets) {
    builder.to(data.list).apply(ParticlesRepresentation3D, {
      type: { name: ParticleTargetRepresentationProvider.name, params: targetProps?.params ?? {} },
      colorTheme: targetProps?.colorTheme ?? { name: colorTheme, params: {} },
    });
  }

  if (hasUntargeted) {
    builder.to(data.list).apply(ParticlesRepresentation3D, {
      type: { name: SpacefillParticlesRepresentationProvider.name, params: { excludeTargets: hasTargets } },
      colorTheme: targetProps?.colorTheme ?? { name: colorTheme, params: {} },
      ...(targetProps?.sizeTheme && { sizeTheme: targetProps.sizeTheme }),
    });

    if (hasFibers) {
      builder.to(data.list).apply(ParticlesRepresentation3D, {
        type: { name: FibersRepresentationProvider.name, params: {} },
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

  builder
    .to(data.list)
    .apply(ParticlesRepresentation3D, { type: { name: SpacefillParticlesRepresentationProvider.name, params: {} } });

  builder.to(data.list).apply(ParticleListUnitcell3D, { attachment: 'center' });

  return builder.commit();
}

/** The particle representations and themes used by `simpleVisuals`. */
export const SimpleParticleVisuals: PluginRegistryEntry = ParticleSpacefill;

const ComplexEntries = [ParticleSpacefill, ParticleFibers, ParticleTarget];

/**
 * The particle representations and themes used by `complexVisuals`: spacefill, fibers, and target with their default
 * themes, and the color themes chosen by what the particle list has (hierarchy, entity, compartment, or particle index).
 */
export const ComplexParticleVisuals: PluginRegistryEntry = {
  particles: {
    representations: ComplexEntries.flatMap((e) => e.particles!.representations!),
    themes: {
      color: [
        ...new Set([
          ...ComplexEntries.flatMap((e) => e.particles!.themes!.color!),
          ParticleIndexColorThemeProvider,
          ParticleEntityColorThemeProvider,
          ParticleCompartmentColorThemeProvider,
          ParticleHierarchyColorThemeProvider,
        ]),
      ],
      size: [...new Set(ComplexEntries.flatMap((e) => e.particles!.themes!.size!))],
    },
  },
};
