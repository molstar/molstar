/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { ParticleSpacefill } from '@molstar/plugin/registry/particles/spacefill';
import { ParticleOrientation } from '@molstar/plugin/registry/particles/orientation';
import { ParticleFibers } from '@molstar/plugin/registry/particles/fibers';
import { ParticleTarget } from '@molstar/plugin/registry/particles/target';
import { BuiltInParticlesFormats } from '@molstar/plugin/state/formats/particles/catalog';
import { ParticleListFromRelionStar } from '@molstar/plugin/state/formats/particles/star';
import { ParticleListFromDynamoTbl } from '@molstar/plugin/state/formats/particles/tbl';
import { ParticleListFromCryoEtDataPortalNdjson } from '@molstar/plugin/state/formats/particles/ndjson';
import { ParticleListFromArtiatomiEm } from '@molstar/plugin/state/formats/particles/em';
import { ParticleListFromMmcifAssembly } from '@molstar/plugin/state/formats/particles/mmcif-assembly';
import { ParticleTrajectoryFromSimularium } from '@molstar/plugin/state/formats/particles/simularium';
import { ParticleListFromTrajectory, ParticleListWithTargets } from '@molstar/plugin/state/transforms/particles/ops';
import { ParticleListUnitcell3D } from '@molstar/plugin/state/transforms/particles/unitcell';
import { ParticlesRepresentation3D } from '@molstar/plugin/state/transforms/particles/representation';
import { AnimateParticleTrajectory } from '@molstar/plugin/state/animation/built-in/particles';
import { ExternalColorThemes } from '@molstar/plugin/themes/external';

/*
 * Particles (cryo-ET picks, subtomogram placements, mesoscale assemblies, Simularium trajectories) as one opt-in
 * entries: `DefaultRegistry` holds none of it, an application composes it (`[...DefaultRegistry, ...ParticlesRegistry]`).
 * Importing this module registers the particle transformers, which a session holding particles needs.
 */

/** The color and size themes made for particles: registered only with `ParticlesRegistry` */
export const ParticleColorThemes = [
  BuiltInColorThemes['particle-attribute'],
  BuiltInColorThemes['particle-compartment'],
  BuiltInColorThemes['particle-entity'],
  BuiltInColorThemes['particle-hierarchy'],
  BuiltInColorThemes['particle-index'],
] as readonly ColorTheme.Provider<any, any>[];
export const ParticleSizeThemes = [BuiltInSizeThemes['particle-size']] as readonly SizeTheme.Provider<any, any>[];

/** The particles scope's themes: every built-in theme (the general ones apply to particles too) and the external ones */
export const ParticleThemes: PluginRegistryEntry = {
  particles: {
    themes: {
      color: [
        ...Object.values(BuiltInColorThemes),
        ...(ExternalColorThemes.particles?.themes?.color ?? []),
      ] as readonly ColorTheme.Provider<any, any>[],
      size: Object.values(BuiltInSizeThemes) as readonly SizeTheme.Provider<any, any>[],
    },
  },
};

export const ParticleRepresentations: PluginRegistryEntry = {
  particles: {
    representations: [ParticleSpacefill, ParticleOrientation, ParticleFibers, ParticleTarget].flatMap(
      (e) => e.particles!.representations!,
    ),
  },
};

export const ParticleActions: PluginRegistryEntry = {
  actions: [
    ParticleListFromRelionStar,
    ParticleListFromDynamoTbl,
    ParticleListFromCryoEtDataPortalNdjson,
    ParticleListFromArtiatomiEm,
    ParticleListFromMmcifAssembly,
    ParticleTrajectoryFromSimularium,
    ParticleListFromTrajectory,
    ParticleListWithTargets,
    ParticleListUnitcell3D,
    ParticlesRepresentation3D,
  ],
};

export const ParticleFormats: PluginRegistryEntry = {
  formats: [...BuiltInParticlesFormats],
};

export const ParticleAnimations: PluginRegistryEntry = {
  animations: [AnimateParticleTrajectory],
};

/** Everything for particles, in registration order (themes, representations, actions, formats, animations) */
export const ParticlesRegistry: readonly PluginRegistryEntry[] = [
  ParticleThemes,
  ParticleRepresentations,
  ParticleActions,
  ParticleFormats,
  ParticleAnimations,
];
