/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getParticleSpacingRadius, Particle } from '@molstar/model/model/particles/particle-list';
import type { Location } from '@molstar/model/model/location';
import type { SizeTheme } from '../size.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';

const Description = `Assigns a size based on the per-particle bounding sphere radius stored in the particle list. Without radius data, half the typical distance between neighboring particles (or the default size, when set).`;

export const ParticleSizeThemeParams = {
  scale: PD.Numeric(1.0, { min: 0, max: 100, step: 0.1 }),
  defaultSize: PD.Numeric(
    0,
    { min: 0, max: 1000, step: 0.1 },
    {
      description: 'The radius of particles without one; 0: half the typical distance between neighbors',
    },
  ),
};
export type ParticleSizeThemeParams = typeof ParticleSizeThemeParams;
export function getParticleSizeThemeParams(ctx: ThemeDataContext) {
  return ParticleSizeThemeParams;
}

export function ParticleSizeTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ParticleSizeThemeParams>,
): SizeTheme<ParticleSizeThemeParams> {
  // Without radii and a size set, the spacing's (picks of one kind about touching), measured once per list
  const fallback =
    props.defaultSize > 0 ? props.defaultSize : ctx.particles ? getParticleSpacingRadius(ctx.particles) : 1;
  function size(location: Location): number {
    let s = fallback;
    if (Particle.isLocation(location)) {
      const { particles, index } = location;
      s = particles.radii?.[index] ?? fallback;
    }
    return s * props.scale;
  }

  return {
    factory: ParticleSizeTheme,
    granularity: 'groupInstance',
    size,
    props,
    description: Description,
  };
}

export const ParticleSizeThemeProvider: SizeTheme.Provider<ParticleSizeThemeParams, 'particle-size'> = {
  name: 'particle-size',
  label: 'Particle Size',
  category: '',
  factory: ParticleSizeTheme,
  getParams: getParticleSizeThemeParams,
  defaultValues: PD.getDefaultValues(ParticleSizeThemeParams),
  isApplicable: (ctx: ThemeDataContext) => !!ctx.particles,
};
