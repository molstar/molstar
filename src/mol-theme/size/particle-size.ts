/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Particle } from '../../mol-model/particles/particle-list';
import { Location } from '../../mol-model/location';
import type { SizeTheme } from '../size';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';

const Description = `Assigns a size based on the per-particle bounding sphere radius stored in the particle list. Falls back to a uniform default size when no radius data is available.`;

export const ParticleSizeThemeParams = {
    scale: PD.Numeric(1.0, { min: 0, max: 100, step: 0.1 }),
    defaultSize: PD.Numeric(1.0, { min: 0, max: 1000, step: 0.1 }),
};
export type ParticleSizeThemeParams = typeof ParticleSizeThemeParams
export function getParticleSizeThemeParams(ctx: ThemeDataContext) {
    return ParticleSizeThemeParams;
}

export function ParticleSizeTheme(ctx: ThemeDataContext, props: PD.Values<ParticleSizeThemeParams>): SizeTheme<ParticleSizeThemeParams> {
    function size(location: Location): number {
        let s = props.defaultSize;
        if (Particle.isLocation(location)) {
            const { particles, index } = location;
            s = particles.radii?.[index] || props.defaultSize;
        }
        return s * props.scale;
    }

    return {
        factory: ParticleSizeTheme,
        granularity: 'groupInstance',
        size,
        props,
        description: Description
    };
}

export const ParticleSizeThemeProvider: SizeTheme.Provider<ParticleSizeThemeParams, 'particle-size'> = {
    name: 'particle-size',
    label: 'Particle Radius',
    category: '',
    factory: ParticleSizeTheme,
    getParams: getParticleSizeThemeParams,
    defaultValues: PD.getDefaultValues(ParticleSizeThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.particles,
};
