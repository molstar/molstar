/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import { getPaletteParams, getPalette } from '@molstar/core/util/color/palette';
import type { TableLegend, ScaleLegend } from '@molstar/core/util/legend';
import { ColorLists, getColorListFromName } from '@molstar/core/util/color/lists';
import { ColorThemeCategory } from './categories.js';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every particle a unique color based on its index in the particle list.';

export const ParticleIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleIndexColorThemeParams = typeof ParticleIndexColorThemeParams

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
    if (ctx.particles) return ctx.particles;
    return undefined;
}

export function getParticleIndexColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ParticleIndexColorThemeParams);
    const particles = getParticleList(ctx);
    if (particles) {
        const { count } = particles;
        if (count > ColorLists[DefaultList].list.length) {
            params.palette.defaultValue.name = 'colors';
            params.palette.defaultValue.params = {
                ...params.palette.defaultValue.params,
                list: { kind: 'interpolate', colors: getColorListFromName(DefaultList).list }
            };
        }
    }
    return params;
}

export function ParticleIndexColorTheme(ctx: ThemeDataContext, props: PD.Values<ParticleIndexColorThemeParams>): ColorTheme<ParticleIndexColorThemeParams> {
    let color: LocationColor;
    let legend: ScaleLegend | TableLegend | undefined;

    const particles = getParticleList(ctx);
    if (particles) {
        const { count } = particles;
        const palette = getPalette(count, props);
        legend = palette.legend;

        const pick = (index: number) => index >= 0 && index < count ? palette.color(index) : DefaultColor;

        color = (location: Location): Color => {
            if (Particle.isLocation(location)) {
                return pick(location.index);
            }
            return DefaultColor;
        };
    } else {
        color = () => DefaultColor;
    }

    return {
        factory: ParticleIndexColorTheme,
        granularity: 'instance',
        color,
        props,
        description: Description,
        legend
    };
}

export const ParticleIndexColorThemeProvider: ColorTheme.Provider<ParticleIndexColorThemeParams, 'particle-index'> = {
    name: 'particle-index',
    label: 'Particle Index',
    category: ColorThemeCategory.Particle,
    factory: ParticleIndexColorTheme,
    getParams: getParticleIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!getParticleList(ctx)
};
