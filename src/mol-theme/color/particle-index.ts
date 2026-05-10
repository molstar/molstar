/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../../mol-theme/theme';
import { getPaletteParams, getPalette } from '../../mol-util/color/palette';
import { TableLegend, ScaleLegend } from '../../mol-util/legend';
import { ColorLists, getColorListFromName } from '../../mol-util/color/lists';
import { ColorThemeCategory } from './categories';
import { Particle } from '../../mol-model/particles/particle-list';

const DefaultList = 'dark-2';
const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives every particle a unique color based on its index in the particle list.';

export const ParticleIndexColorThemeParams = {
    ...getPaletteParams({ type: 'colors', colorList: DefaultList }),
};
export type ParticleIndexColorThemeParams = typeof ParticleIndexColorThemeParams
export function getParticleIndexColorThemeParams(ctx: ThemeDataContext) {
    const params = PD.clone(ParticleIndexColorThemeParams);
    if (ctx.particles) {
        const { count } = ctx.particles;
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

    if (ctx.particles) {
        const { count } = ctx.particles;
        const palette = getPalette(count, props);
        legend = palette.legend;

        color = (location: Location): Color => {
            if (Particle.isLocation(location)) {
                return palette.color(location.index);
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
    category: ColorThemeCategory.Misc,
    factory: ParticleIndexColorTheme,
    getParams: getParticleIndexColorThemeParams,
    defaultValues: PD.getDefaultValues(ParticleIndexColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.particles
};
