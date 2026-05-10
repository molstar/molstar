/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import type { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorThemeCategory } from './categories';
import { Volume } from '../../mol-model/volume/volume';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { Grid } from '../../mol-model/volume/grid';
import { ScaleLegend } from '../../mol-util/legend';
import { clamp, normalize } from '../../mol-math/interpolate';

const BlackColor = Color(0x000000);
const WhiteColor = Color(0xFFFFFF);
const Description = 'Assign grayscale color based on normalized volume values using IMOD-style black/white levels, gamma, and invert.';

export const VolumeWindowingColorThemeParams = {
    levels: PD.Interval([0, 255], { min: 0, max: 255, step: 1 }, { label: 'Black / White', description: 'Display range mapped from black to white, in IMOD-style 8-bit display levels.' }),
    gamma: PD.Numeric(1, { min: 0.1, max: 5, step: 0.01 }, { description: 'Gamma applied after scalar black/white mapping.' }),
    invert: PD.Boolean(false, { description: 'Invert the displayed grayscale after scalar mapping.' }),
};
export type VolumeWindowingColorThemeParams = typeof VolumeWindowingColorThemeParams

export function getVolumeWindowingColorThemeParams(ctx: ThemeDataContext) {
    return PD.clone(VolumeWindowingColorThemeParams);
}

export function applyVolumeWindowing(value: number, levels: [number, number], gamma: number, invert: boolean) {
    const black = Math.min(levels[0], levels[1]) / 255;
    const white = Math.max(levels[0], levels[1]) / 255;
    const intensity = clamp((clamp(value, 0, 1) - black) / Math.max(white - black, 0.0001), 0, 1);
    const corrected = Math.pow(intensity, 1 / Math.max(gamma, 0.0001));
    return invert ? 1 - corrected : corrected;
}

function grayscaleColor(value: number) {
    const channel = Math.round(clamp(value, 0, 1) * 255);
    return Color.fromRgb(channel, channel, channel);
}

export function getVolumeWindowingPalette(levels: [number, number], gamma: number, invert: boolean) {
    const colors = new Array<Color>(256);
    for (let i = 0; i < 256; ++i) {
        colors[i] = grayscaleColor(applyVolumeWindowing(i / 255, levels, gamma, invert));
    }
    return colors;
}

export function VolumeWindowingColorTheme(ctx: ThemeDataContext, props: PD.Values<VolumeWindowingColorThemeParams>): ColorTheme<VolumeWindowingColorThemeParams, any> {
    const legend = ScaleLegend(
        props.levels[0].toString(),
        props.levels[1].toString(),
        props.invert ? [[WhiteColor, 0], [BlackColor, 1]] : [[BlackColor, 0], [WhiteColor, 1]]
    );

    if (ctx.volume) {
        const { min, max } = ctx.volume.grid.stats;
        const palette = {
            colors: getVolumeWindowingPalette(props.levels, props.gamma, props.invert),
            filter: 'linear' as const,
            domain: [0, 1] as [number, number],
            defaultColor: BlackColor,
        };

        if (ctx.locationKinds?.includes('direct-location')) {
            return {
                factory: VolumeWindowingColorTheme as any,
                granularity: 'direct',
                palette,
                props,
                description: Description,
                legend,
            };
        }

        const getTrilinearlyInterpolated = Grid.makeGetTrilinearlyInterpolated(ctx.volume.grid, 'none');

        const color: LocationColor = (location: Location): Color => {
            if (!isPositionLocation(location)) return BlackColor;

            const rawValue = getTrilinearlyInterpolated(location.position);
            if (Number.isNaN(rawValue)) return BlackColor;

            const normalizedValue = normalize(rawValue, min, max);
            return grayscaleColor(applyVolumeWindowing(normalizedValue, props.levels, props.gamma, props.invert));
        };

        return {
            factory: VolumeWindowingColorTheme as any,
            granularity: 'vertex',
            preferSmoothing: true,
            color,
            props,
            description: Description,
            legend,
        };
    }

    return {
        factory: VolumeWindowingColorTheme as any,
        granularity: 'uniform',
        color: () => BlackColor,
        props,
        description: Description,
        legend,
    };
}

export const VolumeWindowingColorThemeProvider: ColorTheme.Provider<VolumeWindowingColorThemeParams, 'volume-windowing'> = {
    name: 'volume-windowing',
    label: 'Volume Windowing',
    category: ColorThemeCategory.Misc,
    factory: VolumeWindowingColorTheme,
    getParams: getVolumeWindowingColorThemeParams,
    defaultValues: PD.getDefaultValues(VolumeWindowingColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !Volume.Segmentation.get(ctx.volume),
};
