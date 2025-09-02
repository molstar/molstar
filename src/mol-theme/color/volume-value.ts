/**
 * Copyright (c) 2021-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorNames } from '../../mol-util/color/names';
import { Location } from '../../mol-model/location';
import { Volume } from '../../mol-model/volume/volume';
import { ColorThemeCategory } from './categories';
import { clamp, normalize } from '../../mol-math/interpolate';
import { Color } from '../../mol-util/color/color';
import { Grid } from '../../mol-model/volume/grid';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';

const Description = 'Assign color based on the given value of a volume cell.';

export const VolumeValueColorThemeParams = {
    colorList: PD.ColorList({
        kind: 'interpolate',
        colors: [
            [ColorNames.white, 0],
            [ColorNames.red, 0.25],
            [ColorNames.white, 0.5],
            [ColorNames.blue, 0.75],
            [ColorNames.white, 1]
        ]
    }, { offsets: true, isEssential: true }),
    domain: PD.MappedStatic('auto', {
        custom: PD.Interval([-1, 1], { step: 0.001 }),
        auto: PD.Group({
            symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
        })
    }),
    isRelative: PD.Boolean(false, { description: 'If true the value is treated as relative to the volume mean and sigma.' }),
    defaultColor: PD.Color(Color(0xcccccc)),
};
export type VolumeValueColorThemeParams = typeof VolumeValueColorThemeParams
export function getVolumeValueColorThemeParams(ctx: ThemeDataContext) {
    return VolumeValueColorThemeParams; // TODO return copy
}

export function VolumeValueColorTheme(ctx: ThemeDataContext, props: PD.Values<VolumeValueColorThemeParams>): ColorTheme<VolumeValueColorThemeParams, any> {
    if (ctx.volume) {
        const { min, max, mean, sigma } = ctx.volume.grid.stats;
        const domain: [number, number] = props.domain.name === 'custom' ? props.domain.params : [min, max];
        const { colorList, defaultColor } = props;

        if (props.domain.name === 'auto' && props.isRelative) {
            domain[0] = (domain[0] - mean) / sigma;
            domain[1] = (domain[1] - mean) / sigma;
        }

        if (props.domain.name === 'auto' && props.domain.params.symmetric) {
            const max = Math.max(Math.abs(domain[0]), Math.abs(domain[1]));
            domain[0] = -max;
            domain[1] = max;
        }

        if (ctx.locationKinds?.includes('direct-location')) {
            // this is for direct-volume rendering where the location is the volume grid cell
            // and we only need to provide the domain and palette here

            const normalizedDomain = [
                normalize(domain[0], min, max),
                normalize(domain[1], min, max)
            ] as [number, number];

            const palette = ColorTheme.Palette(colorList.colors, colorList.kind, normalizedDomain, defaultColor);

            return {
                factory: VolumeValueColorTheme as any,
                granularity: 'direct',
                props,
                description: Description,
                palette,
            };
        } else {
            const getTrilinearlyInterpolated = Grid.makeGetTrilinearlyInterpolated(ctx.volume.grid, 'none');

            const color = (location: Location): Color => {
                if (!isPositionLocation(location)) {
                    return props.defaultColor;
                }

                const value = getTrilinearlyInterpolated(location.position);
                if (isNaN(value)) return props.defaultColor;

                return (clamp((value - domain[0]) / (domain[1] - domain[0]), 0, 1) * ColorTheme.PaletteScale) as Color;
            };

            const palette = ColorTheme.Palette(colorList.colors, colorList.kind, undefined, defaultColor);

            return {
                factory: VolumeValueColorTheme as any,
                granularity: 'vertex',
                preferSmoothing: true,
                color,
                palette,
                props,
                description: Description,
            };
        }
    } else {
        return {
            factory: VolumeValueColorTheme as any,
            granularity: 'uniform',
            color: () => props.defaultColor,
            props,
            description: Description,
        };
    }
}

export const VolumeValueColorThemeProvider: ColorTheme.Provider<VolumeValueColorThemeParams, 'volume-value'> = {
    name: 'volume-value',
    label: 'Volume Value',
    category: ColorThemeCategory.Misc,
    factory: VolumeValueColorTheme,
    getParams: getVolumeValueColorThemeParams,
    defaultValues: PD.getDefaultValues(VolumeValueColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !Volume.Segmentation.get(ctx.volume),
};