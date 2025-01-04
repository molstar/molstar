/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '../../mol-util/color';
import { Location } from '../../mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { Grid, Volume } from '../../mol-model/volume';
import { isPositionLocation } from '../../mol-geo/util/location-iterator';
import { clamp } from '../../mol-math/interpolate';
import { ColorThemeCategory } from './categories';

const Description = `Assigns a color based on volume value at a given vertex.`;

export const VolumeDataColorThemeParams = {
    coloring: PD.MappedStatic('absolute-value', {
        'absolute-value': PD.Group({
            domain: PD.MappedStatic('auto', {
                custom: PD.Interval([-1, 1], { step: 0.001 }),
                auto: PD.Group({
                    symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                })
            }),
            list: PD.ColorList('red-white-blue', { presetKind: 'scale' })
        }),
        'relative-value': PD.Group({
            domain: PD.MappedStatic('auto', {
                custom: PD.Interval([-1, 1], { step: 0.001 }),
                auto: PD.Group({
                    symmetric: PD.Boolean(false, { description: 'If true the automatic range is determined as [-|max|, |max|].' })
                })
            }),
            list: PD.ColorList('red-white-blue', { presetKind: 'scale' })
        })
    }),
    defaultColor: PD.Color(Color(0xcccccc)),
};
export type VolumeDataColorThemeParams = typeof VolumeDataColorThemeParams

export function VolumeDataColorTheme(ctx: ThemeDataContext, props: PD.Values<VolumeDataColorThemeParams>): ColorTheme<VolumeDataColorThemeParams> {
    let color: LocationColor;
    let palette: ColorTheme.Palette | undefined;

    if (ctx.volume) {
        const coloring = props.coloring.params;
        const { stats } = ctx.volume.grid;
        const domain: [number, number] = coloring.domain.name === 'custom' ? coloring.domain.params : [stats.min, stats.max];

        const isRelative = props.coloring.name === 'relative-value';
        if (coloring.domain.name === 'auto' && isRelative) {
            domain[0] = (domain[0] - stats.mean) / stats.sigma;
            domain[1] = (domain[1] - stats.mean) / stats.sigma;
        }

        if (coloring.domain.name === 'auto' && coloring.domain.params.symmetric) {
            const max = Math.max(Math.abs(domain[0]), Math.abs(domain[1]));
            domain[0] = -max;
            domain[1] = max;
        }

        const getTrilinearlyInterpolated = Grid.makeGetTrilinearlyInterpolated(ctx.volume.grid, isRelative ? 'relative' : 'none');

        color = (location: Location): Color => {
            if (!isPositionLocation(location)) {
                return props.defaultColor;
            }

            const value = getTrilinearlyInterpolated(location.position);
            if (isNaN(value)) return props.defaultColor;

            return (clamp((value - domain[0]) / (domain[1] - domain[0]), 0, 1) * ColorTheme.PaletteScale) as Color;
        };

        palette = ColorTheme.Palette(coloring.list.colors, coloring.list.kind);
    } else {
        color = () => props.defaultColor;
    }

    return {
        factory: VolumeDataColorTheme,
        granularity: 'vertex',
        preferSmoothing: true,
        color,
        palette,
        props,
        description: Description,
        // TODO: figure out how to do legend for this
    };
}

export const VolumeDataColorThemeProvider: ColorTheme.Provider<VolumeDataColorThemeParams, 'volume-data'> = {
    name: 'volume-data',
    label: 'Volume Data',
    category: ColorThemeCategory.Misc,
    factory: VolumeDataColorTheme,
    getParams: () => VolumeDataColorThemeParams,
    defaultValues: PD.getDefaultValues(VolumeDataColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume && !Volume.Segmentation.get(ctx.volume),
};