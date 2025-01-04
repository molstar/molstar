/**
 * Copyright (c) 2021-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorNames } from '../../mol-util/color/names';
import { ColorTypeDirect } from '../../mol-geo/geometry/color-data';
import { Volume } from '../../mol-model/volume/volume';
import { ColorThemeCategory } from './categories';
import { normalize } from '../../mol-math/interpolate';

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
};
export type VolumeValueColorThemeParams = typeof VolumeValueColorThemeParams
export function getVolumeValueColorThemeParams(ctx: ThemeDataContext) {
    return VolumeValueColorThemeParams; // TODO return copy
}

export function VolumeValueColorTheme(ctx: ThemeDataContext, props: PD.Values<VolumeValueColorThemeParams>): ColorTheme<VolumeValueColorThemeParams, ColorTypeDirect> {
    let palette: ColorTheme.Palette | undefined;

    if (ctx.volume) {
        const { min, max } = ctx.volume.grid.stats;
        const domain: [number, number] = props.domain.name === 'custom' ? props.domain.params : [min, max];
        const { colorList } = props;

        if (props.domain.name === 'auto' && props.domain.params.symmetric) {
            const max = Math.max(Math.abs(domain[0]), Math.abs(domain[1]));
            domain[0] = -max;
            domain[1] = max;
        }

        const normalizedDomain = [
            normalize(domain[0], min, max),
            normalize(domain[1], min, max)
        ] as [number, number];

        palette = ColorTheme.Palette(colorList.colors, colorList.kind, normalizedDomain);
    }

    return {
        factory: VolumeValueColorTheme,
        granularity: 'direct',
        props,
        description: Description,
        palette,
    };
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