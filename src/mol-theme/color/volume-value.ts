/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { Color, ColorScale } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { ColorNames } from '../../mol-util/color/names';

const DefaultColor = Color(0xCCCCCC);
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
};
export type VolumeValueColorThemeParams = typeof VolumeValueColorThemeParams
export function getVolumeValueColorThemeParams(ctx: ThemeDataContext) {
    return VolumeValueColorThemeParams; // TODO return copy
}

export function VolumeValueColorTheme(ctx: ThemeDataContext, props: PD.Values<VolumeValueColorThemeParams>): ColorTheme<VolumeValueColorThemeParams> {
    const scale = ColorScale.create({ domain: [0, 1], listOrName: props.colorList.colors });

    const colors: Color[] = [];
    for (let i = 0; i < 256; ++i) {
        colors[i] = scale.color(i / 255);
    }

    const palette: ColorTheme.Palette = { colors, filter: 'linear' };

    return {
        factory: VolumeValueColorTheme,
        granularity: 'direct',
        color: () => DefaultColor,
        props: props,
        description: Description,
        legend: scale.legend,
        palette,
    };
}

export const VolumeValueColorThemeProvider: ColorTheme.Provider<VolumeValueColorThemeParams, 'volume-value'> = {
    name: 'volume-value',
    label: 'Volume Value',
    category: ColorTheme.Category.Misc,
    factory: VolumeValueColorTheme,
    getParams: getVolumeValueColorThemeParams,
    defaultValues: PD.getDefaultValues(VolumeValueColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.volume,
};