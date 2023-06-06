/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ColorTheme } from '../color';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';
import { defaults } from '../../mol-util';
import { ColorThemeCategory } from './categories';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives everything the same, uniform color.';

export const UniformColorThemeParams = {
    value: PD.Color(DefaultColor),
    saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
};
export type UniformColorThemeParams = typeof UniformColorThemeParams
export function getUniformColorThemeParams(ctx: ThemeDataContext) {
    return UniformColorThemeParams; // TODO return copy
}

export function UniformColorTheme(ctx: ThemeDataContext, props: PD.Values<UniformColorThemeParams>): ColorTheme<UniformColorThemeParams> {
    let color = defaults(props.value, DefaultColor);
    color = Color.saturate(color, props.saturation);
    color = Color.lighten(color, props.lightness);

    return {
        factory: UniformColorTheme,
        granularity: 'uniform',
        color: () => color,
        props: props,
        description: Description,
        legend: TableLegend([['uniform', color]])
    };
}

export const UniformColorThemeProvider: ColorTheme.Provider<UniformColorThemeParams, 'uniform'> = {
    name: 'uniform',
    label: 'Uniform',
    category: ColorThemeCategory.Misc,
    factory: UniformColorTheme,
    getParams: getUniformColorThemeParams,
    defaultValues: PD.getDefaultValues(UniformColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true
};