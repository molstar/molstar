/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from '../color';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ThemeDataContext } from '../theme';
import { TableLegend } from '../../mol-util/legend';
import { defaults } from '../../mol-util';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives everything the same, uniform color.';

export const UniformColorThemeParams = {
    value: PD.Color(DefaultColor),
};
export type UniformColorThemeParams = typeof UniformColorThemeParams
export function getUniformColorThemeParams(ctx: ThemeDataContext) {
    return UniformColorThemeParams; // TODO return copy
}

export function UniformColorTheme(ctx: ThemeDataContext, props: PD.Values<UniformColorThemeParams>): ColorTheme<UniformColorThemeParams> {
    const color = defaults(props.value, DefaultColor);

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
    category: ColorTheme.Category.Misc,
    factory: UniformColorTheme,
    getParams: getUniformColorThemeParams,
    defaultValues: PD.getDefaultValues(UniformColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true
};