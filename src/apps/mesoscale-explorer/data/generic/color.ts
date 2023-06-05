/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../../../../mol-util/color';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { ThemeDataContext } from '../../../../mol-theme/theme';
import { TableLegend } from '../../../../mol-util/legend';
import { defaults } from '../../../../mol-util';
import { ColorThemeCategory } from '../../../../mol-theme/color/categories';
import { ColorTheme } from '../../../../mol-theme/color';

const DefaultColor = Color(0xCCCCCC);
const Description = 'Gives everything the same, uniform color.';

export const GenericUniformColorThemeParams = {
    value: PD.Color(DefaultColor),
    saturation: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
    lightness: PD.Numeric(0, { min: -6, max: 6, step: 0.1 }),
};
export type GenericUniformColorThemeParams = typeof GenericUniformColorThemeParams
export function getGenericUniformColorThemeParams(ctx: ThemeDataContext) {
    return GenericUniformColorThemeParams; // TODO return copy
}

export function GenericUniformColorTheme(ctx: ThemeDataContext, props: PD.Values<GenericUniformColorThemeParams>): ColorTheme<GenericUniformColorThemeParams> {
    let color = defaults(props.value, DefaultColor);
    color = Color.saturate(color, props.saturation);
    color = Color.lighten(color, props.lightness);

    return {
        factory: GenericUniformColorTheme,
        granularity: 'uniform',
        color: () => color,
        props: props,
        description: Description,
        legend: TableLegend([['uniform', color]])
    };
}

export const GenericUniformColorThemeProvider: ColorTheme.Provider<GenericUniformColorThemeParams, 'generic-uniform'> = {
    name: 'generic-uniform',
    label: 'Generic Uniform',
    category: ColorThemeCategory.Misc,
    factory: GenericUniformColorTheme,
    getParams: getGenericUniformColorThemeParams,
    defaultValues: PD.getDefaultValues(GenericUniformColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => true
};