/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme, TableLegend } from '../color';
import { Color } from 'mol-util/color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from 'mol-theme/theme';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives everything the same, uniform color.'

export const UniformColorThemeParams = {
    value: PD.Color('Color Value', '', DefaultColor),
}
export function getUniformColorThemeParams(ctx: ThemeDataContext) {
    return UniformColorThemeParams // TODO return copy
}
export type UniformColorThemeProps = PD.Values<typeof UniformColorThemeParams>

export function UniformColorTheme(ctx: ThemeDataContext, props: UniformColorThemeProps): ColorTheme<UniformColorThemeProps> {
    const color = props.value || DefaultColor

    return {
        granularity: 'uniform',
        color: () => color,
        props: props,
        description: Description,
        legend: TableLegend([['uniform', color]])
    }
}

export const UniformColorThemeProvider: ColorTheme.Provider<typeof UniformColorThemeParams> = {
    factory: UniformColorTheme, params: getUniformColorThemeParams
}