/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme, ColorThemeProps, TableLegend } from '../color';
import { Color } from 'mol-util/color';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives everything the same, uniform color.'

export function UniformColorTheme(props: ColorThemeProps): ColorTheme {
    const color = props.value || DefaultColor

    return {
        granularity: 'uniform',
        color: () => color,
        description: Description,
        legend: TableLegend([['uniform', color]])
    }
}