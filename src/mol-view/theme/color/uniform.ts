/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme, ColorThemeProps } from '../color';
import { Color } from 'mol-util/color';

const DefaultColor = Color(0xCCCCCC)

export function UniformColorTheme(props: ColorThemeProps): ColorTheme {
    const color = props.value || DefaultColor

    return {
        granularity: 'uniform',
        color: () => color
    }
}