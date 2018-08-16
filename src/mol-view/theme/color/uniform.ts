/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme, ColorThemeProps } from '../color';
import { Color } from 'mol-util/color';

const DefaultColor = 0xCCCCCC as Color

export function UniformColorTheme(props: ColorThemeProps): ColorTheme {
    const color = props.value || DefaultColor

    return {
        kind: 'uniform',
        color: () => color
    }
}