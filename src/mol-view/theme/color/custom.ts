/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { ColorThemeProps, ColorTheme } from '../color';
import { defaults } from 'mol-util';

const DefaultColor = Color(0xCCCCCC)

export function CustomColorTheme(props: ColorThemeProps): ColorTheme {
    const value = defaults(props.value, DefaultColor)
    return {
        granularity: defaults(props.granularity, 'uniform'),
        color: defaults(props.color, () => value)
    }
}