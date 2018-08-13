/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorThemeProps } from '../..';
import { ColorTheme } from '.';

const DefaultColor = 0xCCCCCC;

export function UniformColorTheme(props: ColorThemeProps): ColorTheme {
    const color = props.value || DefaultColor

    return {
        kind: 'uniform',
        color: () => color
    }
}