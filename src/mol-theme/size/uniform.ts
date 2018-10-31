/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SizeTheme, SizeThemeProps } from '../size';

const DefaultSize = 1
const DefaultFactor = 1

export function UniformSizeTheme(props: SizeThemeProps): SizeTheme {
    const value = props.value || DefaultSize
    const factor = props.factor || DefaultFactor
    const size = value * factor

    return {
        granularity: 'uniform',
        size: () => size
    }
}