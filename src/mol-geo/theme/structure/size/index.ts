/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SizeType, LocationSize } from '../../../util/size-data';
import { SizeThemeProps } from '../..';
import { PhysicalSizeTheme } from './physical';
import { UniformSizeTheme } from './uniform';

export interface SizeTheme {
    kind: SizeType
    size: LocationSize
}

export function SizeTheme(props: SizeThemeProps): SizeTheme {
    switch (props.name) {
        case 'physical': return PhysicalSizeTheme(props)
        case 'uniform': return UniformSizeTheme(props)
    }
}