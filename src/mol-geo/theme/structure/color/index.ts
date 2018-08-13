/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorThemeProps } from '../..';

import { ElementIndexColorTheme } from './element-index';
import { CarbohydrateSymbolColorTheme } from './carbohydrate-symbol';
import { ChainIdColorTheme } from './chain-id';
import { ElementSymbolColorTheme } from './element-symbol';
import { UnitIndexColorTheme } from './unit-index';
import { UniformColorTheme } from './uniform';
import { ColorType, LocationColor } from '../../../util/color-data';

export interface ColorTheme {
    kind: ColorType
    color: LocationColor
}

export function ColorTheme(props: ColorThemeProps): ColorTheme {
    switch (props.name) {
        case 'element-index': return ElementIndexColorTheme(props)
        case 'carbohydrate-symbol': return CarbohydrateSymbolColorTheme(props)
        case 'chain-id': return ChainIdColorTheme(props)
        case 'element-symbol': return ElementSymbolColorTheme(props)
        case 'unit-index': return UnitIndexColorTheme(props)
        case 'uniform': return UniformColorTheme(props)
    }
}