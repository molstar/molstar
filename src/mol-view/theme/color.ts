/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Structure } from 'mol-model/structure';
import { ColorType, LocationColor } from 'mol-geo/util/color-data';

import { ElementIndexColorTheme } from './color/element-index';
import { CarbohydrateSymbolColorTheme } from './color/carbohydrate-symbol';
import { ChainIdColorTheme } from './color/chain-id';
import { ElementSymbolColorTheme } from './color/element-symbol';
import { UnitIndexColorTheme } from './color/unit-index';
import { UniformColorTheme } from './color/uniform';

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

export interface ColorThemeProps {
    name: 'element-index' | 'chain-id'| 'unit-index' | 'uniform' | 'carbohydrate-symbol' | 'element-symbol'
    domain?: [number, number]
    value?: Color
    structure?: Structure
}

export const ColorThemeInfo = {
    'element-index': {},
    'carbohydrate-symbol': {},
    'chain-id': {},
    'element-symbol': {},
    'unit-index': {},
    'uniform': {}
}
export type ColorThemeName = keyof typeof ColorThemeInfo
export const ColorThemeNames = Object.keys(ColorThemeInfo)