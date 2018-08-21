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
import { CrossLinkColorTheme } from './color/cross-link';
import { ShapeGroupColorTheme } from './color/shape-group';

export interface ColorTheme {
    kind: ColorType
    color: LocationColor
}

export function ColorTheme(props: ColorThemeProps): ColorTheme {
    switch (props.name) {
        case 'element-index': return ElementIndexColorTheme(props)
        case 'carbohydrate-symbol': return CarbohydrateSymbolColorTheme(props)
        case 'cross-link': return CrossLinkColorTheme(props)
        case 'chain-id': return ChainIdColorTheme(props)
        case 'element-symbol': return ElementSymbolColorTheme(props)
        case 'unit-index': return UnitIndexColorTheme(props)
        case 'uniform': return UniformColorTheme(props)
        case 'shape-group': return ShapeGroupColorTheme(props)
    }
}

export interface ColorThemeProps {
    name: ColorThemeName
    domain?: [number, number]
    value?: Color
    structure?: Structure
}

export const ColorThemeInfo = {
    'element-index': {},
    'carbohydrate-symbol': {},
    'cross-link': {},
    'chain-id': {},
    'element-symbol': {},
    'unit-index': {},
    'uniform': {},
    'shape-group': {}
}
export type ColorThemeName = keyof typeof ColorThemeInfo
export const ColorThemeNames = Object.keys(ColorThemeInfo)