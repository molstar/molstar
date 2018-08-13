/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Structure } from 'mol-model/structure';

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

//

export interface SizeThemeProps {
    name: 'physical' | 'uniform'
    value?: number
    factor?: number
    structure?: Structure
}

export const SizeThemeInfo = {
    'physical': {},
    'uniform': {}
}
export type SizeThemeName = keyof typeof SizeThemeInfo
export const SizeThemeNames = Object.keys(SizeThemeInfo)