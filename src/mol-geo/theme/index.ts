/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';

export interface UniformColorTheme {
    name: 'uniform'
    value: Color
}

export interface ScaleColorTheme {
    name:  'atom-index' | 'chain-id'| 'instance-index'
    domain?: [number, number]
}

export interface TableColorTheme {
    name:  'carbohydrate-symbol' | 'element-symbol'
}

export type ColorTheme = UniformColorTheme | ScaleColorTheme | TableColorTheme

export interface UniformSizeTheme {
    name: 'uniform',
    value: number
}

export interface ScaleSizeTheme {
    name: 'physical' // van-der-Waals for atoms, given radius for coarse spheres
    factor?: number // scaling factor
}

export type SizeTheme = UniformSizeTheme | ScaleSizeTheme