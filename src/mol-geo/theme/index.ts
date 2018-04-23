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
    name:  'atom-index' | 'chain-id' | 'element-symbol' | 'instance-index'
    domain?: [number, number]
}

export type ColorTheme = UniformColorTheme | ScaleColorTheme

export interface UniformSizeTheme {
    name: 'uniform',
    value: number
}

export interface ScaleSizeTheme {
    name: 'vdw'
}

export type SizeTheme = UniformSizeTheme | ScaleSizeTheme