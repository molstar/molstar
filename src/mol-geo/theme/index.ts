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
    name: 'instance-id' | 'element-symbol' | 'atom-id'
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