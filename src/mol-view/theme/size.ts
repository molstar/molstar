/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { SizeType, LocationSize } from 'mol-geo/util/size-data';

import { PhysicalSizeTheme } from './size/physical';
import { UniformSizeTheme } from './size/uniform';

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