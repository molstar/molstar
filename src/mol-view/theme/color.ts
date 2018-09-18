/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Structure } from 'mol-model/structure';
import { Location } from 'mol-model/location';
import { ColorType } from 'mol-geo/geometry/color-data';

import { ElementIndexColorTheme } from './color/element-index';
import { CarbohydrateSymbolColorTheme } from './color/carbohydrate-symbol';
import { ChainIdColorTheme } from './color/chain-id';
import { ElementSymbolColorTheme } from './color/element-symbol';
import { UnitIndexColorTheme } from './color/unit-index';
import { UniformColorTheme } from './color/uniform';
import { CrossLinkColorTheme } from './color/cross-link';
import { ShapeGroupColorTheme } from './color/shape-group';
import { CustomColorTheme } from './color/custom';
import { ResidueNameColorTheme } from './color/residue-name';
import { SequenceIdColorTheme } from './color/sequence-id';
import { SecondaryStructureColorTheme } from './color/secondary-structure';

export type LocationColor = (location: Location, isSecondary: boolean) => Color

export interface ScaleLegend {
    kind: 'scale-legend'
    minLabel: string,
    maxLabel: string,
    colors: Color[]
}
export function ScaleLegend(minLabel: string, maxLabel: string, colors: Color[]): ScaleLegend {
    return { kind: 'scale-legend', minLabel, maxLabel, colors }
}

export interface TableLegend {
    kind: 'table-legend'
    table: [ string, Color ][]
}
export function TableLegend(table: [ string, Color ][]): TableLegend {
    return { kind: 'table-legend', table }
}

export interface ColorTheme {
    granularity: ColorType
    color: LocationColor
    description?: string
    legend?: ScaleLegend | TableLegend
}

export function ColorTheme(props: ColorThemeProps): ColorTheme {
    switch (props.name) {
        case 'carbohydrate-symbol': return CarbohydrateSymbolColorTheme(props)
        case 'chain-id': return ChainIdColorTheme(props)
        case 'cross-link': return CrossLinkColorTheme(props)
        case 'custom': return CustomColorTheme(props)
        case 'element-index': return ElementIndexColorTheme(props)
        case 'element-symbol': return ElementSymbolColorTheme(props)
        case 'residue-name': return ResidueNameColorTheme(props)
        case 'secondary-structure': return SecondaryStructureColorTheme(props)
        case 'sequence-id': return SequenceIdColorTheme(props)
        case 'shape-group': return ShapeGroupColorTheme(props)
        case 'unit-index': return UnitIndexColorTheme(props)
        case 'uniform': return UniformColorTheme(props)
    }
}

export interface ColorThemeProps {
    name: ColorThemeName
    domain?: [number, number]
    value?: Color
    structure?: Structure
    color?: LocationColor
    granularity?: ColorType,
    description?: string,
    legend?: ScaleLegend | TableLegend
}

export const ColorThemeInfo = {
    'carbohydrate-symbol': {},
    'chain-id': {},
    'cross-link': {},
    'custom': {},
    'element-index': {},
    'element-symbol': {},
    'residue-name': {},
    'secondary-structure': {},
    'sequence-id': {},
    'shape-group': {},
    'unit-index': {},
    'uniform': {},
}
export type ColorThemeName = keyof typeof ColorThemeInfo
export const ColorThemeNames = Object.keys(ColorThemeInfo)