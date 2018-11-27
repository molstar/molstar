/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Link, ElementIndex, Unit } from 'mol-model/structure';

import { SaccharideColors, MonosaccharidesColorTable } from 'mol-model/structure/structure/carbohydrates/constants';
import { Location } from 'mol-model/location';
import { ColorTheme, LocationColor } from '../color';
import { Color } from 'mol-util/color';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { ThemeDataContext } from '../theme';
import { TableLegend } from 'mol-util/color/tables';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Assigns colors according to the Symbol Nomenclature for Glycans (SNFG).'

// name: ColorThemeName
// domain?: [number, number]
// value?: Color
// list?: Color[]
// map?: ColorMap<any>

export const CarbohydrateSymbolColorThemeParams = {
    // domain: PD.Interval('Color Domain', '', [0, 1]),
    // value: PD.Color('Color Value', '', DefaultColor),
}
export type CarbohydrateSymbolColorThemeParams = typeof CarbohydrateSymbolColorThemeParams
export function getCarbohydrateSymbolColorThemeParams(ctx: ThemeDataContext) {
    return CarbohydrateSymbolColorThemeParams // TODO return copy
}

export function CarbohydrateSymbolColorTheme(ctx: ThemeDataContext, props: PD.Values<CarbohydrateSymbolColorThemeParams>): ColorTheme<CarbohydrateSymbolColorThemeParams> {
    let color: LocationColor

    if (ctx.structure) {
        const { elements, getElementIndex, getAnomericCarbon } = ctx.structure.carbohydrates

        const getColor = (unit: Unit, index: ElementIndex) => {
            const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[index]
            const anomericCarbon = getAnomericCarbon(unit, residueIndex)
            if (anomericCarbon !== undefined) {
                const idx = getElementIndex(unit, anomericCarbon)
                if (idx !== undefined) return elements[idx].component.color
            }
            return DefaultColor
        }

        color = (location: Location, isSecondary: boolean) => {
            if (isSecondary) {
                return SaccharideColors.Secondary
            } else {
                if (StructureElement.isLocation(location)) {
                    return getColor(location.unit, location.element)
                } else if (Link.isLocation(location)) {
                    return getColor(location.aUnit, location.aUnit.elements[location.aIndex])
                }
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        factory: CarbohydrateSymbolColorTheme,
        granularity: 'group',
        color: color,
        props: props,
        description: Description,
        legend: TableLegend(MonosaccharidesColorTable)
    }
}

export const CarbohydrateSymbolColorThemeProvider: ColorTheme.Provider<CarbohydrateSymbolColorThemeParams> = {
    label: 'Carbohydrate Symbol',
    factory: CarbohydrateSymbolColorTheme,
    getParams: getCarbohydrateSymbolColorThemeParams
}