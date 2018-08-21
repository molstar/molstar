/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Link, ElementIndex, Unit } from 'mol-model/structure';

import { SaccharideColors } from 'mol-model/structure/structure/carbohydrates/constants';
import { Location } from 'mol-model/location';
import { ColorThemeProps, ColorTheme } from '../color';
import { LocationColor } from 'mol-geo/util/color-data';
import { Color } from 'mol-util/color';

const DefaultColor = 0xCCCCCC as Color

export function CarbohydrateSymbolColorTheme(props: ColorThemeProps): ColorTheme {
    let colorFn: LocationColor

    if (props.structure) {
        const { elements, getElementIndex, getAnomericCarbon } = props.structure.carbohydrates

        const getColor = (unit: Unit, index: ElementIndex) => {
            const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[index]
            const anomericCarbon = getAnomericCarbon(unit, residueIndex)
            if (anomericCarbon !== undefined) {
                const idx = getElementIndex(unit, anomericCarbon)
                if (idx !== undefined) return elements[idx].component.color
            }
            return DefaultColor
        }

        colorFn = (location: Location, isSecondary: boolean) => {
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
        colorFn = () => DefaultColor
    }

    return {
        kind: 'group',
        color: colorFn
    }
}