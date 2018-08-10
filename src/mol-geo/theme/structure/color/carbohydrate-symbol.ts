/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Link, ElementIndex, Unit } from 'mol-model/structure';

import { ColorData, createElementColor } from '../../../util/color-data';
import { Color } from 'mol-util/color';
import { LocationIterator, LocationValue } from '../../../representation/structure/visual/util/location-iterator';
import { ColorTheme } from '../..';
import { SaccharideColors } from 'mol-model/structure/structure/carbohydrates/constants';

const DefaultColor = 0xCCCCCC;

export function carbohydrateSymbolColorData(locationIt: LocationIterator, props: ColorTheme, colorData?: ColorData) {
    let colorFn: (locationValue: LocationValue) => Color

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

        colorFn = (locationValue: LocationValue) => {
            const { location: l } = locationValue
            if (locationValue.isSecondary) {
                return SaccharideColors.Secondary
            } else {
                if (StructureElement.isLocation(l)) {
                    return getColor(l.unit, l.element)
                } else if (Link.isLocation(l)) {
                    return getColor(l.aUnit, l.aUnit.elements[l.aIndex])
                }
            }
            return DefaultColor
        }
    } else {
        colorFn = () => DefaultColor
    }

    return createElementColor(locationIt, colorFn, colorData)
}