/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { StructureElement, Link, Unit } from 'mol-model/structure';
import { OrderedSet } from 'mol-data/int';
import { ColorThemeProps, ColorTheme, LocationColor } from '../color';

const DefaultColor = Color(0xCCCCCC)

export function ElementIndexColorTheme(props: ColorThemeProps): ColorTheme {
    let colorFn: LocationColor

    if (props.structure) {
        const { units } = props.structure
        const unitCount = units.length
        const cummulativeElementCount = new Map<number, number>()

        let elementCount = 0
        for (let i = 0; i < unitCount; ++i) {
            cummulativeElementCount.set(i, elementCount)
            elementCount += units[i].elements.length
        }
        const scale = ColorScale.create({ domain: [ 0, elementCount ] })

        colorFn = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const unitIndex = Unit.findUnitById(location.unit.id, units)
                const unitElementIndex = OrderedSet.findPredecessorIndex(location.unit.elements, location.element)
                return scale.color(cummulativeElementCount.get(unitIndex) || 0 + unitElementIndex)
            } else if (Link.isLocation(location)) {
                const unitId = Unit.findUnitById(location.aUnit.id, units)
                return scale.color(cummulativeElementCount.get(unitId) || 0 + location.aIndex)
            }
            return DefaultColor
        }
    } else {
        colorFn = () => DefaultColor
    }

    return {
        kind: 'groupInstance',
        color: colorFn
    }
}