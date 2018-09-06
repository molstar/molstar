/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { StructureElement, Link } from 'mol-model/structure';
import { OrderedSet } from 'mol-data/int';
import { ColorThemeProps, ColorTheme, LocationColor } from '../color';

const DefaultColor = Color(0xCCCCCC)

export function ElementIndexColorTheme(props: ColorThemeProps): ColorTheme {
    let color: LocationColor

    if (props.structure) {
        const { units } = props.structure
        const unitCount = units.length
        const cummulativeElementCount = new Map<number, number>()
        const unitIdIndex = new Map<number, number>()

        let elementCount = 0
        for (let i = 0; i < unitCount; ++i) {
            cummulativeElementCount.set(i, elementCount)
            elementCount += units[i].elements.length
            unitIdIndex.set(units[i].id, i)
        }
        const scale = ColorScale.create({ domain: [ 0, elementCount - 1 ] })

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.unit.id)!
                const unitElementIndex = OrderedSet.findPredecessorIndex(location.unit.elements, location.element)
                return scale.color(cummulativeElementCount.get(unitIndex)! + unitElementIndex)
            } else if (Link.isLocation(location)) {
                const unitIndex = unitIdIndex.get(location.aUnit.id)!
                return scale.color(cummulativeElementCount.get(unitIndex)! + location.aIndex)
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return { granularity: 'groupInstance', color }
}