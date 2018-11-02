/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { StructureElement, Link } from 'mol-model/structure';
import { ColorTheme, ColorThemeProps, LocationColor } from '../color';

const DefaultColor = Color(0xCCCCCC)
const Description = 'Gives every polymer a unique color based on the position (index) of the polymer in the list of polymers in the structure.'

export function PolymerIndexColorTheme(props: ColorThemeProps): ColorTheme {
    let color: LocationColor
    let scale: ColorScale | undefined = undefined

    if (props.structure) {
        const { units } = props.structure
        let polymerCount = 0
        for (let i = 0, il = units.length; i <il; ++i) {
            if (units[i].polymerElements.length > 0) ++polymerCount
        }
        scale = ColorScale.create({ domain: [ 0, polymerCount - 1 ] })
        const unitIdColor = new Map<number, Color>()
        for (let i = 0, j = 0, il = units.length; i <il; ++i) {
            if (units[i].polymerElements.length > 0) {
                unitIdColor.set(units[i].id, scale.color(j))
                ++j
            }
        }

        color = (location: Location): Color => {
            let color: Color | undefined
            if (StructureElement.isLocation(location)) {
                color = unitIdColor.get(location.unit.id)
            } else if (Link.isLocation(location)) {
                color = unitIdColor.get(location.aUnit.id)
            }
            return color !== undefined ? color : DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return {
        granularity: 'instance',
        color,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}