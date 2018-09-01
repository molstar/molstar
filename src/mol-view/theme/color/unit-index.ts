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

export function UnitIndexColorTheme(props: ColorThemeProps): ColorTheme {
    let color: LocationColor

    if (props.structure) {
        const { units } = props.structure
        const scale = ColorScale.create({ domain: [ 0, units.length - 1 ] })
        const unitIdColor = new Map<number, Color>()
        for (let i = 0, il = units.length; i <il; ++i) {
            unitIdColor.set(units[i].id, scale.color(units[i].id))
        }

        color = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                return unitIdColor.get(location.unit.id)!
            } else if (Link.isLocation(location)) {
                return unitIdColor.get(location.aUnit.id)!
            }
            return DefaultColor
        }
    } else {
        color = () => DefaultColor
    }

    return { kind: 'instance', color }
}