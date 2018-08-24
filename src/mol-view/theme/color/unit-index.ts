/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { Unit, StructureElement, Link } from 'mol-model/structure';
import { ColorTheme, ColorThemeProps, LocationColor } from '../color';

const DefaultColor = Color(0xCCCCCC)

export function UnitIndexColorTheme(props: ColorThemeProps): ColorTheme {
    let colorFn: LocationColor

    if (props.structure) {
        const { units } = props.structure
        const unitCount = units.length

        const scale = ColorScale.create({ domain: [ 0, unitCount ] })

        colorFn = (location: Location): Color => {
            if (StructureElement.isLocation(location)) {
                return scale.color(Unit.findUnitById(location.unit.id, units))
            } else if (Link.isLocation(location)) {
                return scale.color(Unit.findUnitById(location.aUnit.id, units))
            }
            return DefaultColor
        }
    } else {
        colorFn = () => DefaultColor
    }

    return {
        kind: 'instance',
        color: colorFn
    }
}