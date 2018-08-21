/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, StructureElement, Link } from 'mol-model/structure';

import { ColorScale, Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorThemeProps, ColorTheme } from '../color';

const DefaultColor = 0xCCCCCC as Color

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
}

export function ChainIdColorTheme(props: ColorThemeProps): ColorTheme {
    const l = StructureElement.create()

    function colorFn(location: Location): Color {
        if (StructureElement.isLocation(location)) {
            const map = location.unit.model.properties.asymIdSerialMap
            const scale = ColorScale.create({ domain: [ 0, map.size - 1 ] })
            const asym_id = getAsymId(location.unit)
            return scale.color(map.get(asym_id(location)) || 0)
        } else if (Link.isLocation(location)) {
            const map = location.aUnit.model.properties.asymIdSerialMap
            const scale = ColorScale.create({ domain: [ 0, map.size - 1 ] })
            const asym_id = getAsymId(location.aUnit)
            l.unit = location.aUnit
            l.element = location.aUnit.elements[location.aIndex]
            return scale.color(map.get(asym_id(l)) || 0)
        }
        return DefaultColor
    }

    return {
        kind: 'group',
        color: colorFn
    }
}