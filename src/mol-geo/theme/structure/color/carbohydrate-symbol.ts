/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureProperties, StructureElement, Link } from 'mol-model/structure';

import { ColorData, createElementColor } from '../../../util/color-data';
import { ColorScale, Color } from 'mol-util/color';
import { LocationIterator, LocationValue } from '../../../representation/structure/visual/util/location-iterator';

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.label_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
}

export function carbohydrateSymbolColorData(locationIt: LocationIterator, colorData?: ColorData) {
    const l = StructureElement.create()

    function colorFn(locationValue: LocationValue): Color {
        const { location } = locationValue
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
        return 0
    }

    return createElementColor(locationIt, colorFn, colorData)
}