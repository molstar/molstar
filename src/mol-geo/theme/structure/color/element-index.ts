/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale } from 'mol-util/color';
import { createElementInstanceColor, ColorData } from '../../../util/color-data';
import { LocationIterator, LocationValue } from '../../../representation/structure/visual/util/location-iterator';

export function elementIndexColorData(locationIt: LocationIterator, colorData?: ColorData) {
    const { elementCount, instanceCount } = locationIt

    const domain = [ 0, instanceCount * elementCount - 1 ]
    const scale = ColorScale.create({ domain })
    return createElementInstanceColor(
        locationIt,
        (value: LocationValue) => scale.color(value.instanceIndex * elementCount + value.elementIndex),
        colorData
    )
}