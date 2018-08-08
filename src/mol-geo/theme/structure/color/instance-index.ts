/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorScale } from 'mol-util/color';
import { createInstanceColor, ColorData } from '../../../util/color-data';
import { LocationIterator, LocationValue } from '../../../representation/structure/visual/util/location-iterator';

export function instanceIndexColorData(locationIt: LocationIterator, colorData?: ColorData) {
    const { instanceCount } = locationIt
    const domain = [ 0, instanceCount - 1 ]
    const scale = ColorScale.create({ domain })

    function colorFn(locationValue: LocationValue) {
        return scale.color(locationValue.instanceIndex)
    }

    return createInstanceColor(locationIt, colorFn, colorData)
}